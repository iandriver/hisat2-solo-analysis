//! External-memory primitives for the doubling loop.
//!
//! The whole-genome graph build needs 625 GiB in C++ because `lateGeneration`
//! keeps three live `PathNode` arrays resident. Nothing about the algorithm
//! requires that: `createNewNodesMaker` is an equi-join
//! (`past_nodes.to = from_table.from`) implemented as a random-access probe, and
//! an equi-join can always be run as sort-merge over two sequential streams.
//!
//! Records are fixed size, so runs are plain files and merging is a k-way heap.
//! Everything here is written so the memory budget is an explicit parameter and
//! can be set absurdly low in tests, which is how the external path gets
//! exercised at 200 bp instead of 3.1 Gb.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

/// A path node, packed to 18 bytes rather than the 32 a four-`u64` struct costs.
///
/// The fields have genuinely different ranges and pretending otherwise is what
/// makes the C++ node 32 bytes: `from`/`to` are reference-graph node ids and sit
/// near 3.05e9 at whole-genome scale, while the ranks reach 5.92e9. So 32 bits
/// for the former and 40 for the latter, which is 18 bytes and takes the
/// three-array peak from 568 GB to ~320 GB before any spilling.
///
/// `to == u32::MAX` marks a node that is already fully sorted, matching the
/// `setSorted()` sentinel. Real `to` values cannot reach it at 3.05e9, but the
/// width is the thing to widen first for a pangenome.
pub const REC: usize = 18;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Rec {
    pub from: u32,
    pub to: u32,
    pub k0: u64, // stored in 5 bytes
    pub k1: u64, // stored in 5 bytes
}

pub const SORTED: u32 = u32::MAX;

impl Rec {
    pub fn is_sorted(&self) -> bool { self.to == SORTED }
    pub fn key(&self) -> (u64, u64) { (self.k0, self.k1) }

    #[inline]
    pub fn write(&self, b: &mut [u8]) {
        debug_assert!(self.k0 < (1u64 << 40) && self.k1 < (1u64 << 40),
                      "rank exceeds the 40-bit field");
        b[0..4].copy_from_slice(&self.from.to_le_bytes());
        b[4..8].copy_from_slice(&self.to.to_le_bytes());
        b[8..13].copy_from_slice(&self.k0.to_le_bytes()[0..5]);
        b[13..18].copy_from_slice(&self.k1.to_le_bytes()[0..5]);
    }
    #[inline]
    pub fn read(b: &[u8]) -> Rec {
        let mut k0 = [0u8; 8]; k0[0..5].copy_from_slice(&b[8..13]);
        let mut k1 = [0u8; 8]; k1[0..5].copy_from_slice(&b[13..18]);
        Rec {
            from: u32::from_le_bytes(b[0..4].try_into().unwrap()),
            to:   u32::from_le_bytes(b[4..8].try_into().unwrap()),
            k0: u64::from_le_bytes(k0),
            k1: u64::from_le_bytes(k1),
        }
    }
}

/// Append-only record file.
///
/// Block-based rather than a `BufWriter` of 18-byte writes: at whole-genome
/// scale something like 10^9 records pass through here per generation, and
/// `write_all` on an 18-byte slice pays its bounds checks and copy on every one
/// of them. Appending into a block and flushing whole blocks turns that into a
/// pointer bump.
pub struct RecWriter { f: File, buf: Vec<u8>, pub n: u64 }

const BLOCK: usize = REC * 4096;   // ~72 KB

impl RecWriter {
    pub fn create(p: &Path) -> std::io::Result<RecWriter> {
        Ok(RecWriter { f: File::create(p)?, buf: Vec::with_capacity(BLOCK), n: 0 })
    }
    pub fn from_file(f: File) -> RecWriter {
        RecWriter { f, buf: Vec::with_capacity(BLOCK), n: 0 }
    }
    #[inline]
    pub fn push(&mut self, r: Rec) -> std::io::Result<()> {
        if self.buf.len() + REC > BLOCK { self.spill()?; }
        let mut rec = [0u8; REC];
        r.write(&mut rec);
        self.buf.extend_from_slice(&rec);
        self.n += 1;
        Ok(())
    }
    #[cold]
    #[inline(never)]
    fn spill(&mut self) -> std::io::Result<()> {
        self.f.write_all(&self.buf)?;
        self.buf.clear();
        Ok(())
    }
    pub fn finish(mut self) -> std::io::Result<u64> {
        if !self.buf.is_empty() { self.f.write_all(&self.buf)?; self.buf.clear(); }
        Ok(self.n)
    }
}

/// Sequential record reader.
///
/// Also block-based: `read_exact` of 18 bytes through a `BufReader` was costing
/// more per record than the sort comparison it feeds.
pub struct RecReader { f: File, buf: Vec<u8>, pos: usize, filled: usize }

impl RecReader {
    pub fn open(p: &Path) -> std::io::Result<RecReader> {
        Self::open_buf(p, 1 << 20)
    }
    /// The k-way merge opens every run at once, so a fixed per-run buffer makes
    /// merge memory proportional to the RUN COUNT -- that is, to n/budget --
    /// while the sort phase stays inside its budget. At 20 Mb that was ~305 runs
    /// at 1 MB each: 305 MB of buffers for a 1.2 MB sort budget.
    pub fn open_buf(p: &Path, cap: usize) -> std::io::Result<RecReader> {
        let cap = (cap / REC).max(64) * REC;
        Ok(RecReader { f: File::open(p)?, buf: vec![0u8; cap], pos: 0, filled: 0 })
    }
    #[inline]
    pub fn next(&mut self) -> std::io::Result<Option<Rec>> {
        if self.pos + REC > self.filled { return self.refill(); }
        let r = Rec::read(&self.buf[self.pos..]);
        self.pos += REC;
        Ok(Some(r))
    }
    #[cold]
    #[inline(never)]
    fn refill(&mut self) -> std::io::Result<Option<Rec>> {
        // `read` may return short at any point, not only at EOF, so a record can
        // straddle two of them -- carry whatever is left to the front
        let rem = self.filled - self.pos;
        self.buf.copy_within(self.pos..self.filled, 0);
        self.filled = rem;
        self.pos = 0;
        while self.filled < REC {
            let n = self.f.read(&mut self.buf[self.filled..])?;
            if n == 0 { return Ok(None); }
            self.filled += n;
        }
        let r = Rec::read(&self.buf[0..]);
        self.pos = REC;
        Ok(Some(r))
    }
}

/// A record file written in numbered segments, so that a reader making one
/// forward pass can free the space behind it.
///
/// Almost every file in the doubling loop is written once and read once. Held as
/// a single file, each is a full copy of the node list sitting on disk for the
/// whole of the pass that consumes it -- three of them alive at the join is what
/// made peak scratch 7.4 copies of an 18-byte node. In segments, a consuming
/// reader hands the space back as it goes and the same pass costs one copy that
/// shrinks while its output grows.
///
/// The segment size is a scratch-vs-inode tradeoff and nothing else: records
/// never cross a boundary, and the concatenation is byte for byte what a single
/// file would have held.
pub const SEG_RECORDS_DEFAULT: usize = 4 << 20;   // 4M records = 72 MB

/// `HT2_SEG` overrides it, so a 200 bp graph can be made to produce dozens of
/// segments and exercise the boundary handling the way `budget` does for runs.
pub fn seg_records() -> usize {
    use std::sync::OnceLock;
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| std::env::var("HT2_SEG").ok()
                        .and_then(|x| x.parse().ok())
                        .filter(|&x| x > 0)
                        .unwrap_or(SEG_RECORDS_DEFAULT))
}

pub fn seg_path(base: &Path, i: usize) -> PathBuf {
    let mut p = base.as_os_str().to_owned();
    p.push(format!(".s{:05}", i));
    PathBuf::from(p)
}

/// Every segment index that still exists for `base`, ascending.
///
/// This lists the directory rather than probing 0, 1, 2, ... until a miss: a
/// consuming reader deletes a PREFIX of the segments, so a base that was read
/// part-way has a hole at the front and probing from zero would stop
/// immediately and leak everything behind it.
fn seg_indices(base: &Path) -> Vec<usize> {
    let dir = match base.parent() { Some(d) if !d.as_os_str().is_empty() => d, _ => Path::new(".") };
    let name = match base.file_name().and_then(|n| n.to_str()) { Some(n) => n, None => return Vec::new() };
    let prefix = format!("{name}.s");
    let mut out: Vec<usize> = Vec::new();
    if let Ok(rd) = std::fs::read_dir(dir) {
        for e in rd.flatten() {
            if let Some(f) = e.file_name().to_str() {
                if let Some(rest) = f.strip_prefix(&prefix) {
                    if let Ok(i) = rest.parse::<usize>() { out.push(i); }
                }
            }
        }
    }
    out.sort_unstable();
    out
}

/// Delete every segment of `base`, including any left behind a hole.
///
/// Costs one directory listing, so it is for the handful of long-lived bases a
/// consuming reader may have left part-way through. Everything else uses
/// `seg_truncate`.
pub fn seg_remove(base: &Path) {
    if seg_path(base, 0).exists() {
        // fast path first: a fully-written base has no hole
        let mut i = 0usize;
        while std::fs::remove_file(seg_path(base, i)).is_ok() { i += 1; }
        if !seg_path(base, i + 1).exists() { return; }
    }
    for i in seg_indices(base) { let _ = std::fs::remove_file(seg_path(base, i)); }
}

/// Delete a base's segments from the front, stopping at the first gap.
///
/// No directory listing, so it is cheap enough to call before every run file --
/// and correct for them, because a run is either untouched or fully consumed.
pub fn seg_truncate(base: &Path) {
    let mut i = 0usize;
    while std::fs::remove_file(seg_path(base, i)).is_ok() { i += 1; }
}

/// Total records across every segment.
pub fn seg_len(base: &Path) -> u64 {
    seg_indices(base).into_iter()
        .filter_map(|i| std::fs::metadata(seg_path(base, i)).ok())
        .map(|m| m.len() / REC as u64)
        .sum()
}

/// Append every segment of `src` after the segments already in `dst`.
/// Metadata only -- no data is copied, and a short segment in the middle is
/// harmless because a reader only cares about the order.
pub fn seg_append(dst: &Path, src: &Path) -> std::io::Result<()> {
    let mut next = seg_indices(dst).last().map(|&i| i + 1).unwrap_or(0);
    for i in seg_indices(src) {
        std::fs::rename(seg_path(src, i), seg_path(dst, next))?;
        next += 1;
    }
    Ok(())
}

/// Move every segment of `from` onto `to`. Metadata only -- no data is copied.
pub fn seg_rename(from: &Path, to: &Path) -> std::io::Result<()> {
    seg_remove(to);
    for i in seg_indices(from) { std::fs::rename(seg_path(from, i), seg_path(to, i))?; }
    Ok(())
}

pub struct SegWriter { base: PathBuf, seg: usize, in_seg: usize, cap: usize, w: Option<RecWriter>, pub n: u64 }

impl SegWriter {
    pub fn create(base: &Path) -> std::io::Result<SegWriter> {
        Self::create_cap(base, seg_records())
    }
    /// A run file is only `budget` records long, so the default segment size can
    /// be larger than the whole run -- which silently turns segmentation off for
    /// exactly the files a merge needs to shrink. Callers that write short files
    /// pass their own cap.
    pub fn create_cap(base: &Path, cap: usize) -> std::io::Result<SegWriter> {
        seg_truncate(base);
        Ok(SegWriter { base: base.to_path_buf(), seg: 0, in_seg: 0, cap: cap.max(1),
                       w: Some(RecWriter::create(&seg_path(base, 0))?), n: 0 })
    }
    /// Continue an existing base, so a second pass can append to what a first
    /// one wrote without copying it.
    pub fn append(base: &Path) -> std::io::Result<SegWriter> {
        let idx = seg_indices(base);
        match idx.last() {
            None => SegWriter::create(base),
            Some(&last) => {
                let p = seg_path(base, last);
                let in_seg = std::fs::metadata(&p)?.len() as usize / REC;
                let f = std::fs::OpenOptions::new().append(true).open(&p)?;
                Ok(SegWriter { base: base.to_path_buf(), seg: last, in_seg,
                               cap: seg_records(), w: Some(RecWriter::from_file(f)), n: 0 })
            }
        }
    }
    pub fn push(&mut self, r: Rec) -> std::io::Result<()> {
        if self.in_seg == self.cap {
            self.w.take().unwrap().finish()?;
            self.seg += 1;
            self.in_seg = 0;
            self.w = Some(RecWriter::create(&seg_path(&self.base, self.seg))?);
        }
        self.w.as_mut().unwrap().push(r)?;
        self.in_seg += 1;
        self.n += 1;
        Ok(())
    }
    pub fn finish(mut self) -> std::io::Result<u64> {
        self.w.take().unwrap().finish()?;
        Ok(self.n)
    }
}

pub struct SegReader { base: PathBuf, seg: usize, r: Option<RecReader>, consume: bool, cap: usize }

impl SegReader {
    /// `consume` unlinks each segment as it is exhausted. The reader is closed
    /// first, so the space comes back immediately rather than at process exit.
    pub fn open(base: &Path, consume: bool) -> std::io::Result<SegReader> {
        Self::open_buf(base, consume, 1 << 20)
    }
    pub fn open_buf(base: &Path, consume: bool, cap: usize) -> std::io::Result<SegReader> {
        let p = seg_path(base, 0);
        let r = if p.exists() { Some(RecReader::open_buf(&p, cap)?) } else { None };
        Ok(SegReader { base: base.to_path_buf(), seg: 0, r, consume, cap })
    }
    /// One record. The common path is a single call into the current segment;
    /// crossing a boundary happens once every few million records and is kept
    /// out of line so it does not weigh on the hot one.
    #[inline]
    pub fn next(&mut self) -> std::io::Result<Option<Rec>> {
        if let Some(r) = self.r.as_mut() {
            if let Some(x) = r.next()? { return Ok(Some(x)); }
        }
        self.advance()
    }

    #[cold]
    #[inline(never)]
    fn advance(&mut self) -> std::io::Result<Option<Rec>> {
        loop {
            if self.r.is_none() { return Ok(None); }
            self.r = None;   // close before unlinking
            if self.consume { let _ = std::fs::remove_file(seg_path(&self.base, self.seg)); }
            self.seg += 1;
            let p = seg_path(&self.base, self.seg);
            if !p.exists() { return Ok(None); }
            self.r = Some(RecReader::open_buf(&p, self.cap)?);
            if let Some(x) = self.r.as_mut().unwrap().next()? { return Ok(Some(x)); }
        }
    }
}

struct HeapItem { key: (u64, u64, u32), idx: usize, rec: Rec }
impl PartialEq for HeapItem { fn eq(&self, o: &Self) -> bool { self.key == o.key } }
impl Eq for HeapItem {}
impl Ord for HeapItem {
    // reversed: BinaryHeap is a max-heap and a merge wants the minimum
    fn cmp(&self, o: &Self) -> Ordering { o.key.cmp(&self.key).then(o.idx.cmp(&self.idx)) }
}
impl PartialOrd for HeapItem { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }

/// Which ordering an external sort should produce.
#[derive(Clone, Copy, PartialEq)]
pub enum By {
    To,
    From,
    Key,
    /// `generateEdges` order: label bucket, then the ranking the edge points at.
    /// Ties inside a bucket resolve by the order the edges were pushed, which is
    /// `sortEdgesFrom` order -- so `(from, to)` finishes the key rather than
    /// leaning on the sort being stable across a join.
    Row,
    /// `nextRow` order: a STABLE sort on ranking alone applied to a list already
    /// in `Row` order. `k1` carries `label * 2 + decremented` so the tie-break
    /// can reproduce that stability exactly, including the one case where two
    /// records reach the same ranking from different ones because the
    /// second-to-last node was dropped.
    Rank,
}

fn keyof(r: &Rec, by: By) -> (u64, u64, u32) {
    match by {
        By::To   => (r.to as u64, 0, 0),
        By::From => (r.from as u64, 0, 0),
        // `from` is part of the key, not decoration: mergeUpdateRank collapses a
        // run of equal-key nodes only when they share a `from`, and keeps the
        // FIRST of the run, so ties have to break deterministically or the
        // surviving node's `to` changes.
        By::Key  => (r.k0, r.k1, r.from),
        By::Row  => ((r.k1 << 40) | r.k0, r.from as u64, r.to),
        By::Rank => ((r.k0 << 4) | r.k1, r.from as u64, r.to),
    }
}

/// External merge sort. `budget` is the number of records held in RAM at once,
/// so a test can set it to 1000 and still exercise runs and a k-way merge.
///
/// Stability matters: `mergeUpdateRank`'s collapse rule and the GFM row
/// tie-break both depend on ties keeping their input order, so runs are sorted
/// with a stable sort and the merge breaks ties by run index.
///
/// `consume` deletes `src` as soon as the run phase has finished reading it.
/// Scratch is the binding resource at whole-genome scale and a source nobody
/// will read again is a full copy of the node file sitting on disk for the
/// length of the merge.
pub fn sort_external(src: &Path, dst: &Path, by: By, budget: usize, tmp: &Path,
                     consume: bool) -> std::io::Result<u64>
{
    // The run phase reads `src` exactly once, so a consuming read hands each
    // segment back as its records land in a run -- the source shrinks at the
    // same rate the runs grow instead of both being resident.
    let mut runs: Vec<PathBuf> = Vec::new();
    {
        let mut rd = SegReader::open(src, consume)?;
        let mut buf: Vec<Rec> = Vec::with_capacity(budget);
        loop {
            buf.clear();
            while buf.len() < budget {
                match rd.next()? { Some(r) => buf.push(r), None => break }
            }
            if buf.is_empty() { break; }
            buf.sort_by_key(|r| keyof(r, by));
            // Runs are segmented as well. Freeing a run only when it is fully
            // exhausted does not help when a merge drains every run in step:
            // eleven half-empty runs still occupy their full size on disk, which
            // is why `runs + dst` measured 1.9 copies rather than one.
            let p = tmp.join(format!("run{}.bin", runs.len()));
            let mut w = SegWriter::create_cap(&p, (budget / 8).max(1))?;
            for r in &buf { w.push(*r)?; }
            w.finish()?;
            runs.push(p);
        }
    }
    if consume { seg_remove(src); }
    if runs.is_empty() { SegWriter::create(dst)?.finish()?; return Ok(0); }

    // Merge in passes of at most MAX_FANIN runs. A whole-genome sort produces
    // thousands of runs -- 5.9e9 records at an 18 MB budget is ~5,900 -- and a
    // single k-way merge would need one open file and one buffer per run, which
    // hits the file-descriptor limit long before it hits the memory budget.
    const MAX_FANIN: usize = 64;
    let mut pass = 0usize;
    while runs.len() > MAX_FANIN {
        let mut next: Vec<PathBuf> = Vec::new();
        for (gi, group) in runs.chunks(MAX_FANIN).enumerate() {
            if group.len() == 1 { next.push(group[0].clone()); continue; }
            let out = tmp.join(format!("merge{pass}_{gi}.bin"));
            merge_runs(group, &out, by, budget)?;
            for p in group { seg_truncate(p); }
            next.push(out);
        }
        runs = next;
        pass += 1;
    }
    // The last pass writes the sorted output in segments, and frees each run as
    // it empties, so `runs + dst` stays at about one copy of the data.
    let n = merge_runs_seg(&runs, dst, by, budget)?;
    for p in runs { seg_truncate(&p); }
    Ok(n)
}

/// The final merge: segmented output, and runs freed as they empty.
fn merge_runs_seg(runs: &[PathBuf], dst: &Path, by: By, budget: usize) -> std::io::Result<u64> {
    let cap = ((budget * REC) / runs.len().max(1)).clamp(8 * 1024, 1 << 20);
    let mut rds: Vec<Option<SegReader>> = runs.iter()
        .map(|p| SegReader::open_buf(p, true, cap).map(Some)).collect::<Result<_,_>>()?;
    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    for i in 0..rds.len() {
        match rds[i].as_mut().unwrap().next()? {
            Some(r) => heap.push(HeapItem { key: keyof(&r, by), idx: i, rec: r }),
            None => { rds[i] = None; seg_truncate(&runs[i]); }
        }
    }
    let mut w = SegWriter::create(dst)?;
    while let Some(it) = heap.pop() {
        w.push(it.rec)?;
        let i = it.idx;
        match rds[i].as_mut().unwrap().next()? {
            Some(r) => heap.push(HeapItem { key: keyof(&r, by), idx: i, rec: r }),
            None => { rds[i] = None; seg_truncate(&runs[i]); }
        }
    }
    w.finish()
}

/// One k-way merge over `runs`, stable: ties break by run index, which is input
/// order because the runs were produced in order.
///
/// Each run is closed and unlinked the moment it is exhausted. Every record read
/// out of a run is written into `dst`, so freeing runs as they empty keeps
/// `runs + dst` at roughly one copy of the data for the whole merge instead of
/// two -- which is the difference between a sort costing 3x its input in scratch
/// and costing 2x.
fn merge_runs(runs: &[PathBuf], dst: &Path, by: By, budget: usize) -> std::io::Result<u64> {
    // Split one buffer budget across the runs instead of giving each its own.
    let cap = ((budget * REC) / runs.len().max(1)).clamp(8 * 1024, 1 << 20);
    let mut rds: Vec<Option<SegReader>> = runs.iter()
        .map(|p| SegReader::open_buf(p, true, cap).map(Some)).collect::<Result<_,_>>()?;
    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    for i in 0..rds.len() {
        match rds[i].as_mut().unwrap().next()? {
            Some(r) => heap.push(HeapItem { key: keyof(&r, by), idx: i, rec: r }),
            None => { rds[i] = None; seg_truncate(&runs[i]); }
        }
    }
    let mut w = SegWriter::create(dst)?;
    while let Some(it) = heap.pop() {
        w.push(it.rec)?;
        let i = it.idx;
        match rds[i].as_mut().unwrap().next()? {
            Some(r) => heap.push(HeapItem { key: keyof(&r, by), idx: i, rec: r }),
            // dropping the reader closes the file, so the unlink frees the space
            // now rather than at process exit
            None => { rds[i] = None; seg_truncate(&runs[i]); }
        }
    }
    w.finish()
}

/// External sort over 8-byte `(u32, u32)` pairs, for the edge list.
///
/// The edge file has to end up in `sortEdgesFrom` order because ties inside a
/// GFM label bucket resolve by edge order, so this cannot be skipped just
/// because the doubling re-sorts its own records.
pub fn sort_pairs_external(src: &Path, dst: &Path, budget: usize, tmp: &Path)
    -> std::io::Result<u64>
{
    let mut runs: Vec<PathBuf> = Vec::new();
    {
        let mut f = BufReader::with_capacity(1 << 20, File::open(src)?);
        let mut b = [0u8; 8];
        let mut v: Vec<(u32, u32)> = Vec::with_capacity(budget);
        loop {
            v.clear();
            while v.len() < budget {
                match f.read_exact(&mut b) {
                    Ok(()) => v.push((u32::from_le_bytes(b[0..4].try_into().unwrap()),
                                      u32::from_le_bytes(b[4..8].try_into().unwrap()))),
                    Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                    Err(e) => return Err(e),
                }
            }
            if v.is_empty() { break; }
            v.sort_unstable();
            let p = tmp.join(format!("erun{}.bin", runs.len()));
            let mut w = BufWriter::with_capacity(1 << 20, File::create(&p)?);
            for &(x, y) in &v { w.write_all(&x.to_le_bytes())?; w.write_all(&y.to_le_bytes())?; }
            w.flush()?;
            runs.push(p);
        }
    }
    let ecap = ((budget * 8) / runs.len().max(1)).clamp(8 * 1024, 1 << 18);
    let mut rds: Vec<BufReader<File>> = runs.iter()
        .map(|p| File::open(p).map(|f| BufReader::with_capacity(ecap, f)))
        .collect::<Result<_, _>>()?;
    // the raw edge file is never read again once the runs exist
    let _ = std::fs::remove_file(src);
    let mut heap: BinaryHeap<(std::cmp::Reverse<(u32, u32)>, usize)> = BinaryHeap::new();
    let mut rd8 = |r: &mut BufReader<File>| -> std::io::Result<Option<(u32, u32)>> {
        let mut b = [0u8; 8];
        match r.read_exact(&mut b) {
            Ok(()) => Ok(Some((u32::from_le_bytes(b[0..4].try_into().unwrap()),
                               u32::from_le_bytes(b[4..8].try_into().unwrap())))),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
            Err(e) => Err(e),
        }
    };
    for i in 0..rds.len() {
        if let Some(p) = rd8(&mut rds[i])? { heap.push((std::cmp::Reverse(p), i)); }
    }
    let mut w = BufWriter::with_capacity(1 << 20, File::create(dst)?);
    let mut n = 0u64;
    while let Some((std::cmp::Reverse((x, y)), i)) = heap.pop() {
        w.write_all(&x.to_le_bytes())?; w.write_all(&y.to_le_bytes())?; n += 1;
        if let Some(p) = rd8(&mut rds[i])? { heap.push((std::cmp::Reverse(p), i)); }
    }
    w.flush()?;
    for p in runs { let _ = std::fs::remove_file(p); }
    Ok(n)
}

/// Worker count: `HT2_THREADS`, else what the machine reports.
pub fn threads() -> usize {
    use std::sync::OnceLock;
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| std::env::var("HT2_THREADS").ok()
                        .and_then(|x| x.parse::<usize>().ok())
                        .filter(|&x| x > 0)
                        .unwrap_or_else(|| std::thread::available_parallelism()
                                            .map(|n| n.get()).unwrap_or(1)))
}
