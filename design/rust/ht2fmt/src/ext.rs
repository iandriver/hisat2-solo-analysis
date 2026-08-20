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

    pub fn write(&self, b: &mut [u8]) {
        debug_assert!(self.k0 < (1u64 << 40) && self.k1 < (1u64 << 40),
                      "rank exceeds the 40-bit field");
        b[0..4].copy_from_slice(&self.from.to_le_bytes());
        b[4..8].copy_from_slice(&self.to.to_le_bytes());
        b[8..13].copy_from_slice(&self.k0.to_le_bytes()[0..5]);
        b[13..18].copy_from_slice(&self.k1.to_le_bytes()[0..5]);
    }
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
pub struct RecWriter { w: BufWriter<File>, buf: [u8; REC], pub n: u64 }

impl RecWriter {
    pub fn create(p: &Path) -> std::io::Result<RecWriter> {
        Ok(RecWriter { w: BufWriter::with_capacity(1 << 20, File::create(p)?), buf: [0; REC], n: 0 })
    }
    pub fn push(&mut self, r: Rec) -> std::io::Result<()> {
        r.write(&mut self.buf);
        self.w.write_all(&self.buf)?;
        self.n += 1;
        Ok(())
    }
    pub fn finish(mut self) -> std::io::Result<u64> { self.w.flush()?; Ok(self.n) }
}

/// Sequential record reader.
pub struct RecReader { r: BufReader<File>, buf: [u8; REC] }

impl RecReader {
    pub fn open(p: &Path) -> std::io::Result<RecReader> {
        Self::open_buf(p, 1 << 20)
    }
    /// The k-way merge opens every run at once, so a fixed per-run buffer makes
    /// merge memory proportional to the RUN COUNT -- that is, to n/budget --
    /// while the sort phase stays inside its budget. At 20 Mb that was ~305 runs
    /// at 1 MB each: 305 MB of buffers for a 1.2 MB sort budget.
    pub fn open_buf(p: &Path, cap: usize) -> std::io::Result<RecReader> {
        Ok(RecReader { r: BufReader::with_capacity(cap, File::open(p)?), buf: [0; REC] })
    }
    pub fn next(&mut self) -> std::io::Result<Option<Rec>> {
        match self.r.read_exact(&mut self.buf) {
            Ok(()) => Ok(Some(Rec::read(&self.buf))),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
            Err(e) => Err(e),
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
pub fn sort_external(src: &Path, dst: &Path, by: By, budget: usize, tmp: &Path)
    -> std::io::Result<u64>
{
    let mut runs: Vec<PathBuf> = Vec::new();
    {
        let mut rd = RecReader::open(src)?;
        let mut buf: Vec<Rec> = Vec::with_capacity(budget);
        loop {
            buf.clear();
            while buf.len() < budget {
                match rd.next()? { Some(r) => buf.push(r), None => break }
            }
            if buf.is_empty() { break; }
            buf.sort_by_key(|r| keyof(r, by));
            let p = tmp.join(format!("run{}.bin", runs.len()));
            let mut w = RecWriter::create(&p)?;
            for r in &buf { w.push(*r)?; }
            w.finish()?;
            runs.push(p);
        }
    }
    if runs.is_empty() { RecWriter::create(dst)?.finish()?; return Ok(0); }

    // Merge in passes of at most MAX_FANIN runs. A whole-genome sort produces
    // thousands of runs -- 5.9e9 records at an 18 MB budget is ~5,900 -- and a
    // single k-way merge would need one open file and one buffer per run, which
    // hits the file-descriptor limit long before it hits the memory budget.
    const MAX_FANIN: usize = 64;
    let mut pass = 0usize;
    while runs.len() > 1 {
        let mut next: Vec<PathBuf> = Vec::new();
        for (gi, group) in runs.chunks(MAX_FANIN).enumerate() {
            if group.len() == 1 { next.push(group[0].clone()); continue; }
            let out = tmp.join(format!("merge{pass}_{gi}.bin"));
            merge_runs(group, &out, by, budget)?;
            for p in group { let _ = std::fs::remove_file(p); }
            next.push(out);
        }
        runs = next;
        pass += 1;
    }
    std::fs::rename(&runs[0], dst)?;
    Ok(std::fs::metadata(dst)?.len() / REC as u64)
}

/// One k-way merge over `runs`, stable: ties break by run index, which is input
/// order because the runs were produced in order.
fn merge_runs(runs: &[PathBuf], dst: &Path, by: By, budget: usize) -> std::io::Result<u64> {
    // Split one buffer budget across the runs instead of giving each its own.
    let cap = ((budget * REC) / runs.len().max(1)).clamp(8 * 1024, 1 << 20);
    let mut rds: Vec<RecReader> = runs.iter()
        .map(|p| RecReader::open_buf(p, cap)).collect::<Result<_,_>>()?;
    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    for (i, rd) in rds.iter_mut().enumerate() {
        if let Some(r) = rd.next()? { heap.push(HeapItem { key: keyof(&r, by), idx: i, rec: r }); }
    }
    let mut w = RecWriter::create(dst)?;
    while let Some(it) = heap.pop() {
        w.push(it.rec)?;
        if let Some(r) = rds[it.idx].next()? {
            heap.push(HeapItem { key: keyof(&r, by), idx: it.idx, rec: r });
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
