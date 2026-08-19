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
        Ok(RecReader { r: BufReader::with_capacity(1 << 20, File::open(p)?), buf: [0; REC] })
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
pub enum By { To, From, Key }

fn keyof(r: &Rec, by: By) -> (u64, u64, u32) {
    match by {
        By::To   => (r.to as u64, 0, 0),
        By::From => (r.from as u64, 0, 0),
        // `from` is part of the key, not decoration: mergeUpdateRank collapses a
        // run of equal-key nodes only when they share a `from`, and keeps the
        // FIRST of the run, so ties have to break deterministically or the
        // surviving node's `to` changes.
        By::Key  => (r.k0, r.k1, r.from),
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

    let mut rds: Vec<RecReader> = runs.iter().map(|p| RecReader::open(p)).collect::<Result<_,_>>()?;
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
    let n = w.finish()?;
    for p in runs { let _ = std::fs::remove_file(p); }
    Ok(n)
}
