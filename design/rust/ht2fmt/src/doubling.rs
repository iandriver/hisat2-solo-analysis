//! The prefix-doubling loop in external memory, shared by `ht2ext` (which
//! checks the generation curve) and `ht2wg` (which goes on to emit an index).
//!
//! `createNewNodesMaker` is an equi-join, `past_nodes.to = from_table.from`,
//! implemented in C++ as a random-access probe into a resident `from_table`.
//! That probe is the only reason all three arrays have to be in RAM, and it is
//! the reason the whole-genome build needs 625 GiB. An equi-join can always be
//! run as sort-merge over two sequential streams, so each generation becomes:
//!
//!   1. sort the current nodes by `to`
//!   2. merge against the same nodes sorted by `from`, emitting the joined nodes
//!   3. sort the output by key
//!   4. stream `mergeUpdateRank` over it
//!
//! Nothing is resident except one sort buffer and one block of equal-`key.first`
//! nodes. The memory budget is an explicit argument, so setting it to a few
//! hundred records forces real runs and a real k-way merge on a 200 bp graph.

use super::ext;
use super::ext::{By, Rec, RecReader, RecWriter, SORTED};
use super::graph;
use std::fs;
use std::path::{Path, PathBuf};

/// Sort-merge join: for each node A, emit one record per node B with
/// `B.from == A.to`. Already-sorted nodes pass through untouched (generations
/// past the first pruning only re-join the unsorted ones).
fn join(a_by_to: &Path, b_by_from: &Path, out: &Path, gen: u32) -> std::io::Result<u64> {
    let mut a = RecReader::open(a_by_to)?;
    let mut b = RecReader::open(b_by_from)?;
    let mut w = RecWriter::create(out)?;
    let mut bcur = b.next()?;
    let mut group: Vec<Rec> = Vec::new();
    let mut group_key: u32 = u32::MAX;
    let mut acur = a.next()?;

    while let Some(an) = acur {
        if gen > 4 && an.is_sorted() { w.push(an)?; acur = a.next()?; continue; }
        if an.to != group_key {
            // advance B to the group for this `to`; both streams are ordered, so
            // this only ever moves forward
            while let Some(bn) = bcur {
                if bn.from < an.to { bcur = b.next()?; } else { break; }
            }
            group.clear();
            while let Some(bn) = bcur {
                if bn.from == an.to { group.push(bn); bcur = b.next()?; } else { break; }
            }
            group_key = an.to;
        }
        for bn in &group {
            let (k0, k1) = if gen <= 3 {
                // keys packed into one integer, 3 bits per character
                let shift = 3u32 * (1u32 << (gen - 1));
                ((an.k0 << shift) + bn.k0, 0)
            } else {
                (an.k0, bn.k0)
            };
            w.push(Rec { from: an.from, to: bn.to, k0, k1 })?;
        }
        acur = a.next()?;
    }
    w.finish()
}

/// `mergeUpdateRank`'s block walk, streamed.
///
/// Only one block of equal `key.first` is ever resident. The rule is the one the
/// in-memory version verified generation-for-generation: a run of equal-key nodes
/// sharing a single `from` is one path node seen several ways and collapses to
/// its first member; a run spanning several `from` values is genuinely several
/// nodes at the same rank and all of them survive.
/// A reader with two records of lookahead, which `mergeUpdateRank` needs.
struct Peek { r: RecReader, buf: Vec<Rec> }

impl Peek {
    fn new(p: &Path) -> std::io::Result<Peek> { Ok(Peek { r: RecReader::open(p)?, buf: Vec::new() }) }
    fn fill(&mut self, n: usize) -> std::io::Result<()> {
        while self.buf.len() < n {
            match self.r.next()? { Some(x) => self.buf.push(x), None => break }
        }
        Ok(())
    }
    fn at(&mut self, i: usize) -> std::io::Result<Option<Rec>> {
        self.fill(i + 1)?;
        Ok(self.buf.get(i).copied())
    }
    fn take(&mut self) -> std::io::Result<Option<Rec>> {
        self.fill(1)?;
        Ok(if self.buf.is_empty() { None } else { Some(self.buf.remove(0)) })
    }
}

/// `mergeUpdateRank`'s block walk, streamed.
///
/// Only one block of equal `key.first` is ever resident. The collapse rule is
/// the one the in-memory version verified generation-for-generation: a run of
/// equal-key nodes sharing a single `from` is one path node seen several ways
/// and collapses to its first member; a run spanning several `from` values is
/// genuinely several nodes at the same rank and all survive.
///
/// The part that is easy to miss — and that a block-at-a-time rewrite silently
/// drops — is the lookahead after a multi-node block: the record immediately
/// following can be DISCARDED outright when it is a block of one, the previous
/// written node is already sorted, and the two share a `from`. Without it the
/// node counts drift the moment pruning starts.
fn merge_update_rank(src: &Path, dst: &Path) -> std::io::Result<(u64, u64)> {
    let mut p = Peek::new(src)?;
    let mut w = RecWriter::create(dst)?;
    let mut ranks: u64 = 0;
    let mut out_n: u64 = 0;
    let mut block: Vec<Rec> = Vec::new();
    let mut prev_written: Option<Rec> = None;

    loop {
        let head = match p.take()? { Some(x) => x, None => break };
        block.clear();
        let k0 = head.k0;
        block.push(head);
        while let Some(x) = p.at(0)? {
            if x.k0 != k0 { break; }
            block.push(p.take()?.unwrap());
        }

        if block.len() == 1 {
            let mut n = block[0];
            n.k0 = ranks; ranks += 1;
            w.push(n)?; out_n += 1;
            prev_written = Some(n);
            continue;
        }

        block.sort_by_key(|x| x.k1);
        let mut i = 0usize;
        while i < block.len() {
            let mut shift = 1usize;
            while i + shift < block.len() && block[i].key() == block[i + shift].key() { shift += 1; }
            let merge = (i..i + shift).all(|j| block[j].from == block[i].from);
            if !merge {
                for j in i..i + shift {
                    let mut n = block[j]; n.k0 = ranks;
                    w.push(n)?; out_n += 1;
                    prev_written = Some(n);
                }
                ranks += 1;
            } else {
                let fold = match prev_written {
                    Some(pw) => pw.is_sorted() && pw.from == block[i].from,
                    None => false,
                };
                if !fold {
                    let mut n = block[i];
                    n.to = SORTED;
                    n.k0 = ranks; ranks += 1;
                    w.push(n)?; out_n += 1;
                    prev_written = Some(n);
                }
            }
            i += shift;
        }

        // The lookahead. `node` in the C++ is the first record of the next
        // block; it is skipped when it is not itself part of a cluster, the
        // previous written node is sorted, and they share a `from`.
        if let (Some(c), Some(pw)) = (p.at(0)?, prev_written) {
            let solo = match p.at(1)? { Some(n) => n.k0 != c.k0, None => true };
            if solo && pw.is_sorted() && c.from == pw.from {
                p.take()?;
            }
        }
    }
    w.finish()?;
    Ok((out_n, ranks))
}

/// `mergeUpdateRank`'s generation-4 body (`gbwt_graph.h:2160`), streamed.
///
/// Generation 4 does NOT use the block walk. It runs `nextMaximalSet`, which
/// collapses each maximal run that shares a single `from` into its first member,
/// then marks every key that occurs exactly once as sorted and renumbers. Using
/// the block walk here instead keeps too many nodes -- 248 against 241 on the
/// 200 bp graph -- because the two rules disagree about which runs may collapse.
///
/// `nextMaximalSet` scans forward from the run's head until `from` changes,
/// remembering the last position at which the key changed, and consumes only up
/// to that boundary. So the buffer is one run of equal `from`, not the file.
fn merge_update_rank_gen4(src: &Path, dst: &Path) -> std::io::Result<(u64, u64)> {
    let collapsed = dst.with_extension("collapse");
    {
        let mut r = RecReader::open(src)?;
        let mut w = RecWriter::create(&collapsed)?;
        let mut buf: Vec<Rec> = Vec::new();
        let mut prev_key: Option<(u64, u64)> = None;
        let mut eof = false;
        loop {
            if buf.is_empty() && !eof {
                match r.next()? { Some(x) => buf.push(x), None => eof = true }
            }
            if buf.is_empty() { break; }
            let r0 = buf[0];
            if prev_key == Some(r0.key()) {
                w.push(r0)?;
                buf.remove(0);
                prev_key = Some(r0.key());
                continue;
            }
            let mut r1 = 1usize;
            let mut i = 1usize;
            loop {
                while buf.len() <= i && !eof {
                    match r.next()? { Some(x) => buf.push(x), None => eof = true }
                }
                if buf.len() <= i { r1 = buf.len(); break; }
                if buf[i - 1].key() != buf[i].key() { r1 = i; }
                if buf[i].from != r0.from { break; }
                i += 1;
            }
            w.push(r0)?;
            prev_key = Some(buf[r1 - 1].key());
            buf.drain(0..r1);
        }
        w.finish()?;
    }

    // A key that occurs exactly once becomes sorted; equal keys share a rank.
    let mut r = RecReader::open(&collapsed)?;
    let mut w = RecWriter::create(dst)?;
    let mut group: Vec<Rec> = Vec::new();
    let mut ranks: u64 = 0;
    let mut out_n: u64 = 0;
    let mut cur = r.next()?;
    while let Some(head) = cur {
        group.clear();
        group.push(head);
        cur = r.next()?;
        while let Some(x) = cur {
            if x.key() != head.key() { break; }
            group.push(x);
            cur = r.next()?;
        }
        let solo = group.len() == 1;
        for g in group.iter() {
            let mut n = *g;
            if solo { n.to = SORTED; }
            n.k0 = ranks; n.k1 = 0;
            w.push(n)?; out_n += 1;
        }
        ranks += 1;
    }
    w.finish()?;
    let _ = fs::remove_file(&collapsed);
    Ok((out_n, ranks))
}

/// The reference graph, on disk, plus the converged path-node file.
pub struct Doubled {
    /// Path nodes after convergence: `from` is a reference-graph node id and
    /// `k0` is the node's final rank. In key order.
    pub cur: PathBuf,
    pub n_path_nodes: u64,
    /// `(generation, temp_nodes, nodes, ranks)`, the four numbers a build log
    /// prints per line.
    pub curve: Vec<(u32, u64, u64, u64)>,
    pub graph: graph::GraphOnDisk,
}

/// Generation 0 plus the loop, entirely on disk.
///
/// The reference graph is always built fragment by fragment straight to record
/// files, so nothing about it stays resident and generation 0 is produced by
/// streaming those files rather than indexing `Vec`s. Edges are in
/// `sortEdgesFrom` order and nodes are in id order, so the label lookup each
/// edge needs is a single forward pass, not a random probe.
pub fn run(fa: &str, snp: &str, hap: &str, wd: &Path, budget: usize, chunk: u32,
           verbose: bool) -> std::io::Result<Doubled>
{
    use std::io::Read;
    fs::create_dir_all(wd)?;
    let cur = wd.join("cur.bin");
    let code = |l: u8| -> u64 { match l { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
                                          b'Y' => 4, _ => 5 } };

    let g = graph::build_fragmented_to_disk(fa, snp, hap, chunk, wd)?;
    if verbose {
        println!("graph on disk: {} nodes, {} edges, chunk {} kb", g.n_nodes, g.n_edges, chunk / 1024);
    }
    {
        let mut nf = std::io::BufReader::with_capacity(1 << 20, fs::File::open(wd.join("nodes.bin"))?);
        let mut ef = std::io::BufReader::with_capacity(1 << 20, fs::File::open(wd.join("edges.bin"))?);
        let mut w = RecWriter::create(&cur)?;
        let mut nbuf = [0u8; 5];
        let mut ebuf = [0u8; 8];
        let mut node_i: u64 = 0;
        let mut cur_label: u8 = 0;
        loop {
            match ef.read_exact(&mut ebuf) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            }
            let f = u32::from_le_bytes(ebuf[0..4].try_into().unwrap());
            let t = u32::from_le_bytes(ebuf[4..8].try_into().unwrap());
            while node_i <= f as u64 {
                nf.read_exact(&mut nbuf)?;
                cur_label = nbuf[0];
                node_i += 1;
            }
            w.push(Rec { from: f, to: t, k0: code(cur_label), k1: 0 })?;
        }
        w.push(Rec { from: g.last_node, to: g.last_node, k0: 5, k1: 0 })?;
        w.finish()?;
    }

    let n0 = fs::metadata(&cur)?.len() / ext::REC as u64;
    let mut curve: Vec<(u32, u64, u64, u64)> = vec![(0, n0, n0, 0)];
    if verbose {
        println!("Generation 0 ({n0} -> {n0} nodes, 0 ranks)   [budget {budget} records = {} KB]",
                 budget * ext::REC / 1024);
    }

    let (by_to, by_from, joined, sorted_k) =
        (wd.join("by_to.bin"), wd.join("by_from.bin"), wd.join("joined.bin"), wd.join("sorted.bin"));
    let mut gen = 0u32;
    let mut n_path_nodes = n0;
    loop {
        gen += 1;
        ext::sort_external(&cur, &by_to, By::To, budget, wd)?;
        ext::sort_external(&cur, &by_from, By::From, budget, wd)?;
        let temp = join(&by_to, &by_from, &joined, gen)?;

        let (nodes, ranks) = if gen <= 3 {
            fs::rename(&joined, &cur)?;
            (temp, 0u64)
        } else {
            ext::sort_external(&joined, &sorted_k, By::Key, budget, wd)?;
            let (n, rk) = if gen == 4 { merge_update_rank_gen4(&sorted_k, &cur)? }
                          else        { merge_update_rank(&sorted_k, &cur)? };
            (n, rk)
        };
        if verbose { println!("Generation {gen} ({temp} -> {nodes} nodes, {ranks} ranks)"); }
        curve.push((gen, temp, nodes, ranks));
        n_path_nodes = nodes;
        if gen > 3 && ranks == nodes { break; }
        if gen > 64 { panic!("doubling did not converge"); }
    }
    for p in [&by_to, &by_from, &joined, &sorted_k] { let _ = fs::remove_file(p); }
    Ok(Doubled { cur, n_path_nodes, curve, graph: g })
}

/// Compare a curve against the `Generation` lines of a `hisat2-build` log.
pub fn check_curve(curve: &[(u32, u64, u64, u64)], log_path: &str) -> bool {
    let log = fs::read_to_string(log_path).expect("log");
    let theirs: Vec<(u32, u64, u64, u64)> = log.lines()
        .filter(|l| l.starts_with("Generation "))
        .map(|l| {
            let g: u32 = l["Generation ".len()..].split(' ').next().unwrap().parse().unwrap();
            let inner = &l[l.find('(').unwrap() + 1..];
            let t: Vec<&str> = inner.split_whitespace().collect();
            (g, t[0].parse().unwrap(), t[2].parse().unwrap(), t[4].parse().unwrap())
        }).collect();
    let mut bad = 0;
    for i in 0..curve.len().max(theirs.len()) {
        let (o, t) = (curve.get(i), theirs.get(i));
        if o != t { bad += 1; println!("  gen {i}: ours {:?} theirs {:?}", o, t); }
    }
    if bad == 0 { println!("\nEXTERNAL CURVE MATCHES HISAT2 ({} generations)", curve.len()); true }
    else { println!("\n{bad} generations differ"); false }
}
