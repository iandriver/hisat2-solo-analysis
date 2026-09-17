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
use super::ext::{By, Rec, SegReader, SegWriter, SORTED};
use super::ckpt;
use super::graph;
use std::fs;
use std::path::{Path, PathBuf};

/// Sort-merge join: for each node A, emit one record per node B with
/// `B.from == A.to`. Already-sorted nodes pass through untouched (generations
/// past the first pruning only re-join the unsorted ones).
fn join(a_by_to: &Path, b_by_from: &Path, out: &Path, gen: u32) -> std::io::Result<u64> {
    // Both inputs are read once, strictly forward, and dead afterwards.
    let mut a = SegReader::open(a_by_to, true)?;
    let mut b = SegReader::open(b_by_from, true)?;
    let mut w = SegWriter::create(out)?;
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
/// The join, split by the value it joins on.
///
/// `a` is ordered by `to` and `b` by `from`, and a record's group is every `b`
/// with `b.from == a.to` -- so a range of that one value cuts both sides at once
/// and no group is ever split. The output is in `a` order, which is `to` order,
/// so concatenating the partitions in cut order is exactly what the serial join
/// would have written.
///
/// Sorted nodes carry `to == SORTED == u32::MAX` and so all land in the last
/// partition, where the passthrough runs as before.
fn join_par(a_by_to: &Path, b_by_from: &Path, out: &Path, gen: u32, budget: usize,
            nthreads: usize) -> std::io::Result<u64>
{
    let ai = ext::SegIndex::open(a_by_to)?;
    let bi = ext::SegIndex::open(b_by_from)?;
    if ai.len == 0 { ext::SegWriter::create(out)?.finish()?; return Ok(0); }

    // cuts taken from `a`, evenly by record position -- the join is linear in
    // `a`, so equal shares of `a` are equal shares of the work
    let want = nthreads.min(ai.len as usize);
    let mut cuts: Vec<u32> = Vec::new();
    for p in 1..want {
        let i = (p as u64 * ai.len) / want as u64;
        cuts.push(ai.at(i)?.to);
    }
    cuts.dedup();

    let nparts = cuts.len() + 1;
    let mut abnd = vec![0u64; nparts + 1];
    let mut bbnd = vec![0u64; nparts + 1];
    abnd[nparts] = ai.len; bbnd[nparts] = bi.len;
    for (ci, &c) in cuts.iter().enumerate() {
        abnd[ci + 1] = ai.lower_bound((c as u64, 0, 0), By::To)?;
        bbnd[ci + 1] = bi.lower_bound((c as u64, 0, 0), By::From)?;
    }

    let cap = ((budget * ext::REC) / (2 * nthreads)).clamp(8 * 1024, 1 << 20);
    let parts: Vec<PathBuf> = (0..nparts).map(|p| {
        let mut q = out.as_os_str().to_owned(); q.push(format!(".j{p}")); PathBuf::from(q)
    }).collect();
    let counts: Vec<std::sync::atomic::AtomicU64> =
        (0..nparts).map(|_| std::sync::atomic::AtomicU64::new(0)).collect();
    let err: std::sync::Mutex<Option<std::io::Error>> = std::sync::Mutex::new(None);
    let next = std::sync::atomic::AtomicUsize::new(0);
    {
        let (ai, bi, abnd, bbnd, parts, counts, err, next) =
            (&ai, &bi, &abnd, &bbnd, &parts, &counts, &err, &next);
        std::thread::scope(|scope| {
            for _ in 0..nthreads.min(nparts) {
                scope.spawn(move || loop {
                    let p = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if p >= nparts { break; }
                    let r = (|| -> std::io::Result<u64> {
                        let mut a = ext::SegSlice::new(ai, abnd[p], abnd[p + 1], cap, false)?;
                        let mut b = ext::SegSlice::new(bi, bbnd[p], bbnd[p + 1], cap, false)?;
                        let mut w = ext::SegWriter::create(&parts[p])?;
                        let mut bcur = b.next(cap)?;
                        let mut group: Vec<Rec> = Vec::new();
                        let mut group_key: u32 = u32::MAX;
                        let mut acur = a.next(cap)?;
                        let mut first = true;
                        while let Some(an) = acur {
                            if gen > 4 && an.is_sorted() { w.push(an)?; acur = a.next(cap)?; continue; }
                            if an.to != group_key || first {
                                while let Some(bn) = bcur {
                                    if bn.from < an.to { bcur = b.next(cap)?; } else { break; }
                                }
                                group.clear();
                                while let Some(bn) = bcur {
                                    if bn.from == an.to { group.push(bn); bcur = b.next(cap)?; } else { break; }
                                }
                                group_key = an.to;
                                first = false;
                            }
                            for bn in &group {
                                let (k0, k1) = if gen <= 3 {
                                    let shift = 3u32 * (1u32 << (gen - 1));
                                    ((an.k0 << shift) + bn.k0, 0)
                                } else {
                                    (an.k0, bn.k0)
                                };
                                w.push(Rec { from: an.from, to: bn.to, k0, k1 })?;
                            }
                            acur = a.next(cap)?;
                        }
                        w.finish()
                    })();
                    match r {
                        Ok(n) => { counts[p].store(n, std::sync::atomic::Ordering::Relaxed); }
                        Err(e) => { *err.lock().unwrap() = Some(e); }
                    }
                });
            }
        });
    }
    if let Some(e) = err.into_inner().unwrap() { return Err(e); }
    ext::seg_remove(out);
    let mut n = 0u64;
    for p in 0..nparts {
        let c = counts[p].load(std::sync::atomic::Ordering::Relaxed);
        if c == 0 { ext::seg_remove(&parts[p]); continue; }
        ext::seg_append(out, &parts[p])?;
        n += c;
    }
    // The inputs are NOT removed here. The caller removes them once the phase
    // that produced this output is checkpointed -- a restart before that point
    // re-runs the join, and it needs them.
    Ok(n)
}

/// A reader with two records of lookahead, which `mergeUpdateRank` needs.
struct Peek { r: SegReader, buf: Vec<Rec> }

impl Peek {
    fn new(p: &Path) -> std::io::Result<Peek> { Ok(Peek { r: SegReader::open(p, true)?, buf: Vec::new() }) }
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
fn merge_update_rank(src: &Path, dst: &Path) -> std::io::Result<(u64, u64, u64)> {
    let mut p = Peek::new(src)?;
    let mut w = SegWriter::create(dst)?;
    let mut ranks: u64 = 0;
    let mut out_n: u64 = 0;
    // How many nodes come out already sorted. The next generation needs it to
    // know whether the shortcut in `join_late` is worth taking, and counting
    // here costs nothing.
    let mut n_sorted: u64 = 0;
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
            if n.is_sorted() { n_sorted += 1; }
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
                    if n.is_sorted() { n_sorted += 1; }
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
                    n_sorted += 1;
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
    Ok((out_n, ranks, n_sorted))
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
fn merge_update_rank_gen4(src: &Path, dst: &Path) -> std::io::Result<(u64, u64, u64)> {
    let collapsed = dst.with_extension("collapse");
    {
        let mut r = SegReader::open(src, true)?;
        let mut w = SegWriter::create(&collapsed)?;
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
    let mut r = SegReader::open(&collapsed, true)?;
    let mut w = SegWriter::create(dst)?;
    let mut group: Vec<Rec> = Vec::new();
    let mut ranks: u64 = 0;
    let mut out_n: u64 = 0;
    let mut n_sorted: u64 = 0;
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
            if n.is_sorted() { n_sorted += 1; }
            w.push(n)?; out_n += 1;
        }
        ranks += 1;
    }
    w.finish()?;
    ext::seg_remove(&collapsed);
    Ok((out_n, ranks, n_sorted))
}


/// The join, when almost every node is already sorted.
///
/// From generation 5 the unsorted set collapses -- on a 20 Mb reference it goes
/// 675,653 -> 116,474 -> 28,919 -> 4,631 -> 544 out of 21.1M nodes -- yet the
/// plain path sorts all 21 million by `to` and again by `from` every generation
/// to move 99.9% of them through untouched. That is 19.6s of the 31.6s the
/// doubling spends sorting.
///
/// The shortcut is exact, not an approximation, and rests on one fact: a sorted
/// node carries `to == SORTED == u32::MAX`, the largest `to` there is. A stable
/// sort by `to` therefore already puts every sorted node at the end, in its
/// original order, so
///
///     by_to  ==  [unsorted, sorted by `to`]  ++  [sorted, in `cur` order]
///
/// and `join` walking that emits its computed head followed by the sorted nodes
/// verbatim. Neither half needs the other to be materialised.
///
/// The from-table shrinks the same way. Only nodes some unsorted node actually
/// looks up can be reached, so one streaming filter against the `to` values
/// replaces a full sort by `from`.
///
/// Three sequential reads of `cur` and one write, against two full external
/// sorts and a merge.
fn join_late(cur: &Path, joined: &Path, wd: &Path, gen: u32, budget: usize, resume: bool)
    -> std::io::Result<u64>
{
    let (un, mini) = (wd.join("un.bin"), wd.join("mini.bin"));
    let (u_by_to, m_by_from) = (wd.join("ubt.bin"), wd.join("mbf.bin"));

    // pass 1 -- the unsorted nodes, and every `to` they will look up
    let mut tos: Vec<u32> = Vec::new();
    {
        let mut r = SegReader::open(cur, false)?;
        let mut w = SegWriter::create(&un)?;
        while let Some(n) = r.next()? {
            if !n.is_sorted() { w.push(n)?; tos.push(n.to); }
        }
        w.finish()?;
    }
    tos.sort_unstable();
    tos.dedup();

    // pass 2 -- the only nodes the join can reach, and the sorted tail, in one
    // read. The tail is written straight out in `cur` order, which is exactly
    // where a stable sort by `to` would have left it.
    let tail = wd.join("tail.bin");
    let mut n_tail = 0u64;
    {
        // Not consumed when resume is on: the caller deletes `cur` only once
        // `joined` is whole, which is what a restart falls back to.
        let mut r = SegReader::open(cur, !resume)?;
        let mut w = SegWriter::create(&mini)?;
        let mut t = SegWriter::create(&tail)?;
        while let Some(n) = r.next()? {
            // BOTH, not either: a sorted node still sits in the from-table, so
            // it can be somebody's lookup target as well as part of the tail
            if tos.binary_search(&n.from).is_ok() { w.push(n)?; }
            if n.is_sorted() { t.push(n)?; n_tail += 1; }
        }
        w.finish()?; t.finish()?;
    }
    drop(tos);

    ext::sort_external(&un, &u_by_to, By::To, budget, wd, true)?;
    ext::sort_external(&mini, &m_by_from, By::From, budget, wd, true)?;
    let n = join(&u_by_to, &m_by_from, joined, gen)?;
    ext::seg_remove(&u_by_to);
    ext::seg_remove(&m_by_from);

    // and the tail goes on the end, by renaming its segments -- no third pass
    // and no copy
    ext::seg_append(joined, &tail)?;
    Ok(n + n_tail)
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
           verbose: bool, ck: &mut ckpt::Ckpt, resume: bool) -> std::io::Result<Doubled>
{
    run_with(fa, snp, hap, "", "", wd, budget, chunk, verbose, ck, resume)
}

pub fn run_with(fa: &str, snp: &str, hap: &str, ss: &str, exon: &str,
                wd: &Path, budget: usize, chunk: u32,
                verbose: bool, ck: &mut ckpt::Ckpt, resume: bool) -> std::io::Result<Doubled>
{
    use std::io::Read;
    fs::create_dir_all(wd)?;
    // The node file alternates between two names. A generation reads one and
    // writes the other, so the state a checkpoint names is never the state the
    // next step is overwriting -- which is the whole of what makes a restart
    // safe, and cannot be had from a single `cur.bin` that each generation
    // consumes in place.
    let cur_of = |g: u32| wd.join(format!("cur{}.bin", g % 2));
    let code = |l: u8| -> u64 { match l { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
                                          b'Y' => 4, _ => 5 } };

    let resumed = !ck.stage.is_empty();
    // Once the loop has converged there is nothing left to re-enter, and the
    // node file is gone -- `generateEdges` consumed it. Coming back in would run
    // generation 11 over an empty file, converge again on nothing, and hand a
    // zero-node graph to everything downstream. It looks like progress.
    if ck.stage == "doubled" || ck.stage == "edges" {
        if verbose {
            println!("doubling: from checkpoint -- {} generations, {} path nodes",
                     ck.curve.len(), ck.path_nodes);
        }
        return Ok(Doubled {
            cur: cur_of(ck.gen), n_path_nodes: ck.path_nodes, curve: ck.curve.clone(),
            graph: graph::GraphOnDisk { n_nodes: ck.g_nodes, n_edges: ck.g_edges,
                                        last_node: ck.g_last, text_len: ck.g_text },
        });
    }
    let g = if resumed {
        graph::GraphOnDisk { n_nodes: ck.g_nodes, n_edges: ck.g_edges,
                             last_node: ck.g_last, text_len: ck.g_text }
    } else {
        let g = graph::build_fragmented_to_disk_with(fa, snp, hap, ss, exon, chunk, wd)?;
        ck.stage = "graph".into();
        ck.g_nodes = g.n_nodes; ck.g_edges = g.n_edges;
        ck.g_last = g.last_node; ck.g_text = g.text_len;
        if resume { ckpt::save(wd, ck)?; }
        ckpt::crash_point("graph", 0);
        g
    };
    if verbose {
        println!("graph on disk: {} nodes, {} edges, chunk {} kb{}", g.n_nodes, g.n_edges,
                 chunk / 1024, if resumed { " (from checkpoint)" } else { "" });
    }
    let cur = cur_of(ck.gen);
    if ck.gen == 0 {
        let mut nf = std::io::BufReader::with_capacity(1 << 20, fs::File::open(wd.join("nodes.bin"))?);
        let mut ef = std::io::BufReader::with_capacity(1 << 20, fs::File::open(wd.join("edges.bin"))?);
        let mut w = SegWriter::create(&cur)?;
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

    let n0 = ext::seg_len(&cur_of(0));
    let mut curve: Vec<(u32, u64, u64, u64)> = if ck.gen > 0 { ck.curve.clone() }
                                               else { vec![(0, n0, n0, 0)] };
    if verbose && ck.gen == 0 {
        println!("Generation 0 ({n0} -> {n0} nodes, 0 ranks)   [budget {budget} records = {} KB]",
                 budget * ext::REC / 1024);
    }

    let (by_to, by_from, joined, sorted_k) =
        (wd.join("by_to.bin"), wd.join("by_from.bin"), wd.join("joined.bin"), wd.join("sorted.bin"));
    let mut gen = ck.gen;
    let mut n_path_nodes = if ck.gen > 0 { ck.path_nodes } else { n0 };
    // Where the doubling's time actually goes, per phase, summed over all
    // generations -- the question of what to thread first is not answerable
    // from the source.
    let (mut t_sort, mut t_join, mut t_rank) = (0f64, 0f64, 0f64);
    let (mut p_sort, mut p_join, mut p_rank) = (0f64, 0f64, 0f64);
    // nodes already marked sorted coming into the next generation
    let mut n_sorted: u64 = ck.sorted;
    loop {
        gen += 1;
        let cur_in = cur_of(gen - 1);
        let cur_out = cur_of(gen);
        let mut t = std::time::Instant::now();
        // `cur` is dead the moment both orderings of it exist -- the join reads
        // `by_to` and `by_from`, and whichever branch follows rewrites `cur`
        // from scratch. Scratch is the binding resource here, so every file gets
        // deleted at the point it stops being readable rather than at the end of
        // the generation.
        // `join` honours the sorted passthrough only past generation 4, and the
        // shortcut only pays while the unsorted set is small enough that its
        // `to` values fit in the sort budget's worth of memory.
        let unsorted = n_path_nodes.saturating_sub(n_sorted);
        let late = gen > 4 && n_sorted > 0
                   && unsorted * 4 < n_path_nodes
                   && unsorted <= budget as u64;
        // ---- phase 1: everything up to `joined` -------------------------
        // `cur_in` is NOT consumed. It is the only state a restart can fall back
        // to until `joined` is whole, and consuming it to save a copy is exactly
        // what makes a crash here unrecoverable.
        // Generations 1-3 have no second phase -- no key sort, no ranking -- so
        // the join writes straight into the generation's output file and there
        // is nothing to hand over afterwards. Renaming `joined` onto `cur_out`
        // was the one step in a generation that could not be repeated: a restart
        // landing after it found `joined` already moved, cleared the
        // destination, and moved nothing onto it.
        let out = if gen <= 3 { cur_out.clone() } else { joined.clone() };
        let mut temp = ck.temp;
        if ck.phase.is_empty() || ck.phase == "sorted" {
            if late {
                temp = join_late(&cur_in, &out, wd, gen, budget, resume)?;
                t_join += t.elapsed().as_secs_f64(); t = std::time::Instant::now();
            } else {
                if ck.phase.is_empty() {
                    ext::sort_external(&cur_in, &by_to, By::To, budget, wd, false)?;
                    // Consumed only when resume is off. Keeping `cur_in` alive
                    // across the two sorts is what costs the extra copy of the
                    // node file -- and it is also the only thing a crash here
                    // could fall back to.
                    ext::sort_external(&cur_in, &by_from, By::From, budget, wd, !resume)?;
                    if resume { ck.phase = "sorted".into(); ckpt::save(wd, ck)?; }
                    ckpt::crash_point("sorted", gen);
                }
                t_sort += t.elapsed().as_secs_f64(); t = std::time::Instant::now();
                let nt = ext::threads();
                temp = if nt > 1 {
                    join_par(&by_to, &by_from, &out, gen, budget, nt)?
                } else {
                    let n = join(&by_to, &by_from, &out, gen)?;
                    ext::seg_remove(&by_to);
                    ext::seg_remove(&by_from);
                    n
                };
                t_join += t.elapsed().as_secs_f64(); t = std::time::Instant::now();
            }
            // Every generation records this, 1-3 included. Their output is
            // `cur_out` rather than `joined`, but the point is the same: once it
            // is whole, the join's inputs can go, and a restart that has not
            // seen this mark will try to redo the join from files that are
            // already gone.
            if resume { ck.phase = "joined".into(); ck.temp = temp; ckpt::save(wd, ck)?; }
            ckpt::crash_point("joined", gen);
            ext::seg_remove(&by_to);
            ext::seg_remove(&by_from);
            ext::seg_remove(&cur_in);
        }

        // ---- phase 2: the key sort, then the ranking --------------------
        let (nodes, ranks) = if gen <= 3 {
            (temp, 0u64)
        } else {
            if ck.phase == "joined" || !resume {
                ext::sort_external(&joined, &sorted_k, By::Key, budget, wd, !resume)?;
                if resume { ck.phase = "keyed".into(); ckpt::save(wd, ck)?; }
                ckpt::crash_point("keyed", gen);
                ext::seg_remove(&joined);
            }
            t_sort += t.elapsed().as_secs_f64(); t = std::time::Instant::now();
            let (n, rk, ns) = if gen == 4 { merge_update_rank_gen4(&sorted_k, &cur_out)? }
                              else        { merge_update_rank(&sorted_k, &cur_out)? };
            n_sorted = ns;
            (n, rk)
        };
        t_rank += t.elapsed().as_secs_f64();
        if verbose {
            let unsorted = nodes.saturating_sub(ranks);
            println!("Generation {gen} ({temp} -> {nodes} nodes, {ranks} ranks)   \
                      [{:.1}s sort, {:.1}s join, {:.1}s rank; {unsorted} unsorted{}]",
                     t_sort - p_sort, t_join - p_join, t_rank - p_rank,
                     if late { ", late" } else { "" });
            p_sort = t_sort; p_join = t_join; p_rank = t_rank;
        }
        if nodes == 0 {
            panic!("generation {gen} produced no nodes -- its input was empty, \
                    which means a restart lost state rather than resumed it");
        }
        curve.push((gen, temp, nodes, ranks));
        n_path_nodes = nodes;
        // The generation is finished only once `cur_out` is whole, so this is
        // the last thing that happens in it and `sorted.bin` outlives it until
        // then -- a restart between the two redoes `mergeUpdateRank` rather than
        // finding both its input and its output half-written.
        // Convergence is recorded by the SAME save that records the
        // generation. Two saves leave a window in which the state says
        // "generation 10 done, keep going" -- and going on runs an eleventh
        // generation over a converged file, which converges again and puts a
        // generation in the curve that the C++ never had.
        let converged = gen > 3 && ranks == nodes;
        ck.stage = if converged { "doubled".into() } else { "gen".to_string() };
        ck.phase = String::new(); ck.gen = gen;
        ck.curve = curve.clone(); ck.path_nodes = nodes; ck.sorted = n_sorted;
        if resume { ckpt::save(wd, ck)?; }
        ckpt::crash_point("gen", gen);
        ext::seg_remove(&sorted_k);
        ext::seg_remove(&joined);
        if converged { break; }
        if gen > 64 { panic!("doubling did not converge"); }
    }
    let cur = cur_of(gen);
    for p in [&by_to, &by_from, &joined, &sorted_k, &cur_of(gen + 1)] { ext::seg_remove(p); }
    if verbose {
        let tot = t_sort + t_join + t_rank;
        println!("  doubling phases: sort {:.1}s ({:.0}%), join {:.1}s ({:.0}%), mergeUpdateRank {:.1}s ({:.0}%)",
                 t_sort, 100.0 * t_sort / tot, t_join, 100.0 * t_join / tot,
                 t_rank, 100.0 * t_rank / tot);
    }
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
