//! Rung 3, step 2 — the prefix doubling, checked against the whole generation
//! curve rather than just its first line.
//!
//! `PathGraph` sorts the reference graph's paths by doubling their length:
//! generation `g` means "sorted by paths of length 2^g", and construction runs
//! `while(!isSorted())` — until every path node has a distinct rank. Each
//! generation prints
//!
//!     Generation g (temp_nodes -> nodes, ranks)
//!
//! so a build log is a per-generation oracle for three numbers at once. Nothing
//! about this can be matched by a plausible-looking implementation: the counts
//! depend on the exact pruning rule, and a wrong rule diverges within a
//! generation or two and never recovers.
//!
//! Structure reproduced from `gbwt_graph.h`. HISAT2 splits the loop into four
//! routines for memory reasons, and the split is visible in the arithmetic:
//!
//!   * `makeFromRef` (:1821)  — one path node per reference-graph EDGE, keyed by
//!     the label of its `from` node, plus one final node for 'Z' keyed 5.
//!   * `generationOne`/`earlyGeneration` (:1866, :1932), generations 1-3 — join
//!     and no pruning at all, so `nodes == temp_nodes` and `ranks == 0`. Keys
//!     are packed into one integer, shifting the left key up by
//!     `3 * 2^(g-1)` bits, which is why this can only run while 3 bits per
//!     character still fit.
//!   * `firstPruneGeneration` (:1946), generation 4 — the same join, then a full
//!     sort by key and the first pruning pass.
//!   * `lateGeneration` (:1979), generations 5+ — already-sorted nodes pass
//!     through untouched and only unsorted ones are re-joined, against a
//!     `from_table` built by sorting a copy by `from`. The output arrives
//!     grouped by `key.first` already, so there is no full sort here — only the
//!     within-block sort that `mergeUpdateRank` does itself.
//!
//! `mergeUpdateRank` (:2158) is where the pruning lives, and it has two
//! different bodies: a `nextMaximalSet` pass at generation 4 and the block walk
//! everywhere else. Both are reproduced.

#[path = "../graph.rs"]
mod graph;
// graph.rs streams the edge list through ext::sort_pairs_external, so every
// binary that includes it needs ext in scope too.
#[allow(dead_code)]
#[path = "../ext.rs"]
mod ext;

use std::collections::HashMap;
use std::convert::TryInto;
use std::{env, fs};

const UNSET: u32 = u32::MAX;

#[derive(Clone, Copy, PartialEq, Eq)]
struct PN { from: u32, to: u32, k0: u32, k1: u32 }

impl PN {
    fn is_sorted(&self) -> bool { self.to == UNSET }
    fn key(&self) -> (u32, u32) { (self.k0, self.k1) }
}

/// Offsets into a `from`-sorted node list, indexed by reference-graph node id.
fn csr(sorted_by_from: &[PN], n_ref_nodes: usize) -> Vec<u32> {
    let mut start = vec![0u32; n_ref_nodes + 1];
    for n in sorted_by_from { start[n.from as usize + 1] += 1; }
    for i in 1..start.len() { start[i] += start[i - 1]; }
    start
}

/// `mergeUpdateRank`, generation-4 body (`gbwt_graph.h:2160`).
fn next_maximal_set(nodes: &[PN], range: (usize, usize)) -> (usize, usize) {
    let n = nodes.len();
    if range.1 >= n { return (0, 0); }
    let r0 = range.1;
    let mut r1 = r0 + 1;
    if r0 > 0 && nodes[r0 - 1].key() == nodes[r0].key() { return (r0, r1); }
    for i in r1..n {
        if nodes[i - 1].key() != nodes[i].key() { r1 = i; }
        if nodes[i].from != nodes[r0].from { return (r0, r1); }
    }
    (r0, n)
}

fn merge_update_rank(nodes: &mut Vec<PN>, generation: u32) -> u32 {
    let ranks;
    if generation == 4 {
        // Collapse each maximal run that shares a `from` into its first member.
        let mut curr = 0usize;
        let mut range = (0usize, 0usize);
        loop {
            range = next_maximal_set(nodes, range);
            if range.0 >= range.1 { break; }
            nodes[curr] = nodes[range.0];
            curr += 1;
        }
        nodes.truncate(curr);

        // A node whose key is unique among its neighbours is now fully sorted.
        let mut candidate = Some(0usize);
        let mut key = nodes[0].key();
        for i in 1..nodes.len() {
            if nodes[i].key() != key {
                if let Some(c) = candidate { nodes[c].to = UNSET; }
                candidate = Some(i);
                key = nodes[i].key();
            } else {
                candidate = None;
            }
        }
        if let Some(c) = candidate { nodes[c].to = UNSET; }

        // Renumber: equal keys share a rank.
        let mut r = 0u32;
        let mut key = nodes[0].key();
        for i in 0..nodes.len() {
            let k = nodes[i].key();
            if k != key { key = k; r += 1; }
            nodes[i].k0 = r;
            nodes[i].k1 = 0;
        }
        ranks = r + 1;
    } else {
        let n = nodes.len();
        let mut block_start = 0usize;
        let mut curr = 0usize;
        let mut node = 0usize;
        let mut r = 0u32;
        loop {
            node += 1;
            if node == n || nodes[node].k0 != nodes[block_start].k0 {
                if node - block_start == 1 {
                    nodes[block_start].k0 = r; r += 1;
                    nodes[curr] = nodes[block_start]; curr += 1;
                } else {
                    nodes[block_start..node].sort_by_key(|x| x.k1);
                    while block_start != node {
                        let mut shift = 1usize;
                        while block_start + shift != node
                              && nodes[block_start].key() == nodes[block_start + shift].key() {
                            shift += 1;
                        }
                        // A run sharing one `from` is one path node seen several
                        // ways, so it can collapse; a run spanning several `from`
                        // values is genuinely several nodes at the same rank.
                        let merge = (block_start..block_start + shift)
                            .all(|i| nodes[i].from == nodes[block_start].from);
                        if !merge {
                            for i in block_start..block_start + shift {
                                nodes[i].k0 = r;
                                nodes[curr] = nodes[i]; curr += 1;
                            }
                            r += 1;
                        } else if curr == 0 || !nodes[curr - 1].is_sorted()
                                  || nodes[curr - 1].from != nodes[block_start].from {
                            nodes[block_start].to = UNSET;
                            nodes[block_start].k0 = r; r += 1;
                            nodes[curr] = nodes[block_start]; curr += 1;
                        }
                        block_start += shift;
                    }
                    if node == n { break; }
                    if node + 1 == n
                       && nodes[curr - 1].is_sorted() && nodes[node].from == nodes[curr - 1].from {
                        break;
                    }
                    // The node just past an unsorted cluster can fold into the
                    // previous one when it is not itself part of a cluster, the
                    // previous node is sorted, and they share a `from`.
                    if node + 1 < n && nodes[node].k0 != nodes[node + 1].k0
                       && nodes[curr - 1].is_sorted() && nodes[node].from == nodes[curr - 1].from {
                        node += 1;
                    }
                }
                block_start = node;
            }
            if node == n { break; }
        }
        nodes.truncate(curr);
        ranks = r;
    }
    ranks
}

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: ht2path <reference.fa> <snp> <haplotype> [<build.log>]");
        std::process::exit(2);
    }
    let b = graph::build(&a[1], &a[2], &a[3]);
    let n_ref_nodes = b.nodes.len();

    // Generation 0: one path node per edge, keyed by the from-node's label.
    let code = |l: u8| -> u32 { match l { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
                                          b'Y' => 4, _ => panic!("bad label") } };
    let mut nodes: Vec<PN> = b.edges.iter()
        .map(|&(f, t)| PN { from: f, to: t, k0: code(b.nodes[f as usize].0), k1: 0 })
        .collect();
    nodes.push(PN { from: b.last_node, to: b.last_node, k0: 5, k1: 0 });

    let mut curve: Vec<(u32, usize, usize, u32)> = vec![(0, nodes.len(), nodes.len(), 0)];
    println!("Generation 0 ({} -> {} nodes, 0 ranks)", nodes.len(), nodes.len());

    // generationOne sorts by `from` with a counting sort; everything after keeps
    // whatever order its own stage leaves behind.
    nodes.sort_by_key(|n| n.from);

    let mut generation = 0u32;
    let mut ranks = 0u32;
    loop {
        generation += 1;
        // For generations up to 4 the join reads the current list, which is in
        // `from` order. From 5 on the list is in rank order and the join reads a
        // separate `from`-sorted copy (`from_table`).
        let from_sorted: Vec<PN> = if generation <= 4 {
            nodes.clone()
        } else {
            let mut v = nodes.clone();
            v.sort_by_key(|n| n.from);
            v
        };
        let start = csr(&from_sorted, n_ref_nodes);

        let mut new: Vec<PN> = Vec::with_capacity(nodes.len() * 2);
        for node in &nodes {
            if generation > 4 && node.is_sorted() { new.push(*node); continue; }
            let (lo, hi) = (start[node.to as usize] as usize, start[node.to as usize + 1] as usize);
            for j in lo..hi {
                let o = &from_sorted[j];
                let (k0, k1) = if generation <= 3 {
                    // keys packed into one integer, 3 bits per character
                    let shift = 3u32 * (1u32 << (generation - 1));
                    ((node.k0 << shift) + o.k0, 0)
                } else {
                    (node.k0, o.k0)
                };
                new.push(PN { from: node.from, to: o.to, k0, k1 });
            }
        }
        let temp_nodes = new.len();
        nodes = new;

        if generation <= 3 {
            ranks = 0; // no pruning at all in the early generations
        } else {
            if generation == 4 { nodes.sort_by_key(|n| (n.k0, n.k1)); }
            ranks = merge_update_rank(&mut nodes, generation);
        }
        println!("Generation {generation} ({temp_nodes} -> {} nodes, {ranks} ranks)", nodes.len());
        curve.push((generation, temp_nodes, nodes.len(), ranks));
        if generation > 3 && ranks as usize == nodes.len() { break; }
        if generation > 64 { eprintln!("did not converge"); std::process::exit(1); }
    }

    // ---- generateEdges (gbwt_graph.h:2367) ------------------------------
    // With the nodes sorted, each reference edge contributes one GFM row per
    // path node whose `from` is that edge's `to`, labelled by the character of
    // the edge's `from` node. Rows are bucketed by label and sorted by the
    // ranking they point at -- which is the BWT order.
    nodes.sort_by_key(|n| n.from);
    for n in nodes.iter_mut() { n.to = b.nodes[n.from as usize].1; }
    let mut start = vec![0u32; n_ref_nodes + 1];
    for n in &nodes { start[n.from as usize + 1] += 1; }
    for i in 1..start.len() { start[i] += start[i - 1]; }

    let mut row_out: Vec<(u8, u8)> = Vec::new();
    let mut m_out: Vec<(u8, u32)> = Vec::new();
    let mut edge_from_count: HashMap<u32, u32> = HashMap::new();
    let mut floc_out: Vec<u32> = Vec::new();
    let mut buckets: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 6]; // (ranking, from)
    for &(f, t) in &b.edges {
        let li = match b.nodes[f as usize].0 {
            b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, b'Y' => 4, b'Z' => 5,
            _ => panic!("bad label"),
        };
        for j in start[t as usize] as usize..start[t as usize + 1] as usize {
            // PathEdge(edge->from, nodes[j].key.first, label) -- the stored
            // `from` is the REFERENCE edge's source, not the target path node's
            // own `from`. Only the merge-join reads it, so using the wrong one
            // leaves the BWT rows correct and the out-degrees wrong.
            buckets[li].push((nodes[j].k0, f));
        }
    }
    // radix by ranking, stable within a bucket
    for bk in buckets.iter_mut() { bk.sort_by_key(|&(r, _)| r); }
    let n_edges: usize = buckets.iter().map(|b| b.len()).sum();
    nodes.sort_by_key(|n| n.k0);

    println!("\ngenerateEdges: {} path nodes, {} path edges (A {} C {} G {} T {} Y {} Z {})",
             nodes.len(), n_edges, buckets[0].len(), buckets[1].len(), buckets[2].len(),
             buckets[3].len(), buckets[4].len(), buckets[5].len());
    // generateEdges drops a node further down, so record the count the header
    // should agree with before that happens.
    let n_path_nodes = nodes.len();
    println!("  numNodes = {} (path nodes - 1), gbwtLen = {} (path edges)",
             n_path_nodes - 1, n_edges);

    // ---- nextRow: the GFM row order (gbwt_graph.h:1609) -------------------
    // Rows are NOT the (label, ranking) order generateEdges leaves behind. They
    // are: path nodes in rank order, and within each node its INCOMING edges.
    // F marks each node's first row. M is a separate cursor over the same nodes
    // advancing by out-degree. So the tail of generateEdges has to run first:
    //   * rewrite edge.from from a reference-graph node id to a path-node index,
    //     and set each node's out-degree while doing it (:2561);
    //   * relabel 'Y' to 'Z' and drop the second-to-last node (:2576);
    //   * re-sort edges by (ranking, from) and build the CSR (:2604).
    {
        let mut ed: Vec<(u32, u32, u8)> = Vec::with_capacity(n_edges);  // (ranking, from, label)
        for (li, bk) in buckets.iter().enumerate() {
            let lab = b"ACGTYZ"[li];
            for &(rank, from) in bk { ed.push((rank, from, lab)); }
        }
        // The merge-join of gbwt_graph.h:2563, literally. Instrumenting a debug
        // hisat2-build on this graph shows the two lists arrive positionally
        // aligned -- edge i's `from` equals node i's `from` -- so the scan pairs
        // them off, advancing the node only on a mismatch and giving a node two
        // rows where two edges land on it.
        let mut outdeg = vec![0u32; nodes.len()];
        {
            let (mut ni, mut ei) = (0usize, 0usize);
            while ni < nodes.len() && ei < ed.len() {
                if ed[ei].1 == nodes[ni].from {
                    ed[ei].1 = ni as u32;
                    ei += 1;
                    outdeg[ni] += 1;
                } else {
                    ni += 1;
                    if ni < nodes.len() { outdeg[ni] = 0; }
                }
            }
            if ei < ed.len() { println!("  NOTE merge-join consumed {ei} of {} edges", ed.len()); }
        }
        // Drop the second-to-last node, keeping its out-degree on the survivor.
        let n = nodes.len();
        outdeg[n - 1] = outdeg[n - 2];
        nodes[n - 2] = nodes[n - 1];
        outdeg[n - 2] = outdeg[n - 1];
        nodes.pop(); outdeg.pop();
        for e in ed.iter_mut() {
            if e.2 == b'Y' { e.2 = b'Z'; }
            else if e.0 as usize >= nodes.len() { e.0 -= 1; }
        }
        // A STABLE sort on ranking alone. PathEdgeToCmp reads as (to, from), but
        // the radix sort is fed the edges in label-bucket order and preserves it
        // within a ranking, so ties resolve by label, not by `from`. Sorting by
        // (ranking, from) instead swaps pairs inside a node's edge group.
        ed.sort_by_key(|&(r, _, _)| r);

        // The M bitvector is a SECOND, independent cursor over the same nodes
        // (gbwt_graph.h:1630): it emits one run per node of length equal to that
        // node's OUT-degree, with M=1 on the first row of each run, and carries
        // the node's genomic position. It is not aligned with the F runs, which
        // are in-degree runs -- both walk the same nodes but at different rates.
        let mut mpos: Vec<(u8, u32)> = Vec::with_capacity(ed.len());
        for (i, nd) in nodes.iter().enumerate() {
            // A node with out-degree 0 still emits ONE row. nextRow sets M
            // before testing whether to advance (gbwt_graph.h:1632), so the run
            // length is max(1, outdegree) -- which is exactly what makes the row
            // count come out at gbwtLen when the dropped node leaves a zero
            // behind: 241 out-degrees over 240 nodes, 242 rows.
            let runs = outdeg[i].max(1);
            for k in 0..runs {
                mpos.push((if k == 0 { 1 } else { 0 }, nd.to));
            }
        }
        // F_loc for node k is the running sum of in-degrees before it
        // (nextFLocation, gbwt_graph.h:1640), sampled once per M==1 row.
        let mut indeg = vec![0u32; nodes.len()];
        for e in ed.iter() { if (e.0 as usize) < indeg.len() { indeg[e.0 as usize] += 1; } }
        let mut floc: Vec<u32> = Vec::with_capacity(nodes.len());
        let mut acc = 0u32;
        for i in 0..nodes.len() { floc.push(acc); acc += indeg[i]; }
        m_out = mpos;
        floc_out = floc;

        // rows in nextRow order, with the F bit
        let mut rows: Vec<(u8, u8)> = Vec::with_capacity(ed.len());
        let mut i = 0usize;
        for node in 0..nodes.len() as u32 {
            let mut first = true;
            while i < ed.len() && ed[i].0 == node {
                rows.push((ed[i].2, if first { 1 } else { 0 }));
                first = false; i += 1;
            }
        }
        println!("  nextRow order: {} rows ({} edges consumed of {})", rows.len(), i, ed.len());
        row_out = rows;
    }

    if a.len() >= 6 {
        // Cross-check against the header of the index HISAT2 actually wrote.
        let idx = fs::read(&a[5]).expect("index");
        let u = |o: usize| u32::from_le_bytes(idx[o..o + 4].try_into().unwrap()) as usize;
        let (their_gbwt, their_nodes) = (u(12), u(16));
        println!("  index header: gbwtLen {their_gbwt}, numNodes {their_nodes}");
        let mut ok = their_nodes == n_path_nodes - 1 && their_gbwt == n_edges;

        // fchr is a content check, not just a total: it says how many rows carry
        // each label, so it fails if the rows are right in number but wrong in
        // kind. Its offset is derived the way the reader derives it -- graph
        // indexes reserve 6 index_t per side and pack 2 rows per byte.
        let (line_rate, n_pat) = (u(20), u(44));
        let mut o = 48 + n_pat * 4;
        let n_frag = u(o);
        o += 4 + n_frag * 12;
        let side_sz = 1usize << line_rate;
        let side_gbwt_sz = side_sz - 24;
        let gbwt_sz = their_gbwt / 2 + 1;
        let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
        let rows_per_side = side_gbwt_sz * 2;
        let off_rate = u(28) as u32;
        o += num_sides * side_sz;
        let n_z = u(o);
        o += 4 + n_z * 4;
        let (zoff_at, fchr_at) = (o - 4 - n_z * 4, o);
        let fchr: Vec<usize> = (0..5).map(|i| u(o + i * 4)).collect();
        let theirs: Vec<usize> = (0..4).map(|i| fchr[i + 1] - fchr[i]).collect();
        let ours: Vec<usize> = (0..4).map(|i| buckets[i].len()).collect();
        println!("  fchr rows per label: ours {ours:?} theirs {theirs:?}");
        if ours != theirs { ok = false; }

        // Unpack the BWT the index actually stores and compare row by row.
        // Graph side layout, from gfm.h:4886: within each side's sideGbwtSz
        // bytes, [0, sz/2) is the BWT at 4 rows/byte low-pair-first,
        // [sz/2, 3sz/4) the F bitvector at 8 rows/byte with
        // F_bpi = bpi + ((sideCur & 1) << 2), then M the same way.
        let gbwt_off = 44 + 4 + n_pat * 4 + 4 + n_frag * 12;
        let sgs = side_gbwt_sz;
        let rows_per_side = sgs * 2;
        let getrow = |r: usize| -> (u8, u8) {
            let side = r / rows_per_side;
            let off = r % rows_per_side;
            let base = gbwt_off + side * side_sz;
            let sc = off >> 2;
            let bpi = off & 3;
            let ch = (idx[base + sc] >> (bpi * 2)) & 3;
            let f_sc = (sgs + sc) >> 1;
            let f_bpi = bpi + ((sc & 1) << 2);
            let f = (idx[base + f_sc] >> f_bpi) & 1;
            (ch, f)
        };
        let mut cbad = 0usize; let mut fbad = 0usize; let mut first_bad = None;
        for r in 0..row_out.len().min(their_gbwt) {
            let (sc, sf) = getrow(r);
            let ours_c = match row_out[r].0 { b'A' => 0u8, b'C' => 1, b'G' => 2, b'T' => 3,
                                             _ => 0 /* Z is stored as A and not counted */ };
            if sc != ours_c {
                cbad += 1;
                if first_bad.is_none() { first_bad = Some(r); }
                if cbad <= 8 {
                    println!("      row {r}: ours {} (F={}) theirs {} (F={})",
                             row_out[r].0 as char, row_out[r].1, b"ACGT"[sc as usize] as char, sf);
                }
            }
            if sf != row_out[r].1 { fbad += 1; }
        }
        let n = row_out.len().min(their_gbwt);
        println!("  BWT rows: {}/{} characters match, {}/{} F bits match{}",
                 n - cbad, n, n - fbad, n,
                 match first_bad { Some(r) => format!("; first char mismatch at row {r}"), None => String::new() });
        if cbad != 0 || fbad != 0 { ok = false; }

        // ---- emit the gbwt block and byte-compare it -----------------------
        // Graph side layout (gfm.h:4886, :4926). Within each side's sideGbwtSz
        // bytes: [0, sz/2) BWT at 4 rows/byte low-pair-first, [sz/2, 3sz/4) the
        // F bitvector at 8 rows/byte with F_bpi = bpi + ((sideCur & 1) << 2),
        // then M the same way. The final 6 index_t are F_locSave, M_occSave and
        // occSave[0..3] -- the values as of the START of the side, not the end.
        {
            let gbwt_tot = num_sides * side_sz;
            let mut blk = vec![0u8; gbwt_tot];
            let mut occ = [0u32; 4];
            let mut m_occ = 0u32;
            let mut f_loc = 0u32;
            let mut z_offs: Vec<u32> = Vec::new();
            let mut fchr_c = [0u32; 4];
            let mut sa_sample: Vec<u32> = Vec::new();
            let off_mask: u32 = u32::MAX << off_rate;
            let mut floc_i = 0usize;

            for si in 0..gbwt_tot / side_sz * rows_per_side {
                let side = si / rows_per_side;
                let off = si % rows_per_side;
                if off == 0 {
                    // The C++ writes these when a side FILLS, using occSave --
                    // the counts as they stood when that side opened. Writing at
                    // the side's start instead, the live counters already hold
                    // exactly those values, so the save variables are not just
                    // unnecessary here but wrong: they lag by a full side.
                    let base = side * side_sz + side_sz - 24;
                    for (k, v) in [f_loc, m_occ, occ[0], occ[1], occ[2], occ[3]].iter().enumerate() {
                        blk[base + k * 4..base + k * 4 + 4].copy_from_slice(&v.to_le_bytes());
                    }
                }
                let (ch, f, m, pos) = if si < row_out.len() {
                    let c = row_out[si].0;
                    let (m, p) = if si < m_out.len() { m_out[si] } else { (0, 0) };
                    (c, row_out[si].1, m, p)
                } else {
                    (b'A', 0, 0, 0)   // padding past the end, counted as 'A'
                };
                let mut count = true;
                let code = match ch {
                    b'A' => 0u8, b'C' => 1, b'G' => 2, b'T' => 3,
                    _ => { count = false; z_offs.push(si as u32); 0 }  // 'Z' is not representable
                };
                if si < row_out.len() && count { fchr_c[code as usize] += 1; }
                if m == 1 {
                    if floc_i < floc_out.len() { f_loc = floc_out[floc_i]; floc_i += 1; }
                    if (m_occ & off_mask) == m_occ { sa_sample.push(pos); }
                }
                if count { occ[code as usize] += 1; }
                if m == 1 { m_occ += 1; }

                let base = side * side_sz;
                let sc = off >> 2;
                let bpi = off & 3;
                blk[base + sc] |= code << (bpi * 2);
                let f_sc = (side_gbwt_sz + sc) >> 1;
                let f_bpi = bpi + ((sc & 1) << 2);
                blk[base + f_sc] |= f << f_bpi;
                let m_sc = f_sc + (side_gbwt_sz >> 2);
                blk[base + m_sc] |= m << f_bpi;
            }

            // The M bitvector is the one piece still open. Total M=1 rows match
            // (one per node), but the RUN LENGTHS do not: HISAT2 gives 238 nodes
            // a run of 1 and exactly two a run of 2, with no zeros, while
            // grouping edges by `from` produces zeros and twos in other places.
            // So out-degree is credited to a specific path node, not shared
            // across the nodes that happen to have the same `from` -- and the
            // merge-join at gbwt_graph.h:2561 cannot be read as a group-by,
            // since replaying it literally on label-ordered edges stalls after
            // exactly one bucket (81 of 242 edges).
            if std::env::var("HT2_MDBG").is_ok() {
                let getm = |r: usize| -> u8 {
                    let side = r / rows_per_side; let off = r % rows_per_side;
                    let sc = off >> 2; let bpi = off & 3;
                    let f_sc = (side_gbwt_sz + sc) >> 1;
                    let m_sc = f_sc + (side_gbwt_sz >> 2);
                    (idx[gbwt_off + side * side_sz + m_sc] >> (bpi + ((sc & 1) << 2))) & 1
                };
                let mut runs: Vec<usize> = Vec::new();
                let mut cur = 0usize;
                for i in 0..row_out.len() {
                    if getm(i) == 1 { if i > 0 { runs.push(cur); } cur = 0; }
                    cur += 1;
                }
                runs.push(cur);
                println!("      their runs: {} nodes, {} with length > 1 {:?}",
                         runs.len(), runs.iter().filter(|&&r| r > 1).count(),
                         runs.iter().enumerate().filter(|(_, &r)| r > 1)
                             .map(|(i, r)| (i, *r)).take(6).collect::<Vec<_>>());
                let mut oruns: Vec<usize> = Vec::new();
                let mut oc = 0usize;
                for i in 0..m_out.len() {
                    if m_out[i].0 == 1 { if i > 0 { oruns.push(oc); } oc = 0; }
                    oc += 1;
                }
                oruns.push(oc);
                println!("      ours  runs: {} nodes, {} with length > 1",
                         oruns.len(), oruns.iter().filter(|&&r| r > 1).count());
            }
            {
                let names = ["F_locSave", "M_occSave", "occA", "occC", "occG", "occT"];
                for side in 0..(gbwt_tot / side_sz).min(3) {
                    let base = side * side_sz + side_sz - 24;
                    let o: Vec<u32> = (0..6).map(|k| u32::from_le_bytes(
                        blk[base + k*4..base + k*4 + 4].try_into().unwrap())).collect();
                    let t: Vec<u32> = (0..6).map(|k| u32::from_le_bytes(
                        idx[gbwt_off + base + k*4..gbwt_off + base + k*4 + 4].try_into().unwrap())).collect();
                    if o != t {
                        println!("      side {side} tallies differ:");
                        for k in 0..6 { if o[k] != t[k] {
                            println!("        {:<10} ours {:>10} theirs {:>10}", names[k], o[k], t[k]); } }
                    }
                }
            }
            let theirs_blk = &idx[gbwt_off..gbwt_off + gbwt_tot];
            let ndiff = blk.iter().zip(theirs_blk.iter()).filter(|(a, b)| a != b).count();
            let firstd = blk.iter().zip(theirs_blk.iter()).position(|(a, b)| a != b);
            println!("  gbwt block: {}/{} bytes match{}", gbwt_tot - ndiff, gbwt_tot,
                     match firstd { Some(o) => format!("; first differing byte at {o} (side {}, offset {})",
                                                       o / side_sz, o % side_sz), None => String::new() });
            if ndiff != 0 { ok = false; }

            // fchr and zOffs, from the same walk
            let mut fchr = [0u32; 5];
            for i in 0..4 { fchr[i + 1] = fchr[i] + fchr_c[i]; }
            let their_fchr: Vec<u32> = (0..5).map(|i| u(fchr_at + i * 4) as u32).collect();
            let fok = (0..5).all(|i| fchr[i] == their_fchr[i]);
            println!("  fchr: ours {:?} theirs {:?} -> {}", fchr, their_fchr,
                     if fok { "match" } else { "DIFFER" });
            if !fok { ok = false; }
            let their_nz = u(zoff_at) ;
            let their_z: Vec<u32> = (0..their_nz).map(|i| u(zoff_at + 4 + i * 4) as u32).collect();
            let zok = z_offs == their_z;
            println!("  zOffs: ours {:?} theirs {:?} -> {}", z_offs, their_z,
                     if zok { "match" } else { "DIFFER" });
            if !zok { ok = false; }

            // .2.ht2 -- the SA sample, one position per 2^offRate M-marked rows
            if let Ok(sabuf) = fs::read(a[5].replace(".1.ht2", ".2.ht2")) {
                let n = (sabuf.len() - 4) / 4;
                let mut sbad = 0usize;
                for k in 0..n.min(sa_sample.len()) {
                    let t = u32::from_le_bytes(sabuf[4 + k * 4..8 + k * 4].try_into().unwrap());
                    if t != sa_sample[k] { sbad += 1; }
                }
                println!("  .2.ht2 SA sample: {}/{} match (ours {} entries, theirs {})",
                         n.min(sa_sample.len()) - sbad, n.min(sa_sample.len()), sa_sample.len(), n);
                if sbad != 0 || sa_sample.len() != n { ok = false; }
            }
        }

        // ---- graph ftab / eftab ---------------------------------------
        // Unlike the linear ftab, which falls out of the suffix-array walk,
        // the graph ftab is built by QUERYING the finished index: for each of
        // the 4^ftabChars prefixes, walk the GFM backward and record the row
        // range (gfm.h:4993). mapGLF is an LF step over the BWT followed by a
        // hop through the node structure -- rank over M to find which node the
        // row belongs to, then select over F to find where that node's incoming
        // edges begin.
        {
            let their_nz = u(zoff_at);
            let ftab_chars = u(32);
            let ftab_len = (1usize << (ftab_chars * 2)) + 1;
            let gl = their_gbwt;                       // gbwtLen
            let mut bwt = vec![0u8; gl];
            let mut fb = vec![0u8; gl];
            let mut mb = vec![0u8; gl];
            for r in 0..gl {
                let side = r / rows_per_side; let off = r % rows_per_side;
                let base = gbwt_off + side * side_sz;
                let sc = off >> 2; let bpi = off & 3;
                bwt[r] = (idx[base + sc] >> (bpi * 2)) & 3;
                let f_sc = (side_gbwt_sz + sc) >> 1;
                let f_bpi = bpi + ((sc & 1) << 2);
                fb[r] = (idx[base + f_sc] >> f_bpi) & 1;
                mb[r] = (idx[base + f_sc + (side_gbwt_sz >> 2)] >> f_bpi) & 1;
            }
            // occ excludes the 'Z' row: it is stored as 'A' but was never counted
            let zset: std::collections::HashSet<usize> =
                (0..their_nz).map(|i| u(zoff_at + 4 + i * 4)).collect();
            let mut occ = vec![[0u32; 4]; gl + 1];
            let mut rankm = vec![0u32; gl + 1];
            let mut self_f: Vec<u32> = Vec::new();
            for r in 0..gl {
                occ[r + 1] = occ[r];
                if !zset.contains(&r) { occ[r + 1][bwt[r] as usize] += 1; }
                rankm[r + 1] = rankm[r] + mb[r] as u32;
                if fb[r] == 1 { self_f.push(r as u32); }
            }
            let fchr_v: Vec<u32> = (0..5).map(|i| u(fchr_at + i * 4) as u32).collect();
            let map_lf = |row: usize, c: usize| -> u32 { fchr_v[c] + occ[row][c] };
            // rankm[k] is the number of M bits in rows [0, k) -- exclusive,
            // which is what rank_M(initFromRow_bit(x)) computes. Row r belongs
            // to node rankm[r+1] - 1.
            let rank_m = |k: usize| -> u32 { rankm[k.min(gl)] };
            let sel_f = |k: u32| -> u32 { *self_f.get(k as usize - 1).unwrap_or(&(gl as u32)) };

            let mut tftab: Vec<(u32, u32)> = Vec::with_capacity(ftab_len - 1);
            for i in 0..ftab_len - 1 {
                let mut q = i;
                let (mut top, mut bot) = (0u32, gl as u32);
                let mut j = 0usize;
                while j < ftab_chars {
                    let c = q & 3; q >>= 2;
                    let nt = map_lf(top as usize, c);
                    let nb = map_lf(bot as usize, c);
                    if nt >= nb { top = nt; bot = nb; break; }
                    let node_top = rank_m((nt + 1) as usize) - 1;
                    let node_bot = rank_m(nb as usize);
                    top = sel_f(node_top + 1);
                    bot = sel_f(node_bot + 1);
                    if top >= bot { break; }
                    j += 1;
                }
                if top >= bot || j < ftab_chars {
                    let v = if i == 0 { 0 } else { tftab[i - 1].1 };
                    tftab.push((v, v));
                } else {
                    tftab.push((top, bot));
                }
            }
            // finalise: ftab[i+1] = tFtab[i].second, with a gap pushed to eftab
            let mut ftab_o = vec![0u32; ftab_len];
            let mut eftab_o: Vec<u32> = Vec::new();
            ftab_o[0] = tftab[0].0; ftab_o[1] = tftab[0].1;
            for i in 1..ftab_len - 1 {
                if ftab_o[i] != tftab[i].0 {
                    let (lo, hi) = (ftab_o[i], tftab[i].0);
                    ftab_o[i] = (eftab_o.len() as u32 / 2) ^ u32::MAX;
                    eftab_o.push(lo); eftab_o.push(hi);
                }
                ftab_o[i + 1] = tftab[i].1;
            }
            let their_eftab_len = u(36);
            // ftab starts right after fchr[5]
            let ftab_at = fchr_at + 20;
            let eftab_at = ftab_at + ftab_len * 4;
            if std::env::var("HT2_FDBG").is_ok() {
                // decode their ftab into ranges: ftab[i] is the start of prefix
                // i's range (= end of prefix i-1's), with a gap pushed to eftab
                let raw = |i: usize| -> u32 { u(ftab_at + i * 4) as u32 };
                let lo_of = |i: usize| -> u32 {
                    let v = raw(i);
                    if v as usize <= gl { v } else { u(eftab_at + ((v ^ u32::MAX) as usize) * 8) as u32 }
                };
                let hi_of = |i: usize| -> u32 {
                    let v = raw(i);
                    if v as usize <= gl { v } else { u(eftab_at + ((v ^ u32::MAX) as usize) * 8 + 4) as u32 }
                };
                let mut shown = 0;
                for i in 0..ftab_len - 1 {
                    let theirs = (hi_of(i), lo_of(i + 1));
                    if theirs.0 != theirs.1 || shown < 4 {
                        if shown < 14 {
                            println!("      prefix {i}: ours {:?} theirs {:?}", tftab[i], theirs);
                            shown += 1;
                        } else { break; }
                    }
                }
            }
            let fbad = (0..ftab_len).filter(|&i| ftab_o[i] != u(ftab_at + i * 4) as u32).count();
            println!("  ftab: {}/{} entries match", ftab_len - fbad, ftab_len);
            let ebad = (0..their_eftab_len.min(eftab_o.len()))
                .filter(|&i| eftab_o[i] != u(eftab_at + i * 4) as u32).count();
            println!("  eftab: {}/{} match (ours {} entries, theirs {})",
                     their_eftab_len.min(eftab_o.len()) - ebad,
                     their_eftab_len.min(eftab_o.len()), eftab_o.len(), their_eftab_len);
            if fbad != 0 || ebad != 0 || eftab_o.len() != their_eftab_len { ok = false; }
        }

        println!("  {}", if ok { "GRAPH GFM MATCHES: sizes, label counts, BWT rows, F bits, gbwt block, fchr, zOffs, SA sample" }
                         else { "GRAPH GFM DIFFERS" });
        if !ok { std::process::exit(1); }
    }

    if a.len() < 5 { return; }
    let log = fs::read_to_string(&a[4]).expect("log");
    let theirs: Vec<(u32, usize, usize, u32)> = log.lines()
        .filter(|l| l.starts_with("Generation "))
        .map(|l| {
            let g: u32 = l["Generation ".len()..].split(' ').next().unwrap().parse().unwrap();
            let inner = &l[l.find('(').unwrap() + 1..];
            let t: Vec<&str> = inner.split_whitespace().collect();
            (g, t[0].parse().unwrap(), t[2].parse().unwrap(), t[4].parse().unwrap())
        }).collect();

    println!("\n  gen |          ours          |         HISAT2         |");
    let mut bad = 0;
    for i in 0..curve.len().max(theirs.len()) {
        let o = curve.get(i);
        let t = theirs.get(i);
        let f = |x: Option<&(u32, usize, usize, u32)>| match x {
            Some(v) => format!("{:>9} {:>9} {:>6}", v.1, v.2, v.3),
            None => "        -         -      -".to_string(),
        };
        let ok = o.is_some() && o == t;
        if !ok { bad += 1; }
        println!("  {:>3} | {} | {} | {}", i, f(o), f(t), if ok { "" } else { "<-- differs" });
    }
    if bad == 0 { println!("\nGENERATION CURVE MATCHES ({} generations)", curve.len()); }
    else { println!("\n{bad} generations differ"); std::process::exit(1); }
}
