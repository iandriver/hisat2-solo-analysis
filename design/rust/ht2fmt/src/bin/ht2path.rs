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
    let mut buckets: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 6]; // (ranking, from)
    for &(f, t) in &b.edges {
        let li = match b.nodes[f as usize].0 {
            b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, b'Y' => 4, b'Z' => 5,
            _ => panic!("bad label"),
        };
        for j in start[t as usize] as usize..start[t as usize + 1] as usize {
            buckets[li].push((nodes[j].k0, nodes[j].from));
        }
    }
    for bk in buckets.iter_mut() { bk.sort_unstable(); }
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
        // All edges leaving reference node X attach to the FIRST path node with
        // that `from`; any later one keeps out-degree 0, which is what the C++
        // merge-join does when several path nodes share a `from`.
        let mut first_of: HashMap<u32, u32> = HashMap::new();
        for (i, nd) in nodes.iter().enumerate() {
            first_of.entry(nd.from).or_insert(i as u32);
        }
        let mut outdeg = vec![0u32; nodes.len()];
        for e in ed.iter_mut() {
            if let Some(&ni) = first_of.get(&e.1) { e.1 = ni; outdeg[ni as usize] += 1; }
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
        o += ((gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz) * side_sz;
        let n_z = u(o);
        o += 4 + n_z * 4;
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

        println!("  {}", if ok { "GRAPH GFM MATCHES: sizes, label counts, BWT rows and F bits" }
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
