//! `RefGraph` — the reference graph HISAT2 builds before any path doubling.
//!
//! Shared by `ht2graph` (which checks its size) and `ht2path` (which doubles
//! over it). See `gbwt_graph.h:606` for the walk and `:1015` for the subset
//! construction; the README records what each rule is for and how it was found.

use std::collections::{HashMap, VecDeque};
use std::fs;

pub const ALT_SGL: u8 = 0;
pub const ALT_DEL: u8 = 1;
pub const ALT_INS: u8 = 2;

pub struct Alt { pub pos: u32, pub len: u32, pub seq: u64, pub typ: u8 }
pub struct Hap { pub left: u32, pub right: u32, pub alts: Vec<u32> }

pub struct Built {
    pub text: Vec<u8>,
    pub nodes: Vec<(u8, u32)>,
    pub edges: Vec<(u32, u32)>,
    pub last_node: u32,
    pub n_alts: usize,
    pub n_haps: usize,
    pub dropped_snps: usize,
    pub dropped_haps: usize,
    pub out_of_order_haps: usize,
    pub pre_edges: usize,
    pub base_edges: usize,
    pub pre_nodes: usize,
}

/// One maximal unambiguous stretch, used to map chromosome coordinates onto the
/// joined text the graph is built over.
struct Run { chrom_off: u32, joined_off: u32, len: u32 }


/// Everything read from the input files, before any graph is built. Split out
/// so the fragmented builder can construct one range at a time from it without
/// re-reading or re-parsing.
pub struct Parsed {
    pub text: Vec<u8>,
    pub alts: Vec<Alt>,
    pub haps: Vec<Hap>,
    pub dropped_snps: usize,
    pub dropped_haps: usize,
    pub out_of_order_haps: usize,
}

pub fn parse(fa: &str, snp: &str, hap: &str) -> Parsed {
    let mut text: Vec<u8> = Vec::new();
    let mut runs: HashMap<String, Vec<Run>> = HashMap::new();
    let mut cur = String::new();
    let mut chrom_off: u32 = 0;
    let mut in_run = false;
    for line in fs::read_to_string(fa).expect("fa").lines() {
        if let Some(rest) = line.strip_prefix('>') {
            cur = rest.split_whitespace().next().unwrap_or("").to_string();
            runs.entry(cur.clone()).or_default();
            chrom_off = 0;
            in_run = false;
            continue;
        }
        for c in line.bytes() {
            let v = match c.to_ascii_uppercase() {
                b'A' => Some(0u8), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None,
            };
            match v {
                Some(v) => {
                    let r = runs.get_mut(&cur).expect("sequence before any header");
                    if !in_run {
                        r.push(Run { chrom_off, joined_off: text.len() as u32, len: 0 });
                        in_run = true;
                    }
                    r.last_mut().unwrap().len += 1;
                    text.push(v);
                }
                None => in_run = false,
            }
            chrom_off += 1;
        }
    }
    let len = text.len() as u32;
    let joined = |chrom: &str, pos: u32| -> Option<u32> {
        let rs = runs.get(chrom)?;
        let i = rs.partition_point(|r| r.chrom_off <= pos);
        if i == 0 { return None; }
        let r = &rs[i - 1];
        if pos - r.chrom_off >= r.len { None } else { Some(r.joined_off + (pos - r.chrom_off)) }
    };

    let mut alts: Vec<Alt> = Vec::new();
    let mut alt_id: HashMap<String, u32> = HashMap::new();
    let mut dropped_snps = 0usize;
    for line in fs::read_to_string(snp).expect("snp").lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 { continue; }
        let pos = match joined(f[2], f[3].parse::<u32>().expect("snp pos")) {
            Some(p) => p, None => { dropped_snps += 1; continue; }
        };
        let (typ, len_, seq) = match f[1] {
            "single" => (ALT_SGL, 1u32, "ACGT".find(f[4].as_bytes()[0] as char).expect("allele") as u64),
            "deletion" => (ALT_DEL, f[4].parse::<u32>().expect("del len"), 0u64),
            "insertion" => {
                let mut s = 0u64;
                for c in f[4].bytes() { s = (s << 2) | "ACGT".find(c as char).expect("ins base") as u64; }
                (ALT_INS, f[4].len() as u32, s)
            }
            other => panic!("unknown snp type '{other}'"),
        };
        alt_id.insert(f[0].to_string(), alts.len() as u32);
        alts.push(Alt { pos, len: len_, seq, typ });
    }
    let mut haps: Vec<Hap> = Vec::new();
    let mut dropped_haps = 0usize;
    for line in fs::read_to_string(hap).expect("haplotype").lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 { continue; }
        let ids: Option<Vec<u32>> = f[4].split(',').map(|s| alt_id.get(s).copied()).collect();
        let (l, r) = (joined(f[1], f[2].parse::<u32>().expect("hap left")),
                      joined(f[1], f[3].parse::<u32>().expect("hap right")));
        match (ids, l, r) {
            (Some(ids), Some(left), Some(right)) if !ids.is_empty() =>
                haps.push(Hap { left, right, alts: ids }),
            _ => dropped_haps += 1,
        }
    }
    Parsed { text, alts, haps, dropped_snps, dropped_haps, out_of_order_haps: 0 }
}

/// Build the nodes and edges for reference range `[a, b)` with LOCAL indices:
/// node 0 is the anchor standing for position `a - 1`, positions `a..b` occupy
/// `1..=b-a`, and `b - a + 1` is the tail. Mirrors `RefGraph`'s per-fragment
/// indexing (`j - curr_pos + 2`, `gbwt_graph.h:1092`).
///
/// Returns (nodes, edges, count of backbone edges).
pub fn build_range(p: &Parsed, a: u32, b: u32, first: bool, last: bool)
    -> (Vec<(u8, u32)>, Vec<(u32, u32)>, usize)
{
    build_range_with(&p.text, &p.alts, &p.haps, a, b, first, last)
}

/// As `build_range`, but over an explicit variant set -- which the local-index
/// retry needs, since a graph that explodes is rebuilt from a thinned list.
pub fn build_range_with(text: &[u8], alts: &[Alt], in_haps: &[Hap],
                        a: u32, b: u32, first: bool, last: bool)
    -> (Vec<(u8, u32)>, Vec<(u32, u32)>, usize)
{
    let span = (b - a) as usize;
    let mut nodes: Vec<(u8, u32)> = Vec::with_capacity(span + 2);
    let mut edges: Vec<(u32, u32)> = Vec::with_capacity(span + 2);
    // local 0: the real head on the first fragment, else the previous
    // fragment's last backbone node
    // The anchor is always labelled 'Y'. reverse_determinize locates the head
    // by that label, so a fragment without one has no start for its forward
    // numbering pass. Using 'Y' rather than the real base is safe because cuts
    // sit at least RELAX bases from any variant, so the anchor has exactly one
    // successor and no merge can depend on its label — and for fi > 0 the node
    // is never emitted anyway, the previous fragment owns it.
    let _ = first;
    nodes.push((b'Y', a.saturating_sub(1)));
    for i in a..b {
        nodes.push((b"ACGT"[text[i as usize] as usize], i));
        edges.push((nodes.len() as u32 - 2, nodes.len() as u32 - 1));
    }
    // local b-a+1: the real tail on the last fragment, else a stand-in for the
    // next fragment's first node, which that fragment will own
    // Both cases push 'Z': on the last fragment it is the real tail, otherwise
    // it is a stand-in that reverse_determinize can walk back from and that the
    // stitch drops, since the next fragment owns that position.
    let _ = last;
    nodes.push((b'Z', b));
    edges.push((nodes.len() as u32 - 2, nodes.len() as u32 - 1));
    let base_edges = edges.len();

    // haplotypes wholly inside this range, with positions shifted local
    let off = a;
    let haps: Vec<&Hap> = in_haps.iter().filter(|h| h.left >= a && h.right < b).collect();
    let mut out_of_order_haps = 0usize;
    for h in haps.into_iter() {
        let mut pass = true;
        for w in 0..h.alts.len().saturating_sub(1) {
            let s1 = &alts[h.alts[w] as usize];
            let s2 = &alts[h.alts[w + 1] as usize];
            let bad = match s1.typ {
                ALT_INS => s1.pos > s2.pos,
                ALT_DEL if s2.typ == ALT_DEL => s1.pos + s1.len >= s2.pos,
                ALT_DEL => s1.pos + s1.len - 1 >= s2.pos,
                _ => s1.pos >= s2.pos,
            };
            if bad { pass = false; break; }
        }
        if !pass { out_of_order_haps += 1; continue; }

        let mut prev: Option<u8> = None;
        let mut id_i = 0usize;
        let mut j = h.left;
        while j <= h.right {
            let altp = h.alts.get(id_i).map(|&x| &alts[x as usize]);
            match altp {
                Some(alt) if alt.pos == j => {
                    match alt.typ {
                        ALT_SGL => {
                            nodes.push((b"ACGT"[alt.seq as usize], alt.pos));
                            if prev != Some(ALT_DEL) {
                                let from = if j == h.left { alt.pos - off } else { nodes.len() as u32 - 2 };
                                edges.push((from, nodes.len() as u32 - 1));
                            }
                            if j == h.right { edges.push((nodes.len() as u32 - 1, alt.pos - off + 2)); }
                        }
                        ALT_DEL => {
                            let from = if j == h.left { alt.pos - off } else { nodes.len() as u32 - 1 };
                            j += alt.len - 1;
                            let to = if j == h.right { alt.pos - off + alt.len + 1 } else { nodes.len() as u32 };
                            edges.push((from, to));
                        }
                        _ => {
                            for k in 0..alt.len {
                                let bp = ((alt.seq >> ((alt.len - k - 1) * 2)) & 3) as usize;
                                nodes.push((b"ACGT"[bp], u32::MAX));
                                if prev == Some(ALT_DEL) && k == 0 { continue; }
                                let from = if k == 0 && j == h.left { alt.pos - off } else { nodes.len() as u32 - 2 };
                                edges.push((from, nodes.len() as u32 - 1));
                            }
                            if j == h.right { edges.push((nodes.len() as u32 - 1, alt.pos - off + 1)); }
                        }
                    }
                    id_i += 1;
                    prev = Some(alt.typ);
                    if alt.typ == ALT_INS { continue; }
                }
                _ => {
                    nodes.push((b"ACGT"[text[j as usize] as usize], j));
                    if prev != Some(ALT_DEL) {
                        let from = if j == h.left && prev.is_none() { j - off } else { nodes.len() as u32 - 2 };
                        edges.push((from, nodes.len() as u32 - 1));
                    }
                    if j == h.right { edges.push((nodes.len() as u32 - 1, j - off + 2)); }
                    prev = Some(ALT_SGL);
                }
            }
            j += 1;
        }
    }

    (nodes, edges, base_edges)
}

pub fn build(fa: &str, snp: &str, hap: &str) -> Built {
    let p = parse(fa, snp, hap);
    let len = p.text.len() as u32;
    let (nodes, edges, base_edges) = build_range(&p, 0, len, true, true);
    let (pre_nodes, pre_edges) = (nodes.len(), edges.len());
    let (nodes, edges, last_node) = if std::env::var("HT2_NO_DET").is_ok() {
        (nodes, edges, len + 1)
    } else {
        reverse_determinize(&nodes, &edges, len + 1)
    };
    Built {
        text: p.text, nodes, edges, last_node,
        n_alts: p.alts.len(), n_haps: p.haps.len(),
        dropped_snps: p.dropped_snps, dropped_haps: p.dropped_haps,
        out_of_order_haps: p.out_of_order_haps,
        pre_edges, base_edges, pre_nodes,
    }
}


/// `RefGraph::reverseDeterminize` (`gbwt_graph.h:1015`).
///
/// A GCSA-style index requires the graph to be reverse deterministic: no node
/// may have two incoming edges from nodes carrying the same label. The naive
/// haplotype walk violates that on purpose — an insertion always emits a
/// duplicate of the reference node it lands on — so this subset construction is
/// not an optional tidy-up, it is part of the definition of the graph.
///
/// Runs BACKWARD from the tail 'Z': repeatedly take a composite node, collect
/// the predecessors of all its members, group them by label, and make each
/// group a composite node. Identical member sets are shared, which is what
/// collapses the duplicates.
pub fn reverse_determinize(nodes: &[(u8, u32)], edges: &[(u32, u32)], last_node: u32)
    -> (Vec<(u8, u32)>, Vec<(u32, u32)>, u32)
{
    let mut by_to: Vec<(u32, u32)> = edges.iter().map(|&(f, t)| (t, f)).collect();
    by_to.sort_unstable();
    let preds_of = |n: u32, out: &mut Vec<u32>| {
        let i = by_to.partition_point(|&(t, _)| t < n);
        for &(t, f) in &by_to[i..] { if t != n { break; } out.push(f); }
    };

    // Member lists live in one arena rather than a Vec per composite node.
    // Measured on a 900 kb graph, determinisation merged 268 nodes out of
    // 903,755 -- 99.97% of composite nodes carry exactly one member, so a
    // per-node heap allocation plus a hash entry was spending ~100 bytes to
    // hold a single u32.
    struct CNode { label: u8, value: u32, mstart: u32, mlen: u32, id: u32 }
    let mut arena: Vec<u32> = Vec::with_capacity(nodes.len() + 16);
    let mut cnodes: Vec<CNode> = Vec::with_capacity(nodes.len());

    // Singletons resolve through a flat table indexed by node id; only genuine
    // multi-member sets reach the map.
    let mut single: Vec<u32> = vec![u32::MAX; nodes.len()];
    let mut multi: HashMap<Box<[u32]>, u32> = HashMap::new();

    arena.push(last_node);
    cnodes.push(CNode {
        label: nodes[last_node as usize].0,
        value: nodes[last_node as usize].1,
        mstart: 0, mlen: 1, id: 0,
    });
    single[last_node as usize] = 0;

    let mut active: VecDeque<u32> = VecDeque::from(vec![0]);
    let mut cedges: Vec<(u32, u32)> = Vec::new();
    let mut first_node = 0u32;
    let mut preds: Vec<u32> = Vec::new();

    while let Some(cid) = active.pop_front() {
        preds.clear();
        let (ms, ml) = (cnodes[cid as usize].mstart as usize, cnodes[cid as usize].mlen as usize);
        for k in 0..ml {
            let m = arena[ms + k];
            preds_of(m, &mut preds);
        }
        if preds.len() >= 2 {
            preds.sort_unstable();
            preds.dedup();
            // stable, so members stay in increasing id order -- the member list
            // is the map key, so its order has to be fixed
            preds.sort_by_key(|&n| nodes[n as usize].0);
        }
        let mut i = 0usize;
        while i < preds.len() {
            let n0 = preds[i];
            let (label, mut value) = nodes[n0 as usize];
            let mstart = arena.len() as u32;
            arena.push(n0);
            i += 1;
            if label == b'Y' && first_node == 0 { first_node = cnodes.len() as u32; }
            while i < preds.len() {
                let n = preds[i];
                let (l, v) = nodes[n as usize];
                if l != label { break; }
                arena.push(n);
                if v != u32::MAX {
                    value = if value == u32::MAX { v } else { value.max(v) };
                }
                i += 1;
            }
            let mlen = arena.len() as u32 - mstart;

            let existing = if mlen == 1 {
                let e = single[n0 as usize];
                if e == u32::MAX { None } else { Some(e) }
            } else {
                multi.get(&arena[mstart as usize..]).copied()
            };
            match existing {
                None => {
                    let new_id = cnodes.len() as u32;
                    if mlen == 1 {
                        single[n0 as usize] = new_id;
                    } else {
                        multi.insert(arena[mstart as usize..].into(), new_id);
                    }
                    cnodes.push(CNode { label, value, mstart, mlen, id: 0 });
                    active.push_back(new_id);
                    cedges.push((new_id, cid));
                }
                Some(e) => {
                    arena.truncate(mstart as usize); // candidate not kept
                    cedges.push((e, cid));
                }
            }
            cnodes[cid as usize].id += 1; // in-degree, consumed by the ordering pass
        }
    }

    let mut fwd: Vec<(u32, u32)> = cedges.clone();
    fwd.sort_unstable();
    let mut out_nodes: Vec<(u8, u32)> = Vec::with_capacity(cnodes.len());
    let mut last_out = 0u32;
    cnodes[first_node as usize].id = 0;
    out_nodes.push((cnodes[first_node as usize].label, cnodes[first_node as usize].value));
    let mut act: VecDeque<u32> = VecDeque::from(vec![first_node]);
    while let Some(cid) = act.pop_front() {
        let mut i = fwd.partition_point(|&(f, _)| f < cid);
        while i < fwd.len() && fwd[i].0 == cid {
            let succ = fwd[i].1;
            cnodes[succ as usize].id -= 1;
            if cnodes[succ as usize].id == 0 {
                act.push_back(succ);
                cnodes[succ as usize].id = out_nodes.len() as u32;
                out_nodes.push((cnodes[succ as usize].label, cnodes[succ as usize].value));
                if cnodes[succ as usize].label == b'Z' { last_out = out_nodes.len() as u32 - 1; }
            }
            i += 1;
        }
    }
    drop(fwd);
    let mut out_edges: Vec<(u32, u32)> = cedges.iter()
        .map(|&(f, t)| (cnodes[f as usize].id, cnodes[t as usize].id)).collect();
    // `reverseDeterminize` ends with sortEdgesFrom (gbwt_graph.h:2409). The
    // doubling re-sorts anyway, so the global path never needed it -- but the
    // fragmented path emits edges per fragment, and without this the two differ
    // in order while agreeing as multisets.
    out_edges.sort_unstable();
    (out_nodes, out_edges, last_out)
}


/// Fragment boundaries for step 3.5.
///
/// `RefGraph` switches to per-fragment automata at `jlen >= 1 << 16`
/// (`gbwt_graph.h:383`) and determinises each chunk separately. That is not a
/// parallelisation detail: it is what bounds the determinisation working set to
/// one chunk instead of one genome, and reproducing only its output while
/// ignoring its structure is how an 812 GB projection got into the plan.
///
/// Cuts must land where no composite node can span them, so a boundary is
/// pushed forward until it is at least `RELAX` bases clear of every variant —
/// the same 128 the C++ uses when computing its alt-free ranges.
pub const RELAX: u32 = 128;

pub fn fragment_bounds(len: u32, alts: &[Alt], chunk: u32) -> Vec<(u32, u32)> {
    // positions that must not be cut through, as sorted inclusive ranges
    let mut bad: Vec<(u32, u32)> = alts.iter().map(|a| {
        let lo = a.pos.saturating_sub(RELAX + 1);
        let hi = match a.typ {
            ALT_DEL => a.pos + a.len + RELAX,
            _       => a.pos + 1 + RELAX,
        };
        (lo, hi.min(len))
    }).collect();
    bad.sort_unstable();
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(bad.len());
    for r in bad {
        match merged.last_mut() {
            Some(l) if r.0 <= l.1 => { if r.1 > l.1 { l.1 = r.1; } }
            _ => merged.push(r),
        }
    }
    let blocked = |p: u32| -> Option<u32> {
        let i = merged.partition_point(|&(lo, _)| lo <= p);
        if i == 0 { return None; }
        let (lo, hi) = merged[i - 1];
        if p >= lo && p <= hi { Some(hi + 1) } else { None }
    };

    let mut out = Vec::new();
    let mut a = 0u32;
    while a < len {
        let mut b = a.saturating_add(chunk).min(len);
        if b < len {
            while let Some(next) = blocked(b) {
                b = next.min(len);
                if b >= len { break; }
            }
        }
        if b <= a { b = len; }
        out.push((a, b));
        a = b;
    }
    out
}


/// Step 3.5 — build the graph one fragment at a time.
///
/// Each fragment is determinised standalone with everything resident, its nodes
/// and edges are appended to the running output, and the fragment is dropped.
/// Peak becomes one chunk instead of one genome.
///
/// Two things make the stitch exact. Cuts sit at least `RELAX` bases from every
/// variant, so the graph is a plain chain there and no composite node can span a
/// boundary. And the shared boundary node is located by its genomic `value`,
/// which determinisation preserves and which — in a variant-free stretch — only
/// the backbone node carries.
pub fn build_fragmented(fa: &str, snp: &str, hap: &str, chunk: u32) -> Built {
    let p = parse(fa, snp, hap);
    let len = p.text.len() as u32;
    let bounds = fragment_bounds(len, &p.alts, chunk);

    let mut g_nodes: Vec<(u8, u32)> = Vec::new();
    let mut g_edges: Vec<(u32, u32)> = Vec::new();
    let mut pre_nodes = 0usize;
    let mut pre_edges = 0usize;
    let mut base_edges = 0usize;
    let mut anchor_global: u32 = u32::MAX; // global id of the node for position a-1
    let mut pending: Vec<u32> = Vec::new(); // edges into a stand-in tail, awaiting the next fragment
    let mut last_node_global: u32 = 0;

    for (fi, &(a, b)) in bounds.iter().enumerate() {
        let first = fi == 0;
        let last = fi + 1 == bounds.len();
        let (mut n, mut e, base_e) = build_range(&p, a, b, first, last);
        pre_nodes += n.len();
        pre_edges += e.len();
        base_edges += base_e;

        // tail to walk back from: the real 'Z' on the last fragment, otherwise
        // the synthetic anchor standing in for the next fragment's first node
        let tail_local = (b - a + 1) as u32;
        let (dn, de, dlast) = reverse_determinize(&n, &e, tail_local);
        n.clear(); e.clear();

        // local -> global. Node 0 of a determinised fragment is its head, which
        // for fi>0 is the boundary node the previous fragment already emitted.
        let mut map: Vec<u32> = vec![u32::MAX; dn.len()];
        let mut emitted_boundary: u32 = u32::MAX;
        let mut tail_local: Vec<u32> = Vec::new();
        for (li, &(lab, val)) in dn.iter().enumerate() {
            if lab == b'Y' {
                map[li] = if first { let g = g_nodes.len() as u32; g_nodes.push((lab, val)); g }
                          else { anchor_global };
                continue;
            }
            if lab == b'Z' && !last {
                tail_local.push(li as u32);   // stand-in; the next fragment owns it
                continue;
            }
            let g = g_nodes.len() as u32;
            // The first real node of this fragment resolves the previous
            // fragment's edge into its stand-in tail.
            if val == a && !pending.is_empty() {
                for &src in pending.iter() { g_edges.push((src, g)); }
                pending.clear();
            }
            map[li] = g;
            g_nodes.push((lab, val));
            if lab == b'Z' { last_node_global = g; }
            if !last && val == b - 1 { emitted_boundary = g; }
        }
        let tail_set: std::collections::HashSet<u32> = tail_local.into_iter().collect();
        // The anchor's outgoing backbone edge is the SAME edge the previous
        // fragment handed over as pending -- both run from the boundary node to
        // this fragment's first base. Emitting the fragment's own copy as well
        // double-counts it, one edge per boundary.
        let anchor_local: u32 = dn.iter().position(|&(l, _)| l == b'Y').unwrap_or(usize::MAX) as u32;
        for &(f, t) in de.iter() {
            if !first && f == anchor_local { continue; }
            if tail_set.contains(&t) {
                // backbone edge into the stand-in tail: hold it until the next
                // fragment emits the node it really points at
                if map[f as usize] != u32::MAX { pending.push(map[f as usize]); }
                continue;
            }
            let (gf, gt) = (map[f as usize], map[t as usize]);
            if gf == u32::MAX || gt == u32::MAX { continue; }
            g_edges.push((gf, gt));
        }
        let _ = dlast;
        anchor_global = emitted_boundary;
    }

    g_edges.sort_unstable();
    Built {
        text: p.text, nodes: g_nodes, edges: g_edges, last_node: last_node_global,
        n_alts: p.alts.len(), n_haps: p.haps.len(),
        dropped_snps: p.dropped_snps, dropped_haps: p.dropped_haps,
        out_of_order_haps: p.out_of_order_haps,
        pre_edges, base_edges, pre_nodes,
    }
}


/// Step 1 — build the graph fragment by fragment and stream it to disk.
///
/// Same construction as `build_fragmented`, but nodes and edges are appended to
/// record files instead of Vecs, so nothing accumulates. What stays resident is
/// one fragment's working set, the reference text, and the variant tables —
/// none of which grow with how much of the genome has already been processed.
///
/// Writes `nodes.bin` (label `u8` + value `u32`, 5 bytes) and `edges.bin`
/// (`from`, `to` as `u32`, 8 bytes, in `sortEdgesFrom` order).
pub struct GraphOnDisk { pub n_nodes: u64, pub n_edges: u64, pub last_node: u32, pub text_len: u32 }

pub fn build_fragmented_to_disk(fa: &str, snp: &str, hap: &str, chunk: u32, dir: &std::path::Path)
    -> std::io::Result<GraphOnDisk>
{
    use std::io::Write;
    let p = parse(fa, snp, hap);
    let len = p.text.len() as u32;
    let bounds = fragment_bounds(len, &p.alts, chunk);

    let mut nw = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(dir.join("nodes.bin"))?);
    let raw_edges = dir.join("edges.raw");
    let mut ew = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&raw_edges)?);
    let (mut n_nodes, mut n_edges) = (0u64, 0u64);
    let mut anchor_global: u32 = u32::MAX;
    let mut pending: Vec<u32> = Vec::new();
    let mut last_node_global: u32 = 0;

    for (fi, &(a, b)) in bounds.iter().enumerate() {
        let first = fi == 0;
        let last = fi + 1 == bounds.len();
        let (n, e, _) = build_range(&p, a, b, first, last);
        let tail_local = (b - a + 1) as u32;
        let (dn, de, _) = reverse_determinize(&n, &e, tail_local);
        drop(n); drop(e);

        let mut map: Vec<u32> = vec![u32::MAX; dn.len()];
        let mut emitted_boundary: u32 = u32::MAX;
        let mut tail_local_ids: Vec<u32> = Vec::new();
        for (li, &(lab, val)) in dn.iter().enumerate() {
            if lab == b'Y' {
                if first {
                    map[li] = n_nodes as u32;
                    nw.write_all(&[lab])?; nw.write_all(&val.to_le_bytes())?; n_nodes += 1;
                } else { map[li] = anchor_global; }
                continue;
            }
            if lab == b'Z' && !last { tail_local_ids.push(li as u32); continue; }
            let g = n_nodes as u32;
            if val == a && !pending.is_empty() {
                for &src in pending.iter() {
                    ew.write_all(&src.to_le_bytes())?; ew.write_all(&g.to_le_bytes())?; n_edges += 1;
                }
                pending.clear();
            }
            map[li] = g;
            nw.write_all(&[lab])?; nw.write_all(&val.to_le_bytes())?; n_nodes += 1;
            if lab == b'Z' { last_node_global = g; }
            if !last && val == b - 1 { emitted_boundary = g; }
        }
        let tail_set: std::collections::HashSet<u32> = tail_local_ids.into_iter().collect();
        let anchor_local: u32 = dn.iter().position(|&(l, _)| l == b'Y').unwrap_or(usize::MAX) as u32;
        for &(f, t) in de.iter() {
            if !first && f == anchor_local { continue; }
            if tail_set.contains(&t) {
                if map[f as usize] != u32::MAX { pending.push(map[f as usize]); }
                continue;
            }
            let (gf, gt) = (map[f as usize], map[t as usize]);
            if gf == u32::MAX || gt == u32::MAX { continue; }
            ew.write_all(&gf.to_le_bytes())?; ew.write_all(&gt.to_le_bytes())?; n_edges += 1;
        }
        anchor_global = emitted_boundary;
    }
    nw.flush()?; ew.flush()?; drop(nw); drop(ew);

    // sortEdgesFrom, externally
    super::ext::sort_pairs_external(&raw_edges, &dir.join("edges.bin"), 1 << 20, dir)?;
    let _ = std::fs::remove_file(&raw_edges);
    Ok(GraphOnDisk { n_nodes, n_edges, last_node: last_node_global, text_len: len })
}

/// The post-condition the construction exists to establish
/// (`assert(isReverseDeterministic(...))`, `gbwt_graph.h:796`).
pub fn reverse_deterministic(nodes: &[(u8, u32)], edges: &[(u32, u32)]) -> usize {
    let mut by_to: Vec<(u32, u8)> = edges.iter().map(|&(f, t)| (t, nodes[f as usize].0)).collect();
    by_to.sort_unstable();
    by_to.windows(2).filter(|w| w[0] == w[1]).count()
}
