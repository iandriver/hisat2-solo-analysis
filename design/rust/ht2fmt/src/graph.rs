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

pub fn build(fa: &str, snp: &str, hap: &str) -> Built {
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

    let mut nodes: Vec<(u8, u32)> = Vec::with_capacity(len as usize + 2);
    let mut edges: Vec<(u32, u32)> = Vec::with_capacity(len as usize + 2);
    nodes.push((b'Y', 0));
    for i in 0..len {
        nodes.push((b"ACGT"[text[i as usize] as usize], i));
        edges.push((nodes.len() as u32 - 2, nodes.len() as u32 - 1));
    }
    nodes.push((b'Z', len));
    edges.push((nodes.len() as u32 - 2, nodes.len() as u32 - 1));
    let base_edges = edges.len();

    let mut out_of_order_haps = 0usize;
    for h in &haps {
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
                                let from = if j == h.left { alt.pos } else { nodes.len() as u32 - 2 };
                                edges.push((from, nodes.len() as u32 - 1));
                            }
                            if j == h.right { edges.push((nodes.len() as u32 - 1, alt.pos + 2)); }
                        }
                        ALT_DEL => {
                            let from = if j == h.left { alt.pos } else { nodes.len() as u32 - 1 };
                            j += alt.len - 1;
                            let to = if j == h.right { alt.pos + alt.len + 1 } else { nodes.len() as u32 };
                            edges.push((from, to));
                        }
                        _ => {
                            for k in 0..alt.len {
                                let bp = ((alt.seq >> ((alt.len - k - 1) * 2)) & 3) as usize;
                                nodes.push((b"ACGT"[bp], u32::MAX));
                                if prev == Some(ALT_DEL) && k == 0 { continue; }
                                let from = if k == 0 && j == h.left { alt.pos } else { nodes.len() as u32 - 2 };
                                edges.push((from, nodes.len() as u32 - 1));
                            }
                            if j == h.right { edges.push((nodes.len() as u32 - 1, alt.pos + 1)); }
                        }
                    }
                    id_i += 1;
                    prev = Some(alt.typ);
                    if alt.typ == ALT_INS { continue; }
                }
                _ => {
                    nodes.push((b"ACGT"[text[j as usize] as usize], j));
                    if prev != Some(ALT_DEL) {
                        let from = if j == h.left && prev.is_none() { j } else { nodes.len() as u32 - 2 };
                        edges.push((from, nodes.len() as u32 - 1));
                    }
                    if j == h.right { edges.push((nodes.len() as u32 - 1, j + 2)); }
                    prev = Some(ALT_SGL);
                }
            }
            j += 1;
        }
    }

    let (pre_nodes, pre_edges) = (nodes.len(), edges.len());
    let (nodes, edges, last_node) = reverse_determinize(&nodes, &edges, len + 1);
    Built {
        text, nodes, edges, last_node,
        n_alts: alts.len(), n_haps: haps.len(),
        dropped_snps, dropped_haps, out_of_order_haps,
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

    struct CNode { label: u8, value: u32, members: Vec<u32>, id: u32 }
    let mut cnodes: Vec<CNode> = vec![CNode {
        label: nodes[last_node as usize].0,
        value: nodes[last_node as usize].1,
        members: vec![last_node],
        id: 0,
    }];
    let mut cnode_map: HashMap<Vec<u32>, u32> = HashMap::new();
    cnode_map.insert(vec![last_node], 0);
    let mut active: VecDeque<u32> = VecDeque::from(vec![0]);
    let mut cedges: Vec<(u32, u32)> = Vec::new();
    let mut first_node = 0u32;
    let mut preds: Vec<u32> = Vec::new();

    while let Some(cid) = active.pop_front() {
        preds.clear();
        for k in 0..cnodes[cid as usize].members.len() {
            let m = cnodes[cid as usize].members[k];
            preds_of(m, &mut preds);
        }
        if preds.len() >= 2 {
            preds.sort_unstable();
            preds.dedup();
            // stable, so members stay in increasing id order — the member list
            // is the map key, so its order has to be fixed
            preds.sort_by_key(|&n| nodes[n as usize].0);
        }
        let mut i = 0usize;
        while i < preds.len() {
            let n0 = preds[i];
            let (label, mut value) = nodes[n0 as usize];
            let mut members = vec![n0];
            i += 1;
            if label == b'Y' && first_node == 0 { first_node = cnodes.len() as u32; }
            while i < preds.len() {
                let n = preds[i];
                let (l, v) = nodes[n as usize];
                if l != label { break; }
                members.push(n);
                if v != u32::MAX {
                    value = if value == u32::MAX { v } else { value.max(v) };
                }
                i += 1;
            }
            match cnode_map.get(&members) {
                None => {
                    let new_id = cnodes.len() as u32;
                    cnode_map.insert(members.clone(), new_id);
                    cnodes.push(CNode { label, value, members, id: 0 });
                    active.push_back(new_id);
                    cedges.push((new_id, cid));
                }
                Some(&existing) => cedges.push((existing, cid)),
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
    let out_edges: Vec<(u32, u32)> = cedges.iter()
        .map(|&(f, t)| (cnodes[f as usize].id, cnodes[t as usize].id)).collect();
    (out_nodes, out_edges, last_out)
}

/// The post-condition the construction exists to establish
/// (`assert(isReverseDeterministic(...))`, `gbwt_graph.h:796`).
pub fn reverse_deterministic(nodes: &[(u8, u32)], edges: &[(u32, u32)]) -> usize {
    let mut by_to: Vec<(u32, u8)> = edges.iter().map(|&(f, t)| (t, nodes[f as usize].0)).collect();
    by_to.sort_unstable();
    by_to.windows(2).filter(|w| w[0] == w[1]).count()
}
