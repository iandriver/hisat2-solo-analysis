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
/// Annotation ALTs. For these two, `Alt::pos` is `left` and `Alt::len` is
/// `right` -- an absolute coordinate, not a length -- because HISAT2 unions
/// the fields (`alt.h:56`). `Alt::seq` is the packed `excluded << 8 | fw`.
pub const ALT_SS:   u8 = 3;
pub const ALT_EXON: u8 = 4;

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
    /// per sequence: name, full length including ambiguous, and the maximal
    /// unambiguous runs as (chrom_off, joined_off, len). Local index windows are
    /// laid out in CHROM coordinates, so the mapping is needed to find which
    /// joined bases a window covers.
    pub seqs: Vec<(String, u32, Vec<(u32, u32, u32)>)>,
    pub text: Vec<u8>,
    pub alts: Vec<Alt>,
    pub haps: Vec<Hap>,
    /// The variant IDs, in the same order as `alts`, as one blob of
    /// newline-terminated names plus a start offset each. `.8.ht2` is exactly
    /// this blob behind a sentinel and a count, and 12.3M separate `String`s
    /// would cost more in headers than the names do in bytes.
    pub alt_names: Vec<u8>,
    pub alt_name_at: Vec<u32>,
    pub dropped_snps: usize,
    pub dropped_haps: usize,
    pub out_of_order_haps: usize,
}

/// Iterate a file's lines without holding it. `read_to_string` on a 3.2 GB
/// FASTA doubles peak memory for the length of the parse, which is the largest
/// single allocation in a whole-genome build and buys nothing.
fn for_each_line<F: FnMut(&str)>(path: &str, what: &str, mut f: F) {
    use std::io::BufRead;
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("{what}: {e}"));
    let mut r = std::io::BufReader::with_capacity(1 << 20, file);
    let mut buf: Vec<u8> = Vec::new();
    loop {
        buf.clear();
        if r.read_until(b'\n', &mut buf).expect(what) == 0 { break; }
        while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') { buf.pop(); }
        f(std::str::from_utf8(&buf).expect("utf8"));
    }
}

pub fn parse(fa: &str, snp: &str, hap: &str) -> Parsed {
    parse_with(fa, snp, hap, "", "")
}

/// As `parse`, plus the `--ss` and `--exon` files. Either may be `""`.
pub fn parse_with(fa: &str, snp: &str, hap: &str, ss: &str, exon: &str) -> Parsed {
    let mut text: Vec<u8> = Vec::new();
    let mut runs: HashMap<String, Vec<Run>> = HashMap::new();
    let mut cur = String::new();
    let mut chrom_off: u32 = 0;
    let mut in_run = false;
    // the true length of each sequence, ambiguity included -- runs only record
    // the unambiguous stretches, so a sequence ending in N would look short
    let mut plen: Vec<(String, u32)> = Vec::new();
    for_each_line(fa, "fa", |line| {
        if let Some(rest) = line.strip_prefix('>') {
            // close the previous sequence before adopting the new name
            if !cur.is_empty() { plen.push((cur.clone(), chrom_off)); }
            cur = rest.split_whitespace().next().unwrap_or("").to_string();
            runs.entry(cur.clone()).or_default();
            chrom_off = 0;
            in_run = false;
            return;
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
    });
    let len = text.len() as u32;
    let joined = |chrom: &str, pos: u32| -> Option<u32> {
        let rs = runs.get(chrom)?;
        let i = rs.partition_point(|r| r.chrom_off <= pos);
        if i == 0 { return None; }
        let r = &rs[i - 1];
        if pos - r.chrom_off >= r.len { None } else { Some(r.joined_off + (pos - r.chrom_off)) }
    };

    let mut alts: Vec<Alt> = Vec::new();
    let mut alt_names: Vec<u8> = Vec::new();
    let mut alt_name_at: Vec<u32> = Vec::new();
    let mut alt_id: HashMap<String, u32> = HashMap::new();
    let mut dropped_snps = 0usize;
    for_each_line(snp, "snp", |line| {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 { return; }
        let pos = match joined(f[2], f[3].parse::<u32>().expect("snp pos")) {
            Some(p) => p, None => { dropped_snps += 1; return; }
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
        alt_name_at.push(alt_names.len() as u32);
        alt_names.extend_from_slice(f[0].as_bytes());
        alt_names.push(b'\n');
        alts.push(Alt { pos, len: len_, seq, typ });
    });
    let mut haps: Vec<Hap> = Vec::new();
    let mut dropped_haps = 0usize;
    for_each_line(hap, "haplotype", |line| {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 { return; }
        let ids: Option<Vec<u32>> = f[4].split(',').map(|s| alt_id.get(s).copied()).collect();
        let (lc, rc) = (f[2].parse::<u32>().expect("hap left"),
                        f[3].parse::<u32>().expect("hap right"));
        let (l, r) = (joined(f[1], lc), joined(f[1], rc));
        // gfm.h's `inside_Ns`: a haplotype is dropped when either end falls in
        // an N run OR when the span crosses a fragment boundary. `joined`
        // already rejects an end inside an N run; the span test is the second
        // half, and it is what a single isolated N in the middle of a
        // haplotype trips. Joined coordinates close up the N run, so an
        // unbroken span is exactly one whose joined width still equals its
        // chromosome width.
        let same_frag = match (l, r) { (Some(a), Some(b)) => b - a == rc - lc, _ => false };
        match (ids, l, r) {
            (Some(ids), Some(left), Some(right)) if !ids.is_empty() && same_frag =>
                haps.push(Hap { left, right, alts: ids }),
            _ => dropped_haps += 1,
        }
    });

    // ---- --ss and --exon (gfm.h:1661 and :1787) --------------------------
    // Both files carry EXONIC coordinates, which gfm.h converts to intronic by
    // moving each end one base inward, and both are rejected unless the two
    // ends land in the same unambiguous run -- the `inside_Ns` walk, which here
    // is the same "joined width still equals chromosome width" test the
    // haplotypes use. `joined` standing in for `checkPosToSzs` is exact: both
    // ask whether a position is inside an unambiguous run of that sequence.
    let mut ss_seq: HashMap<u64, u32> = HashMap::new();
    if !ss.is_empty() {
        for_each_line(ss, "ss", |line| {
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 4 || f[0].starts_with('#') { return; }
            let (fl, fr) = (f[1].parse::<u32>().expect("ss left"),
                            f[2].parse::<u32>().expect("ss right"));
            if fl == 0 || fr == 0 { return; }
            let (cl, cr) = (fl + 1, fr - 1);
            if cl >= cr { return; }
            // gfm.h checks the file's own ends are in a run before converting.
            if joined(f[0], fl).is_none() || joined(f[0], fr).is_none() { return; }
            let (left, right) = match (joined(f[0], cl), joined(f[0], cr)) {
                (Some(l), Some(r)) if r - l == cr - cl => (l, r),
                _ => return,
            };
            // The repeat guard. gfm.h counts each junction's 16 bp + 16 bp
            // flank, and its `continue` sits INSIDE this window check, so a
            // junction repeating the one before it is dropped outright rather
            // than merely left uncounted -- but only when the window is in
            // range at all.
            const SEQLEN: u32 = 16;
            if left >= SEQLEN && (right + 1 + SEQLEN) as usize <= text.len() {
                let mut key = 0u64;
                for i in (left - SEQLEN)..left { key = key << 2 | text[i as usize] as u64; }
                for i in (right + 1)..(right + 1 + SEQLEN) { key = key << 2 | text[i as usize] as u64; }
                if let Some(prev) = alts.last() {
                    if prev.pos == left && prev.len == right { return; }
                }
                *ss_seq.entry(key).or_insert(0) += 1;
            }
            alt_name_at.push(alt_names.len() as u32);
            alt_names.extend_from_slice(b"ss\n");
            alts.push(Alt { pos: left, len: right, seq: (f[3] == "+") as u64, typ: ALT_SS });
        });
        // Second pass: exclude every junction whose flank pair is not unique.
        // Assigned only where the window is in range, so one that is not keeps
        // the `false` it was pushed with.
        const SEQLEN: u32 = 16;
        for a in alts.iter_mut() {
            if a.typ != ALT_SS { continue; }
            let (left, right) = (a.pos, a.len);
            if left >= SEQLEN && (right + 1 + SEQLEN) as usize <= text.len() {
                let mut key = 0u64;
                for i in (left - SEQLEN)..left { key = key << 2 | text[i as usize] as u64; }
                for i in (right + 1)..(right + 1 + SEQLEN) { key = key << 2 | text[i as usize] as u64; }
                if ss_seq.get(&key).copied().unwrap_or(0) > 1 { a.seq |= 1 << 8; }
            }
        }
    }
    if !exon.is_empty() {
        for_each_line(exon, "exon", |line| {
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 4 || f[0].starts_with('#') { return; }
            let (fl, fr) = (f[1].parse::<u32>().expect("exon left"),
                            f[2].parse::<u32>().expect("exon right"));
            if fl == 0 || fr == 0 { return; }
            let (cl, cr) = (fl + 1, fr - 1);
            if cl >= cr { return; }
            // No `checkPosToSzs` here: gfm.h:1787 runs the `inside_Ns` walk and
            // nothing else, and exons never reach the graph.
            let (left, right) = match (joined(f[0], cl), joined(f[0], cr)) {
                (Some(l), Some(r)) if r - l == cr - cl => (l, r),
                _ => return,
            };
            alt_name_at.push(alt_names.len() as u32);
            alt_names.extend_from_slice(b"exon\n");
            alts.push(Alt { pos: left, len: right, seq: (f[3] == "+") as u64, typ: ALT_EXON });
        });
    }
    alt_name_at.push(alt_names.len() as u32);

    if !cur.is_empty() { plen.push((cur.clone(), chrom_off)); }
    let mut seqs: Vec<(String, u32, Vec<(u32, u32, u32)>)> = Vec::new();
    for (name, full) in plen.iter() {
        let rs: Vec<(u32, u32, u32)> = runs.get(name)
            .map(|v| v.iter().map(|r| (r.chrom_off, r.joined_off, r.len)).collect())
            .unwrap_or_default();
        seqs.push((name.clone(), *full, rs));
    }
    // `_alts` is sorted before anything is built with it (gfm.h:1863), and the
    // names and every haplotype's alt indices are permuted to match. The
    // fixtures' variant files happen to arrive sorted, so this is a no-op on
    // them -- but a variant file that is not position-sorted would otherwise
    // build a different index here than `hisat2-build` does.
    {
        // HISAT2 orders by `type <` over its own enum -- SGL 1, INS 2, DEL 3,
        // SPLICESITE 5, EXON 6 -- with insertions forced ahead of everything at
        // the same position (`alt.h:91`). Only the order matters.
        let rank = |t: u8| -> u8 {
            match t {
                ALT_INS => 1, ALT_SGL => 2, ALT_DEL => 4, ALT_SS => 5, ALT_EXON => 6,
                other => panic!("unknown ALT type code {other}"),
            }
        };
        let mut ord: Vec<u32> = (0..alts.len() as u32).collect();
        // Stable, which is what pairing each ALT with its index achieves in the
        // C++: equal variants keep their file order.
        ord.sort_by_key(|&i| {
            let a = &alts[i as usize];
            (a.pos, rank(a.typ), a.len, a.seq)
        });
        if ord.iter().enumerate().any(|(i, &o)| i as u32 != o) {
            let mut new_alts: Vec<Alt> = Vec::with_capacity(alts.len());
            let mut new_names: Vec<u8> = Vec::with_capacity(alt_names.len());
            let mut new_at: Vec<u32> = Vec::with_capacity(alt_name_at.len());
            let mut back = vec![0u32; alts.len()];
            for (i, &o) in ord.iter().enumerate() {
                back[o as usize] = i as u32;
                let a = &alts[o as usize];
                new_alts.push(Alt { pos: a.pos, len: a.len, seq: a.seq, typ: a.typ });
                new_at.push(new_names.len() as u32);
                new_names.extend_from_slice(
                    &alt_names[alt_name_at[o as usize] as usize..alt_name_at[o as usize + 1] as usize]);
            }
            new_at.push(new_names.len() as u32);
            alts = new_alts; alt_names = new_names; alt_name_at = new_at;
            for h in haps.iter_mut() {
                for a in h.alts.iter_mut() { *a = back[*a as usize]; }
            }
        }
        // NOT a stable sort: see `gnusort` at the end of this file. `HT2_HAPSORT=stable`
        // restores the old behaviour for comparison.
        if std::env::var("HT2_HAPSORT").map(|v| v == "stable").unwrap_or(false) {
            haps.sort_by_key(|h| (h.left, h.right));
        } else {
            gnusort::sort(&mut haps, &|a: &Hap, b: &Hap| (a.left, a.right) < (b.left, b.right));
        }
    }
    Parsed { seqs, text, alts, haps, alt_names, alt_name_at,
             dropped_snps, dropped_haps, out_of_order_haps: 0 }
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
    // `curr_pos`/`curr_len` in the C++ are a reference FRAGMENT, not a chunk of
    // our choosing, and its `break` on a haplotype running past the window end
    // is unreachable once haplotypes spanning a fragment boundary are dropped
    // at parse time (below). So a plain containment filter is exact here.
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
                    if alt.typ == ALT_INS {
                        // The C++ is `for(j = left; j <= right; j++) { if(prev == INS) j--; ... }`,
                        // so the increment lands BEFORE the bound test and the decrement after it.
                        // An insertion on the haplotype's last position therefore ends the walk,
                        // and a variant co-located with it is never applied. A bare `continue`
                        // here re-tests the unchanged j and applies that variant -- two extra
                        // edges per such site, 94,265 of them genome-wide.
                        if j >= h.right { break; }
                        continue;
                    }
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
    let mut bad: Vec<(u32, u32)> = alts.iter().filter_map(|a| {
        let lo = a.pos.saturating_sub(RELAX + 1);
        let hi = match a.typ {
            ALT_DEL  => a.pos + a.len + RELAX,
            // `len` is `right`, so this is the whole intron. It has to be: the
            // splice-site edge runs from `left` to `right + 2`, and
            // `build_fragmented` can only hand an edge to the fragment directly
            // after the one that made it. This is where the annotation costs
            // something -- on GENCODE v32 it takes the longest uncuttable block
            // from 41 kb to 2.47 Mb.
            ALT_SS   => a.len + 1 + RELAX,
            // Exons are not in the graph at all, so they block nothing
            // (`gbwt_graph.h:404` skips them in exactly this loop).
            ALT_EXON => return None,
            // SGL is `pos + 1`; INS is `pos` in the C++ and one base wider
            // here. Deliberate: a wider range only ever pushes a cut later, and
            // where the cuts land is free (see `check_bounds`).
            _        => a.pos + 1 + RELAX,
        };
        Some((lo, hi.min(len)))
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


/// Every cut must be clear of every variant that contributes graph structure.
/// That is the only thing the stitched graph depends on: where the cuts land
/// beyond it is free, which is why ht2wg's 256 kb chunking reproduces
/// hisat2-build's 1 Mb chunking byte for byte. Cheap enough to run over a whole
/// genome, so it is the check for the annotation case rather than a diff
/// against boundaries the C++ happens to pick.
pub fn check_bounds(bounds: &[(u32, u32)], alts: &[Alt]) -> Result<(), String> {
    let cuts: Vec<u32> = bounds.iter().skip(1).map(|&(a, _)| a).collect();
    if cuts.windows(2).any(|w| w[0] >= w[1]) {
        return Err("fragment bounds are not strictly increasing".into());
    }
    for a in alts {
        let (lo, hi) = match a.typ {
            ALT_EXON => continue,
            ALT_DEL  => (a.pos.saturating_sub(RELAX + 1), a.pos + a.len + RELAX),
            ALT_SS   => (a.pos.saturating_sub(RELAX + 1), a.len + 1 + RELAX),
            _        => (a.pos.saturating_sub(RELAX + 1), a.pos + 1 + RELAX),
        };
        let i = cuts.partition_point(|&c| c < lo);
        if let Some(&c) = cuts.get(i) {
            if c <= hi {
                return Err(format!(
                    "cut at {c} falls inside a type-{} variant spanning [{lo}, {hi}]", a.typ));
            }
        }
    }
    Ok(())
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

/// libstdc++'s `std::sort`, reproduced exactly.
///
/// `_haplotypes` is ordered by `EList::sort` -> `std::sort` (`ds.h:830`) under a
/// comparator that returns false for full `(left, right)` ties (`alt.h:223`).
/// `std::sort` is NOT stable, so the order among tied haplotypes is unspecified
/// and differs by standard library -- hisat2-build's own `.7` is not
/// reproducible across platforms. The E3 index was built on Linux, so matching
/// it byte for byte means matching libstdc++ specifically; a stable sort agrees
/// with neither libstdc++ nor libc++ once ties exist.
///
/// This is introsort as libstdc++ implements it: quicksort with a median-of-3
/// pivot moved to `first`, a recursion depth limit of `2*floor(log2(n))` after
/// which it falls back to heapsort, a threshold of 16 below which partitioning
/// stops, and a final insertion-sort pass over the whole range. Every step has
/// to match, not just the algorithm class -- the tie order is a function of the
/// exact swap sequence.
pub mod gnusort {
    const THRESHOLD: usize = 16;

    /// Set by `sort` if the depth limit was ever hit. The heapsort fallback is
    /// the least-tested branch here, so it is worth knowing whether it ran.
    pub static HEAP_FALLBACKS: std::sync::atomic::AtomicUsize =
        std::sync::atomic::AtomicUsize::new(0);

    pub fn sort<T, F: Fn(&T, &T) -> bool>(v: &mut [T], lt: &F) {
        if v.is_empty() { return; }
        let n = v.len();
        let depth = 2 * (usize::BITS - 1 - n.leading_zeros());
        introsort_loop(v, 0, n, depth, lt);
        final_insertion_sort(v, lt);
    }

    fn introsort_loop<T, F: Fn(&T, &T) -> bool>(
        v: &mut [T], first: usize, mut last: usize, mut depth: u32, lt: &F)
    {
        while last - first > THRESHOLD {
            if depth == 0 {
                HEAP_FALLBACKS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                // __partial_sort(first, last, last) == make_heap then sort_heap
                make_heap(v, first, last, lt);
                sort_heap(v, first, last, lt);
                return;
            }
            depth -= 1;
            let cut = unguarded_partition_pivot(v, first, last, lt);
            introsort_loop(v, cut, last, depth, lt);
            last = cut;
        }
    }

    fn unguarded_partition_pivot<T, F: Fn(&T, &T) -> bool>(
        v: &mut [T], first: usize, last: usize, lt: &F) -> usize
    {
        let mid = first + (last - first) / 2;
        move_median_to_first(v, first, first + 1, mid, last - 1, lt);
        unguarded_partition(v, first + 1, last, first, lt)
    }

    fn move_median_to_first<T, F: Fn(&T, &T) -> bool>(
        v: &mut [T], result: usize, a: usize, b: usize, c: usize, lt: &F)
    {
        if lt(&v[a], &v[b]) {
            if lt(&v[b], &v[c]) { v.swap(result, b); }
            else if lt(&v[a], &v[c]) { v.swap(result, c); }
            else { v.swap(result, a); }
        } else if lt(&v[a], &v[c]) { v.swap(result, a); }
        else if lt(&v[b], &v[c]) { v.swap(result, c); }
        else { v.swap(result, b); }
    }

    /// The pivot sits at `pivot` (the range's original `first`) and the scan
    /// starts one past it, so no swap here ever moves the pivot element.
    fn unguarded_partition<T, F: Fn(&T, &T) -> bool>(
        v: &mut [T], mut first: usize, mut last: usize, pivot: usize, lt: &F) -> usize
    {
        loop {
            while lt(&v[first], &v[pivot]) { first += 1; }
            last -= 1;
            while lt(&v[pivot], &v[last]) { last -= 1; }
            if first >= last { return first; }
            v.swap(first, last);
            first += 1;
        }
    }

    fn final_insertion_sort<T, F: Fn(&T, &T) -> bool>(v: &mut [T], lt: &F) {
        let n = v.len();
        if n > THRESHOLD {
            insertion_sort(v, 0, THRESHOLD, lt);
            for i in THRESHOLD..n { linear_insert(v, i, lt); }
        } else {
            insertion_sort(v, 0, n, lt);
        }
    }

    fn insertion_sort<T, F: Fn(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, lt: &F) {
        if first == last { return; }
        for i in (first + 1)..last {
            if lt(&v[i], &v[first]) {
                // __move_backward3(first, i, i+1) then *first = val
                v[first..=i].rotate_right(1);
            } else {
                linear_insert(v, i, lt);
            }
        }
    }

    /// `__unguarded_linear_insert`: shift down while the held value is smaller.
    /// The swap form produces the same permutation as the C++ hole-and-move.
    fn linear_insert<T, F: Fn(&T, &T) -> bool>(v: &mut [T], last: usize, lt: &F) {
        let mut l = last;
        while l > 0 && lt(&v[l], &v[l - 1]) { v.swap(l, l - 1); l -= 1; }
    }

    // ---- the heapsort fallback -----------------------------------------
    fn adjust_heap<T, F: Fn(&T, &T) -> bool>(
        v: &mut [T], first: usize, mut hole: usize, len: usize, lt: &F)
    {
        let top = hole;
        let mut second = hole;
        while second < (len - 1) / 2 {
            second = 2 * (second + 1);
            if lt(&v[first + second], &v[first + second - 1]) { second -= 1; }
            v.swap(first + hole, first + second);
            hole = second;
        }
        if len & 1 == 0 && second == (len - 2) / 2 {
            second = 2 * (second + 1);
            v.swap(first + hole, first + second - 1);
            hole = second - 1;
        }
        // __push_heap
        while hole > top {
            let parent = (hole - 1) / 2;
            if !lt(&v[first + parent], &v[first + hole]) { break; }
            v.swap(first + hole, first + parent);
            hole = parent;
        }
    }

    fn make_heap<T, F: Fn(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, lt: &F) {
        let len = last - first;
        if len < 2 { return; }
        let mut parent = (len - 2) / 2;
        loop {
            adjust_heap(v, first, parent, len, lt);
            if parent == 0 { return; }
            parent -= 1;
        }
    }

    fn sort_heap<T, F: Fn(&T, &T) -> bool>(v: &mut [T], first: usize, mut last: usize, lt: &F) {
        while last - first > 1 {
            last -= 1;
            v.swap(first, last);              // __pop_heap into the vacated slot
            adjust_heap(v, first, 0, last - first, lt);
        }
    }
}
