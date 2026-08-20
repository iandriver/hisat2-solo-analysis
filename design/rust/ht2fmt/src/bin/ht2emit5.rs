//! Emit `.5.ht2` / `.6.ht2` — the hierarchical local indexes — and byte-compare.
//!
//! Each local index is a GFM over a 57,344 bp window at `u16` width with
//! `offRate` 3 and `ftabChars` 6, so this reuses the construction verified for
//! the global index rather than introducing any new format.

#[path = "../graph.rs"]
mod graph;
#[allow(dead_code)]
#[path = "../ext.rs"]
mod ext;
#[path = "../gfmbuild.rs"]
mod gfmbuild;

use std::convert::TryInto;
use std::{env, fs};

const LOCAL_SIZE: u32 = (1 << 16) - (1 << 13);      // 57,344
const LOCAL_OVERLAP: u32 = 1024;
const LOCAL_INTERVAL: u32 = LOCAL_SIZE - LOCAL_OVERLAP;
const LOCAL_LINE_RATE: u32 = 7;
const LOCAL_OFF_RATE: u32 = 3;
const LOCAL_FTAB_CHARS: u32 = 6;
const LOCAL_MAX_GBWT: u32 = (1 << 16) - (1 << 11);   // 63,488

/// `HGFM::selectAlts` -- keep `k` variants evenly spaced through the original
/// list and rebuild one single-variant haplotype per kept variant.
fn select_alts(orig: &[graph::Alt], k: usize) -> (Vec<graph::Alt>, Vec<graph::Hap>) {
    let n = orig.len();
    let k = k.min(n);
    let mut alts = Vec::with_capacity(k);
    let mut haps = Vec::with_capacity(k);
    if k == 0 { return (alts, haps); }
    for i in 0..k {
        let o = &orig[(i as u64 * n as u64 / k as u64) as usize];
        alts.push(graph::Alt { pos: o.pos, len: o.len, seq: o.seq, typ: o.typ });
    }
    for (a, alt) in alts.iter().enumerate() {
        let right = if alt.typ == graph::ALT_DEL { alt.pos + alt.len - 1 } else { alt.pos };
        haps.push(graph::Hap { left: alt.pos, right, alts: vec![a as u32] });
    }
    (alts, haps)
}

fn put16(v: &mut Vec<u8>, x: u32) { v.extend_from_slice(&(x as u16).to_le_bytes()); }
fn put32(v: &mut Vec<u8>, x: u32) { v.extend_from_slice(&x.to_le_bytes()); }

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 5 {
        eprintln!("usage: ht2emit5 <reference.fa> <snp> <haplotype> <index_prefix>");
        std::process::exit(2);
    }
    let p = graph::parse(&a[1], &a[2], &a[3]);
    let len = p.text.len() as u32;
    // Windows are laid out in CHROMOSOME coordinates, ambiguity included --
    // hgfm.h accrues rec.off + rec.len per RefRecord and asserts the count at
    // :2172. A 1,000,000 bp sequence with a 100,000 bp N run therefore gets 18
    // windows, not the 16 its 900,000 joined bases would suggest.
    let mut plan: Vec<(u32, u32, u32, u32, u32)> = Vec::new(); // tidx, cstart, joined_start, wlen(joined), chrom_len
    for (ti, (_name, full, runs)) in p.seqs.iter().enumerate() {
        let nw = ((*full + LOCAL_INTERVAL - 1) / LOCAL_INTERVAL).max(1);
        for w in 0..nw {
            let cs = w * LOCAL_INTERVAL;
            let ce = (cs + LOCAL_SIZE).min(*full);
            // joined bases covered by chrom range [cs, ce)
            let mut jstart = u32::MAX;
            let mut cnt = 0u32;
            for &(co, jo, rl) in runs.iter() {
                let lo = co.max(cs);
                let hi = (co + rl).min(ce);
                if lo < hi {
                    if jstart == u32::MAX { jstart = jo + (lo - co); }
                    cnt += hi - lo;
                }
            }
            if jstart == u32::MAX { jstart = 0; }
            plan.push((ti as u32, cs, jstart, cnt, *full));
        }
    }
    let n_windows = plan.len() as u32;
    println!("{len} joined bp over {} sequence(s) -> {n_windows} local index(es)", p.seqs.len());

    let ftab_len = (1usize << (2 * LOCAL_FTAB_CHARS)) + 1;
    let side_sz = 1usize << LOCAL_LINE_RATE;

    let mut o5: Vec<u8> = Vec::new();
    let mut o6: Vec<u8> = Vec::new();
    put32(&mut o5, 1);
    put32(&mut o5, n_windows);
    put32(&mut o5, LOCAL_LINE_RATE);
    put32(&mut o5, 2);
    put32(&mut o5, LOCAL_OFF_RATE);
    put32(&mut o5, LOCAL_FTAB_CHARS);
    put32(&mut o5, (-1i32) as u32);
    put32(&mut o6, 1);

    for (w, &(tidx, cstart, a0, wlen, _full)) in plan.iter().enumerate() {
        let w = w as u32;
        let b0 = a0 + wlen;
        let wtext: Vec<u8> = p.text[a0 as usize..b0 as usize].to_vec();
        if wlen == 0 {
            // a window entirely inside an N run: header only, which is what the
            // `len == 0` early return in LocalGFM::readIntoMemory is for
            put32(&mut o5, tidx); put32(&mut o5, cstart); put32(&mut o5, a0);
            put16(&mut o5, 0); put16(&mut o5, 0); put16(&mut o5, 0); put16(&mut o5, 0);
            println!("  window {w}: chrom {cstart}, empty (inside an N run)");
            continue;
        }
        // graph over just this window; the same construction as the global one
        // A local graph that grows past local_max_gbwt nodes is thrown away and
        // rebuilt from a thinned variant set; hgfm.h:1954 binary-searches for the
        // largest set that fits. tiny.log shows zero explosions, which is why a
        // single window matched without any of this; clean.log shows 13.
        // Haplotypes index into the FULL alt list, so restricting the alts to a
        // window means remapping every haplotype's indices -- `alt_map` in the
        // C++ (hgfm.h:2337). Without it a haplotype in window 1 still points at
        // a global alt index and reads off the end.
        let mut alt_map: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
        let mut orig: Vec<graph::Alt> = Vec::new();
        for (gi, x) in p.alts.iter().enumerate() {
            let keep = match x.typ {
                graph::ALT_DEL => x.pos >= a0 && x.pos + x.len <= b0,
                _              => x.pos >= a0 && x.pos < b0,
            };
            if keep {
                alt_map.insert(gi as u32, orig.len() as u32);
                // LOCAL coordinates: the C++ does `alts.back().pos -= curr_sztot`
                // (hgfm.h:2355), so a local index is built over a window that
                // starts at 0. Window 0 cannot tell the difference.
                orig.push(graph::Alt { pos: x.pos - a0, len: x.len, seq: x.seq, typ: x.typ });
            }
        }
        let mut alts: Vec<graph::Alt> = orig.iter()
            .map(|x| graph::Alt { pos: x.pos, len: x.len, seq: x.seq, typ: x.typ }).collect();
        let mut haps: Vec<graph::Hap> = p.haps.iter()
            .filter(|h| h.left >= a0 && h.right < b0)
            .filter_map(|h| {
                let m: Option<Vec<u32>> = h.alts.iter().map(|a| alt_map.get(a).copied()).collect();
                m.map(|alts| graph::Hap { left: h.left - a0, right: h.right - a0, alts })
            })
            .collect();
        let (mut bs_active, mut bs_lo, mut bs_hi) = (false, 0usize, 0usize);
        let mut bs_cur = alts.len();
        let mut dropped = 0usize;
        let g;
        loop {
            let (n, e, _) = graph::build_range_with(&wtext, &alts, &haps, 0, wlen, true, true);
            let (gn, ge, last) = graph::reverse_determinize(&n, &e, wlen + 1);
            let built = gfmbuild::build_gfm(&gn, &ge, last, wlen, LOCAL_OFF_RATE, LOCAL_MAX_GBWT);
            let exploded = match &built {
                None => true,
                Some(b) => b.gbwt_len > LOCAL_MAX_GBWT,
            };
            if !exploded && bs_active && bs_hi > bs_lo + 1 {
                bs_lo = bs_cur;
                bs_cur = bs_lo + (bs_hi - bs_lo) / 2;
                let (na, nh) = select_alts(&orig, bs_cur);
                alts = na; haps = nh;
                continue;
            }
            if exploded {
                if !bs_active { bs_active = true; bs_lo = 0; }
                bs_hi = bs_cur;
                bs_cur = if bs_hi <= bs_lo + 1 { bs_lo } else { bs_lo + (bs_hi - bs_lo) / 2 };
                let (na, nh) = select_alts(&orig, bs_cur);
                alts = na; haps = nh;
                continue;
            }
            dropped = orig.len() - alts.len();
            g = built.unwrap();
            break;
        }
        let _ = dropped;

        let linear = g.len + 1 == g.gbwt_len;
        let gbwt_sz = if linear { g.gbwt_len / 4 + 1 } else { g.gbwt_len / 2 + 1 } as usize;
        let side_gbwt_sz = side_sz - 2 * if linear { 4 } else { 6 };
        let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
        let rows_per_side = if linear { side_gbwt_sz * 4 } else { side_gbwt_sz * 2 };
        let gbwt_tot = num_sides * side_sz;

        let sect_at = o5.len();
        put32(&mut o5, tidx);
        put32(&mut o5, cstart);     // localOffset -- chromosome coordinates
        put32(&mut o5, a0);         // joinedOffset
        put16(&mut o5, g.len);
        put16(&mut o5, g.gbwt_len);
        put16(&mut o5, g.num_nodes);
        // eftabLen is discovered below, so remember where to patch it
        let eftab_len_at = o5.len();
        put16(&mut o5, 0);
        // plen is the window's CHROMOSOME span, ambiguity included -- the same
        // relationship the global index has, where plen[0] is 1,000,000 while
        // len is 900,000. Only a reference with N runs can tell them apart.
        let cspan = (cstart + LOCAL_SIZE).min(_full) - cstart;
        put16(&mut o5, 1); put16(&mut o5, cspan);           // nPat, plen
        // rstarts, local: one record per unambiguous stretch inside the window,
        // as (offset in the window's joined text, sequence, offset in the
        // window's pattern). A window that opens inside an N run has a non-zero
        // third field -- 46,231 for the one that straddles this reference's gap.
        let ce = cstart + cspan;
        let mut frags: Vec<(u32, u32, u32)> = Vec::new();
        let mut jacc = 0u32;
        for &(co, _jo, rl) in p.seqs[tidx as usize].2.iter() {
            let lo = co.max(cstart);
            let hi = (co + rl).min(ce);
            if lo < hi { frags.push((jacc, 0, lo - cstart)); jacc += hi - lo; }
        }
        if frags.is_empty() { frags.push((0, 0, 0)); }
        put16(&mut o5, frags.len() as u32);
        for &(j, t, c) in &frags { put16(&mut o5, j); put16(&mut o5, t); put16(&mut o5, c); }

        // ---- gbwt block ----
        let mut blk = vec![0u8; gbwt_tot];
        let (mut occ, mut m_occ, mut f_loc) = ([0u32; 4], 0u32, 0u32);
        let mut indeg = vec![0u32; g.num_nodes as usize + 1];
        {
            let mut node = 0usize;
            for r in 0..g.rows.len() {
                if g.rows[r].1 == 1 && r > 0 { node += 1; }
                if node < indeg.len() { indeg[node] += 1; }
            }
        }
        let mut floc_of_node: Vec<u32> = Vec::with_capacity(indeg.len());
        { let mut acc = 0u32; for i in 0..indeg.len() { floc_of_node.push(acc); acc += indeg[i]; } }
        let mut floc_i = 0usize;
        for si in 0..num_sides * rows_per_side {
            let side = si / rows_per_side;
            let off = si % rows_per_side;
            if off == 0 {
                let base = side * side_sz + side_sz - 12;
                for (k, v) in [f_loc, m_occ, occ[0], occ[1], occ[2], occ[3]].iter().enumerate() {
                    blk[base + k * 2..base + k * 2 + 2].copy_from_slice(&(*v as u16).to_le_bytes());
                }
            }
            let (ch, f, m) = if si < g.rows.len() { (g.rows[si].0, g.rows[si].1, g.mrun[si].0) }
                             else { (b'A', 0, 0) };
            let mut count = true;
            let cc = match ch { b'A' => 0u8, b'C' => 1, b'G' => 2, b'T' => 3,
                                _ => { count = false; 0 } };
            if m == 1 { if floc_i < floc_of_node.len() { f_loc = floc_of_node[floc_i]; floc_i += 1; } }
            if count { occ[cc as usize] += 1; }
            if m == 1 { m_occ += 1; }
            let base = side * side_sz;
            let sc = off >> 2; let bpi = off & 3;
            blk[base + sc] |= cc << (bpi * 2);
            let f_sc = (side_gbwt_sz + sc) >> 1;
            let f_bpi = bpi + ((sc & 1) << 2);
            blk[base + f_sc] |= f << f_bpi;
            blk[base + f_sc + (side_gbwt_sz >> 2)] |= m << f_bpi;
        }
        o5.extend_from_slice(&blk);

        put16(&mut o5, g.z_offs.len() as u32);
        for &z in &g.z_offs { put16(&mut o5, z); }
        for i in 0..5 { put16(&mut o5, g.fchr[i]); }

        // ---- ftab / eftab, by querying the block just written ----
        let gl = g.gbwt_len as usize;
        let (mut bwt, mut fb, mut mb) = (vec![0u8; gl], vec![0u8; gl], vec![0u8; gl]);
        for r in 0..gl {
            let side = r / rows_per_side; let off = r % rows_per_side;
            let base = side * side_sz; let sc = off >> 2; let bpi = off & 3;
            bwt[r] = (blk[base + sc] >> (bpi * 2)) & 3;
            let f_sc = (side_gbwt_sz + sc) >> 1;
            let f_bpi = bpi + ((sc & 1) << 2);
            fb[r] = (blk[base + f_sc] >> f_bpi) & 1;
            mb[r] = (blk[base + f_sc + (side_gbwt_sz >> 2)] >> f_bpi) & 1;
        }
        let zset: std::collections::HashSet<usize> = g.z_offs.iter().map(|&z| z as usize).collect();
        let mut occp = vec![[0u32; 4]; gl + 1];
        let mut rankm = vec![0u32; gl + 1];
        let mut fsel: Vec<u32> = Vec::new();
        for r in 0..gl {
            occp[r + 1] = occp[r];
            if !zset.contains(&r) { occp[r + 1][bwt[r] as usize] += 1; }
            rankm[r + 1] = rankm[r] + mb[r] as u32;
            if fb[r] == 1 { fsel.push(r as u32); }
        }
        let mut tftab: Vec<(u32, u32)> = Vec::with_capacity(ftab_len - 1);
        for i in 0..ftab_len - 1 {
            let mut q = i;
            let (mut top, mut bot) = (0u32, gl as u32);
            let mut j = 0usize;
            while j < LOCAL_FTAB_CHARS as usize {
                let c = q & 3; q >>= 2;
                let nt = g.fchr[c] + occp[top as usize][c];
                let nb = g.fchr[c] + occp[bot as usize][c];
                if nt >= nb { top = nt; bot = nb; break; }
                let node_top = rankm[(nt as usize + 1).min(gl)] - 1;
                let node_bot = rankm[(nb as usize).min(gl)];
                top = *fsel.get(node_top as usize).unwrap_or(&(gl as u32));
                bot = *fsel.get(node_bot as usize).unwrap_or(&(gl as u32));
                if top >= bot { break; }
                j += 1;
            }
            if top >= bot || j < LOCAL_FTAB_CHARS as usize {
                let v = if i == 0 { 0 } else { tftab[i - 1].1 };
                tftab.push((v, v));
            } else { tftab.push((top, bot)); }
        }
        let mut ftab_o = vec![0u32; ftab_len];
        let mut eftab_o: Vec<u32> = Vec::new();
        ftab_o[0] = tftab[0].0; ftab_o[1] = tftab[0].1;
        for i in 1..ftab_len - 1 {
            if ftab_o[i] != tftab[i].0 {
                let (lo, hi) = (ftab_o[i], tftab[i].0);
                ftab_o[i] = (eftab_o.len() as u32 / 2) ^ 0xFFFF;   // u16 marker
                eftab_o.push(lo); eftab_o.push(hi);
            }
            ftab_o[i + 1] = tftab[i].1;
        }
        for i in 0..ftab_len { put16(&mut o5, ftab_o[i]); }
        for &v in &eftab_o { put16(&mut o5, v); }
        o5[eftab_len_at..eftab_len_at + 2].copy_from_slice(&(eftab_o.len() as u16).to_le_bytes());

        for &s in &g.sa_sample { put16(&mut o6, s); }
        println!("  window {w}: @{sect_at} chrom {cstart} joined [{a0},{b0}) len {} gbwtLen {} numNodes {} eftabLen {} sides {} offs {}",
                 g.len, g.gbwt_len, g.num_nodes, eftab_o.len(), num_sides, g.sa_sample.len());
    }
    o5.push(0);

    let t5 = fs::read(format!("{}.5.ht2", a[4])).expect(".5.ht2");
    let t6 = fs::read(format!("{}.6.ht2", a[4])).expect(".6.ht2");
    for (name, ours, theirs) in [(".5", &o5, &t5), (".6", &o6, &t6)] {
        let n = ours.len().min(theirs.len());
        let d = (0..n).filter(|&i| ours[i] != theirs[i]).count();
        let first = (0..n).find(|&i| ours[i] != theirs[i]);
        if d == 0 && ours.len() == theirs.len() {
            println!("{name}.ht2: BYTE-IDENTICAL ({} bytes)", ours.len());
        } else {
            println!("{name}.ht2: {} of {n} bytes match (ours {} bytes, theirs {}){}",
                     n - d, ours.len(), theirs.len(),
                     match first { Some(f) => format!("; first differing byte at {f}"), None => String::new() });
        }
    }
    let _ = TryInto::<u32>::try_into(0u32);
}
