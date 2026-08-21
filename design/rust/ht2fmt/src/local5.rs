//! Emit `.5.ht2` / `.6.ht2` — the hierarchical local indexes.
//!
//! Each local index is a GFM over a 57,344 bp window at `u16` width with
//! `offRate` 3 and `ftabChars` 6, so this reuses the construction verified for
//! the global index rather than introducing any new format.

use super::gfmbuild;
use super::graph;
use std::convert::TryInto;
use std::io::Write;

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

/// Emit both files for `p`, one window's section at a time.
///
/// Nothing bigger than a single local index is ever buffered: `eftabLen` is only
/// known once that window's ftab is built, so the section is held long enough to
/// patch it and then written. A whole-genome `.5.ht2` is several GB and cannot
/// be assembled in memory.
pub fn emit<W5: Write, W6: Write>(p: &graph::Parsed, w5: &mut W5, w6: &mut W6,
                                  width: usize, verbose: bool)
    -> std::io::Result<()>
{
    // A local index is `local_index_t` = u16 throughout whatever the global
    // index width is, so only four fields change between `.5.ht2` and
    // `.5.ht2l`: the local-index count, and each window's tidx, localOffset and
    // joinedOffset.
    let put_idx = |v: &mut Vec<u8>, x: u32| v.extend_from_slice(&(x as u64).to_le_bytes()[0..width]);
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
    if verbose {
        println!("{len} joined bp over {} sequence(s) -> {n_windows} local index(es)", p.seqs.len());
    }

    let mut hdr: Vec<u8> = Vec::new();
    put32(&mut hdr, 1);
    put_idx(&mut hdr, n_windows);
    put32(&mut hdr, LOCAL_LINE_RATE);
    put32(&mut hdr, 2);
    put32(&mut hdr, LOCAL_OFF_RATE);
    put32(&mut hdr, LOCAL_FTAB_CHARS);
    put32(&mut hdr, (-1i32) as u32);
    w5.write_all(&hdr)?;
    w6.write_all(&1u32.to_le_bytes())?;
    let mut off5 = hdr.len();


/// One local index: its `.5` section and its `.6` samples.
///
/// Depends on nothing but the parsed reference and this window's own plan entry
/// -- no shared state, no ordering constraint -- which is what lets the whole
/// stage run one window per core. `hgfm.h` threads exactly this loop.
struct WinOut { sec: Vec<u8>, sa: Vec<u8>, note: String }

fn build_window(p: &graph::Parsed, (tidx, cstart, a0, wlen, _full): (u32, u32, u32, u32, u32),
                width: usize) -> WinOut
{
    let put_idx = |v: &mut Vec<u8>, x: u32| v.extend_from_slice(&(x as u64).to_le_bytes()[0..width]);
    let ftab_len = (1usize << (2 * LOCAL_FTAB_CHARS)) + 1;
    let side_sz = 1usize << LOCAL_LINE_RATE;
    let mut sec: Vec<u8> = Vec::new();
    let mut sa: Vec<u8> = Vec::new();
        let b0 = a0 + wlen;
        let wtext: Vec<u8> = p.text[a0 as usize..b0 as usize].to_vec();
        if wlen == 0 {
            // a window entirely inside an N run: header only, which is what the
            // `len == 0` early return in LocalGFM::readIntoMemory is for
            put_idx(&mut sec, tidx); put_idx(&mut sec, cstart); put_idx(&mut sec, a0);
            put16(&mut sec, 0); put16(&mut sec, 0); put16(&mut sec, 0); put16(&mut sec, 0);
            return WinOut { sec, sa, note: format!("chrom {cstart}, empty (inside an N run)") };
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

        put_idx(&mut sec, tidx);
        put_idx(&mut sec, cstart);     // localOffset -- chromosome coordinates
        put_idx(&mut sec, a0);         // joinedOffset
        put16(&mut sec, g.len);
        put16(&mut sec, g.gbwt_len);
        put16(&mut sec, g.num_nodes);
        // eftabLen is discovered below, so remember where to patch it
        let eftab_len_at = sec.len();
        put16(&mut sec, 0);
        // plen is the window's CHROMOSOME span, ambiguity included -- the same
        // relationship the global index has, where plen[0] is 1,000,000 while
        // len is 900,000. Only a reference with N runs can tell them apart.
        let cspan = (cstart + LOCAL_SIZE).min(_full) - cstart;
        put16(&mut sec, 1); put16(&mut sec, cspan);           // nPat, plen
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
        put16(&mut sec, frags.len() as u32);
        for &(j, t, c) in &frags { put16(&mut sec, j); put16(&mut sec, t); put16(&mut sec, c); }

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
                // A linear side reserves 4 `local_index_t` for the occ tallies;
                // a graph side reserves 6, the extra two being F_locSave and
                // M_occSave. Both hold the counts as of the side's start.
                let saves: &[u32] = if linear { &[occ[0], occ[1], occ[2], occ[3]] }
                                    else { &[f_loc, m_occ, occ[0], occ[1], occ[2], occ[3]] };
                let base = side * side_sz + side_sz - 2 * saves.len();
                for (k, v) in saves.iter().enumerate() {
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
            // A linear side is BWT for its whole length -- there are no F and M
            // bitvectors, because every row is its own node. Writing them here
            // ran off the end of the block: the linear side packs 4 rows/byte
            // across all 120 bytes, so the graph formula's F offset lands past
            // the side entirely.
            if !linear {
                let f_sc = (side_gbwt_sz + sc) >> 1;
                let f_bpi = bpi + ((sc & 1) << 2);
                blk[base + f_sc] |= f << f_bpi;
                blk[base + f_sc + (side_gbwt_sz >> 2)] |= m << f_bpi;
            }
        }
        sec.extend_from_slice(&blk);

        put16(&mut sec, g.z_offs.len() as u32);
        for &z in &g.z_offs { put16(&mut sec, z); }
        for i in 0..5 { put16(&mut sec, g.fchr[i]); }

        // ---- ftab / eftab, by querying the block just written ----
        let gl = g.gbwt_len as usize;
        let (mut bwt, mut fb, mut mb) = (vec![0u8; gl], vec![0u8; gl], vec![0u8; gl]);
        for r in 0..gl {
            let side = r / rows_per_side; let off = r % rows_per_side;
            let base = side * side_sz; let sc = off >> 2; let bpi = off & 3;
            bwt[r] = (blk[base + sc] >> (bpi * 2)) & 3;
            if linear {
                // Every row is its own node, so rank_M is the identity and
                // select_F is too -- which collapses the graph ftab walk below
                // into the ordinary FM-index one, no separate code path needed.
                fb[r] = 1; mb[r] = 1;
            } else {
                let f_sc = (side_gbwt_sz + sc) >> 1;
                let f_bpi = bpi + ((sc & 1) << 2);
                fb[r] = (blk[base + f_sc] >> f_bpi) & 1;
                mb[r] = (blk[base + f_sc + (side_gbwt_sz >> 2)] >> f_bpi) & 1;
            }
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
        // The graph writer stores only the eftab entries it used; the linear
        // one reserves `ftabChars*2` and zero-fills the rest, which is why
        // every linear local index reports eftabLen 12 exactly.
        // The linear writer absorbs the rows no `ftabChars`-long prefix can
        // reach -- suffixes shorter than 6 characters, and the $ row -- into
        // the final bucket, so its last entry is an eftab pair ending at
        // gbwtLen rather than a bare boundary. The graph writer does not: its
        // last entry stays a plain `gbwtLen - 1` (checked against every graph
        // window in these fixtures), which is why this is linear-only.
        if linear && ftab_o[ftab_len - 1] != g.gbwt_len {
            let (lo, hi) = (ftab_o[ftab_len - 1], g.gbwt_len);
            ftab_o[ftab_len - 1] = (eftab_o.len() as u32 / 2) ^ 0xFFFF;
            eftab_o.push(lo); eftab_o.push(hi);
        }
        // The graph writer stores only the eftab entries it used; the linear
        // one reserves `ftabChars*2` and zero-fills the rest, which is why
        // every linear local index reports eftabLen 12 exactly.
        if linear { eftab_o.resize(2 * LOCAL_FTAB_CHARS as usize, 0); }
        for i in 0..ftab_len { put16(&mut sec, ftab_o[i]); }
        for &v in &eftab_o { put16(&mut sec, v); }
        sec[eftab_len_at..eftab_len_at + 2].copy_from_slice(&(eftab_o.len() as u16).to_le_bytes());

        for &s in &g.sa_sample { sa.extend_from_slice(&(s as u16).to_le_bytes()); }
        let note = format!("chrom {cstart} joined [{a0},{b0}) len {} gbwtLen {} numNodes {} eftabLen {} sides {} offs {}",
                           g.len, g.gbwt_len, g.num_nodes, eftab_o.len(), num_sides, g.sa_sample.len());
        WinOut { sec, sa, note }
}

    // One window per core. Results are handed to the writer in window order,
    // and dispatch is held within `cap` windows of whatever is being written
    // next, so at most that many sections are ever resident -- a few MB -- and
    // the window the writer is waiting for is always already in flight.
    let threads = super::ext::threads().min(plan.len().max(1));
    if threads <= 1 {
        for (w, &e) in plan.iter().enumerate() {
            let o = build_window(p, e, width);
            if verbose { println!("  window {w}: @{off5} {}", o.note); }
            w5.write_all(&o.sec)?; off5 += o.sec.len();
            w6.write_all(&o.sa)?;
        }
    } else {
        use std::collections::BTreeMap;
        use std::sync::{Condvar, Mutex};
        struct State { dispatch: usize, write: usize, done: BTreeMap<usize, WinOut> }
        let cap = threads * 2;
        let n = plan.len();
        let st = Mutex::new(State { dispatch: 0, write: 0, done: BTreeMap::new() });
        let cv = Condvar::new();
        let (st, cv, plan_ref) = (&st, &cv, &plan);
        let mut err: Option<std::io::Error> = None;
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(move || loop {
                    let i = {
                        let mut g = st.lock().unwrap();
                        while g.dispatch >= g.write + cap && g.dispatch < n {
                            g = cv.wait(g).unwrap();
                        }
                        if g.dispatch >= n { return; }
                        let i = g.dispatch; g.dispatch += 1; i
                    };
                    let o = build_window(p, plan_ref[i], width);
                    let mut g = st.lock().unwrap();
                    g.done.insert(i, o);
                    cv.notify_all();
                });
            }
            // the writer, on this thread, so `w5`/`w6` need not be Send
            for w in 0..n {
                let o = {
                    let mut g = st.lock().unwrap();
                    loop {
                        if let Some(o) = g.done.remove(&w) { break o; }
                        g = cv.wait(g).unwrap();
                    }
                };
                {
                    let mut g = st.lock().unwrap();
                    g.write = w + 1;
                    cv.notify_all();
                }
                if verbose { println!("  window {w}: @{off5} {}", o.note); }
                if let Err(e) = w5.write_all(&o.sec) { err = Some(e); break; }
                off5 += o.sec.len();
                if let Err(e) = w6.write_all(&o.sa) { err = Some(e); break; }
            }
            // let any parked worker out if the writer bailed
            let mut g = st.lock().unwrap();
            g.write = n; g.dispatch = n;
            cv.notify_all();
        });
        if let Some(e) = err { return Err(e); }
    }
    // the file ends with a single '\0' (`fout5 << '\0'`)
    w5.write_all(&[0])?;
    Ok(())
}
