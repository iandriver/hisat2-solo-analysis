//! Rung 2, final step — emit a complete `.1.ht2` + `.2.ht2` from a FASTA and
//! byte-compare against `hisat2-build`'s own output.
//!
//! Every earlier binary checked one piece in isolation. This one writes the
//! whole file, which is the only test that cannot be passed by accident: a
//! misunderstanding anywhere — a section order, a padding rule, an off-by-one
//! in the occ tallies, a stray byte in `refnames` — shows up as a differing
//! offset rather than as a number that happens to agree.
//!
//! Writer being reproduced: `GFM::writeFromMemory(justHeader=true)` (gfm.h:6604)
//! for the header, `joinToDisk` for the fragment records, `GFM::buildToDisk`
//! (gfm.h:5148) for everything from the gbwt block onward, and the refnames
//! trailer at gfm.h:2383. `.3.ht2`/`.4.ht2` come from a different writer —
//! `BitPairReference::szsToDisk` (`reference.cpp:653`) — and are covered here
//! too, since they are what the aligner reads the actual sequence from.
//!
//! `.5`-`.8` are deliberately out of scope at this rung: `.5`/`.6` hold the
//! hierarchical local indexes (`hgfm.h:2068`), which only exist once there is a
//! graph to localise, and belong with the SNP/haplotype rung.
//!
//! Two header fields are back-patched by `buildToDisk` rather than written in
//! order, and both are reproduced here as final values:
//!   * `gbwtLen`/`numNodes` at offset 12/16 (gfm.h:5166), written as 0 up front;
//!   * `eftabLen` at offset 36 (gfm.h:5471), likewise, and always `ftabChars*2`
//!     regardless of how many entries actually absorb.
//!
//! Defaults come from `hisat2_build.cpp`: `offRate` 4, `ftabChars` 10, and
//! `lineRate` = `default_lineRate_fm` = 6 for a 32-bit index with no variants
//! (gfm.h:4308). They are asserted against the reference header rather than
//! assumed silently.

#[path = "../refin.rs"]
mod refin;

use refin::read_fasta;
use std::convert::TryInto;
use std::{env, fs};

const HI: u32 = u32::MAX;

/// HISAT2's suffix array: the text is treated as padded past its end with a
/// character LARGER than any real one, so the empty suffix sorts last and any
/// proper prefix sorts after its extension (see `sa_rule_check.rs`).
fn hisat2_suffix_array(t: &[u8]) -> Vec<u32> {
    let n = t.len();
    let m = n + 1;
    let mut sa: Vec<u32> = (0..m as u32).collect();
    let mut rank: Vec<u32> = (0..m).map(|i| if i == n { HI } else { t[i] as u32 }).collect();
    let mut tmp = vec![0u32; m];
    let mut k = 1usize;
    loop {
        let key = |rk: &Vec<u32>, i: u32, k: usize| -> (u32, u32) {
            let a = rk[i as usize];
            let j = i as usize + k;
            (a, if j < m { rk[j] } else { HI })
        };
        sa.sort_unstable_by_key(|&i| key(&rank, i, k));
        tmp[sa[0] as usize] = 0;
        let mut cl = 0u32;
        for w in 1..m {
            if key(&rank, sa[w], k) != key(&rank, sa[w - 1], k) { cl += 1; }
            tmp[sa[w] as usize] = cl;
        }
        std::mem::swap(&mut rank, &mut tmp);
        if cl as usize == m - 1 { break; }
        k <<= 1;
        if k >= m { break; }
    }
    sa
}

fn put(v: &mut Vec<u8>, x: u32) { v.extend_from_slice(&x.to_le_bytes()); }

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: ht2emit <reference.fa> <out_prefix> [<reference.1.ht2>]");
        eprintln!("  writes <out_prefix>.1.ht2 and .2.ht2; byte-compares if a reference is given");
        std::process::exit(2);
    }
    let (fa, out_prefix) = (&a[1], &a[2]);

    let r = read_fasta(fa);
    let t = &r.text;
    let len = t.len() as u32;

    // ---- parameters ----------------------------------------------------
    let line_rate: u32 = 6;      // default_lineRate_fm, 32-bit, no variants
    let off_rate: u32 = 4;
    let ftab_chars: u32 = 10;
    let version: u32 = (2 << 24) | (2 << 16) | (3 << 8); // 2.2.3, no extra tag

    let gbwt_len = len + 1;      // linear FM
    let num_nodes = len + 1;
    let ftab_len = (1usize << (ftab_chars * 2)) + 1;
    let eftab_len = (ftab_chars * 2) as usize;
    let side_sz = 1u32 << line_rate;
    let side_gbwt_sz = side_sz - 4 * 4;
    let side_gbwt_len = side_gbwt_sz << 2;
    let gbwt_sz = gbwt_len / 4 + 1;
    let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
    let offs_len = ((num_nodes as u64 + (1u64 << off_rate) - 1) >> off_rate) as usize;

    // ---- header --------------------------------------------------------
    let mut o: Vec<u8> = Vec::new();
    put(&mut o, 1);                 // endianness sentinel
    put(&mut o, version);
    put(&mut o, len);
    put(&mut o, gbwt_len);          // back-patched by buildToDisk
    put(&mut o, num_nodes);         // back-patched by buildToDisk
    put(&mut o, line_rate);
    put(&mut o, 2);                 // linesPerSide, written as a literal 2
    put(&mut o, off_rate);
    put(&mut o, ftab_chars);
    put(&mut o, eftab_len as u32);  // back-patched by buildToDisk
    put(&mut o, (-1i32) as u32);    // -flags, flags = 1 (entireReverse off)

    // ---- front end -----------------------------------------------------
    put(&mut o, r.names.len() as u32);
    for &p in &r.plen { put(&mut o, p); }
    put(&mut o, r.frags.len() as u32);
    for f in &r.frags { put(&mut o, f.joined_off); put(&mut o, f.text_id); put(&mut o, f.text_off); }

    // ---- the BWT walk (gfm.h:5240 onward) -------------------------------
    let sa = hisat2_suffix_array(t);
    let mut fchr = [0u32; 5];
    let mut ftab = vec![0u32; ftab_len];
    let mut absorb = vec![0u32; ftab_len];
    let mut absorb_cnt = 0u32;
    let mut z_offs: Vec<u32> = Vec::new();
    let mut sa_sample: Vec<u32> = Vec::new();

    let mut occ = [0u32; 4];
    let mut occ_save = [0u32; 4];
    let mut side_buf = vec![0u8; side_sz as usize];
    let mut side_cur = 0usize;
    let mut si: u32 = 0;
    let off_mask: u32 = u32::MAX << off_rate;
    let gbwt_start = o.len();

    while (o.len() - gbwt_start) < (num_sides * side_sz) as usize {
        side_buf[side_cur] = 0;
        for bpi in 0..4 {
            let bwt_char;
            let mut count = true;
            if si <= len {
                let sa_elt = sa[si as usize];
                if sa_elt == 0 {
                    // '$' is not representable and must not be counted, or LF breaks
                    bwt_char = 0u8; count = false;
                    z_offs.push(si);
                } else {
                    bwt_char = t[(sa_elt - 1) as usize];
                    fchr[bwt_char as usize] += 1;
                }
                if (len - sa_elt) >= ftab_chars {
                    let mut suf_int = 0u32;
                    for i in 0..ftab_chars as usize { suf_int = (suf_int << 2) | t[sa_elt as usize + i] as u32; }
                    ftab[(suf_int + 1) as usize] += 1;
                    if absorb_cnt > 0 { absorb[suf_int as usize] = absorb_cnt; absorb_cnt = 0; }
                } else {
                    absorb_cnt += 1;
                }
                if (si & off_mask) == si { sa_sample.push(sa_elt); }
            } else {
                // Past the end of the SA: pad with 'A', and the padding IS counted.
                bwt_char = 0;
            }
            if count { occ[bwt_char as usize] += 1; }
            side_buf[side_cur] |= bwt_char << (bpi * 2);   // pack_2b_in_8b, low pair first
            si += 1;
        }
        side_cur += 1;
        if side_cur == side_gbwt_sz as usize {
            side_cur = 0;
            // The four tallies are occSave: the counts as of the START of this side.
            for c in 0..4 {
                let at = side_sz as usize - 16 + c * 4;
                side_buf[at..at + 4].copy_from_slice(&occ_save[c].to_le_bytes());
            }
            occ_save = occ;
            o.extend_from_slice(&side_buf);
            side_buf.iter_mut().for_each(|b| *b = 0);
        }
    }
    if absorb_cnt > 0 { absorb[ftab_len - 1] = absorb_cnt; }

    // ---- zOffs, fchr ----------------------------------------------------
    put(&mut o, z_offs.len() as u32);
    for &z in &z_offs { put(&mut o, z); }
    // exclusive prefix sum, then shift up by one
    for i in 1..4 { fchr[i] += fchr[i - 1]; }
    for i in (1..5).rev() { fchr[i] = fchr[i - 1]; }
    fchr[0] = 0;
    for i in 0..5 { put(&mut o, fchr[i]); }

    // ---- ftab / eftab ---------------------------------------------------
    let mut eftab = vec![0u32; eftab_len];
    let mut ftab_out = vec![0u32; ftab_len];
    let mut ecur = 0usize;
    let hi_of = |i: usize, ftab_out: &Vec<u32>, eftab: &Vec<u32>| -> u32 {
        let v = ftab_out[i];
        if v as u64 <= gbwt_len as u64 { v } else { eftab[((v ^ HI) as usize) * 2 + 1] }
    };
    for i in 1..ftab_len {
        let lo = ftab[i] + hi_of(i - 1, &ftab_out, &eftab);
        if absorb[i] > 0 {
            eftab[ecur * 2] = lo;
            eftab[ecur * 2 + 1] = lo + absorb[i];
            ftab_out[i] = (ecur as u32) ^ HI;
            ecur += 1;
        } else {
            ftab_out[i] = lo;
        }
    }
    for i in 0..ftab_len { put(&mut o, ftab_out[i]); }
    for i in 0..eftab_len { put(&mut o, eftab[i]); }

    // ---- refnames (gfm.h:2383): each name + '\n', then a single '\0' -----
    for n in &r.names { o.extend_from_slice(n.as_bytes()); o.push(b'\n'); }
    o.push(0);

    // ---- .2.ht2 ---------------------------------------------------------
    let mut o2: Vec<u8> = Vec::new();
    put(&mut o2, 1); // endian hint for the secondary stream
    for &v in &sa_sample { put(&mut o2, v); }

    // ---- .3.ht2 / .4.ht2 (reference.cpp:665) ----------------------------
    // `.3` is a sentinel, a record count, and one 9-byte RefRecord per maximal
    // unambiguous stretch: the number of ambiguous characters immediately
    // before it, its length, and whether it opens a new sequence.
    let mut o3: Vec<u8> = Vec::new();
    put(&mut o3, 1);
    put(&mut o3, r.recs.len() as u32);
    for rec in &r.recs {
        put(&mut o3, rec.off);
        put(&mut o3, rec.len);
        o3.push(if rec.first { 1 } else { 0 });
    }
    // `.4` is the unambiguous text, 2 bits per base, low pair first
    // (`BitpairOutFileBuf::write`, filebuf.h:589) -- the same packing as the BWT.
    let mut o4: Vec<u8> = vec![0u8; (t.len() + 3) / 4];
    for (i, &c) in t.iter().enumerate() { o4[i >> 2] |= c << ((i & 3) * 2); }

    fs::write(format!("{out_prefix}.3.ht2"), &o3).expect("write .3.ht2");
    fs::write(format!("{out_prefix}.4.ht2"), &o4).expect("write .4.ht2");
    fs::write(format!("{out_prefix}.1.ht2"), &o).expect("write .1.ht2");
    fs::write(format!("{out_prefix}.2.ht2"), &o2).expect("write .2.ht2");
    println!("wrote {out_prefix}.1.ht2 ({} bytes) and .2.ht2 ({} bytes)", o.len(), o2.len());
    println!("  len {len}  nPat {}  nFrag {}  numSides {num_sides}  ftabLen {ftab_len}  offsLen {offs_len}",
             r.names.len(), r.frags.len());
    assert_eq!(sa_sample.len(), offs_len, "SA sample count vs derived offsLen");

    // ---- compare --------------------------------------------------------
    if a.len() < 4 { return; }
    let mut bad = false;
    for (suffix, ours) in [(".1.ht2", &o), (".2.ht2", &o2), (".3.ht2", &o3), (".4.ht2", &o4)] {
        let path = a[3].replace(".1.ht2", suffix);
        let theirs = match fs::read(&path) { Ok(b) => b, Err(e) => { println!("  {path}: {e}"); bad = true; continue; } };
        if theirs.len() != ours.len() {
            println!("{suffix}: LENGTH {} vs theirs {}", ours.len(), theirs.len());
            bad = true;
        }
        let n = ours.len().min(theirs.len());
        let mut first = None;
        let mut ndiff = 0usize;
        for i in 0..n {
            if ours[i] != theirs[i] { ndiff += 1; if first.is_none() { first = Some(i); } }
        }
        if ndiff == 0 && ours.len() == theirs.len() {
            println!("{suffix}: BYTE-IDENTICAL ({} bytes)", ours.len());
        } else {
            bad = true;
            println!("{suffix}: {ndiff} differing bytes of {n}; first at offset {:?}", first);
            if let Some(f) = first {
                let lo = f.saturating_sub(8);
                let hi = (f + 24).min(n);
                println!("    ours   {:02x?}", &ours[lo..hi]);
                println!("    theirs {:02x?}", &theirs[lo..hi]);
            }
        }
        // header sanity: fail loudly if our assumed defaults differ from theirs
        if suffix == ".1.ht2" && theirs.len() >= 44 {
            let u = |o: usize| u32::from_le_bytes(theirs[o..o + 4].try_into().unwrap());
            for (name, got, want) in [("lineRate", u(20), line_rate), ("offRate", u(28), off_rate),
                                      ("ftabChars", u(32), ftab_chars), ("version", u(4), version)] {
                if got != want { println!("    NOTE assumed {name} {want}, index has {got}"); }
            }
        }
    }
    if bad { std::process::exit(1); }
    println!("\nRUNG 2 COMPLETE: .1/.2/.3/.4 reproduced byte for byte");
}
