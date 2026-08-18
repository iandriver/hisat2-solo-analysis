//! Build ftab + eftab and check them against a reference index.
//!
//! From GFM::buildToDisk (gfm.h:5148, and the finalisation at :5430):
//!   * walking the SA in row order, a suffix with (len - saElt) >= ftabChars
//!     contributes ftab[sufInt + 1]++, where sufInt packs the first ftabChars
//!     characters two bits each, leftmost most significant;
//!   * a shorter suffix instead bumps absorbCnt, which is stored into
//!     absorbFtab[sufInt] at the next long suffix (and into absorbFtab[ftabLen-1]
//!     if the walk ends while absorbCnt > 0);
//!   * finalisation turns ftab into a running prefix sum. Entries that absorb
//!     short suffixes store [lo, hi] in eftab and keep a marker
//!     `eftabCur ^ INDEX_MAX` in ftab, resolved by ftabLo/ftabHi (gfm.h:2617).
//!
//! The suffix order uses HISAT2's rule: pad past the end of the text with a
//! character LARGER than any real one, so a suffix that is a prefix of another
//! sorts after it and the empty suffix sorts last.
use std::convert::TryInto;
use std::{env, fs};

const HI: u32 = u32::MAX;

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

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 3 { eprintln!("usage: ht2ftab <reference.fa> <index.1.ht2>"); std::process::exit(2); }
    let mut t = Vec::new();
    for ln in fs::read_to_string(&a[1]).expect("fa").lines() {
        if ln.starts_with('>') { continue; }
        for c in ln.bytes() { match c.to_ascii_uppercase() {
            b'A'=>t.push(0u8), b'C'=>t.push(1), b'G'=>t.push(2), b'T'=>t.push(3), _=>{} } }
    }
    let n = t.len() as u64;
    let buf = fs::read(&a[2]).expect("index");
    let u = |o: usize| -> u64 { u32::from_le_bytes(buf[o..o+4].try_into().unwrap()) as u64 };
    let gbwt_len = u(12);
    let ftab_chars = u(32) as usize;
    let eftab_len = u(36) as usize;
    let ftab_len = (1usize << (ftab_chars * 2)) + 1;

    let sa = hisat2_suffix_array(&t);

    // pass 1: raw counts and absorbFtab
    let mut ftab = vec![0u64; ftab_len];
    let mut absorb = vec![0u64; ftab_len];
    let mut absorb_cnt = 0u64;
    for &e in sa.iter() {
        let se = e as u64;
        // The empty suffix (position n) is NOT skipped: buildToDisk sees
        // len - saElt == 0, which is < ftabChars, so it counts as a short
        // suffix and bumps absorbCnt like any other.
        if n - se >= ftab_chars as u64 {
            let mut si = 0u64;
            for i in 0..ftab_chars { si = (si << 2) | t[(se as usize) + i] as u64; }
            ftab[(si + 1) as usize] += 1;
            if absorb_cnt > 0 { absorb[si as usize] = absorb_cnt; absorb_cnt = 0; }
        } else {
            absorb_cnt += 1;
        }
    }
    if absorb_cnt > 0 { absorb[ftab_len - 1] = absorb_cnt; }

    // pass 2: prefix sum, with absorbed entries pushed into eftab
    let mut eftab = vec![0u64; eftab_len];
    let mut ftab_out = vec![0u64; ftab_len];
    let mut ecur = 0usize;
    let hi_of = |i: usize, ftab_out: &Vec<u64>, eftab: &Vec<u64>| -> u64 {
        let v = ftab_out[i];
        if v <= gbwt_len { v } else { eftab[((v ^ u32::MAX as u64) as usize) * 2 + 1] }
    };
    for i in 1..ftab_len {
        let lo = ftab[i] + hi_of(i - 1, &ftab_out, &eftab);
        if absorb[i] > 0 {
            let hi = lo + absorb[i];
            eftab[ecur * 2] = lo;
            eftab[ecur * 2 + 1] = hi;
            ftab_out[i] = (ecur as u64) ^ (u32::MAX as u64);
            ecur += 1;
        } else {
            ftab_out[i] = lo;
        }
    }

    // compare
    let n_pat = u(44);
    let n_frag = u(44 + 4 + n_pat as usize * 4);
    let line_rate = u(20) as u32;
    let side_sz = 1u64 << line_rate;
    let side_gbwt_sz = side_sz - 16;
    let gbwt_sz = gbwt_len / 4 + 1;
    let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
    let gbwt_off = 44 + 4 + n_pat as usize * 4 + 4 + n_frag as usize * 12;
    let zoff_at = gbwt_off + (num_sides * side_sz) as usize;
    let ftab_at = zoff_at + 4 + 4 + 20;
    let eftab_at = ftab_at + ftab_len * 4;

    let mut bad = 0usize; let mut shown = 0;
    for i in 0..ftab_len {
        let stored = u(ftab_at + i * 4);
        if stored != ftab_out[i] {
            bad += 1;
            if shown < 5 { println!("  ftab[{i}] ours {} stored {}", ftab_out[i], stored); shown += 1; }
        }
    }
    println!("ftab:  {}/{} entries match", ftab_len - bad, ftab_len);
    let mut ebad = 0usize;
    for i in 0..eftab_len {
        let stored = u(eftab_at + i * 4);
        if stored != eftab[i] { ebad += 1; println!("  eftab[{i}] ours {} stored {}", eftab[i], stored); }
    }
    println!("eftab: {}/{} entries match  (eftabLen {}, used {} slots)", eftab_len - ebad, eftab_len, eftab_len, ecur * 2);
    if bad == 0 && ebad == 0 { println!("\nFTAB + EFTAB EXACT"); } else { std::process::exit(1); }
}
