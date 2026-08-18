//! Recover HISAT2's FULL suffix array by LF-stepping the stored BWT.
//!
//! `.2.ht2` only samples one row in 2^offRate, which is too coarse to localise
//! where our SA and theirs part company. Walking the BWT gives every row.
//!
//! From blockwise_sa.h the empty suffix is appended LAST, so row `len` is the
//! `$` row. Starting there and applying LF repeatedly visits the row of suffix
//! len-1, then len-2, and so on -- which both reconstructs the text backwards
//! and labels every row with its SA value. Two built-in checks: the
//! reconstruction must equal the reference, and the final row must equal the
//! stored zOffs (the row of suffix 0).
use std::convert::TryInto;
use std::{env, fs};

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 3 { eprintln!("usage: ht2sa <index.1.ht2> <reference.fa>"); std::process::exit(2); }
    let buf = fs::read(&a[1]).expect("index");
    let u = |o: usize| -> u64 { u32::from_le_bytes(buf[o..o+4].try_into().unwrap()) as u64 };
    let len = u(8); let gbwt_len = u(12); let line_rate = u(20) as u32;
    let n_pat = u(44); let n_frag = u(44 + 4 + n_pat as usize * 4);
    let side_sz = 1u64 << line_rate;
    let side_gbwt_sz = side_sz - 16;
    let side_gbwt_len = side_gbwt_sz * 4;
    let gbwt_sz = gbwt_len / 4 + 1;
    let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
    let gbwt_off = 44 + 4 + n_pat as usize * 4 + 4 + n_frag as usize * 12;
    let zoff_at = gbwt_off + (num_sides * side_sz) as usize;
    let z = u(zoff_at + 4);
    let fchr: Vec<u64> = (0..5).map(|i| u(zoff_at + 8 + i * 4)).collect();

    let mut bwt = vec![0u8; gbwt_len as usize];
    for r in 0..gbwt_len as usize {
        let side = r as u64 / side_gbwt_len; let off = r as u64 % side_gbwt_len;
        bwt[r] = ((buf[gbwt_off + (side * side_sz + (off >> 2)) as usize] >> (2 * (off & 3))) & 3) as u8;
    }
    let mut occ = vec![[0u32; 4]; gbwt_len as usize + 1];
    for r in 0..gbwt_len as usize {
        occ[r + 1] = occ[r];
        if r as u64 != z { occ[r + 1][bwt[r] as usize] += 1; }
    }

    let mut t = Vec::new();
    for ln in fs::read_to_string(&a[2]).expect("fa").lines() {
        if ln.starts_with('>') { continue; }
        for c in ln.bytes() { match c.to_ascii_uppercase() {
            b'A'=>t.push(0u8), b'C'=>t.push(1), b'G'=>t.push(2), b'T'=>t.push(3), _=>{} } }
    }

    // $ is last: row `len` is the empty suffix.
    let mut sa = vec![u32::MAX; gbwt_len as usize];
    let mut r = len;
    sa[r as usize] = len as u32;
    let mut out = Vec::with_capacity(t.len());
    for k in 0..t.len() {
        let c = bwt[r as usize] as usize;
        out.push(c as u8);
        r = fchr[c] + occ[r as usize][c] as u64;
        sa[r as usize] = (t.len() - 1 - k) as u32;
    }
    out.reverse();
    println!("text reconstruction: {}", if out == t { "MATCHES reference" } else { "DIFFERS" });
    println!("final row after {} LF steps: {}  (stored zOffs {})  {}",
             t.len(), r, z, if r == z { "match" } else { "MISMATCH" });
    let unlabelled = sa.iter().filter(|&&v| v == u32::MAX).count();
    println!("rows labelled: {}/{} (unlabelled {})", gbwt_len as usize - unlabelled, gbwt_len, unlabelled);
    if out == t && r == z {
        fs::write("their_sa.bin", sa.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
        println!("wrote their_sa.bin ({} u32 entries)", sa.len());
    }
}
