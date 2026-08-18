//! Is the SA stored in .2.ht2 actually sorted as suffixes?
//!
//! A suffix array is unique, so if HISAT2's sampled values are in true
//! lexicographic suffix order then ours cannot legitimately differ, and the
//! discrepancy must be in how we line the rows up. If they are NOT sorted,
//! then whatever nextSuffix() emits is not a plain suffix array, which is a
//! finding in itself.
use std::{env, fs};

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 3 { eprintln!("usage: ht2sacheck <index.2.ht2> <reference.fa>"); std::process::exit(2); }
    let sab = fs::read(&a[1]).expect("read .2.ht2");
    let mut t = Vec::new();
    for ln in fs::read_to_string(&a[2]).expect("fa").lines() {
        if ln.starts_with('>') { continue; }
        for c in ln.bytes() { match c.to_ascii_uppercase() {
            b'A'=>t.push(0u8), b'C'=>t.push(1), b'G'=>t.push(2), b'T'=>t.push(3), _=>{} } }
    }
    let n = t.len();
    let get = |k: usize| u32::from_le_bytes(sab[4+k*4..8+k*4].try_into().unwrap()) as usize;
    let cnt = (sab.len()-4)/4;
    // compare suffixes t[i..] vs t[j..]; shorter-is-prefix sorts first ($ ends it)
    let cmp = |i: usize, j: usize| -> std::cmp::Ordering {
        let (mut x, mut y) = (i, j);
        while x < n && y < n {
            if t[x] != t[y] { return t[x].cmp(&t[y]); }
            x += 1; y += 1;
        }
        (n - i).cmp(&(n - j))   // the one that ran out first is shorter, so smaller
    };
    let mut bad = 0usize; let mut shown = 0;
    for k in 1..cnt {
        let (p, q) = (get(k-1), get(k));
        if p >= n || q >= n { continue; }          // the $ entry (value == n)
        if cmp(p, q) != std::cmp::Ordering::Less {
            bad += 1;
            if shown < 6 {
                println!("  NOT SORTED at sampled rows {}->{}: SA {} -> {}", (k-1)*16, k*16, p, q);
                shown += 1;
            }
        }
    }
    println!("\nsampled adjacent pairs out of order: {bad} of {}", cnt-1);
    println!("(sampled rows are 16 apart, so a true suffix array must still be strictly increasing here)");
    // where does the $ entry sit?
    for k in 0..cnt { if get(k) >= n { println!("$ (value {}) appears at sampled row {}", get(k), k*16); } }
}
