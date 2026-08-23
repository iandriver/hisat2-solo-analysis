//! Predict a `.7` variant database against an existing index, without building.
//!
//! `.7` is a straight serialisation of the alt and haplotype lists, so it is a
//! pure function of the three input files -- no graph, no doubling. That makes
//! it checkable in minutes against an index that already exists, which is the
//! only practical way to iterate on the haplotype sort order: the ordering
//! among tied `(left, right)` haplotypes comes from libstdc++'s `std::sort`,
//! and the whole-genome artifact is the only oracle for it on this machine.
//!
//! usage: ht2predict <ref.fa> <snp> <haplotype> <index.7.ht2l> [--large]

#[path = "../graph.rs"]
mod graph;
#[path = "../ext.rs"]
mod ext;

use std::{env, fs, process};

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 5 {
        eprintln!("usage: ht2predict <ref.fa> <snp> <haplotype> <index.7.ht2[l]>");
        process::exit(2);
    }
    let large = a[4].ends_with("l");
    let w = if large { 8usize } else { 4 };
    let b = fs::read(&a[4]).expect(".7");
    let mut o = 0usize;
    let mut rd = |n: usize, o: &mut usize| -> u64 {
        let mut v = [0u8; 8];
        v[..n].copy_from_slice(&b[*o..*o + n]);
        *o += n;
        u64::from_le_bytes(v)
    };
    assert_eq!(rd(4, &mut o), 1, "endianness sentinel");
    let n_alt = rd(w, &mut o) as usize;
    o += n_alt * (w + 4 + w + 8);            // pos, type, len, seq
    let n_hap = rd(w, &mut o) as usize;
    println!("artifact: {n_alt} alts, {n_hap} haplotypes");

    let p = graph::parse(&a[1], &a[2], &a[3]);
    println!("ours:     {} alts, {} haplotypes", p.alts.len(), p.haps.len());
    let fb = graph::gnusort::HEAP_FALLBACKS.load(std::sync::atomic::Ordering::Relaxed);
    println!("heapsort fallbacks during the haplotype sort: {fb}");

    if p.alts.len() != n_alt || p.haps.len() != n_hap {
        println!("\nCOUNT MISMATCH -- stopping before the record comparison");
        process::exit(1);
    }

    let (mut same, mut diff) = (0usize, 0usize);
    let mut first: Option<String> = None;
    for i in 0..n_hap {
        let l = rd(w, &mut o);
        let r = rd(w, &mut o);
        let n = rd(w, &mut o) as usize;
        let ids: Vec<u64> = (0..n).map(|_| rd(w, &mut o)).collect();
        let h = &p.haps[i];
        let ours: Vec<u64> = h.alts.iter().map(|&x| x as u64).collect();
        if h.left as u64 == l && h.right as u64 == r && ours == ids { same += 1; }
        else {
            diff += 1;
            if first.is_none() {
                first = Some(format!("  at {i}\n    ours     ({}, {}, {:?})\n    artifact ({l}, {r}, {ids:?})",
                                     h.left, h.right, ours));
            }
        }
    }
    println!("\nrecords matching: {same} / {n_hap}   differing: {diff}");
    if let Some(f) = first { println!("first difference:\n{f}"); }
    if o != b.len() { println!("WARNING: consumed {o} of {} bytes", b.len()); }
    if diff == 0 { println!("\n.7 WILL BE BYTE-IDENTICAL"); } else { process::exit(1); }
}
