//! Rebuild individual local-index windows and compare their headers against an
//! existing `.5`. One window takes seconds once the inputs are parsed, so a
//! change to the variant-thinning search can be tested without the 17-hour
//! whole-genome build it would otherwise need.
//!
//! usage: ht2win <ref.fa> <snp> <haplotype> <index.5.ht2[l]> <win>[,<win>...]

#[path = "../graph.rs"]
mod graph;
#[path = "../ext.rs"]
mod ext;
#[path = "../gfmbuild.rs"]
mod gfmbuild;
#[path = "../local5.rs"]
mod local5;

use std::{env, fs, process};

/// Walk a `.5` and return (len, gbwtLen, numNodes, eftabLen) per window.
fn scan(path: &str, w: usize) -> Vec<(u16, u16, u16, u16)> {
    let b = fs::read(path).expect(".5");
    let u32a = |o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
    let u16a = |o: usize| u16::from_le_bytes(b[o..o + 2].try_into().unwrap());
    let idxa = |o: usize| -> u64 {
        let mut v = [0u8; 8]; v[..w].copy_from_slice(&b[o..o + w]); u64::from_le_bytes(v)
    };
    let mut p = 4usize;
    let n = idxa(p) as usize; p += w;
    let lr = u32a(p) as usize; p += 4; p += 4; p += 4;
    let fc = u32a(p) as usize; p += 4; p += 4;
    let ftab_len = (1usize << (2 * fc)) + 1;
    let side = 1usize << lr;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        p += 3 * w;
        let (ln, gb, nn, ef) = (u16a(p), u16a(p + 2), u16a(p + 4), u16a(p + 6));
        p += 8;
        out.push((ln, gb, nn, ef));
        if ln == 0 { continue; }
        let linear = ln + 1 == gb || gb == 0;
        let gl = if gb == 0 { ln as usize + 1 } else { gb as usize };
        let gs = if linear { gl / 4 + 1 } else { gl / 2 + 1 };
        let sgs = side - 2 * if linear { 4 } else { 6 };
        let tot = ((gs + sgs - 1) / sgs) * side;
        let npat = u16a(p) as usize; p += 2 + npat * 2;
        let nfrag = u16a(p) as usize; p += 2 + nfrag * 6;
        p += tot;
        let nz = u16a(p) as usize; p += 2 + nz * 2;
        p += 10 + ftab_len * 2 + ef as usize * 2;
    }
    out
}

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 6 {
        eprintln!("usage: ht2win <ref.fa> <snp> <hap> <index.5.ht2[l]> <win>[,<win>...]");
        process::exit(2);
    }
    let w = if a[4].ends_with("l") { 8usize } else { 4 };
    let want: Vec<usize> = a[5].split(',').map(|x| x.parse().expect("window index")).collect();
    let theirs = scan(&a[4], w);
    println!("reference index: {} windows", theirs.len());

    let p = graph::parse(&a[1], &a[2], &a[3]);
    let plan = local5::window_plan(&p);
    println!("our plan: {} windows\n", plan.len());

    let mut bad = 0;
    for &i in &want {
        let e = plan[i];
        let o = local5::build_window(&p, e, w);
        let h = &o.sec[3 * w..3 * w + 8];
        let g = |k: usize| u16::from_le_bytes(h[k * 2..k * 2 + 2].try_into().unwrap());
        let (ours, t) = ((g(0), g(1), g(2), g(3)), theirs[i]);
        let ok = ours == t;
        if !ok { bad += 1; }
        println!("window {i} (joined @{}, len {}): {}", e.2, e.3, if ok { "MATCH" } else { "DIFFERS" });
        println!("   ours len={} gbwtLen={} numNodes={} eftabLen={}", ours.0, ours.1, ours.2, ours.3);
        if !ok { println!("   ref  len={} gbwtLen={} numNodes={} eftabLen={}", t.0, t.1, t.2, t.3); }
    }
    println!("\n{} of {} windows match", want.len() - bad, want.len());
    if bad > 0 { process::exit(1); }
}
