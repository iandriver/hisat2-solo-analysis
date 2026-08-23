//! Emit `.5.ht2` / `.6.ht2` and byte-compare against an index `hisat2-build`
//! wrote.
//!
//! The emitter itself lives in `local5.rs`, shared with `ht2wg`, so what this
//! checks is the same code that a whole-genome build runs.

#[path = "../graph.rs"]
mod graph;
#[allow(dead_code)]
#[path = "../ext.rs"]
mod ext;
#[path = "../gfmbuild.rs"]
mod gfmbuild;
#[path = "../local5.rs"]
mod local5;

use std::{env, fs};

fn main() -> std::io::Result<()> {
    let a: Vec<String> = env::args().collect();
    if a.len() < 5 {
        eprintln!("usage: ht2emit5 <reference.fa> <snp> <haplotype> <index_prefix>");
        std::process::exit(2);
    }
    let large = env::var("HT2_LARGE").is_ok();
    let (w, ext) = if large { (8usize, "ht2l") } else { (4usize, "ht2") };
    let p = graph::parse(&a[1], &a[2], &a[3]);
    let (mut o5, mut o6): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    local5::emit(&p, &mut o5, &mut o6, w, true)?;

    let t5 = fs::read(format!("{}.5.{ext}", a[4])).expect(".5");
    let t6 = fs::read(format!("{}.6.{ext}", a[4])).expect(".6");
    let mut bad = false;
    for (name, ours, theirs) in [(".5", &o5, &t5), (".6", &o6, &t6)] {
        let n = ours.len().min(theirs.len());
        let d = (0..n).filter(|&i| ours[i] != theirs[i]).count();
        let first = (0..n).find(|&i| ours[i] != theirs[i]);
        if d == 0 && ours.len() == theirs.len() {
            println!("{name}.{ext}: BYTE-IDENTICAL ({} bytes)", ours.len());
        } else {
            bad = true;
            // Keep the mismatching output so the differing bytes can be read
            // back with the same parser that reads theirs.
            let _ = fs::write(format!("{}.ours{}.{ext}", a[4], name), ours.as_slice());
            println!("{name}.{ext}: {} of {n} bytes match (ours {} bytes, theirs {}){}",
                     n - d, ours.len(), theirs.len(),
                     match first { Some(f) => format!("; first differing byte at {f}"), None => String::new() });
        }
    }
    if bad { std::process::exit(1); }
    Ok(())
}
