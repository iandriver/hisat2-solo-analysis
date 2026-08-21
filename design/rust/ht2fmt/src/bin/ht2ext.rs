//! Step 3 — the doubling loop in external memory, checked against the whole
//! generation curve rather than just its first line.
//!
//! The loop itself lives in `doubling.rs`, shared with `ht2wg`, so the curve
//! this binary verifies is the same code that emits the index. See that module
//! for why the join is a sort-merge and what stays resident.
//!
//! The oracle is `ht2path`: the in-memory loop already reproduces HISAT2's
//! generation curve exactly on three graphs, so "identical curve" is a hard test
//! rather than a plausibility check. The reference graph is always built
//! fragment by fragment to disk here — `ht2frag` is what checks that path
//! against the whole-graph construction.

#[path = "../graph.rs"]
mod graph;
#[path = "../ext.rs"]
mod ext;
#[path = "../ckpt.rs"]
mod ckpt;
#[path = "../doubling.rs"]
mod doubling;

use std::env;
use std::path::PathBuf;

fn main() -> std::io::Result<()> {
    let a: Vec<String> = env::args().collect();
    if a.len() < 5 {
        eprintln!("usage: ht2ext <reference.fa> <snp> <haplotype> <workdir> [<build.log>] [budget_records]");
        eprintln!("  budget_records defaults to 4096 -- deliberately tiny, to force real spilling");
        eprintln!("  HT2_CHUNK sets the graph fragment size in bases (default 262144)");
        std::process::exit(2);
    }
    let (fa, snp, hap, wd) = (&a[1], &a[2], &a[3], PathBuf::from(&a[4]));
    let budget: usize = a.get(6).and_then(|x| x.parse().ok()).unwrap_or(4096);
    let chunk: u32 = env::var("HT2_CHUNK").ok().and_then(|x| x.parse().ok()).unwrap_or(1 << 18);

    let mut ck = ckpt::Ckpt::default();
    let d = doubling::run(fa, snp, hap, &wd, budget, chunk, true, &mut ck, false)?;
    println!("\n  resident ceiling during the run: one sort buffer of {} KB plus one key block",
             budget * ext::REC / 1024);

    if a.len() < 6 || a[5] == "-" { return Ok(()); }
    if !doubling::check_curve(&d.curve, &a[5]) { std::process::exit(1); }
    Ok(())
}
