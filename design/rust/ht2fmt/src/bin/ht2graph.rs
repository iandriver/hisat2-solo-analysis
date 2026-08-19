//! Rung 3, step 1 — build the reference graph and check its size against
//! HISAT2's own `Generation 0` line.
//!
//! `PathGraph::makeFromRef` (`gbwt_graph.h:1821`) sets
//! `temp_nodes = base.edges.size() + 1`, and `printInfo` reports `temp_nodes`
//! as the left-hand number. So the first line HISAT2 prints,
//!
//!     Generation 0 (236 -> 236 nodes, 0 ranks)
//!
//! is an exact, free oracle for the edge count of the reference graph. That
//! makes the graph testable before any of the prefix-doubling machinery exists.
//!
//! The construction itself lives in `graph.rs`; this binary is the check.

#[path = "../graph.rs"]
mod graph;
// graph.rs streams the edge list through ext::sort_pairs_external, so every
// binary that includes it needs ext in scope too.
#[allow(dead_code)]
#[path = "../ext.rs"]
mod ext;

use std::{env, fs};

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: ht2graph <reference.fa> <snp> <haplotype> [<build.log>]");
        eprintln!("  reports edges+1; compares against the log's Generation 0 line if given");
        std::process::exit(2);
    }
    let b = graph::build(&a[1], &a[2], &a[3]);

    println!("reference {} bp, {} alts, {} haplotypes ({} snps and {} haplotypes had no \
joined position, {} haplotypes out of order)",
             b.text.len(), b.n_alts, b.n_haps, b.dropped_snps, b.dropped_haps, b.out_of_order_haps);
    println!("  before determinizing: nodes {} edges {} ({} reference + {} haplotype)",
             b.pre_nodes, b.pre_edges, b.base_edges, b.pre_edges - b.base_edges);
    println!("  after  determinizing: nodes {} edges {} (removed {} nodes, {} edges)",
             b.nodes.len(), b.edges.len(),
             b.pre_nodes as i64 - b.nodes.len() as i64, b.pre_edges as i64 - b.edges.len() as i64);

    let bad = graph::reverse_deterministic(&b.nodes, &b.edges);
    println!("  reverse-deterministic: {}",
             if bad == 0 { "yes".to_string() }
             else { format!("NO -- {bad} nodes have two same-label predecessors") });
    if bad > 0 { std::process::exit(1); }

    let temp_nodes = b.edges.len() + 1;
    println!("  temp_nodes = edges + 1 = {temp_nodes}");

    if a.len() < 5 { return; }
    let log = fs::read_to_string(&a[4]).expect("log");
    let g0 = log.lines().find(|l| l.starts_with("Generation 0 "))
        .expect("no 'Generation 0' line in the log");
    // "Generation 0 (236 -> 236 nodes, 0 ranks)" -- the number after '(' is
    // temp_nodes; the generation index itself parses too, so anchor on the paren.
    let theirs: usize = g0.split('(').nth(1).and_then(|t| t.split_whitespace().next())
        .and_then(|t| t.parse().ok()).expect("parse Generation 0");
    println!("  HISAT2 Generation 0 = {theirs}");
    if theirs == temp_nodes {
        println!("\nGRAPH SIZE MATCHES");
    } else {
        println!("\nGRAPH SIZE MISMATCH: ours {temp_nodes}, theirs {theirs} (off by {})",
                 temp_nodes as i64 - theirs as i64);
        std::process::exit(1);
    }
}
