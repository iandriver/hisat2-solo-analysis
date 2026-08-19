//! Step 3.5 — check that building the graph fragment by fragment gives the
//! identical graph, then measure what it costs.
//!
//! Fragmentation is only worth anything if the stitched result is the same
//! graph, so the test is a direct diff of nodes and edges against the global
//! build, not a comparison of counts.

#[path = "../graph.rs"]
mod graph;

use std::env;

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: ht2frag <reference.fa> <snp> <haplotype> [chunk_bp]");
        std::process::exit(2);
    }
    let chunk: u32 = a.get(4).and_then(|x| x.parse().ok()).unwrap_or(1 << 20);

    let p = graph::parse(&a[1], &a[2], &a[3]);
    let len = p.text.len() as u32;
    let bounds = graph::fragment_bounds(len, &p.alts, chunk);
    let widths: Vec<u32> = bounds.iter().map(|&(x, y)| y - x).collect();
    println!("{len} bp, {} alts -> {} fragments (chunk {chunk}); widths min {} max {}",
             p.alts.len(), bounds.len(),
             widths.iter().min().copied().unwrap_or(0),
             widths.iter().max().copied().unwrap_or(0));
    drop(p);

    // HT2_FRAG_ONLY=1 builds only the fragmented path, so /usr/bin/time -l
    // measures its peak rather than the global build's.
    if std::env::var("HT2_FRAG_ONLY").is_ok() {
        let f = graph::build_fragmented(&a[1], &a[2], &a[3], chunk);
        println!("  fragmented: {} nodes, {} edges", f.nodes.len(), f.edges.len());
        return;
    }
    let g = graph::build(&a[1], &a[2], &a[3]);
    let f = graph::build_fragmented(&a[1], &a[2], &a[3], chunk);

    println!("  global    : {} nodes, {} edges, last {}", g.nodes.len(), g.edges.len(), g.last_node);
    println!("  fragmented: {} nodes, {} edges, last {}", f.nodes.len(), f.edges.len(), f.last_node);

    let mut bad = 0usize;
    if g.nodes.len() != f.nodes.len() { println!("  NODE COUNT differs"); bad += 1; }
    else {
        let d = g.nodes.iter().zip(f.nodes.iter()).enumerate()
            .filter(|(_, (x, y))| x != y).take(5).collect::<Vec<_>>();
        if !d.is_empty() {
            bad += 1;
            println!("  {} node entries differ; first few:", 
                     g.nodes.iter().zip(f.nodes.iter()).filter(|(x, y)| x != y).count());
            for (i, (x, y)) in d { println!("    node {i}: global {:?} frag {:?}", x, y); }
        }
    }
    // edges are order-sensitive downstream, so compare as sequences
    if g.edges.len() != f.edges.len() { println!("  EDGE COUNT differs"); bad += 1; }
    else {
        let n = g.edges.iter().zip(f.edges.iter()).filter(|(x, y)| x != y).count();
        if n > 0 {
            bad += 1;
            println!("  {n} edge entries differ in order; first few:");
            for (i, (x, y)) in g.edges.iter().zip(f.edges.iter()).enumerate()
                                .filter(|(_, (x, y))| x != y).take(5) {
                println!("    edge {i}: global {:?} frag {:?}", x, y);
            }
            let mut gs = g.edges.clone(); gs.sort_unstable();
            let mut fs = f.edges.clone(); fs.sort_unstable();
            println!("  as multisets: {}", if gs == fs { "IDENTICAL (ordering only)" } else { "DIFFERENT" });
        }
    }
    if bad == 0 { println!("\nFRAGMENTED GRAPH IDENTICAL TO GLOBAL"); }
    else { println!("\nFRAGMENTED GRAPH DIFFERS"); std::process::exit(1); }
}
