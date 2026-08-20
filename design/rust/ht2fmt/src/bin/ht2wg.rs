//! Step 4 — build a graph index end to end on the external path, and write it.
//!
//! `ht2path` and `ht2emit5` reproduce every section of a HISAT2 graph index byte
//! for byte, but both hold the whole path graph in RAM, which is the 625 GiB the
//! C++ needs. This binary runs the same construction with the fragmented graph
//! builder and the external doubling loop underneath it, and actually WRITES the
//! files rather than comparing against an existing index.
//!
//! Every stage has an in-memory oracle, so `--verify` against a `hisat2-build`
//! index is a byte-level test of the whole chain rather than of its output alone.

#[path = "../graph.rs"]
mod graph;
#[path = "../ext.rs"]
mod ext;
#[path = "../doubling.rs"]
mod doubling;
#[path = "../wgemit.rs"]
mod wgemit;
#[path = "../refin.rs"]
mod refin;
#[path = "../gfmbuild.rs"]
mod gfmbuild;
#[path = "../local5.rs"]
mod local5;

use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::os::unix::fs::FileExt;
use std::path::PathBuf;
use std::{env, process};

const OFF_RATE: u32 = 4;
const FTAB_CHARS: u32 = 10;
const VERSION: u32 = (2 << 24) | (2 << 16) | (3 << 8);   // 2.2.3

/// `writeI32` — always four bytes, whatever `index_t` is.
fn put(v: &mut Vec<u8>, x: u32) { v.extend_from_slice(&x.to_le_bytes()); }

fn tail_of(p: &str) -> &str { p.rsplit('/').next().unwrap_or(p) }

fn compare(ours: &str, theirs: &str) -> bool {
    let a = match fs::read(ours) { Ok(b) => b, Err(e) => { println!("  {ours}: {e}"); return false } };
    let b = match fs::read(theirs) { Ok(b) => b, Err(e) => { println!("  {theirs}: {e}"); return false } };
    let n = a.len().min(b.len());
    let first = (0..n).find(|&i| a[i] != b[i]);
    if a.len() == b.len() && first.is_none() {
        println!("  {:>12}: BYTE-IDENTICAL ({} bytes)", tail_of(ours), a.len());
        true
    } else {
        let d = (0..n).filter(|&i| a[i] != b[i]).count();
        println!("  {:>12}: {d} differing bytes of {n} (ours {}, theirs {}){}",
                 tail_of(ours), a.len(), b.len(),
                 match first { Some(f) => format!("; first at {f}"), None => String::new() });
        false
    }
}

fn main() -> std::io::Result<()> {
    let a: Vec<String> = env::args().collect();
    if a.len() < 6 {
        eprintln!("usage: ht2wg <reference.fa> <snp> <haplotype> <workdir> <out_prefix> \\");
        eprintln!("             [budget_records] [verify_prefix]");
        eprintln!("  HT2_CHUNK sets the graph fragment size in bases (default 262144)");
        eprintln!("  HT2_KEEP=1 leaves the workdir in place");
        process::exit(2);
    }
    // A large index is `index_t = uint64_t`, which changes every `writeIndex`
    // field, moves the default line rate from 7 to 8 (six 8-byte tallies per
    // side instead of six 4-byte ones), and renames the files to `.ht2l`.
    let large = env::var("HT2_LARGE").is_ok();
    let w: usize = if large { 8 } else { 4 };
    let ext = if large { "ht2l" } else { "ht2" };
    let line_rate: u32 = if large { 8 } else { 7 };
    let put_idx = |v: &mut Vec<u8>, x: u64| v.extend_from_slice(&x.to_le_bytes()[0..w]);
    let (fa, snp, hap) = (&a[1], &a[2], &a[3]);
    let wd = PathBuf::from(&a[4]);
    let out = a[5].clone();
    let budget: usize = a.get(6).and_then(|x| x.parse().ok()).unwrap_or(1 << 20);
    let verify = a.get(7).cloned();
    let chunk: u32 = env::var("HT2_CHUNK").ok().and_then(|x| x.parse().ok()).unwrap_or(1 << 18);
    fs::create_dir_all(&wd)?;

    // ---- the reference front end, .3 and .4 -----------------------------
    // Done first and dropped before the graph stage, so the joined text is not
    // resident twice.
    let (names, plen, frags, len) = {
        let r = refin::read_fasta(fa);
        let mut o3: Vec<u8> = Vec::new();
        put(&mut o3, 1);
        put_idx(&mut o3, r.recs.len() as u64);
        for rec in &r.recs {
            put_idx(&mut o3, rec.off as u64); put_idx(&mut o3, rec.len as u64);
            o3.push(rec.first as u8);
        }
        fs::write(format!("{out}.3.{ext}"), &o3)?;
        let mut o4 = vec![0u8; (r.text.len() + 3) / 4];
        for (i, &c) in r.text.iter().enumerate() { o4[i >> 2] |= c << ((i & 3) * 2); }
        fs::write(format!("{out}.4.{ext}"), &o4)?;
        (r.names, r.plen, r.frags, r.text.len() as u32)
    };
    println!("reference: {len} joined bp over {} sequence(s), {} fragment(s)",
             names.len(), frags.len());

    // ---- the path graph, fragmented and external ------------------------
    let d = doubling::run(fa, snp, hap, &wd, budget, chunk, true)?;
    let rs = wgemit::generate_edges(&wd, &d.cur, budget, true)?;
    let _ = fs::remove_file(&d.cur);
    let g = wgemit::geom(rs.gbwt_len, line_rate, w);

    // ---- .1.ht2 through the gbwt block, and .2.ht2 ----------------------
    let mut head: Vec<u8> = Vec::new();
    put(&mut head, 1);                    // endianness sentinel
    put(&mut head, VERSION);
    put_idx(&mut head, len as u64);
    put_idx(&mut head, rs.gbwt_len);
    put_idx(&mut head, rs.n_nodes);
    put(&mut head, line_rate);
    put(&mut head, 2);                    // linesPerSide, written as a literal 2
    put(&mut head, OFF_RATE);
    put(&mut head, FTAB_CHARS);
    let eftab_len_at = head.len() as u64;
    put_idx(&mut head, 0);                // eftabLen, back-patched below
    put(&mut head, (-1i32) as u32);       // -flags
    put_idx(&mut head, names.len() as u64);
    for &p in &plen { put_idx(&mut head, p as u64); }
    put_idx(&mut head, frags.len() as u64);
    for f in &frags {
        put_idx(&mut head, f.joined_off as u64);
        put_idx(&mut head, f.text_id as u64);
        put_idx(&mut head, f.text_off as u64);
    }
    let gbwt_off = head.len() as u64;

    let path1 = format!("{out}.1.{ext}");
    let blk = {
        let mut w1 = BufWriter::with_capacity(1 << 20, File::create(&path1)?);
        let mut w2 = BufWriter::with_capacity(1 << 20, File::create(format!("{out}.2.{ext}"))?);
        w1.write_all(&head)?;
        w2.write_all(&1u32.to_le_bytes())?;
        let blk = wgemit::write_block(&rs, g, OFF_RATE, &mut w1, &mut w2)?;
        // zOffs and fchr fall out of the same walk
        let mut tail: Vec<u8> = Vec::new();
        put_idx(&mut tail, blk.z_offs.len() as u64);
        for &z in &blk.z_offs { put_idx(&mut tail, z); }
        for i in 0..5 { put_idx(&mut tail, blk.fchr[i]); }
        w1.write_all(&tail)?;
        w1.flush()?; w2.flush()?;
        blk
    };
    println!("  gbwt block: {} sides of {} bytes, {} SA samples, zOffs {:?}",
             g.num_sides, g.side_sz, blk.n_sa, blk.z_offs);
    let bad_fchr = (0..4).any(|i| blk.fchr[i + 1] - blk.fchr[i] != rs.bucket[i]);
    if bad_fchr { println!("  WARNING fchr disagrees with the label buckets"); }

    // ---- the ftab, queried against the block just written ---------------
    let (ftab, eftab) = {
        let f = File::open(&path1)?;
        let nav = wgemit::Nav::new(f, gbwt_off, g, rs.gbwt_len, blk.fchr,
                                   blk.z_offs.clone(), blk.f_rank_save);
        wgemit::build_ftab(&nav, FTAB_CHARS, w, true)?
    };

    {
        let mut w1 = OpenOptions::new().append(true).open(&path1)?;
        let mut tail: Vec<u8> = Vec::with_capacity((ftab.len() + eftab.len()) * 4);
        for &v in &ftab { put_idx(&mut tail, v); }
        for &v in &eftab { put_idx(&mut tail, v); }
        // refnames (gfm.h:2383): each name + '\n', then a single '\0'
        for n in &names { tail.extend_from_slice(n.as_bytes()); tail.push(b'\n'); }
        tail.push(0);
        w1.write_all(&tail)?;
        w1.flush()?;
    }
    // eftabLen is only known once the ftab is built, so the header carries a
    // placeholder until here -- which is exactly what buildToDisk does.
    // `out1.seekp(24 + sizeof(index_t) * 3)` (gfm.h:5111)
    OpenOptions::new().write(true).open(&path1)?
        .write_all_at(&(eftab.len() as u64).to_le_bytes()[0..w], eftab_len_at)?;

    println!("wrote {out}.1.{ext} ({} bytes), .2, .3, .4",
             fs::metadata(&path1)?.len());

    // ---- .5.ht2 / .6.ht2, one 57,344 bp window at a time ----------------
    // The reference is parsed a second time here rather than kept across the
    // doubling: holding the joined text plus the variant list through a
    // whole-genome build costs more than re-reading them costs.
    {
        let p = graph::parse(fa, snp, hap);
        let mut w5 = BufWriter::with_capacity(1 << 20, File::create(format!("{out}.5.{ext}"))?);
        let mut w6 = BufWriter::with_capacity(1 << 20, File::create(format!("{out}.6.{ext}"))?);
        local5::emit(&p, &mut w5, &mut w6, w, true)?;
        w5.flush()?; w6.flush()?;
    }
    println!("wrote {out}.5.{ext} ({} bytes) and .6.{ext} ({} bytes)",
             fs::metadata(format!("{out}.5.{ext}"))?.len(),
             fs::metadata(format!("{out}.6.{ext}"))?.len());

    if env::var("HT2_KEEP").is_err() {
        for f in ["nodes.bin", "edges.bin", "rows.bin", "nodeinfo.bin", "floc.bin"] {
            let _ = fs::remove_file(wd.join(f));
        }
    }

    // ---- verify ---------------------------------------------------------
    if let Some(v) = verify {
        println!("\nagainst {v}:");
        let mut ok = doubling::check_curve(&d.curve, &format!("{v}.log"));
        for n in 1..=6 {
            ok &= compare(&format!("{out}.{n}.{ext}"), &format!("{v}.{n}.{ext}"));
        }
        if !ok { process::exit(1); }
        println!("\nEXTERNAL BUILD MATCHES hisat2-build");
    }
    Ok(())
}
