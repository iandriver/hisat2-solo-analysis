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
#[path = "../ckpt.rs"]
mod ckpt;
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
use std::time::Instant;
use std::{env, process};

/// Per-stage wall time. The whole-genome question is which stage to thread
/// first, and that is not answerable from the source.
struct Timer { t0: Instant, last: Instant, marks: Vec<(String, f64)> }
impl Timer {
    fn new() -> Timer { let n = Instant::now(); Timer { t0: n, last: n, marks: Vec::new() } }
    fn mark(&mut self, what: &str) {
        let now = Instant::now();
        self.marks.push((what.to_string(), now.duration_since(self.last).as_secs_f64()));
        self.last = now;
    }
    fn report(&self) {
        let total = self.last.duration_since(self.t0).as_secs_f64();
        println!("\nstage timing (total {total:.1}s):");
        for (what, secs) in &self.marks {
            println!("  {:<22} {:>7.1}s  {:>5.1}%", what, secs, 100.0 * secs / total);
        }
    }
}

const OFF_RATE: u32 = 4;
const FTAB_CHARS: u32 = 10;
const VERSION: u32 = (2 << 24) | (2 << 16) | (3 << 8);   // 2.2.3

/// `writeI32` — always four bytes, whatever `index_t` is.
fn put(v: &mut Vec<u8>, x: u32) { v.extend_from_slice(&x.to_le_bytes()); }

/// `.7.ht2` / `.8.ht2`: the variant database (`gfm.h:1912`). Not part of the
/// graph -- a straight serialisation of the alt and haplotype lists, and where
/// `Zs:Z` gets its rsIDs. The repeat block that follows is empty unless
/// `--repeat-ref` was given, which is why the fixtures' `.7` ends at the
/// haplotypes.
fn alt_db(p: &graph::Parsed, w: usize) -> (Vec<u8>, Vec<u8>) {
    let put_idx = |v: &mut Vec<u8>, x: u64| v.extend_from_slice(&x.to_le_bytes()[0..w]);
    let mut o7: Vec<u8> = Vec::new();
    put(&mut o7, 1);
    put_idx(&mut o7, p.alts.len() as u64);
    for x in &p.alts {
        put_idx(&mut o7, x.pos as u64);
        // our own type codes are not HISAT2's ALT_TYPE enum
        put(&mut o7, match x.typ {
            graph::ALT_SGL  => 1,
            graph::ALT_INS  => 2,
            graph::ALT_DEL  => 3,
            graph::ALT_SS   => 5,
            graph::ALT_EXON => 6,
            other => panic!("unknown ALT type code {other}"),
        });
        put_idx(&mut o7, x.len as u64);
        o7.extend_from_slice(&x.seq.to_le_bytes());
    }
    put_idx(&mut o7, p.haps.len() as u64);
    for h in &p.haps {
        put_idx(&mut o7, h.left as u64);
        put_idx(&mut o7, h.right as u64);
        put_idx(&mut o7, h.alts.len() as u64);
        for &i in &h.alts { put_idx(&mut o7, i as u64); }
    }
    let mut o8: Vec<u8> = Vec::new();
    put(&mut o8, 1);
    put_idx(&mut o8, p.alts.len() as u64);
    o8.extend_from_slice(&p.alt_names);
    (o7, o8)
}

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
    // The annotation goes in by environment rather than argv so the five
    // positionals and every existing run script keep working. It is still an
    // input, so it is in the resume fingerprint.
    let ss = env::var("HT2_SS").unwrap_or_default();
    let exon = env::var("HT2_EXON").unwrap_or_default();
    // `local5::emit` does not know about splice sites yet, so an annotated run
    // would write local indexes built as though there were none while `.7` said
    // otherwise -- an index whose ALT table and whose local graphs disagree,
    // which nothing downstream would flag. HT2_NO_LOCAL skips `.5`/`.6`
    // entirely, which is what makes the rest checkable in the meantime.
    let no_local = env::var("HT2_NO_LOCAL").is_ok();
    let annotated = !ss.is_empty() || !exon.is_empty();
    if annotated && !no_local && env::var("HT2_ALTS_ONLY").is_err() {
        eprintln!("HT2_SS/HT2_EXON are not wired into the local indexes yet.");
        eprintln!("Run with HT2_NO_LOCAL=1 (skips .5/.6) or HT2_ALTS_ONLY=1, or unset them.");
        process::exit(2);
    }
    fs::create_dir_all(&wd)?;
    let mut timer = Timer::new();

    // ---- resume ---------------------------------------------------------
    // A whole-genome build runs for hours. It resumes by being re-run with the
    // same arguments; the checkpoint is trusted only when its fingerprint
    // matches, and that fingerprint covers this binary's own mtime and size, so
    // a rebuild invalidates it rather than silently finishing a build that two
    // different versions of the construction started. `HT2_FRESH=1` starts over.
    // Resume costs one extra copy of the node file on disk -- 82 bytes per path
    // node against 57 -- because each phase's input has to outlive the step that
    // replaces it. `HT2_NO_RESUME=1` buys that back and gives up the ability to
    // restart.
    let resume = env::var("HT2_NO_RESUME").is_err();
    let fp = ckpt::fingerprint(fa, snp, hap, &ss, &exon, large, chunk);
    let mut ck = if env::var("HT2_FRESH").is_ok() || !resume { None }
                 else { ckpt::load(&wd, &fp) }
        .unwrap_or_default();
    if !ck.stage.is_empty() {
        println!("resuming: {} complete{}",
                 match ck.stage.as_str() {
                     "graph" => "reference graph".to_string(),
                     "gen"   => format!("generation {}", ck.gen),
                     "edges" => "generateEdges".to_string(),
                     other   => other.to_string(),
                 },
                 if ck.phase.is_empty() { String::new() }
                 else { format!(", phase {}", ck.phase) });
    } else if ckpt::path(&wd).exists() {
        println!("checkpoint does not match these inputs or this binary -- starting over");
    }
    ck.fingerprint = fp;
    // Whatever a crash left that the checkpoint does not vouch for has to go. A
    // half-written run or merge partition does not announce itself: the next
    // file of the same name adopts its segments and the records turn up as
    // though they belonged.
    {
        let mut keep: Vec<String> = vec!["nodes.bin".into(), "edges.bin".into()];
        // Both node files, always. Whichever one a phase does not vouch for is
        // rewritten from scratch by the step that follows, so keeping a partial
        // one costs nothing -- and working out which is which per phase is the
        // sort of reasoning that gets one case wrong.
        keep.push("cur0.bin".into());
        keep.push("cur1.bin".into());
        match ck.phase.as_str() {
            "sorted" => { keep.push("by_to.bin".into()); keep.push("by_from.bin".into()); }
            "joined" => { keep.push("joined.bin".into()); }
            "keyed"  => { keep.push("sorted.bin".into()); }
            _ => {}
        }
        if ck.stage == "edges" {
            for f in ["rows.bin", "nodeinfo.bin", "floc.bin"] { keep.push(f.into()); }
        }
        let refs: Vec<&str> = keep.iter().map(|s| s.as_str()).collect();
        ckpt::clean(&wd, &refs);
    }

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
    timer.mark("reference, .3/.4");

    // HT2_REF_ONLY stops here. .3 and .4 are a pure function of the FASTA --
    // no variants, no graph -- so diffing them against a known-good index is a
    // minutes-long proof that the reference going in is the reference that
    // came out, before committing hours to the graph stages.
    if env::var("HT2_REF_ONLY").is_ok() {
        println!("\nHT2_REF_ONLY: stopping after the reference front end");
        println!("wrote {out}.3.{ext} ({} bytes), {out}.4.{ext} ({} bytes)",
                 fs::metadata(format!("{out}.3.{ext}"))?.len(),
                 fs::metadata(format!("{out}.4.{ext}"))?.len());
        timer.report();
        if let Some(v) = verify {
            println!("\nagainst {v}:");
            let mut ok = true;
            for n in [3u32, 4] {
                ok &= compare(&format!("{out}.{n}.{ext}"), &format!("{v}.{n}.{ext}"));
            }
            if !ok { process::exit(1); }
            println!("\nREFERENCE MATCHES ({v}.3/.4 byte-identical)");
        }
        return Ok(());
    }

    // HT2_ALTS_ONLY stops after the variant database. `.7` and `.8` are a pure
    // function of the FASTA and the variant files -- no graph, no doubling --
    // so this checks the parse, the ordering and the annotation ALTs against a
    // known-good index in seconds rather than at the end of a full build.
    if env::var("HT2_ALTS_ONLY").is_ok() {
        let p = graph::parse_with(fa, snp, hap, &ss, &exon);
        let n_ss = p.alts.iter().filter(|a| a.typ == graph::ALT_SS).count();
        let n_ex = p.alts.iter().filter(|a| a.typ == graph::ALT_EXON).count();
        let n_x  = p.alts.iter().filter(|a| a.typ == graph::ALT_SS && a.seq & (1 << 8) != 0).count();
        let (o7, o8) = alt_db(&p, w);
        fs::write(format!("{out}.7.{ext}"), &o7)?;
        fs::write(format!("{out}.8.{ext}"), &o8)?;
        println!("\nHT2_ALTS_ONLY: {} alts ({} splice sites, {} of them excluded; {} exons), \
{} haplotypes", p.alts.len(), n_ss, n_x, n_ex, p.haps.len());
        println!("wrote {out}.7.{ext} ({} bytes) and {out}.8.{ext} ({} bytes)", o7.len(), o8.len());
        if let Some(v) = verify {
            println!("\nagainst {v}:");
            let mut ok = true;
            for n in [7u32, 8] {
                ok &= compare(&format!("{out}.{n}.{ext}"), &format!("{v}.{n}.{ext}"));
            }
            if !ok { process::exit(1); }
            println!("\nVARIANT DATABASE MATCHES ({v}.7/.8 byte-identical)");
        }
        return Ok(());
    }

    // ---- the path graph, fragmented and external ------------------------
    let d = doubling::run_with(fa, snp, hap, &ss, &exon, &wd, budget, chunk, true, &mut ck, resume)?;
    timer.mark("graph + doubling");
    // generate_edges consumes `cur` -- it reads it once, to sort by `from` --
    // so its own checkpoint has to be the thing that lets a restart skip it
    let rs = if ck.stage == "edges" {
        println!("\ngenerateEdges: from checkpoint -- {} path nodes, {} path edges",
                 ck.e_nodes + 1, ck.e_gbwt);
        wgemit::Rows { rows: wd.join("rows.bin"), nodeinfo: wd.join("nodeinfo.bin"),
                       floc: wd.join("floc.bin"), n_nodes: ck.e_nodes,
                       gbwt_len: ck.e_gbwt, bucket: ck.e_bucket }
    } else {
        let rs = wgemit::generate_edges(&wd, &d.cur, budget, true, resume)?;
        ck.stage = "edges".into();
        ck.e_nodes = rs.n_nodes; ck.e_gbwt = rs.gbwt_len; ck.e_bucket = rs.bucket;
        // Checkpoint BEFORE dropping `cur`, not after. In between, the state on
        // disk says the doubling is done and the node file is still needed --
        // which is true. The other order says the node file is expendable
        // before anything has recorded what replaced it.
        if resume { ckpt::save(&wd, &ck)?; }
        ckpt::crash_point("edges", 0);
        ext::seg_remove(&d.cur);
        rs
    };
    timer.mark("generateEdges");
    // gfm.h:149 -- no variants means one node per base plus the terminator,
    // which HISAT2 encodes as a plain FM index rather than a graph.
    let linear = rs.gbwt_len == len as u64 + 1;
    let g = wgemit::geom(rs.gbwt_len, line_rate, w, linear);

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
    timer.mark("gbwt block + .2");
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

    timer.mark("ftab");
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

    // ---- .5 through .8: the local indexes and the variant database ------
    // The reference is parsed a second time here rather than kept across the
    // doubling: holding the joined text plus the variant list through a
    // whole-genome build costs more than re-reading them costs.
    {
        let p = graph::parse_with(fa, snp, hap, &ss, &exon);
        if !no_local {
            let mut w5 = BufWriter::with_capacity(1 << 20, File::create(format!("{out}.5.{ext}"))?);
            let mut w6 = BufWriter::with_capacity(1 << 20, File::create(format!("{out}.6.{ext}"))?);
            local5::emit(&p, &mut w5, &mut w6, w, true)?;
            w5.flush()?; w6.flush()?;
        }

        let (o7, o8) = alt_db(&p, w);
        fs::write(format!("{out}.7.{ext}"), &o7)?;
        fs::write(format!("{out}.8.{ext}"), &o8)?;
        println!("wrote {out}.7.{ext} ({} bytes) and .8.{ext} ({} bytes) -- {} variants, {} haplotypes",
                 o7.len(), o8.len(), p.alts.len(), p.haps.len());
    }
    timer.mark("local indexes .5/.6");
    if no_local {
        println!("HT2_NO_LOCAL: skipped {out}.5.{ext} and .6.{ext}");
    } else {
        println!("wrote {out}.5.{ext} ({} bytes) and .6.{ext} ({} bytes)",
                 fs::metadata(format!("{out}.5.{ext}"))?.len(),
                 fs::metadata(format!("{out}.6.{ext}"))?.len());
    }

    timer.report();

    if env::var("HT2_KEEP").is_err() {
        for f in ["nodes.bin", "edges.bin", "rows.bin", "nodeinfo.bin", "floc.bin"] {
            let _ = fs::remove_file(wd.join(f));
        }
        let _ = fs::remove_file(ckpt::path(&wd));
    }

    // ---- verify ---------------------------------------------------------
    if let Some(v) = verify {
        println!("\nagainst {v}:");
        let mut ok = doubling::check_curve(&d.curve, &format!("{v}.log"));
        for n in 1..=8 {
            if no_local && (n == 5 || n == 6) { continue; }
            ok &= compare(&format!("{out}.{n}.{ext}"), &format!("{v}.{n}.{ext}"));
        }
        if !ok { process::exit(1); }
        println!("\nEXTERNAL BUILD MATCHES hisat2-build");
    }
    Ok(())
}
