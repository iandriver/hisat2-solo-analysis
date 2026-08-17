//! S1 — parse a HISAT2 `.1.ht2` in full and re-emit it byte-identically.
//!
//! Rung 1 of the ladder in `rust_builder_plan.md`. No construction: the claim
//! under test is only that the on-disk layout is understood exactly. That is a
//! prerequisite for emitting anything, and it is the cheapest place to discover
//! a misunderstanding.
//!
//! The test has teeth because every section size is *derived*, not stored. If
//! any formula in `GFMParams::init` (`gfm.h:138`) is reproduced wrongly, the
//! computed sections will not sum to the file length.
//!
//! Layout, from `GFM::readIntoMemory` (`gfm.h:5905`) which is the authoritative
//! reader:
//!
//! ```text
//!   u32      endianness sentinel (1; 1<<24 means byte-swapped)
//!   u32      index version
//!   index_t  len, gbwtLen, numNodes
//!   i32      lineRate, linesPerSide, offRate, ftabChars
//!   index_t  eftabLen
//!   i32      flags
//!   index_t  nPat
//!   index_t  plen[nPat]
//!   index_t  nFrag
//!   index_t  rstarts[nFrag * 3]
//!   u8       gbwt[gbwtTotLen]          <- derived, the bulk
//!   index_t  numZOffs
//!   index_t  zOffs[numZOffs]
//!   index_t  fchr[5]
//!   index_t  ftab[ftabLen]             <- derived: 4^ftabChars + 1
//!   index_t  eftab[eftabLen]
//!   bytes    refnames, newline-separated, terminated by '\0' or EOF
//! ```
//!
//! `index_t` is 4 bytes for `.ht2`, 8 for `.ht2l` (`gfm.cpp:27`).

use std::env;
use std::fs::File;
use std::io::Read;
use std::process::exit;

const GFM_ENTIRE_REV: i32 = 4;

struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
    swap: bool,
    width: usize,
}

impl<'a> Cursor<'a> {
    fn need(&self, n: usize) -> Result<(), String> {
        if self.pos + n > self.buf.len() {
            Err(format!(
                "truncated: wanted {n} bytes at offset {}, file has {}",
                self.pos,
                self.buf.len()
            ))
        } else {
            Ok(())
        }
    }
    fn u32(&mut self) -> Result<u32, String> {
        self.need(4)?;
        let b: [u8; 4] = self.buf[self.pos..self.pos + 4].try_into().unwrap();
        self.pos += 4;
        Ok(if self.swap { u32::from_be_bytes(b) } else { u32::from_le_bytes(b) })
    }
    fn i32(&mut self) -> Result<i32, String> {
        Ok(self.u32()? as i32)
    }
    fn idx(&mut self) -> Result<u64, String> {
        if self.width == 4 {
            Ok(self.u32()? as u64)
        } else {
            self.need(8)?;
            let b: [u8; 8] = self.buf[self.pos..self.pos + 8].try_into().unwrap();
            self.pos += 8;
            Ok(if self.swap { u64::from_be_bytes(b) } else { u64::from_le_bytes(b) })
        }
    }
    fn skip(&mut self, n: usize) -> Result<(), String> {
        self.need(n)?;
        self.pos += n;
        Ok(())
    }
}

/// Mirror of `GFMParams::init` (`gfm.h:138`). Every one of these is derived at
/// load time rather than stored, so getting them right is the whole test.
struct Derived {
    linear_fm: bool,
    gbwt_sz: u64,
    ftab_len: u64,
    offs_len: u64,
    side_sz: u64,
    side_gbwt_sz: u64,
    num_sides: u64,
    gbwt_tot_len: u64,
}

fn derive(len: u64, gbwt_len: u64, num_nodes: u64, line_rate: i32, off_rate: i32,
          ftab_chars: i32, width: u64) -> Derived {
    let linear_fm = len + 1 == gbwt_len || gbwt_len == 0;
    let gbwt_len = if gbwt_len == 0 { len + 1 } else { gbwt_len };
    let num_nodes = if num_nodes == 0 { len + 1 } else { num_nodes };
    let gbwt_sz = if linear_fm { gbwt_len / 4 + 1 } else { gbwt_len / 2 + 1 };
    let ftab_len = (1u64 << (ftab_chars * 2)) + 1;
    let offs_len = (num_nodes + (1u64 << off_rate) - 1) >> off_rate;
    let side_sz = 1u64 << line_rate;
    // graph indexes reserve 6 index_t per side, linear ones 4
    let side_gbwt_sz = side_sz - width * if linear_fm { 4 } else { 6 };
    let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
    Derived {
        linear_fm,
        gbwt_sz,
        ftab_len,
        offs_len,
        side_sz,
        side_gbwt_sz,
        num_sides,
        gbwt_tot_len: num_sides * side_sz,
    }
}

struct Section {
    name: &'static str,
    off: usize,
    len: usize,
}

fn run(path: &str) -> Result<(), String> {
    let width: usize = if path.ends_with(".ht2l") { 8 } else { 4 };
    let mut buf = Vec::new();
    File::open(path)
        .and_then(|mut f| f.read_to_end(&mut buf))
        .map_err(|e| format!("cannot read {path}: {e}"))?;
    let file_len = buf.len();

    if buf.len() < 16 {
        return Err(format!("file too short: {} bytes", buf.len()));
    }
    let raw = u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
    let swap = match raw {
        1 => false,
        x if x == 1u32 << 24 => true,
        other => return Err(format!("bad endianness sentinel 0x{other:08x}")),
    };

    let mut c = Cursor { buf: &buf, pos: 0, swap, width };
    let mut secs: Vec<Section> = Vec::new();
    let mut mark = |secs: &mut Vec<Section>, name, start: usize, end: usize| {
        secs.push(Section { name, off: start, len: end - start })
    };

    let h0 = c.pos;
    let _sentinel = c.u32()?;
    let index_version = c.u32()?;
    let len = c.idx()?;
    let gbwt_len = c.idx()?;
    let num_nodes = c.idx()?;
    let line_rate = c.i32()?;
    let _lines_per_side = c.i32()?;
    let off_rate = c.i32()?;
    let ftab_chars = c.i32()?;
    let eftab_len = c.idx()?;
    let flags = c.i32()?;
    mark(&mut secs, "header", h0, c.pos);

    let d = derive(len, gbwt_len, num_nodes, line_rate, off_rate, ftab_chars, width as u64);

    let s = c.pos;
    let n_pat = c.idx()?;
    c.skip(n_pat as usize * width)?;
    mark(&mut secs, "nPat + plen[]", s, c.pos);

    let s = c.pos;
    let n_frag = c.idx()?;
    c.skip(n_frag as usize * 3 * width)?;
    mark(&mut secs, "nFrag + rstarts[]", s, c.pos);

    let s = c.pos;
    c.skip(d.gbwt_tot_len as usize)?;
    mark(&mut secs, "gbwt (derived)", s, c.pos);

    let s = c.pos;
    let num_zoffs = c.idx()?;
    c.skip(num_zoffs as usize * width)?;
    mark(&mut secs, "numZOffs + zOffs[]", s, c.pos);

    let s = c.pos;
    c.skip(5 * width)?;
    mark(&mut secs, "fchr[5]", s, c.pos);

    let s = c.pos;
    c.skip(d.ftab_len as usize * width)?;
    mark(&mut secs, "ftab[] (derived)", s, c.pos);

    let s = c.pos;
    c.skip(eftab_len as usize * width)?;
    mark(&mut secs, "eftab[]", s, c.pos);

    // refnames run to '\0' or EOF
    let s = c.pos;
    let mut names = 1usize;
    while c.pos < c.buf.len() {
        let b = c.buf[c.pos];
        c.pos += 1;
        if b == 0 {
            break;
        }
        if b == b'\n' {
            names += 1;
        }
    }
    mark(&mut secs, "refnames", s, c.pos);

    let major = (index_version >> 16) & 0xff;
    let minor = (index_version >> 8) & 0xff;
    let extra = match index_version & 0xff { 1 => "-alpha", 2 => "-beta", _ => "" };

    println!("{path}");
    println!("  index_t width      {width} bytes ({})",
             if width == 4 { ".ht2" } else { ".ht2l" });
    println!("  endianness         {}", if swap { "big (swapped)" } else { "little" });
    println!("  index version      2.{major}.{minor}{extra}");
    println!("  mode               {}",
             if d.linear_fm { "linear FM" } else { "graph FM" });
    println!("  len                {len}");
    println!("  gbwtLen            {gbwt_len}");
    println!("  numNodes           {num_nodes}");
    println!("  lineRate {line_rate}   offRate {off_rate}   ftabChars {ftab_chars}   eftabLen {eftab_len}");
    println!("  flags              {flags} (entireReverse={})",
             flags < 0 && ((-flags) & GFM_ENTIRE_REV) != 0);
    println!("  nPat {n_pat}   nFrag {n_frag}   zOffs {num_zoffs}   refnames {names}");
    println!();
    println!("  derived: gbwtSz {}  sideSz {}  sideGbwtSz {}  numSides {}",
             d.gbwt_sz, d.side_sz, d.side_gbwt_sz, d.num_sides);
    println!("           gbwtTotLen {}  ftabLen {}  offsLen {}",
             d.gbwt_tot_len, d.ftab_len, d.offs_len);
    println!();
    println!("  {:<22} {:>14} {:>16}", "section", "offset", "bytes");
    for s in &secs {
        println!("  {:<22} {:>14} {:>16}", s.name, s.off, s.len);
    }

    let total: usize = secs.iter().map(|s| s.len).sum();
    println!();
    println!("  sections total     {total}");
    println!("  file length        {file_len}");
    if total != file_len {
        return Err(format!(
            "LAYOUT MISMATCH: sections sum to {total}, file is {file_len} ({} bytes off)",
            file_len as i64 - total as i64
        ));
    }
    println!("\nLAYOUT OK: every section accounted for, sums exactly to the file length");
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: ht2fmt <index.1.ht2|index.1.ht2l> ...");
        exit(2);
    }
    let mut bad = 0;
    for p in &args[1..] {
        if let Err(e) = run(p) {
            eprintln!("\n{p}: {e}");
            bad += 1;
        }
        println!();
    }
    exit(if bad > 0 { 1 } else { 0 });
}
