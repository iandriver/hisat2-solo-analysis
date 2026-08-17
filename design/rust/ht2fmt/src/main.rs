//! S1 — read a HISAT2 `.ht2` header and re-emit it byte-identically.
//!
//! Rung 1 of the validation ladder in `rust_builder_plan.md`. No construction,
//! no cleverness: the only claim being tested is that we know the on-disk
//! layout well enough to reproduce it exactly. Everything later depends on
//! that, and it is the cheapest possible place to find out we are wrong.
//!
//! Field order is taken from `GFM::readIntoMemory` (`gfm.h:5905` onward), which
//! is the authoritative reader:
//!
//! ```text
//!   u32     endianness sentinel (1, or 1<<24 if the file is big-endian)
//!   u32     index version
//!   index_t len
//!   index_t gbwtLen
//!   index_t numNodes
//!   i32     lineRate
//!   i32     linesPerSide
//!   i32     offRate
//!   i32     ftabChars
//!   index_t eftabLen
//!   i32     flags
//! ```
//!
//! `index_t` is 4 bytes for `.ht2` and 8 for `.ht2l`; the width is carried by
//! the file extension, not by the header (`gfm.cpp:27`).

use std::env;
use std::fs::File;
use std::io::{Read, Write};
use std::process::exit;

const GFM_ENTIRE_REV: i32 = 4;

struct Reader {
    buf: Vec<u8>,
    pos: usize,
    swap: bool,
    width: usize,
}

impl Reader {
    fn u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        b.copy_from_slice(&self.buf[self.pos..self.pos + 4]);
        self.pos += 4;
        if self.swap { u32::from_be_bytes(b) } else { u32::from_le_bytes(b) }
    }
    fn i32(&mut self) -> i32 {
        self.u32() as i32
    }
    /// index_t: width depends on the file, not the header.
    fn idx(&mut self) -> u64 {
        if self.width == 4 {
            self.u32() as u64
        } else {
            let mut b = [0u8; 8];
            b.copy_from_slice(&self.buf[self.pos..self.pos + 8]);
            self.pos += 8;
            if self.swap { u64::from_be_bytes(b) } else { u64::from_le_bytes(b) }
        }
    }
}

#[derive(Debug)]
struct Header {
    sentinel: u32,
    index_version: u32,
    len: u64,
    gbwt_len: u64,
    num_nodes: u64,
    line_rate: i32,
    lines_per_side: i32,
    off_rate: i32,
    ftab_chars: i32,
    eftab_len: u64,
    flags: i32,
    width: usize,
    swap: bool,
    header_bytes: usize,
}

impl Header {
    fn parse(buf: Vec<u8>, width: usize) -> Result<Header, String> {
        if buf.len() < 64 {
            return Err(format!("file too short: {} bytes", buf.len()));
        }
        // The sentinel decides endianness: written as 1, so reading 1<<24
        // means the file was produced on the other endianness.
        let raw = u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
        let swap = match raw {
            1 => false,
            x if x == 1u32 << 24 => true,
            other => return Err(format!("bad endianness sentinel: 0x{other:08x}")),
        };
        let mut r = Reader { buf, pos: 0, swap, width };
        let h = Header {
            sentinel: r.u32(),
            index_version: r.u32(),
            len: r.idx(),
            gbwt_len: r.idx(),
            num_nodes: r.idx(),
            line_rate: r.i32(),
            lines_per_side: r.i32(),
            off_rate: r.i32(),
            ftab_chars: r.i32(),
            eftab_len: r.idx(),
            flags: r.i32(),
            width,
            swap,
            header_bytes: 0,
        };
        let n = r.pos;
        Ok(Header { header_bytes: n, ..h })
    }

    /// Mirror of `GFM::readIndexVersion` (`gfm.h:2815`).
    fn version_string(&self) -> String {
        let major = (self.index_version >> 16) & 0xff;
        let minor = (self.index_version >> 8) & 0xff;
        let extra = match self.index_version & 0xff {
            1 => "-alpha",
            2 => "-beta",
            _ => "",
        };
        format!("2.{major}.{minor}{extra}")
    }

    fn entire_reverse(&self) -> bool {
        // gfm.h: flags is negative and the ENTIRE_REV bit is tested on -flags
        self.flags < 0 && ((-self.flags) & GFM_ENTIRE_REV) != 0
    }

    /// Re-emit exactly the bytes we claim to understand.
    fn emit(&self) -> Vec<u8> {
        let mut o: Vec<u8> = Vec::with_capacity(self.header_bytes);
        let put_u32 = |o: &mut Vec<u8>, v: u32, swap: bool| {
            o.extend_from_slice(&if swap { v.to_be_bytes() } else { v.to_le_bytes() });
        };
        let put_idx = |o: &mut Vec<u8>, v: u64, swap: bool, width: usize| {
            if width == 4 {
                let v = v as u32;
                o.extend_from_slice(&if swap { v.to_be_bytes() } else { v.to_le_bytes() });
            } else {
                o.extend_from_slice(&if swap { v.to_be_bytes() } else { v.to_le_bytes() });
            }
        };
        put_u32(&mut o, self.sentinel, self.swap);
        put_u32(&mut o, self.index_version, self.swap);
        put_idx(&mut o, self.len, self.swap, self.width);
        put_idx(&mut o, self.gbwt_len, self.swap, self.width);
        put_idx(&mut o, self.num_nodes, self.swap, self.width);
        put_u32(&mut o, self.line_rate as u32, self.swap);
        put_u32(&mut o, self.lines_per_side as u32, self.swap);
        put_u32(&mut o, self.off_rate as u32, self.swap);
        put_u32(&mut o, self.ftab_chars as u32, self.swap);
        put_idx(&mut o, self.eftab_len, self.swap, self.width);
        put_u32(&mut o, self.flags as u32, self.swap);
        o
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: ht2fmt <index.1.ht2|index.1.ht2l>");
        exit(2);
    }
    let path = &args[1];
    let width = if path.ends_with(".ht2l") { 8 } else { 4 };

    let mut buf = Vec::new();
    match File::open(path).and_then(|mut f| f.read_to_end(&mut buf)) {
        Ok(_) => {}
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            exit(1);
        }
    }
    let file_len = buf.len();
    let original = buf.clone();

    let h = match Header::parse(buf, width) {
        Ok(h) => h,
        Err(e) => {
            eprintln!("parse failed: {e}");
            exit(1);
        }
    };

    println!("{path}");
    println!("  index_t width      {} bytes ({})", h.width,
             if h.width == 4 { "small, .ht2" } else { "large, .ht2l" });
    println!("  endianness         {}", if h.swap { "big (swapped)" } else { "little" });
    println!("  index version      {} (raw 0x{:08x})", h.version_string(), h.index_version);
    println!("  len                {}", h.len);
    println!("  gbwtLen            {}", h.gbwt_len);
    println!("  numNodes           {}", h.num_nodes);
    println!("  lineRate           {}", h.line_rate);
    println!("  linesPerSide       {}", h.lines_per_side);
    println!("  offRate            {}", h.off_rate);
    println!("  ftabChars          {}", h.ftab_chars);
    println!("  eftabLen           {}", h.eftab_len);
    println!("  flags              {} (entireReverse={})", h.flags, h.entire_reverse());
    println!("  header bytes       {}", h.header_bytes);
    println!("  file bytes         {file_len}");

    // The actual test: our bytes must equal the file's bytes.
    let emitted = h.emit();
    if emitted.len() != h.header_bytes {
        println!("\nROUND-TRIP FAIL: emitted {} bytes, parsed {}", emitted.len(), h.header_bytes);
        exit(1);
    }
    if emitted[..] != original[..h.header_bytes] {
        println!("\nROUND-TRIP FAIL: bytes differ");
        for i in 0..h.header_bytes {
            if emitted[i] != original[i] {
                println!("  first difference at offset {i}: emitted 0x{:02x}, file 0x{:02x}",
                         emitted[i], original[i]);
                break;
            }
        }
        exit(1);
    }
    println!("\nROUND-TRIP OK: {} header bytes reproduced exactly", h.header_bytes);

    if let Ok(mut f) = File::create("/dev/null") {
        let _ = f.write_all(&emitted);
    }
}
