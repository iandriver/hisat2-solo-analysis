//! Parse `.5.ht2` / `.6.ht2` — the hierarchical local indexes — and account for
//! every byte.
//!
//! Same test with the same teeth as rung 1: every section size is *derived*, so
//! if any formula is wrong the sections will not sum to the file length.
//!
//! Layout, from `HGFM::readIntoMemory` (`hgfm.h:2680`) and
//! `LocalGFM::readIntoMemory` (`hgfm.h:1109`):
//!
//! ```text
//! .5:  u32 endianness sentinel
//!      index_t nLocalGFMs
//!      i32 lineRate, i32 unused, i32 offRate, i32 ftabChars, i32 flags
//!      then per local index:
//!        full_index_t tidx, localOffset, joinedOffset     (4 bytes each)
//!        local_index_t len, gbwtLen, numNodes, eftabLen   (2 bytes each)
//!        local_index_t nPat, plen[nPat]
//!        local_index_t nFrag, rstarts[nFrag*3]
//!        bytes gbwt[gbwtTotLen]
//!        local_index_t numZOffs, zOffs[]
//!        local_index_t fchr[5], ftab[ftabLen], eftab[eftabLen]
//! .6:  u32 endianness sentinel, then each local index's offs[] at local width
//! ```
//!
//! `local_index_t` is 2 bytes — `local_max_gbwt` is `(1<<16) - (1<<11)`, so a
//! local index is sized so its gbwt fits in 16 bits.

use std::convert::TryInto;
use std::{env, fs};

const LOCAL_FTAB_CHARS: usize = 6;

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 { eprintln!("usage: ht2local <index_prefix>"); std::process::exit(2); }
    let b5 = fs::read(format!("{}.5.ht2", a[1])).expect(".5.ht2");
    let b6 = fs::read(format!("{}.6.ht2", a[1])).expect(".6.ht2");
    let u32a = |b: &[u8], o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
    let u16a = |b: &[u8], o: usize| u16::from_le_bytes(b[o..o + 2].try_into().unwrap()) as usize;

    let mut p = 0usize;
    assert_eq!(u32a(&b5, p), 1, "endianness sentinel"); p += 4;
    let n_local = u32a(&b5, p) as usize; p += 4;
    let line_rate = u32a(&b5, p) as usize; p += 4;
    p += 4;
    let off_rate = u32a(&b5, p) as usize; p += 4;
    let ftab_chars = u32a(&b5, p) as usize; p += 4;
    let _flags = u32a(&b5, p); p += 4;
    println!("{}.5.ht2: {} local indexes, lineRate {}, offRate {}, ftabChars {}",
             a[1], n_local, line_rate, off_rate, ftab_chars);
    assert_eq!(ftab_chars, LOCAL_FTAB_CHARS);

    let ftab_len = (1usize << (2 * ftab_chars)) + 1;
    let side_sz = 1usize << line_rate;
    let mut p6 = 4usize;                       // .6 sentinel
    let (mut tot_len, mut tot_gbwt, mut graphs, mut linears) = (0usize, 0usize, 0usize, 0usize);

    for i in 0..n_local {
        let _tidx = u32a(&b5, p); p += 4;
        let _loff = u32a(&b5, p); p += 4;
        let _joff = u32a(&b5, p); p += 4;
        let len = u16a(&b5, p); p += 2;
        let gbwt_len = u16a(&b5, p); p += 2;
        let num_nodes = u16a(&b5, p); p += 2;
        let eftab_len = u16a(&b5, p); p += 2;
        if len == 0 { continue; }               // readIntoMemory returns early

        // GFMParams::init at local width
        let linear = len + 1 == gbwt_len || gbwt_len == 0;
        let gl = if gbwt_len == 0 { len + 1 } else { gbwt_len };
        let nn = if num_nodes == 0 { len + 1 } else { num_nodes };
        let gbwt_sz = if linear { gl / 4 + 1 } else { gl / 2 + 1 };
        let side_gbwt_sz = side_sz - 2 * if linear { 4 } else { 6 };
        let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
        let gbwt_tot = num_sides * side_sz;
        let offs_len = (nn + (1 << off_rate) - 1) >> off_rate;
        if linear { linears += 1; } else { graphs += 1; }
        tot_len += len; tot_gbwt += gl;

        let p_pat = p;
        let n_pat = u16a(&b5, p); p += 2 + n_pat * 2;
        let p_frag = p;
        let n_frag = u16a(&b5, p); p += 2 + n_frag * 6;
        p += gbwt_tot;
        let p_z = p;
        let n_z = u16a(&b5, p); p += 2 + n_z * 2;
        let p_fchr = p;
        p += 5 * 2;                              // fchr
        p += ftab_len * 2;
        p += eftab_len * 2;
        p6 += offs_len * 2;
        if i < 64 {
            println!("  local {i}: tidx {_tidx} localOff {_loff} joinedOff {_joff} len {len} \
gbwtLen {gl} numNodes {nn} eftabLen {eftab_len} nPat {n_pat} nFrag {n_frag} numZOffs {n_z} \
sides {num_sides} offsLen {offs_len} {}", if linear { "linear" } else { "graph" });
            if std::env::var("HT2_LDBG").is_ok() {
                let pl = u16a(&b5, p_pat + 2);
                let rs: Vec<usize> = (0..n_frag * 3).map(|k| u16a(&b5, p_frag + 2 + k * 2)).collect();
                let fc: Vec<usize> = (0..5).map(|k| u16a(&b5, p_fchr + k * 2)).collect();
                let zo: Vec<usize> = (0..n_z).map(|k| u16a(&b5, p_z + 2 + k * 2)).collect();
                println!("        plen[0] {pl}  rstarts {rs:?}  fchr {fc:?}  zOffs {zo:?}");
            }
        }
    }

    // a single trailing NUL closes .5 (`fout5 << '\0'`, hgfm.h:2519), the same
    // terminator refnames gets in .1.ht2
    p += 1;
    println!("  {graphs} graph + {linears} linear local indexes; total len {tot_len}, total gbwt {tot_gbwt}");
    println!("  .5 consumed {p} of {} bytes", b5.len());
    println!("  .6 consumed {p6} of {} bytes", b6.len());
    if p == b5.len() && p6 == b6.len() {
        println!("\nLOCAL INDEX LAYOUT OK: every byte of .5 and .6 accounted for");
    } else {
        println!("\nLAYOUT MISMATCH: .5 off by {}, .6 off by {}",
                 b5.len() as i64 - p as i64, b6.len() as i64 - p6 as i64);
        std::process::exit(1);
    }
}
