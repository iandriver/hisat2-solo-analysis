//! Invert the BWT stored in a linear .1.ht2 and compare with the reference.
//!
//! This settles whether a mismatch lives in our unpacking or in our suffix
//! array, without relying on either. If the stored BWT, fchr and zOffs
//! reconstruct the text, the unpacking is right by construction.
//!
//! Layout per GFM::buildToDisk (gfm.h:5148): rows are laid out in order, 4 per
//! byte, low bit-pair first (pack_2b_in_8b: eight |= two << (off*2)), filling
//! the first sideGbwtSz bytes of each sideSz side; the four per-side occ
//! counts occupy the final 4*sizeof(index_t) bytes and hold the tallies as of
//! the START of that side (occSave, not occ).
use std::env;
use std::fs;

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 3 { eprintln!("usage: ht2invert <index.1.ht2> <reference.fa>"); std::process::exit(2); }
    let buf = fs::read(&a[1]).expect("read index");
    let u = |o: usize| -> u64 { u32::from_le_bytes(buf[o..o+4].try_into().unwrap()) as u64 };

    let len = u(8); let gbwt_len = u(12);
    let line_rate = u(20) as u32; let n_pat = u(44);
    let n_frag = u(44 + 4 + n_pat as usize * 4);
    let side_sz = 1u64 << line_rate;
    let side_gbwt_sz = side_sz - 4 * 4;          // linear FM reserves 4 index_t
    let side_gbwt_len = side_gbwt_sz * 4;
    let gbwt_sz = gbwt_len / 4 + 1;
    let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
    let gbwt_off = 44 + 4 + n_pat as usize * 4 + 4 + n_frag as usize * 12;
    let zoff_at = gbwt_off + (num_sides * side_sz) as usize;
    let z = u(zoff_at + 4);
    let fchr_at = zoff_at + 4 + 4;
    let fchr: Vec<u64> = (0..5).map(|i| u(fchr_at + i * 4)).collect();
    eprintln!("len={len} gbwtLen={gbwt_len} zOffs={z} fchr={fchr:?}");

    // unpack every row
    let mut bwt = vec![0u8; gbwt_len as usize];
    for r in 0..gbwt_len as usize {
        let side = r as u64 / side_gbwt_len;
        let off = r as u64 % side_gbwt_len;
        let byte = buf[gbwt_off + (side * side_sz + (off >> 2)) as usize];
        bwt[r] = ((byte >> (2 * (off & 3))) & 3) as u8;
    }
    // occ[c][r] = occurrences of c in bwt[0..r), skipping the zOffs row
    let mut occ = vec![[0u32; 4]; gbwt_len as usize + 1];
    for r in 0..gbwt_len as usize {
        occ[r + 1] = occ[r];
        if r as u64 != z { occ[r + 1][bwt[r] as usize] += 1; }
    }
    // The empty suffix occupies a row whose F entry is '$', so the character
    // blocks may start one row later. Try both.
    let lf = |r: u64, bias: u64| -> u64 {
        let c = bwt[r as usize] as usize;
        bias + fchr[c] + occ[r as usize][c] as u64
    };

    // Decisive check of the UNPACKING alone: each side stores the occ tallies
    // as of its start (occSave). Recount them from our unpacked rows; padding
    // rows count as 'A', the zOffs row does not count at all.
    println!("\nper-side occ check (unpacking only, no LF involved):");
    let mut agree = 0usize; let mut shown = 0;
    for sd in 0..num_sides as usize {
        let base = gbwt_off + sd * side_sz as usize + side_gbwt_sz as usize;
        let stored: Vec<u64> = (0..4).map(|i| u(base + i * 4)).collect();
        let row0 = sd as u64 * side_gbwt_len;
        let mine: Vec<u64> = (0..4).map(|c| occ[row0.min(gbwt_len) as usize][c] as u64).collect();
        if stored == mine { agree += 1; }
        else if shown < 4 {
            println!("  side {sd:>6} rows[0,{row0}): stored {stored:?}  ours {mine:?}");
            shown += 1;
        }
    }
    println!("  {agree}/{} sides agree", num_sides);

    // reference text
    let mut t = Vec::new();
    for ln in fs::read_to_string(&a[2]).expect("read fa").lines() {
        if ln.starts_with('>') { continue; }
        for ch in ln.bytes() {
            match ch.to_ascii_uppercase() { b'A'=>t.push(0u8), b'C'=>t.push(1), b'G'=>t.push(2), b'T'=>t.push(3), _=>{} }
        }
    }
    eprintln!("reference text: {} bp", t.len());

    // Walk from the row of the empty suffix, emitting the text backwards.
    // Try every plausible starting row rather than assuming one.
    for &(start, bias) in &[(0u64,0u64),(0,1),(z,0),(z,1)] {
        let mut r = start;
        let mut out = Vec::with_capacity(t.len());
        for _ in 0..t.len() { out.push(bwt[r as usize]); r = lf(r, bias); }
        out.reverse();
        let ok = out == t;
        let mut first = usize::MAX;
        if !ok { for i in 0..t.len().min(out.len()) { if out[i] != t[i] { first = i; break; } } }
        println!("start {start:>9} bias {bias}: reconstruction {}  (final row {r}){}",
                 if ok { "MATCHES the reference" } else { "differs" },
                 if ok { String::new() } else { format!(", first difference at {first}") });
    }
}

// (appended) verify_occ is invoked from main via a second pass below.
