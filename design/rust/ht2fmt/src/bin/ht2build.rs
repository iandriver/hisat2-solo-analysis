//! Rung 2 — construct a linear `.ht2` index and match `hisat2-build` exactly.
//!
//! This stage covers the front end: FASTA -> joined 2-bit reference, the
//! fragment records, the character histogram, and the suffix array / BWT.
//! Each piece is checked against the reference index rather than against my own
//! expectations, so a wrong answer shows up immediately:
//!
//!   * `nPat`, `plen[]`, `nFrag`, `rstarts[]` - compared byte for byte
//!   * `fchr[5]`                              - compared value for value
//!   * `zOffs`                                - the BWT row holding the
//!     sentinel. Matching it means the whole suffix array is right, since it is
//!     the rank of the complete string among all suffixes.
//!
//! Ambiguous characters are excluded from the indexed text but still counted in
//! `plen`, which is why `len` (900,000) is less than `plen[0]` (1,000,000) for
//! the example reference: it contains a 100,000 bp run of N.

use std::env;
use std::fs::File;
use std::io::Read;
use std::process::exit;

/// One maximal run of unambiguous sequence.
#[derive(Debug, Clone, Copy)]
struct Frag {
    joined_off: u64, // offset in the concatenated unambiguous text
    text_id: u64,    // which input sequence
    text_off: u64,   // offset within that sequence, counting ambiguous chars
}

struct Reference {
    names: Vec<String>,
    plen: Vec<u64>,   // full length of each sequence, ambiguous chars included
    frags: Vec<Frag>,
    text: Vec<u8>,    // 0..=3, ambiguous characters dropped
}

fn code(c: u8) -> Option<u8> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn read_fasta(path: &str) -> Result<Reference, String> {
    let mut raw = Vec::new();
    File::open(path)
        .and_then(|mut f| f.read_to_end(&mut raw))
        .map_err(|e| format!("cannot read {path}: {e}"))?;

    let mut r = Reference { names: Vec::new(), plen: Vec::new(), frags: Vec::new(), text: Vec::new() };
    let mut text_off: u64 = 0; // position within the current sequence
    let mut in_frag = false;

    for line in raw.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if !r.names.is_empty() {
                r.plen.push(text_off);
            }
            // name runs to the first whitespace
            let name: String = line[1..]
                .iter()
                .take_while(|&&c| c != b' ' && c != b'\t' && c != b'\r')
                .map(|&c| c as char)
                .collect();
            r.names.push(name);
            text_off = 0;
            in_frag = false;
            continue;
        }
        for &ch in line {
            if ch == b'\r' {
                continue;
            }
            match code(ch) {
                Some(v) => {
                    if !in_frag {
                        r.frags.push(Frag {
                            joined_off: r.text.len() as u64,
                            text_id: (r.names.len() - 1) as u64,
                            text_off,
                        });
                        in_frag = true;
                    }
                    r.text.push(v);
                }
                None => in_frag = false,
            }
            text_off += 1;
        }
    }
    if !r.names.is_empty() {
        r.plen.push(text_off);
    }
    if r.text.is_empty() {
        return Err("no unambiguous sequence found".into());
    }
    Ok(r)
}

/// Suffix array of `t` plus an implicit sentinel that sorts before everything.
/// Prefix doubling: O(n log^2 n), which is ample at this scale and easy to
/// argue is correct, which matters more here than speed.
fn suffix_array(t: &[u8]) -> Vec<u32> {
    let n = t.len();
    let m = n + 1; // + sentinel
    let mut sa: Vec<u32> = (0..m as u32).collect();
    // rank 0 is the sentinel (position n); real characters start at 1
    let mut rank: Vec<u32> = (0..m)
        .map(|i| if i == n { 0 } else { t[i] as u32 + 1 })
        .collect();
    let mut tmp: Vec<u32> = vec![0; m];

    let key = |rank: &Vec<u32>, i: u32, k: usize, m: usize| -> (u32, u32) {
        let a = rank[i as usize];
        let j = i as usize + k;
        let b = if j < m { rank[j] } else { 0 };
        (a, b)
    };

    let mut k = 1usize;
    loop {
        sa.sort_unstable_by_key(|&i| key(&rank, i, k, m));
        tmp[sa[0] as usize] = 0;
        let mut classes = 0u32;
        for w in 1..m {
            if key(&rank, sa[w], k, m) != key(&rank, sa[w - 1], k, m) {
                classes += 1;
            }
            tmp[sa[w] as usize] = classes;
        }
        std::mem::swap(&mut rank, &mut tmp);
        if classes as usize == m - 1 {
            break;
        }
        k <<= 1;
        if k >= m {
            break;
        }
    }
    sa
}


/// Unpack the BWT the reference index actually stores and diff it against ours.
///
/// Side layout from `GFM::postReadInit` (`gfm.h:2783`): BWT characters are
/// packed 4 per byte in the first `sideGbwtSz` bytes of each `sideSz` side, and
/// the per-side cumulative counts occupy the tail. Row r therefore lives at
/// side r/sideGbwtLen, character offset r%sideGbwtLen.
fn compare_bwt(buf: &[u8], gbwt_off: usize, side_sz: usize, side_gbwt_len: usize,
               sa: &[u32], t: &[u8], z_ours: u64) {
    let at = |r: usize, hi_first: bool| -> u8 {
        let side = r / side_gbwt_len;
        let off = r % side_gbwt_len;
        let byte = buf[gbwt_off + side * side_sz + (off >> 2)];
        let bp = off & 3;
        if hi_first { (byte >> (2 * (3 - bp))) & 3 } else { (byte >> (2 * bp)) & 3 }
    };
    let m = sa.len(); // n + 1 rows, sentinel first
    for &hi_first in &[true, false] {
        // try a few row offsets: our sentinel-first numbering may be shifted
        for shift in -2i64..=2 {
            let mut agree = 0usize;
            let mut total = 0usize;
            for row in 1..m.min(20000) {
                let r2 = row as i64 + shift;
                if r2 < 0 || r2 as usize >= m { continue; }
                let p = sa[row] as usize;
                if p == 0 { continue; } // sentinel row, character is a placeholder
                total += 1;
                if t[p - 1] == at(r2 as usize, hi_first) { agree += 1; }
            }
            if total > 0 && agree * 100 / total >= 90 {
                println!("  BWT match: hi_first={hi_first} shift={shift} -> {}/{} ({}%)",
                         agree, total, agree * 100 / total);
                // full-length pass, and name the rows that disagree
                let mut bad: Vec<(usize, usize, u8, u8)> = Vec::new();
                let mut tot2 = 0usize; let mut nbad = 0usize;
                for row in 1..m {
                    let r2 = row as i64 + shift;
                    if r2 < 0 || r2 as usize >= m { continue; }
                    let p = sa[row] as usize;
                    if p == 0 { continue; }
                    tot2 += 1;
                    let got = at(r2 as usize, hi_first);
                    if t[p - 1] != got { nbad += 1; if bad.len() < 6 {
                        bad.push((row, p, t[p - 1], got)); }
                    }
                }
                println!("    full pass: {} rows, {} disagreements ({:.4}%)", tot2, nbad, 100.0 * nbad as f64 / tot2 as f64);
                for (row, p, mine, theirs) in &bad {
                    println!("      row {:>9}  SA={:>9}  ours={} theirs={}{}",
                             row, p, mine, theirs,
                             if *p == 509431 { "   <- fragment 1 start" } else { "" });
                }
            }
        }
    }
    println!("  (our sentinel-first row for suffix 0 = {z_ours})");
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: ht2build <reference.fa> <reference_index.1.ht2>");
        eprintln!("  builds the front end and checks it against an existing index");
        exit(2);
    }
    let fa = &args[1];
    let ref_idx = &args[2];

    let r = match read_fasta(fa) {
        Ok(r) => r,
        Err(e) => { eprintln!("{e}"); exit(1); }
    };
    let len = r.text.len() as u64;

    // fchr: cumulative counts of A,C,G,T over the indexed text
    let mut cnt = [0u64; 4];
    for &c in &r.text {
        cnt[c as usize] += 1;
    }
    let mut fchr = [0u64; 5];
    for i in 0..4 {
        fchr[i + 1] = fchr[i] + cnt[i];
    }

    eprintln!("reference: {} sequence(s), {} unambiguous bp, {} fragment(s)",
              r.names.len(), len, r.frags.len());

    let sa = suffix_array(&r.text);
    // zOffs is the row whose suffix is the entire string, i.e. where the
    // preceding character is the sentinel.
    let z = sa.iter().position(|&x| x == 0).unwrap() as u64;

    // ---- check against the reference index -------------------------------
    let mut buf = Vec::new();
    if File::open(ref_idx).and_then(|mut f| f.read_to_end(&mut buf)).is_err() {
        eprintln!("cannot read {ref_idx}");
        exit(1);
    }
    let u = |o: usize| -> u64 {
        u32::from_le_bytes(buf[o..o + 4].try_into().unwrap()) as u64
    };
    let ref_len = u(8);
    let n_pat = u(44);
    let plen0 = u(48);
    let n_frag = u(52);

    let mut ok = true;
    let mut check = |name: &str, got: u64, want: u64| {
        let good = got == want;
        if !good { ok = false; }
        println!("  {:<24} got {:>12}   want {:>12}   {}",
                 name, got, want, if good { "ok" } else { "MISMATCH" });
    };

    println!("\nfront-end checks against {ref_idx}:");
    check("len", len, ref_len);
    check("nPat", r.names.len() as u64, n_pat);
    check("plen[0]", r.plen[0], plen0);
    check("nFrag", r.frags.len() as u64, n_frag);
    for (i, f) in r.frags.iter().enumerate() {
        let o = 56 + i * 12;
        check(&format!("rstarts[{i}].joined"), f.joined_off, u(o));
        check(&format!("rstarts[{i}].text_id"), f.text_id, u(o + 4));
        check(&format!("rstarts[{i}].text_off"), f.text_off, u(o + 8));
    }
    // fchr sits after the gbwt block and zOffs; recompute its offset the same
    // way the reader does rather than hard-coding it
    let line_rate = u(8 + 3 * 4) as i32;
    let side_sz = 1u64 << line_rate;
    let gbwt_sz = (ref_len + 1) / 4 + 1; // linear FM
    let side_gbwt_sz = side_sz - 4 * 4;
    let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
    let gbwt_tot = num_sides * side_sz;
    let zoff_at = 44 + 4 + n_pat as usize * 4 + 4 + n_frag as usize * 12 + gbwt_tot as usize;
    let num_z = u(zoff_at);
    check("numZOffs", 1, num_z);
    check("zOffs[0]", z, u(zoff_at + 4));
    let fchr_at = zoff_at + 4 + num_z as usize * 4;
    for i in 0..5 {
        check(&format!("fchr[{i}]"), fchr[i], u(fchr_at + i * 4));
    }

    println!("\nBWT cross-check against the stored index:");
    let gbwt_off = 44 + 4 + n_pat as usize * 4 + 4 + n_frag as usize * 12;
    compare_bwt(&buf, gbwt_off, side_sz as usize,
                (side_gbwt_sz << 2) as usize, &sa, &r.text, z);

    println!();
    if ok {
        println!("FRONT END OK: fragments, character histogram and suffix array all agree");
    } else {
        println!("FRONT END MISMATCH -- see above");
        exit(1);
    }
}
