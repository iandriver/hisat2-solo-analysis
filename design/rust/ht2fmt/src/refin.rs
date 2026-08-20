//! Reading a FASTA the way `RefReadInParams`/`joinToDisk` read one, streamed.
//!
//! Shared by `ht2emit` (linear) and `ht2wg` (graph) so the front end, `.3.ht2`
//! and `.4.ht2` come from one place. Nothing here holds the file: a 3.1 Gb
//! genome is already 3.1 GB as 2-bit-per-base codes, and slurping the FASTA to
//! find that out doubles it for no reason.

use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Clone, Copy)]
pub struct Frag { pub joined_off: u32, pub text_id: u32, pub text_off: u32 }

/// One `RefRecord` as `.3.ht2` stores it: the number of ambiguous characters
/// immediately preceding an unambiguous stretch, the stretch's length, and
/// whether it opens a new sequence. A sequence that ENDS in ambiguity
/// contributes a further record with `len == 0` carrying that trailing count,
/// which is why `.3.ht2` can hold more records than `.1.ht2` has `rstarts`
/// entries -- `joinToDisk` keeps only the non-empty ones.
#[derive(Clone, Copy)]
pub struct RefRec { pub off: u32, pub len: u32, pub first: bool }

pub struct Reference {
    pub names: Vec<String>,
    pub plen: Vec<u32>,
    pub frags: Vec<Frag>,
    pub recs: Vec<RefRec>,
    pub text: Vec<u8>,
}

pub fn code(c: u8) -> Option<u8> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None,
    }
}

pub fn read_fasta(path: &str) -> Reference {
    let f = BufReader::with_capacity(1 << 20, File::open(path).expect("cannot read reference"));
    let mut r = Reference { names: Vec::new(), plen: Vec::new(), frags: Vec::new(),
                            recs: Vec::new(), text: Vec::new() };
    let mut text_off: u32 = 0;
    let mut in_frag = false;
    let mut amb: u32 = 0;        // ambiguous characters seen since the last record
    let mut seq_started = false; // has this sequence produced a record yet?
    let mut line: Vec<u8> = Vec::new();
    let mut f = f;
    loop {
        line.clear();
        if f.read_until(b'\n', &mut line).expect("read") == 0 { break; }
        while line.last() == Some(&b'\n') || line.last() == Some(&b'\r') { line.pop(); }
        if line.is_empty() { continue; }
        if line[0] == b'>' {
            if !r.names.is_empty() {
                r.plen.push(text_off);
                if amb > 0 { r.recs.push(RefRec { off: amb, len: 0, first: !seq_started }); }
            }
            amb = 0;
            seq_started = false;
            // `_refnames` keeps the WHOLE header line. `_refnames_nospace`
            // (gfm.h:1400) is a separate, whitespace-truncated copy used only to
            // match chromosome names against SNP/splice-site files -- it is not
            // what gets written. A single-sequence reference whose header has no
            // description cannot tell the two apart; multi.fa can.
            r.names.push(line[1..].iter().map(|&c| c as char).collect());
            text_off = 0;
            in_frag = false;
            continue;
        }
        for &ch in line.iter() {
            match code(ch) {
                Some(v) => {
                    if !in_frag {
                        r.frags.push(Frag { joined_off: r.text.len() as u32,
                                            text_id: (r.names.len() - 1) as u32, text_off });
                        r.recs.push(RefRec { off: amb, len: 0, first: !seq_started });
                        amb = 0;
                        seq_started = true;
                        in_frag = true;
                    }
                    r.text.push(v);
                    r.recs.last_mut().unwrap().len += 1;
                }
                None => { in_frag = false; amb += 1; }
            }
            text_off += 1;
        }
    }
    if !r.names.is_empty() {
        r.plen.push(text_off);
        if amb > 0 { r.recs.push(RefRec { off: amb, len: 0, first: !seq_started }); }
    }
    r
}
