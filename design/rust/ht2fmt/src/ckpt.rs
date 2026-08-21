//! Resume state, for a build long enough that losing it matters.
//!
//! The doubling loop is already checkpointed by its own shape: at the end of a
//! generation `cur.bin` is a complete, self-describing state and nothing else is
//! live. All this adds is a note of which generation that is, plus the handful
//! of numbers that cannot be recovered from the files -- the generation curve so
//! far, and the counts `generateEdges` prints.
//!
//! A checkpoint is only trusted when its fingerprint matches, and the
//! fingerprint includes the BUILDER'S OWN mtime and size. Reusing half a build
//! across a code change is the same failure as building fixtures with a stale
//! binary: everything looks consistent and nothing is.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

#[derive(Default, Clone)]
pub struct Ckpt {
    pub fingerprint: String,
    /// "", "graph", "gen", "edges" -- the last stage that finished
    pub stage: String,
    /// Progress inside the generation after `gen`: "", "sorted", "joined",
    /// "keyed". Each names an artifact on disk that a restart can pick up from,
    /// which is why the step that would overwrite its input has to be the step
    /// AFTER the checkpoint, never the same one.
    pub phase: String,
    pub g_nodes: u64,
    pub g_edges: u64,
    pub g_last: u32,
    pub g_text: u32,
    pub gen: u32,
    pub curve: Vec<(u32, u64, u64, u64)>,
    pub path_nodes: u64,
    pub sorted: u64,
    pub e_nodes: u64,
    pub e_gbwt: u64,
    pub e_bucket: [u64; 6],
    /// `temp_nodes` for the generation in progress
    pub temp: u64,
}

pub fn path(wd: &Path) -> PathBuf { wd.join("checkpoint.txt") }

/// Everything that changes the output. Deliberately NOT the thread count or the
/// sort budget: neither affects a single byte, so a run interrupted on eighteen
/// threads can be finished on four.
pub fn fingerprint(fa: &str, snp: &str, hap: &str, large: bool, chunk: u32) -> String {
    let mut s = String::new();
    for f in [fa, snp, hap] {
        let m = std::fs::metadata(f).ok();
        let (len, mt) = match &m {
            Some(m) => (m.len(), m.modified().ok()
                        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
                        .map(|d| d.as_secs()).unwrap_or(0)),
            None => (0, 0),
        };
        let _ = write!(s, "{f}:{len}:{mt} ");
    }
    if let Ok(exe) = std::env::current_exe() {
        if let Ok(m) = std::fs::metadata(&exe) {
            let mt = m.modified().ok()
                .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
                .map(|d| d.as_secs()).unwrap_or(0);
            let _ = write!(s, "exe:{}:{mt} ", m.len());
        }
    }
    let _ = write!(s, "large={} chunk={}", large as u8, chunk);
    s
}

pub fn load(wd: &Path, fingerprint: &str) -> Option<Ckpt> {
    let text = std::fs::read_to_string(path(wd)).ok()?;
    let mut c = Ckpt::default();
    for line in text.lines() {
        let mut it = line.splitn(2, ' ');
        let (k, v) = (it.next().unwrap_or(""), it.next().unwrap_or(""));
        match k {
            "fingerprint" => c.fingerprint = v.to_string(),
            "stage" => c.stage = v.to_string(),
            "phase" => c.phase = v.to_string(),
            "graph" => {
                let n: Vec<u64> = v.split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if n.len() == 4 {
                    c.g_nodes = n[0]; c.g_edges = n[1];
                    c.g_last = n[2] as u32; c.g_text = n[3] as u32;
                }
            }
            "gen" => c.gen = v.parse().unwrap_or(0),
            "nodes" => c.path_nodes = v.parse().unwrap_or(0),
            "sorted" => c.sorted = v.parse().unwrap_or(0),
            "temp" => c.temp = v.parse().unwrap_or(0),
            "curve" => {
                let n: Vec<u64> = v.split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if n.len() == 4 { c.curve.push((n[0] as u32, n[1], n[2], n[3])); }
            }
            "edges" => {
                let n: Vec<u64> = v.split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if n.len() == 8 {
                    c.e_nodes = n[0]; c.e_gbwt = n[1];
                    for i in 0..6 { c.e_bucket[i] = n[2 + i]; }
                }
            }
            _ => {}
        }
    }
    if c.fingerprint != fingerprint { return None; }
    Some(c)
}

pub fn save(wd: &Path, c: &Ckpt) -> std::io::Result<()> {
    let mut s = String::new();
    let _ = writeln!(s, "fingerprint {}", c.fingerprint);
    let _ = writeln!(s, "stage {}", c.stage);
    let _ = writeln!(s, "phase {}", c.phase);
    let _ = writeln!(s, "graph {} {} {} {}", c.g_nodes, c.g_edges, c.g_last, c.g_text);
    let _ = writeln!(s, "gen {}", c.gen);
    let _ = writeln!(s, "nodes {}", c.path_nodes);
    let _ = writeln!(s, "sorted {}", c.sorted);
    let _ = writeln!(s, "temp {}", c.temp);
    for (g, t, n, r) in &c.curve { let _ = writeln!(s, "curve {g} {t} {n} {r}"); }
    let _ = writeln!(s, "edges {} {} {} {} {} {} {} {}", c.e_nodes, c.e_gbwt,
                     c.e_bucket[0], c.e_bucket[1], c.e_bucket[2],
                     c.e_bucket[3], c.e_bucket[4], c.e_bucket[5]);
    // rename, so a crash mid-write cannot leave a checkpoint that half-describes
    // a state nobody was ever in
    let tmp = wd.join("checkpoint.tmp");
    std::fs::write(&tmp, s)?;
    std::fs::rename(&tmp, path(wd))
}

/// Delete everything in the workdir that a resume must not inherit.
///
/// A crash can land anywhere, including part-way through writing a run or a
/// merge partition, and a stale segment does not announce itself -- it gets
/// adopted by the next file of the same name and shows up as records from
/// another generation merged in as though they belonged. So the rule is a
/// whitelist of what a checkpoint actually vouches for; everything else goes.
pub fn clean(wd: &Path, keep: &[&str]) {
    if let Ok(rd) = std::fs::read_dir(wd) {
        for e in rd.flatten() {
            let name = match e.file_name().into_string() { Ok(n) => n, Err(_) => continue };
            if name == "checkpoint.txt" { continue; }
            // segments are `<base>.sNNNNN`
            let base = match name.find(".s") {
                Some(i) if name[i + 2..].chars().all(|c| c.is_ascii_digit())
                           && name.len() > i + 2 => &name[..i],
                _ => name.as_str(),
            };
            if keep.contains(&base) { continue; }
            let _ = std::fs::remove_file(e.path());
        }
    }
}

/// A deterministic crash, for testing resume.
///
/// Racing a `kill -9` against a build finds bugs but cannot say which point it
/// found them at, and cannot prove a point was ever reached. `HT2_CRASH_AT`
/// names a checkpoint boundary and `HT2_CRASH_GEN` narrows it to one
/// generation, so every place a restart can land gets visited on purpose.
pub fn crash_point(what: &str, gen: u32) {
    if let Ok(at) = std::env::var("HT2_CRASH_AT") {
        if at != what { return; }
        if let Ok(g) = std::env::var("HT2_CRASH_GEN") {
            if g.parse::<u32>().ok() != Some(gen) { return; }
        }
        eprintln!("HT2_CRASH_AT={what} gen={gen}");
        std::process::exit(70);
    }
}
