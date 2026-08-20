//! Emission from the external path: `generateEdges`, the GFM rows, the gbwt
//! block, and the graph `ftab` — all in bounded memory.
//!
//! `ht2path` does the same work with everything resident, and is the oracle for
//! this module. The differences are all about what is allowed to be an array:
//!
//!   * the CSR `start[]` over reference-graph node ids is `n_ref_nodes + 1`
//!     entries, 12 GB at whole-genome scale, so the "which path nodes have
//!     `from == t`" lookup becomes a sort-merge join instead;
//!   * `outdeg`, `indeg` and `floc` are one entry per path node — 5.9e9 of them
//!     — so each becomes a stream consumed in node order;
//!   * the `ftab` cannot precompute `occ`/`rank_M`/`select_F` over every row.
//!     Those become side-local queries against the block that was just written,
//!     which is what the per-side tallies exist for.
//!
//! The one thing kept in RAM is an F-bit rank per side: `select_F` is the only
//! primitive with no per-side tally to stand on, and one `u64` per side is
//! 228 MB for a human graph.

use super::ext::{self, By, Rec, SegReader, SegWriter};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

pub const LABELS: &[u8; 6] = b"ACGTYZ";

pub fn label_index(l: u8) -> u64 {
    match l { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, b'Y' => 4, b'Z' => 5,
              _ => panic!("bad label {}", l as char) }
}

/// Side geometry, derived exactly as `GFMParams::init` derives it (gfm.h:138).
///
/// `w` is `sizeof(index_t)`: a large index reserves 48 bytes per side for its
/// six tallies rather than 24, which is why `hisat2-build-l` also moves the
/// default line rate from 7 to 8.
#[derive(Clone, Copy)]
pub struct Geom {
    pub side_sz: usize,
    pub side_gbwt_sz: usize,
    pub rows_per_side: usize,
    pub num_sides: usize,
    pub gbwt_tot: usize,
    pub w: usize,
}

pub fn geom(gbwt_len: u64, line_rate: u32, w: usize) -> Geom {
    let side_sz = 1usize << line_rate;
    let side_gbwt_sz = side_sz - 6 * w;
    let gbwt_sz = (gbwt_len / 2 + 1) as usize;
    let num_sides = (gbwt_sz + side_gbwt_sz - 1) / side_gbwt_sz;
    Geom { side_sz, side_gbwt_sz, rows_per_side: side_gbwt_sz * 2, num_sides,
           gbwt_tot: num_sides * side_sz, w }
}

/// `writeIndex<index_t>` — 4 or 8 bytes little-endian.
pub fn put_idx(v: &mut Vec<u8>, x: u64, w: usize) {
    v.extend_from_slice(&x.to_le_bytes()[0..w]);
}

// ---------------------------------------------------------------------------
// small fixed-width stream helpers

struct U40Writer { w: BufWriter<File> }
impl U40Writer {
    fn create(p: &Path) -> std::io::Result<U40Writer> {
        Ok(U40Writer { w: BufWriter::with_capacity(1 << 20, File::create(p)?) })
    }
    fn push(&mut self, v: u64) -> std::io::Result<()> { self.w.write_all(&v.to_le_bytes()[0..5]) }
    fn finish(mut self) -> std::io::Result<()> { self.w.flush() }
}
struct U40Reader { r: BufReader<File> }
impl U40Reader {
    fn open(p: &Path) -> std::io::Result<U40Reader> {
        Ok(U40Reader { r: BufReader::with_capacity(1 << 20, File::open(p)?) })
    }
    fn next(&mut self) -> std::io::Result<Option<u64>> {
        let mut b = [0u8; 8];
        match self.r.read_exact(&mut b[0..5]) {
            Ok(()) => Ok(Some(u64::from_le_bytes(b))),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
            Err(e) => Err(e),
        }
    }
}

/// One `(label, F)` per GFM row, packed into a byte.
struct RowWriter { w: BufWriter<File> }
impl RowWriter {
    fn create(p: &Path) -> std::io::Result<RowWriter> {
        Ok(RowWriter { w: BufWriter::with_capacity(1 << 20, File::create(p)?) })
    }
    fn push(&mut self, label: u8, f: u8) -> std::io::Result<()> {
        self.w.write_all(&[label_index(label) as u8 | (f << 3)])
    }
    fn finish(mut self) -> std::io::Result<()> { self.w.flush() }
}
pub struct RowReader { r: BufReader<File> }
impl RowReader {
    pub fn open(p: &Path) -> std::io::Result<RowReader> {
        Ok(RowReader { r: BufReader::with_capacity(1 << 20, File::open(p)?) })
    }
    /// `(label byte, F bit)`
    pub fn next(&mut self) -> std::io::Result<Option<(u8, u8)>> {
        let mut b = [0u8; 1];
        match self.r.read_exact(&mut b) {
            Ok(()) => Ok(Some((LABELS[(b[0] & 7) as usize], b[0] >> 3))),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
            Err(e) => Err(e),
        }
    }
}

/// `(out-degree, genomic position)` per path node, in node order.
struct NodeInfoWriter { w: BufWriter<File> }
impl NodeInfoWriter {
    fn create(p: &Path) -> std::io::Result<NodeInfoWriter> {
        Ok(NodeInfoWriter { w: BufWriter::with_capacity(1 << 20, File::create(p)?) })
    }
    fn push(&mut self, outdeg: u32, pos: u32) -> std::io::Result<()> {
        self.w.write_all(&outdeg.to_le_bytes())?; self.w.write_all(&pos.to_le_bytes())
    }
    fn finish(mut self) -> std::io::Result<()> { self.w.flush() }
}
struct NodeInfoReader { r: BufReader<File> }
impl NodeInfoReader {
    fn open(p: &Path) -> std::io::Result<NodeInfoReader> {
        Ok(NodeInfoReader { r: BufReader::with_capacity(1 << 20, File::open(p)?) })
    }
    fn next(&mut self) -> std::io::Result<Option<(u32, u32)>> {
        let mut b = [0u8; 8];
        match self.r.read_exact(&mut b) {
            Ok(()) => Ok(Some((u32::from_le_bytes(b[0..4].try_into().unwrap()),
                               u32::from_le_bytes(b[4..8].try_into().unwrap())))),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
            Err(e) => Err(e),
        }
    }
}

// ---------------------------------------------------------------------------
// generateEdges

/// What the block writer consumes: three streams, all in node or row order.
pub struct Rows {
    pub rows: PathBuf,        // one byte per GFM row: label index | F << 3
    pub nodeinfo: PathBuf,    // (out-degree, position) per path node
    pub floc: PathBuf,        // F-run start row per path node, 5 bytes each
    pub n_nodes: u64,         // path nodes AFTER the second-to-last is dropped
    pub gbwt_len: u64,
    pub bucket: [u64; 6],     // rows per label, which is what fchr encodes
}

/// `generateEdges` (gbwt_graph.h:2367) plus the tail of it that `nextRow` needs,
/// as a chain of external sorts and sequential merges.
///
/// The in-memory version indexes `start[]` by reference-graph node id to find
/// the path nodes a given edge lands on. Here that is a join between the edges
/// keyed by their `to` and the path nodes keyed by their `from`, which is the
/// same equi-join the doubling loop already runs this way.
pub fn generate_edges(wd: &Path, cur: &Path, budget: usize, verbose: bool)
    -> std::io::Result<Rows>
{
    let nodes_bin = wd.join("nodes.bin");
    let edges_bin = wd.join("edges.bin");
    let (nbf0, nbf) = (wd.join("nbf0.bin"), wd.join("nbf.bin"));
    let (elab, ebt) = (wd.join("elab.bin"), wd.join("ebt.bin"));
    let (trip, ed, ed2, edr) =
        (wd.join("trip.bin"), wd.join("ed.bin"), wd.join("ed2.bin"), wd.join("edr.bin"));
    let nbr = wd.join("nbr.bin");

    // (a) path nodes by `from`, which is what the `start[]` CSR indexes: for a
    //     given reference node, the path nodes that hang off it.
    ext::sort_external(cur, &nbf0, By::From, budget, wd, true)?;

    // (b) give each path node the genomic position of its reference node, the
    //     way `n.to = b.nodes[n.from].1` does. `nodes.bin` is in id order and
    //     `nbf0` is sorted by id, so this is one forward pass.
    let n_nodes_raw = {
        let mut r = SegReader::open(&nbf0, true)?;
        let mut nf = BufReader::with_capacity(1 << 20, File::open(&nodes_bin)?);
        let mut w = SegWriter::create(&nbf)?;
        let mut nbuf = [0u8; 5];
        let (mut node_i, mut val) = (0u64, 0u32);
        let mut idx = 0u64;
        while let Some(mut n) = r.next()? {
            while node_i <= n.from as u64 {
                nf.read_exact(&mut nbuf)?;
                val = u32::from_le_bytes(nbuf[1..5].try_into().unwrap());
                node_i += 1;
            }
            n.to = val;
            n.k1 = idx; idx += 1;
            w.push(n)?;
            let _ = &idx;
        }
        w.finish()?
    };
    ext::seg_remove(&nbf0);

    // (c) each reference edge carries the LABEL of its `from` node into the
    //     join, because after sorting by `to` the label is no longer a forward
    //     scan away.
    {
        let mut ef = BufReader::with_capacity(1 << 20, File::open(&edges_bin)?);
        let mut nf = BufReader::with_capacity(1 << 20, File::open(&nodes_bin)?);
        let mut w = SegWriter::create(&elab)?;
        let (mut nbuf, mut ebuf) = ([0u8; 5], [0u8; 8]);
        let (mut node_i, mut lab) = (0u64, 0u8);
        loop {
            match ef.read_exact(&mut ebuf) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            }
            let f = u32::from_le_bytes(ebuf[0..4].try_into().unwrap());
            let t = u32::from_le_bytes(ebuf[4..8].try_into().unwrap());
            while node_i <= f as u64 { nf.read_exact(&mut nbuf)?; lab = nbuf[0]; node_i += 1; }
            w.push(Rec { from: f, to: t, k0: 0, k1: label_index(lab) })?;
        }
        w.finish()?;
    }
    ext::sort_external(&elab, &ebt, By::To, budget, wd, true)?;

    // (d) the join itself: for each edge, one row per path node whose `from` is
    //     the edge's `to`, labelled by the edge's `from` node.
    let mut bucket = [0u64; 6];
    {
        // `ebt` is dead after this; `nbf` is not -- the rank-ordered copy is
        // built from it below, and that sort is what frees it.
        let mut e = SegReader::open(&ebt, true)?;
        let mut n = SegReader::open(&nbf, false)?;
        let mut w = SegWriter::create(&trip)?;
        let mut ncur = n.next()?;
        let mut group: Vec<Rec> = Vec::new();
        let mut group_key = u32::MAX;
        let mut ecur = e.next()?;
        while let Some(ed_) = ecur {
            if ed_.to != group_key {
                while let Some(nn) = ncur { if nn.from < ed_.to { ncur = n.next()?; } else { break; } }
                group.clear();
                while let Some(nn) = ncur {
                    if nn.from == ed_.to { group.push(nn); ncur = n.next()?; } else { break; }
                }
                group_key = ed_.to;
            }
            for g in &group {
                bucket[ed_.k1 as usize] += 1;
                w.push(Rec { from: ed_.from, to: ed_.to, k0: g.k0, k1: ed_.k1 })?;
            }
            ecur = e.next()?;
        }
        w.finish()?;
    }
    ext::seg_remove(&ebt);

    // The merge-join further down walks the nodes in RANK order, not `from`
    // order -- `generateEdges` re-sorts them (`nodes.sort_by_key(|n| n.k0)`)
    // between building the buckets and pairing them off, and missing that makes
    // the scan stall on its first mismatch. Ranks are unique at convergence, so
    // this is just `cur.bin`'s own order with the positions attached.
    //
    // Built HERE rather than beside `nbf`, so the two orderings of the node list
    // are never alive at the same time as the join's inputs. `nbf` was only ever
    // the `start[]` side of that join and is dead the moment it finishes.
    ext::sort_external(&nbf, &nbr, By::Key, budget, wd, true)?;
    ext::sort_external(&trip, &ed, By::Row, budget, wd, true)?;

    let gbwt_len: u64 = bucket.iter().sum();
    if verbose {
        println!("\ngenerateEdges: {} path nodes, {} path edges (A {} C {} G {} T {} Y {} Z {})",
                 n_nodes_raw, gbwt_len, bucket[0], bucket[1], bucket[2], bucket[3],
                 bucket[4], bucket[5]);
        println!("  numNodes = {} (path nodes - 1), gbwtLen = {} (path edges)",
                 n_nodes_raw - 1, gbwt_len);
    }

    // (e) the merge-join of gbwt_graph.h:2563. Instrumenting a debug
    //     hisat2-build shows the two lists arrive positionally aligned -- edge
    //     i's `from` equals node i's `from` -- so the scan pairs them off,
    //     advancing the node only on a mismatch and giving a node two rows where
    //     two edges land on it. Only the out-degrees survive it; the rewritten
    //     `from` is never read again.
    let nodeinfo_raw = wd.join("nodeinfo_raw.bin");
    {
        let mut e = SegReader::open(&ed, false)?;
        let mut n = SegReader::open(&nbr, true)?;
        let mut w = NodeInfoWriter::create(&nodeinfo_raw)?;
        let mut ncur = n.next()?;
        let mut ecur = e.next()?;
        let mut outdeg = 0u32;
        let mut consumed = 0u64;
        while let (Some(nn), Some(ee)) = (ncur, ecur) {
            if ee.from == nn.from { outdeg += 1; ecur = e.next()?; consumed += 1; }
            else { w.push(outdeg, nn.to)?; outdeg = 0; ncur = n.next()?; }
        }
        // whichever stream ran out, the remaining nodes still need an entry
        while let Some(nn) = ncur { w.push(outdeg, nn.to)?; outdeg = 0; ncur = n.next()?; }
        w.finish()?;
        if verbose && consumed < gbwt_len {
            println!("  NOTE merge-join consumed {consumed} of {gbwt_len} edges");
        }
    }

    // (f) drop the second-to-last node, keeping its out-degree on the survivor.
    let n_nodes = n_nodes_raw - 1;
    let nodeinfo = wd.join("nodeinfo.bin");
    {
        let mut r = NodeInfoReader::open(&nodeinfo_raw)?;
        let mut w = NodeInfoWriter::create(&nodeinfo)?;
        let mut hold: Vec<(u32, u32)> = Vec::new();
        while let Some(x) = r.next()? {
            hold.push(x);
            if hold.len() > 2 { let h = hold.remove(0); w.push(h.0, h.1)?; }
        }
        // hold = [second-to-last, last]; the survivor takes the last node's
        // position and the second-to-last node's out-degree
        if hold.len() == 2 { w.push(hold[0].0, hold[1].1)?; }
        else { for h in &hold { w.push(h.0, h.1)?; } }
        w.finish()?;
    }
    let _ = fs::remove_file(&nodeinfo_raw);

    // (g) relabel 'Y' to 'Z' and pull the rankings past the dropped node down.
    //     `k1` becomes `label * 2 + decremented` so the re-sort can reproduce a
    //     stable sort on ranking alone: two records can otherwise reach the same
    //     ranking from different ones and their original order is not recoverable
    //     from the new key.
    {
        let mut r = SegReader::open(&ed, true)?;
        let mut w = SegWriter::create(&ed2)?;
        while let Some(mut x) = r.next()? {
            let mut dec = 0u64;
            if x.k1 == 4 { x.k1 = 5; }
            else if x.k0 >= n_nodes { x.k0 -= 1; dec = 1; }
            x.k1 = x.k1 * 2 + dec;
            w.push(x)?;
        }
        w.finish()?;
    }
    ext::seg_remove(&ed);
    ext::sort_external(&ed2, &edr, By::Rank, budget, wd, true)?;

    // (h) rows in nextRow order, with the F bit, plus the F-run start of every
    //     node. `floc` has an entry for EVERY node, including those with no
    //     incoming edges, because the M walk consumes it one per node.
    let rows = wd.join("rows.bin");
    let floc = wd.join("floc.bin");
    {
        let mut r = SegReader::open(&edr, true)?;
        let mut w = RowWriter::create(&rows)?;
        let mut fw = U40Writer::create(&floc)?;
        let mut emitted = 0u64;
        let mut node = 0u64;
        let mut cur = r.next()?;
        while node < n_nodes {
            fw.push(emitted)?;
            let mut first = true;
            while let Some(x) = cur {
                if x.k0 != node { break; }
                w.push(LABELS[(x.k1 / 2) as usize], if first { 1 } else { 0 })?;
                first = false; emitted += 1;
                cur = r.next()?;
            }
            node += 1;
        }
        w.finish()?; fw.finish()?;
        if verbose { println!("  nextRow order: {emitted} rows of {gbwt_len}"); }
    }
    ext::seg_remove(&edr);
    ext::seg_remove(&nbr);
    Ok(Rows { rows, nodeinfo, floc, n_nodes, gbwt_len, bucket })
}

// ---------------------------------------------------------------------------
// the gbwt block

pub struct BlockOut {
    pub fchr: [u64; 5],
    pub z_offs: Vec<u64>,
    /// F bits before each side — the one rank structure with no per-side tally
    /// to stand on, and the only thing this module keeps in RAM.
    pub f_rank_save: Vec<u64>,
    pub n_sa: u64,
}

/// Write the gbwt block, the `.2.ht2` SA sample, `fchr` and `zOffs` in one pass.
///
/// Graph side layout (gfm.h:4886, :4926). Within each side's `sideGbwtSz` bytes:
/// `[0, sz/2)` BWT at 4 rows/byte low-pair-first, `[sz/2, 3sz/4)` the F
/// bitvector at 8 rows/byte with `F_bpi = bpi + ((sideCur & 1) << 2)`, then M
/// the same way. The final 6 `index_t` are `F_locSave`, `M_occSave` and
/// `occSave[0..3]` — the values as of the START of the side.
pub fn write_block(rs: &Rows, g: Geom, off_rate: u32,
                   w1: &mut impl Write, w2: &mut impl Write) -> std::io::Result<BlockOut>
{
    let mut rr = RowReader::open(&rs.rows)?;
    let mut ni = NodeInfoReader::open(&rs.nodeinfo)?;
    let mut fl = U40Reader::open(&rs.floc)?;

    let mut occ = [0u64; 4];
    let (mut m_occ, mut f_loc) = (0u64, 0u64);
    let mut fchr_c = [0u64; 4];
    let mut z_offs: Vec<u64> = Vec::new();
    let mut f_rank_save: Vec<u64> = Vec::with_capacity(g.num_sides);
    let mut f_count = 0u64;
    let mut n_sa = 0u64;
    let off_mask: u64 = u64::MAX << off_rate;

    // the M walk: a second cursor over the same nodes, advancing by out-degree
    let mut m_left = 0u32;
    let mut m_pos = 0u32;
    let mut m_done = false;

    let mut side = vec![0u8; g.side_sz];
    let mut row_i: u64 = 0;
    for s in 0..g.num_sides {
        for b in side.iter_mut() { *b = 0; }
        {
            let base = g.side_sz - 6 * g.w;
            for (k, v) in [f_loc, m_occ, occ[0], occ[1], occ[2], occ[3]].iter().enumerate() {
                side[base + k * g.w..base + (k + 1) * g.w]
                    .copy_from_slice(&v.to_le_bytes()[0..g.w]);
            }
        }
        f_rank_save.push(f_count);
        for off in 0..g.rows_per_side {
            let (ch, f) = match if row_i < rs.gbwt_len { rr.next()? } else { None } {
                Some(x) => x,
                None => (b'A', 0),      // padding past the end, counted as 'A'
            };
            let in_range = row_i < rs.gbwt_len;
            // the M cursor
            let mut m = 0u8;
            if !m_done {
                if m_left == 0 {
                    match ni.next()? {
                        Some((outdeg, pos)) => { m_left = outdeg.max(1); m_pos = pos; m = 1; }
                        None => { m_done = true; }
                    }
                }
                if !m_done { m_left -= 1; }
            }
            let mut count = true;
            let code = match ch { b'A' => 0u8, b'C' => 1, b'G' => 2, b'T' => 3,
                                  _ => { count = false; z_offs.push(row_i); 0 } };
            if in_range && count { fchr_c[code as usize] += 1; }
            if m == 1 {
                if let Some(v) = fl.next()? { f_loc = v; }
                if (m_occ & off_mask) == m_occ {
                    // A node with no genomic position carries `(index_t)INDEX_MAX`.
                    // The graph is built with 32-bit ids, so the sentinel arrives
                    // here as `u32::MAX` and has to be WIDENED, not zero-extended
                    // -- a large index writes 0xFFFF_FFFF_FFFF_FFFF. Invisible at
                    // 32 bits, where the two are the same number.
                    let v = if m_pos == u32::MAX { u64::MAX } else { m_pos as u64 };
                    w2.write_all(&v.to_le_bytes()[0..g.w])?;
                    n_sa += 1;
                }
            }
            if count { occ[code as usize] += 1; }
            if m == 1 { m_occ += 1; }
            if f == 1 { f_count += 1; }

            let sc = off >> 2;
            let bpi = off & 3;
            side[sc] |= code << (bpi * 2);
            let f_sc = (g.side_gbwt_sz + sc) >> 1;
            let f_bpi = bpi + ((sc & 1) << 2);
            side[f_sc] |= f << f_bpi;
            side[f_sc + (g.side_gbwt_sz >> 2)] |= m << f_bpi;
            row_i += 1;
        }
        w1.write_all(&side)?;
        let _ = s;
    }
    let mut fchr = [0u64; 5];
    for i in 0..4 { fchr[i + 1] = fchr[i] + fchr_c[i]; }
    Ok(BlockOut { fchr, z_offs, f_rank_save, n_sa })
}

// ---------------------------------------------------------------------------
// GFM navigation over the written block, for the ftab

/// Random access to the gbwt block of a partially written `.1.ht2`.
///
/// Every query is answered from ONE side: the per-side tallies hold `occ`,
/// `M_occSave` and `F_locSave` as of that side's start, so a query costs one
/// 128-byte read plus a scan of at most `rows_per_side` rows. `select_F` is the
/// exception — nothing in the file counts F bits — so it binary-searches the
/// in-RAM per-side F rank and then reads one side.
pub struct Nav {
    f: File,
    off: u64,
    g: Geom,
    pub gbwt_len: u64,
    fchr: [u64; 5],
    zset: Vec<u64>,
    f_rank_save: Vec<u64>,
}

impl Nav {
    pub fn new(f: File, off: u64, g: Geom, gbwt_len: u64, fchr: [u64; 5],
               zset: Vec<u64>, f_rank_save: Vec<u64>) -> Nav {
        Nav { f, off, g, gbwt_len, fchr, zset, f_rank_save }
    }

    fn side(&self, s: usize, buf: &mut [u8]) -> std::io::Result<()> {
        self.f.read_exact_at(buf, self.off + (s * self.g.side_sz) as u64)
    }
    fn tally(buf: &[u8], side_sz: usize, w: usize, k: usize) -> u64 {
        let mut b = [0u8; 8];
        let at = side_sz - 6 * w + k * w;
        b[0..w].copy_from_slice(&buf[at..at + w]);
        u64::from_le_bytes(b)
    }
    fn bwt_at(buf: &[u8], off: usize) -> u8 { (buf[off >> 2] >> ((off & 3) * 2)) & 3 }
    fn fbit_at(buf: &[u8], g: &Geom, off: usize) -> u8 {
        let sc = off >> 2;
        (buf[(g.side_gbwt_sz + sc) >> 1] >> ((off & 3) + ((sc & 1) << 2))) & 1
    }
    fn mbit_at(buf: &[u8], g: &Geom, off: usize) -> u8 {
        let sc = off >> 2;
        (buf[((g.side_gbwt_sz + sc) >> 1) + (g.side_gbwt_sz >> 2)] >> ((off & 3) + ((sc & 1) << 2))) & 1
    }

    /// Number of `c` characters in rows `[0, row)`. The `'Z'` row is stored as
    /// `'A'` but was never counted when the block was written, so it has to be
    /// discounted here too or every LF step after it shifts.
    pub fn occ(&self, row: u64, c: usize, buf: &mut [u8]) -> std::io::Result<u64> {
        let s = (row / self.g.rows_per_side as u64) as usize;
        let start = s as u64 * self.g.rows_per_side as u64;
        self.side(s, buf)?;
        let mut n = Self::tally(buf, self.g.side_sz, self.g.w, 2 + c);
        for off in 0..(row - start) as usize {
            if Self::bwt_at(buf, off) as usize == c { n += 1; }
        }
        if c == 0 {
            for &z in &self.zset { if z >= start && z < row { n -= 1; } }
        }
        Ok(n)
    }

    /// M bits in rows `[0, x)` — the exclusive prefix sum
    /// `rank_M(initFromRow_bit(x))` computes.
    pub fn rank_m(&self, x: u64, buf: &mut [u8]) -> std::io::Result<u64> {
        let x = x.min(self.gbwt_len);
        let s = (x / self.g.rows_per_side as u64) as usize;
        let start = s as u64 * self.g.rows_per_side as u64;
        self.side(s, buf)?;
        let mut n = Self::tally(buf, self.g.side_sz, self.g.w, 1);
        for off in 0..(x - start) as usize {
            n += Self::mbit_at(buf, &self.g, off) as u64;
        }
        Ok(n)
    }

    /// Row of the `k`-th F bit, one-based, matching `self_f[k - 1]`.
    pub fn select_f(&self, k: u64, buf: &mut [u8]) -> std::io::Result<u64> {
        if k == 0 { return Ok(self.gbwt_len); }
        // last side whose start has fewer than k F bits before it
        let mut lo = 0usize;
        let mut hi = self.f_rank_save.len();
        while lo + 1 < hi {
            let mid = (lo + hi) / 2;
            if self.f_rank_save[mid] < k { lo = mid; } else { hi = mid; }
        }
        if self.f_rank_save[lo] >= k { return Ok(self.gbwt_len); }
        self.side(lo, buf)?;
        let mut n = self.f_rank_save[lo];
        for off in 0..self.g.rows_per_side {
            if Self::fbit_at(buf, &self.g, off) == 1 {
                n += 1;
                if n == k {
                    let row = lo as u64 * self.g.rows_per_side as u64 + off as u64;
                    return Ok(if row < self.gbwt_len { row } else { self.gbwt_len });
                }
            }
        }
        Ok(self.gbwt_len)
    }

    /// One `mapGLF` step: an LF over the BWT, then a hop through the node
    /// structure — rank over M to find which node the row belongs to, then
    /// select over F to find where that node's incoming edges begin.
    /// `None` means the range went empty, which ends the walk.
    fn step(&self, top: u64, bot: u64, c: usize, buf: &mut [u8])
        -> std::io::Result<Option<(u64, u64)>>
    {
        let nt = self.fchr[c] + self.occ(top, c, buf)?;
        let nb = self.fchr[c] + self.occ(bot, c, buf)?;
        if nt >= nb { return Ok(None); }
        let node_top = self.rank_m(nt + 1, buf)? - 1;
        let node_bot = self.rank_m(nb, buf)?;
        let t = self.select_f(node_top + 1, buf)?;
        let b = self.select_f(node_bot + 1, buf)?;
        if t >= b { return Ok(None); }
        Ok(Some((t, b)))
    }
}

/// Build the graph `ftab` and `eftab` by querying the block just written.
///
/// The straight transcription walks all `4^ftabChars` prefixes independently,
/// which is 10.5M LF steps for `ftabChars` 10 and every one of them a random
/// read into an 11 GB file. But `mapGLF` consumes the prefix from its LAST
/// character first, so prefixes sharing low-order bits share a walk: doing it as
/// a depth-first traversal of that trie costs at most 1.4M steps, and prunes
/// whole subtrees the moment a range goes empty.
pub fn build_ftab(nav: &Nav, ftab_chars: u32, w: usize, verbose: bool)
    -> std::io::Result<(Vec<u64>, Vec<u64>)>
{
    let ftab_len = (1usize << (ftab_chars * 2)) + 1;
    let mut tftab: Vec<(u64, u64)> = vec![(u64::MAX, u64::MAX); ftab_len - 1];
    let mut buf = vec![0u8; nav.g.side_sz];
    let mut steps = 0u64;

    // (top, bot, depth, low bits of i fixed so far)
    let mut stack: Vec<(u64, u64, u32, usize)> = vec![(0, nav.gbwt_len, 0, 0)];
    while let Some((top, bot, depth, bits)) = stack.pop() {
        for c in 0..4usize {
            steps += 1;
            match nav.step(top, bot, c, &mut buf)? {
                None => {}                       // whole subtree is degenerate
                Some((t, b)) => {
                    let nbits = bits | (c << (2 * depth));
                    if depth + 1 == ftab_chars { tftab[nbits] = (t, b); }
                    else { stack.push((t, b, depth + 1, nbits)); }
                }
            }
        }
    }
    if verbose { println!("  ftab: {steps} LF steps over the trie (flat walk would be {})",
                          (ftab_len - 1) as u64 * ftab_chars as u64); }

    // A prefix whose walk broke takes the previous entry's upper bound, so the
    // fill is sequential even though the search was not.
    let mut prev = 0u64;
    for i in 0..ftab_len - 1 {
        if tftab[i].0 == u64::MAX { tftab[i] = (prev, prev); }
        prev = tftab[i].1;
    }

    let mut ftab_o = vec![0u64; ftab_len];
    let mut eftab_o: Vec<u64> = Vec::new();
    ftab_o[0] = tftab[0].0; ftab_o[1] = tftab[0].1;
    for i in 1..ftab_len - 1 {
        if ftab_o[i] != tftab[i].0 {
            let (lo, hi) = (ftab_o[i], tftab[i].0);
            // `ftab[i] = eftabCur ^ (index_t)INDEX_MAX` (gfm.h:5088)
            let index_max = if w == 4 { u32::MAX as u64 } else { u64::MAX };
            ftab_o[i] = (eftab_o.len() as u64 / 2) ^ index_max;
            eftab_o.push(lo); eftab_o.push(hi);
        }
        ftab_o[i + 1] = tftab[i].1;
    }
    Ok((ftab_o, eftab_o))
}
