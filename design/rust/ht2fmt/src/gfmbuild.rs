//! Build one GFM from a graph and produce everything an emitter needs.
//!
//! Factored out so the global index and the ~57 kb local indexes share it: they
//! are the same construction at different widths and parameters, which is the
//! whole reason `.5`/`.6` needed no new format work.

#[allow(dead_code)]
pub struct BuiltGfm {
    pub len: u32,
    pub gbwt_len: u32,
    pub num_nodes: u32,
    pub rows: Vec<(u8, u8)>,   // (label, F bit) in nextRow order
    pub mrun: Vec<(u8, u32)>,  // (M bit, genomic position)
    pub fchr: [u32; 5],
    pub z_offs: Vec<u32>,
    pub sa_sample: Vec<u32>,
}

const UNSET: u32 = u32::MAX;

#[derive(Clone, Copy)]
struct PN { from: u32, to: u32, k0: u32, k1: u32 }
impl PN {
    fn is_sorted(&self) -> bool { self.to == UNSET }
    fn key(&self) -> (u32, u32) { (self.k0, self.k1) }
}

fn next_maximal_set(nodes: &[PN], range: (usize, usize)) -> (usize, usize) {
    let n = nodes.len();
    if range.1 >= n { return (0, 0); }
    let r0 = range.1;
    let mut r1 = r0 + 1;
    if r0 > 0 && nodes[r0 - 1].key() == nodes[r0].key() { return (r0, r1); }
    for i in r1..n {
        if nodes[i - 1].key() != nodes[i].key() { r1 = i; }
        if nodes[i].from != nodes[r0].from { return (r0, r1); }
    }
    (r0, n)
}

fn merge_update_rank(nodes: &mut Vec<PN>, generation: u32) -> u32 {
    if generation == 4 {
        let mut curr = 0usize;
        let mut range = (0usize, 0usize);
        loop {
            range = next_maximal_set(nodes, range);
            if range.0 >= range.1 { break; }
            nodes[curr] = nodes[range.0];
            curr += 1;
        }
        nodes.truncate(curr);
        let mut candidate = Some(0usize);
        let mut key = nodes[0].key();
        for i in 1..nodes.len() {
            if nodes[i].key() != key {
                if let Some(c) = candidate { nodes[c].to = UNSET; }
                candidate = Some(i); key = nodes[i].key();
            } else { candidate = None; }
        }
        if let Some(c) = candidate { nodes[c].to = UNSET; }
        let mut r = 0u32;
        let mut key = nodes[0].key();
        for i in 0..nodes.len() {
            let k = nodes[i].key();
            if k != key { key = k; r += 1; }
            nodes[i].k0 = r; nodes[i].k1 = 0;
        }
        r + 1
    } else {
        let n = nodes.len();
        let (mut block_start, mut curr, mut node) = (0usize, 0usize, 0usize);
        let mut r = 0u32;
        loop {
            node += 1;
            if node == n || nodes[node].k0 != nodes[block_start].k0 {
                if node - block_start == 1 {
                    nodes[block_start].k0 = r; r += 1;
                    nodes[curr] = nodes[block_start]; curr += 1;
                } else {
                    nodes[block_start..node].sort_by_key(|x| x.k1);
                    while block_start != node {
                        let mut shift = 1usize;
                        while block_start + shift != node
                              && nodes[block_start].key() == nodes[block_start + shift].key() { shift += 1; }
                        let merge = (block_start..block_start + shift)
                            .all(|i| nodes[i].from == nodes[block_start].from);
                        if !merge {
                            for i in block_start..block_start + shift {
                                nodes[i].k0 = r;
                                nodes[curr] = nodes[i]; curr += 1;
                            }
                            r += 1;
                        } else if curr == 0 || !nodes[curr - 1].is_sorted()
                                  || nodes[curr - 1].from != nodes[block_start].from {
                            nodes[block_start].to = UNSET;
                            nodes[block_start].k0 = r; r += 1;
                            nodes[curr] = nodes[block_start]; curr += 1;
                        }
                        block_start += shift;
                    }
                    if node == n { break; }
                    if node + 1 == n && nodes[curr - 1].is_sorted()
                       && nodes[node].from == nodes[curr - 1].from { break; }
                    if node + 1 < n && nodes[node].k0 != nodes[node + 1].k0
                       && nodes[curr - 1].is_sorted() && nodes[node].from == nodes[curr - 1].from {
                        node += 1;
                    }
                }
                block_start = node;
            }
            if node == n { break; }
        }
        nodes.truncate(curr);
        r
    }
}

/// Full pipeline: prefix doubling to convergence, `generateEdges`, then the row
/// order `nextRow` produces.
pub fn build_gfm(gnodes: &[(u8, u32)], gedges: &[(u32, u32)], last_node: u32,
                 text_len: u32, off_rate: u32) -> BuiltGfm
{
    let n_ref_nodes = gnodes.len();
    let code = |l: u8| -> u32 { match l { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
                                          b'Y' => 4, _ => 5 } };
    let mut nodes: Vec<PN> = gedges.iter()
        .map(|&(f, t)| PN { from: f, to: t, k0: code(gnodes[f as usize].0), k1: 0 })
        .collect();
    nodes.push(PN { from: last_node, to: last_node, k0: 5, k1: 0 });
    nodes.sort_by_key(|n| n.from);

    let mut generation = 0u32;
    loop {
        generation += 1;
        let from_sorted: Vec<PN> = if generation <= 4 { nodes.clone() }
                                   else { let mut v = nodes.clone(); v.sort_by_key(|n| n.from); v };
        let mut start = vec![0u32; n_ref_nodes + 1];
        for n in &from_sorted { start[n.from as usize + 1] += 1; }
        for i in 1..start.len() { start[i] += start[i - 1]; }

        let mut new: Vec<PN> = Vec::with_capacity(nodes.len() * 2);
        for node in &nodes {
            if generation > 4 && node.is_sorted() { new.push(*node); continue; }
            let (lo, hi) = (start[node.to as usize] as usize, start[node.to as usize + 1] as usize);
            for j in lo..hi {
                let o = &from_sorted[j];
                let (k0, k1) = if generation <= 3 {
                    let shift = 3u32 * (1u32 << (generation - 1));
                    ((node.k0 << shift) + o.k0, 0)
                } else { (node.k0, o.k0) };
                new.push(PN { from: node.from, to: o.to, k0, k1 });
            }
        }
        nodes = new;
        let ranks = if generation <= 3 { 0 } else {
            if generation == 4 { nodes.sort_by_key(|n| (n.k0, n.k1)); }
            merge_update_rank(&mut nodes, generation)
        };
        if generation > 3 && ranks as usize == nodes.len() { break; }
        if generation > 64 { panic!("doubling did not converge"); }
    }

    // generateEdges
    nodes.sort_by_key(|n| n.from);
    for n in nodes.iter_mut() { n.to = gnodes[n.from as usize].1; }
    let mut start = vec![0u32; n_ref_nodes + 1];
    for n in &nodes { start[n.from as usize + 1] += 1; }
    for i in 1..start.len() { start[i] += start[i - 1]; }
    let mut buckets: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 6];
    for &(f, t) in gedges {
        let li = match gnodes[f as usize].0 {
            b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, b'Y' => 4, b'Z' => 5, _ => panic!() };
        for j in start[t as usize] as usize..start[t as usize + 1] as usize {
            buckets[li].push((nodes[j].k0, f));
        }
    }
    for bk in buckets.iter_mut() { bk.sort_by_key(|&(r, _)| r); }
    nodes.sort_by_key(|n| n.k0);

    let mut ed: Vec<(u32, u32, u8)> = Vec::new();
    for (li, bk) in buckets.iter().enumerate() {
        let lab = b"ACGTYZ"[li];
        for &(rank, from) in bk { ed.push((rank, from, lab)); }
    }
    let mut outdeg = vec![0u32; nodes.len()];
    {
        let (mut ni, mut ei) = (0usize, 0usize);
        while ni < nodes.len() && ei < ed.len() {
            if ed[ei].1 == nodes[ni].from { ed[ei].1 = ni as u32; ei += 1; outdeg[ni] += 1; }
            else { ni += 1; if ni < nodes.len() { outdeg[ni] = 0; } }
        }
    }
    let n = nodes.len();
    outdeg[n - 1] = outdeg[n - 2];
    nodes[n - 2] = nodes[n - 1];
    outdeg[n - 2] = outdeg[n - 1];
    nodes.pop(); outdeg.pop();
    for e in ed.iter_mut() {
        if e.2 == b'Y' { e.2 = b'Z'; }
        else if e.0 as usize >= nodes.len() { e.0 -= 1; }
    }
    ed.sort_by_key(|&(r, _, _)| r);

    let mut mrun: Vec<(u8, u32)> = Vec::with_capacity(ed.len());
    for (i, nd) in nodes.iter().enumerate() {
        for k in 0..outdeg[i].max(1) { mrun.push((if k == 0 { 1 } else { 0 }, nd.to)); }
    }
    let mut rows: Vec<(u8, u8)> = Vec::with_capacity(ed.len());
    let mut i = 0usize;
    for node in 0..nodes.len() as u32 {
        let mut first = true;
        while i < ed.len() && ed[i].0 == node {
            rows.push((ed[i].2, if first { 1 } else { 0 }));
            first = false; i += 1;
        }
    }

    // fchr, zOffs and the SA sample fall out of one walk over the rows
    let mut fchr_c = [0u32; 4];
    let mut z_offs = Vec::new();
    let mut sa_sample = Vec::new();
    let off_mask: u32 = u32::MAX << off_rate;
    let mut m_occ = 0u32;
    for r in 0..rows.len() {
        let (m, pos) = mrun[r];
        match rows[r].0 {
            b'A' => fchr_c[0] += 1, b'C' => fchr_c[1] += 1,
            b'G' => fchr_c[2] += 1, b'T' => fchr_c[3] += 1,
            _ => z_offs.push(r as u32),
        }
        if m == 1 {
            if (m_occ & off_mask) == m_occ { sa_sample.push(pos); }
            m_occ += 1;
        }
    }
    let mut fchr = [0u32; 5];
    for i in 0..4 { fchr[i + 1] = fchr[i] + fchr_c[i]; }

    BuiltGfm {
        len: text_len,
        gbwt_len: rows.len() as u32,
        num_nodes: nodes.len() as u32,
        rows, mrun, fchr, z_offs, sa_sample,
    }
}
