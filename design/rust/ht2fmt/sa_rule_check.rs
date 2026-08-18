use std::convert::TryInto;
use std::{env, fs};
fn main(){
    let a:Vec<String>=env::args().collect();
    let tb=fs::read(&a[1]).unwrap();
    let th:Vec<u32>=tb.chunks(4).map(|c| u32::from_le_bytes(c.try_into().unwrap())).collect();
    let mut t=Vec::new();
    for ln in fs::read_to_string(&a[2]).unwrap().lines(){
        if ln.starts_with('>'){continue}
        for c in ln.bytes(){ match c.to_ascii_uppercase(){b'A'=>t.push(0u8),b'C'=>t.push(1),b'G'=>t.push(2),b'T'=>t.push(3),_=>{}} } }
    let n=t.len(); let m=n+1;
    // HISAT2 rule: pad past the end with a character LARGER than any real one,
    // so a suffix that is a prefix of another sorts AFTER it, and the empty
    // suffix sorts last of all.
    const HI: u32 = u32::MAX;
    let mut sa:Vec<u32>=(0..m as u32).collect();
    let mut rank:Vec<u32>=(0..m).map(|i| if i==n {HI} else {t[i] as u32}).collect();
    let mut tmp=vec![0u32;m]; let mut k=1usize;
    loop{
        let key=|rk:&Vec<u32>,i:u32,k:usize|->(u32,u32){
            let a=rk[i as usize]; let j=i as usize+k;
            (a, if j<m {rk[j]} else {HI}) };
        sa.sort_unstable_by_key(|&i| key(&rank,i,k));
        tmp[sa[0] as usize]=0; let mut cl=0u32;
        for w in 1..m { if key(&rank,sa[w],k)!=key(&rank,sa[w-1],k){cl+=1;} tmp[sa[w] as usize]=cl; }
        std::mem::swap(&mut rank,&mut tmp);
        if cl as usize==m-1 {break} k<<=1; if k>=m {break}
    }
    let mut bad=0usize; let mut first=None;
    for r in 0..m { if sa[r]!=th[r] { bad+=1; if first.is_none(){first=Some(r);} } }
    println!("our SA (larger-than-all padding) vs HISAT2's full SA:");
    println!("  mismatches: {} of {}", bad, m);
    if bad==0 { println!("\n  EXACT MATCH -- HISAT2's suffix-array order is reproduced"); }
    else if let Some(r)=first {
        println!("  first at row {}: ours {} theirs {}", r, sa[r], th[r]);
        for i in r.saturating_sub(2)..(r+4).min(m) {
            println!("    row {:>7} ours {:>7} theirs {:>7}", i, sa[i], th[i]); }
    }
}
