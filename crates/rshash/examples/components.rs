//! Micro-benchmark of the streaming lookup's components.
//! usage: cargo run --release --example components -- <index> <reads>
use rshash::{AnyIndex, KmerWord, RsHash, alphabet, input, word::Windows};
use std::time::Instant;

fn run<W: KmerWord>(idx: &RsHash<W>, qs: &[Vec<u8>]) {
    let l = idx.window();
    let lv = &idx.levels()[0];
    let t = Instant::now();
    let mut acc = W::ZERO;
    let mut n = 0u64;
    for q in qs {
        for (v, rc) in Windows::<W, _>::new(q.iter().map(|&c| alphabet::rank(c)), l) {
            acc = acc ^ v ^ rc;
            n += 1;
        }
    }
    let base = t.elapsed().as_nanos() as f64 / n as f64;
    println!("windows only:        {base:.2} ns/kmer ({acc:?})");

    let t = Instant::now();
    let mut s = 0u64;
    for q in qs {
        let mut st = None;
        for (v, rc) in Windows::<W, _>::new(q.iter().map(|&c| alphabet::rank(c)), l) {
            let m = match st {
                None => lv.mp.find(v, rc),
                Some(mut m) => {
                    lv.mp.update(v, rc, &mut m);
                    m
                }
            };
            st = Some(m);
            s = s.wrapping_add(m.value);
        }
    }
    println!("+ rolling minimiser: {:.2} ns/kmer ({s})", t.elapsed().as_nanos() as f64 / n as f64);

    let t = Instant::now();
    let mut s = 0u64;
    let mut last = u64::MAX;
    for q in qs {
        let mut st = None;
        for (v, rc) in Windows::<W, _>::new(q.iter().map(|&c| alphabet::rank(c)), l) {
            let m = match st {
                None => lv.mp.find(v, rc),
                Some(mut m) => {
                    lv.mp.update(v, rc, &mut m);
                    m
                }
            };
            st = Some(m);
            if m.value != last {
                s += lv.r.contains(m.value).is_some() as u64;
                last = m.value;
            }
        }
    }
    println!("+ R1.contains:       {:.2} ns/kmer ({s})", t.elapsed().as_nanos() as f64 / n as f64);

    let t = Instant::now();
    let mut s = 0u64;
    for q in qs {
        for (v, rc) in Windows::<W, _>::new(q.iter().map(|&c| alphabet::rank(c)), l) {
            let m = lv.mp.find(v, rc);
            s = s.wrapping_add(m.value);
        }
    }
    println!("find every window:   {:.2} ns/kmer ({s})", t.elapsed().as_nanos() as f64 / n as f64);

    let mut e = idx.streaming_lookup();
    let t = Instant::now();
    let mut s = 0;
    for q in qs {
        s += e.query(q);
    }
    println!("streaming lookup:    {:.2} ns/kmer ({s})", t.elapsed().as_nanos() as f64 / n as f64);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = AnyIndex::load_file(&a[1]).unwrap();
    let qs = input::read_sequences(&a[2]).unwrap();
    match &idx {
        AnyIndex::W64(i) => run(i, &qs),
        AnyIndex::W128(i) => run(i, &qs),
    }
}
