//! Point-lookup micro-benchmark with pre-generated queries.
//! usage: cargo run --release --example pointbench -- <index> [n]
use rshash::{AnyIndex, KmerWord, RsHash};
use std::time::Instant;

fn run<W: KmerWord>(idx: &RsHash<W>, n: usize) {
    let pos = idx.random_text_windows(n, 7);
    let neg = idx.random_windows(n, 7);
    for (name, ks) in [("pos", &pos), ("neg", &neg)] {
        for _ in 0..3 {
            let t = Instant::now();
            let f = idx.lookup_many(ks);
            println!("{name}: {:.1} ns/kmer ({f})", t.elapsed().as_nanos() as f64 / ks.len() as f64);
        }
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = AnyIndex::load_file(&a[1]).unwrap();
    let n = a.get(2).map_or(2_000_000, |x| x.parse().unwrap());
    match &idx {
        AnyIndex::W64(i) => run(i, n),
        AnyIndex::W128(i) => run(i, n),
    }
}
