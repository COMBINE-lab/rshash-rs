//! End-to-end exactness: every query answer (point lookup, streaming lookup,
//! point locate, streaming locate) must equal a brute-force oracle, for all
//! combinations of levels, last-level variants, window sizes and shapes.

use rshash::query::Hit;
use rshash::text::Text;
use rshash::{KmerWord, Params, RsHash};
use std::collections::{BTreeSet, HashMap};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // splitmix64
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn random_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
}

fn rc_seq(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// A text mixing random sequences, repeats (to create frequent minimisers
/// and push k-mers to deeper levels), palindromes and tiny/empty sequences.
fn make_text(seed: u64) -> Vec<Vec<u8>> {
    let mut rng = Rng(seed);
    let mut seqs = Vec::new();
    for _ in 0..6 {
        let n = 50 + rng.below(400) as usize;
        seqs.push(random_seq(&mut rng, n));
    }
    // highly repetitive: a motif repeated with rare mutations
    let motif = random_seq(&mut rng, 23);
    for _ in 0..3 {
        let mut s = Vec::new();
        for _ in 0..30 {
            s.extend_from_slice(&motif);
        }
        for _ in 0..5 {
            let i = rng.below(s.len() as u64) as usize;
            s[i] = b"ACGT"[rng.below(4) as usize];
        }
        seqs.push(s);
    }
    // low-complexity
    seqs.push(b"A".repeat(120));
    seqs.push(b"AC".repeat(70));
    // reverse-complemented copy of a sequence (both strands present)
    let c = rc_seq(&seqs[0]);
    seqs.push(c);
    seqs.push(Vec::new());
    seqs.push(b"ACGTTGCA".to_vec());
    seqs
}

fn make_queries(seqs: &[Vec<u8>], seed: u64) -> Vec<Vec<u8>> {
    let mut rng = Rng(seed);
    let mut qs = Vec::new();
    for s in seqs {
        qs.push(s.clone());
        qs.push(rc_seq(s));
        // mutated copy
        let mut m = s.clone();
        for _ in 0..(s.len() / 40 + 1) {
            if !m.is_empty() {
                let i = rng.below(m.len() as u64) as usize;
                m[i] = b"ACGT"[rng.below(4) as usize];
            }
        }
        qs.push(m);
    }
    // chimera across sequence boundaries (must not match across them)
    let mut chim = seqs[0][seqs[0].len() - 20..].to_vec();
    chim.extend_from_slice(&seqs[1][..20]);
    qs.push(chim);
    for _ in 0..5 {
        qs.push(random_seq(&mut rng, 200));
    }
    qs
}

/// Brute-force oracle: text windows keyed per shape (or plain window value).
struct Oracle<W: KmerWord> {
    l: u32,
    /// (w_mask, start, end) per shape; empty = contiguous windows
    shapes: Vec<(W, u32, u32)>,
    maps: Vec<HashMap<W, Vec<u64>>>,
}

impl<W: KmerWord> Oracle<W> {
    fn new(idx: &RsHash<W>) -> Self {
        let text = idx.text();
        let g = idx.geometry();
        let l = g.window;
        let shapes: Vec<(W, u32, u32)> = match &g.shapes {
            None => vec![],
            Some(s) => s.shapes.iter().map(|a| (W::from_u128_trunc(a.w_mask), a.start, a.end)).collect(),
        };
        let mut maps = Vec::new();
        let n = text.num_sequences();
        if shapes.is_empty() {
            let mut m: HashMap<W, Vec<u64>> = HashMap::new();
            for i in 0..n {
                let (s, e) = text.sequence_bounds(i);
                if e - s < l as u64 {
                    continue;
                }
                for p in s..=e - l as u64 {
                    m.entry(text.window(p, l)).or_default().push(p);
                }
            }
            maps.push(m);
        } else {
            for &(wm, a, b) in &shapes {
                let mut m: HashMap<W, Vec<u64>> = HashMap::new();
                for i in 0..n {
                    let (s, e) = text.sequence_bounds(i);
                    if e - s < (b - a) as u64 {
                        continue;
                    }
                    // p + a >= s, p + b <= e
                    for sp in s..=e - (b - a) as u64 {
                        let p = sp - a as u64;
                        let win: W = text.window(p, l);
                        m.entry(win.pext(wm)).or_default().push(p);
                    }
                }
                maps.push(m);
            }
        }
        Self { l, shapes, maps }
    }

    fn hits(&self, v: W) -> BTreeSet<Hit> {
        let rc = v.rc(self.l);
        let mut out = BTreeSet::new();
        if self.shapes.is_empty() {
            for &p in self.maps[0].get(&v).into_iter().flatten() {
                out.insert(Hit { pos: p, forward: true });
            }
            for &p in self.maps[0].get(&rc).into_iter().flatten() {
                out.insert(Hit { pos: p, forward: false });
            }
        } else {
            for (i, &(wm, _, _)) in self.shapes.iter().enumerate() {
                for &p in self.maps[i].get(&v.pext(wm)).into_iter().flatten() {
                    out.insert(Hit { pos: p, forward: true });
                }
                for &p in self.maps[i].get(&rc.pext(wm)).into_iter().flatten() {
                    out.insert(Hit { pos: p, forward: false });
                }
            }
        }
        out
    }
}

fn windows_of<W: KmerWord>(q: &[u8], l: u32) -> Vec<W> {
    let ranks: Vec<u8> = q.iter().map(|&c| rshash::alphabet::rank(c)).collect();
    rshash::word::windows::<W>(&ranks, l).map(|x| x.0).collect()
}

fn check_index<W: KmerWord>(seqs: &[Vec<u8>], params: &Params, queries: &[Vec<u8>]) -> (u64, u64) {
    let text = Text::from_sequences(seqs);
    let idx = RsHash::<W>::build(text, params).unwrap_or_else(|e| panic!("{params:?}: {e}"));
    let oracle = Oracle::new(&idx);
    let l = idx.window();
    if std::env::var("RSHASH_TEST_VERBOSE").is_ok() {
        let occ: Vec<usize> = idx.levels().iter().map(|lv| lv.num_occurrences()).collect();
        let last: usize = idx.last_level().iter().map(|t| t.len()).sum();
        eprintln!(
            "L={} K={} levels={} t={:?} ht={} loc={} shapes={} occ={:?} last_kmers={} last_keys={}",
            l,
            idx.kernel(),
            params.levels,
            params.t,
            params.use_ht,
            params.loc,
            params.shapes.len(),
            occ,
            idx.stats().last_level_kmers,
            last
        );
    }
    let mut eng = idx.streaming_lookup();
    let mut leng = idx.streaming_locate();
    let (mut total, mut positive) = (0u64, 0u64);
    for q in queries {
        let wins = windows_of::<W>(q, l);
        let expect: Vec<BTreeSet<Hit>> = wins.iter().map(|&v| oracle.hits(v)).collect();
        // point lookup
        for (i, &v) in wins.iter().enumerate() {
            assert_eq!(
                idx.lookup(v),
                !expect[i].is_empty(),
                "point lookup {params:?} window {i} of {:?}",
                String::from_utf8_lossy(q)
            );
        }
        // streaming lookup
        let mut got = vec![false; wins.len()];
        let n = eng.query_with(q, |i, f| got[i] = f);
        for i in 0..wins.len() {
            assert_eq!(
                got[i],
                !expect[i].is_empty(),
                "streaming lookup {params:?} window {i} of {:?}",
                String::from_utf8_lossy(q)
            );
        }
        assert_eq!(n, expect.iter().filter(|e| !e.is_empty()).count() as u64);
        total += wins.len() as u64;
        positive += n;
        if params.loc {
            for (i, &v) in wins.iter().enumerate() {
                let mut out = Vec::new();
                idx.locate(v, &mut out);
                let got: BTreeSet<Hit> = out.iter().copied().collect();
                assert_eq!(got.len(), out.len(), "duplicate hits");
                assert_eq!(got, expect[i], "point locate {params:?} window {i}");
            }
            leng.query(q, |i, hits| {
                let got: BTreeSet<Hit> = hits.iter().copied().collect();
                assert_eq!(got, expect[i], "streaming locate {params:?} window {i}");
            });
        }
    }
    // serialisation round trip: identical bytes and answers
    let mut bytes = Vec::new();
    idx.write(&mut bytes).unwrap();
    let idx2 = RsHash::<W>::read(&mut &bytes[..]).unwrap();
    let mut bytes2 = Vec::new();
    idx2.write(&mut bytes2).unwrap();
    assert!(bytes == bytes2, "re-serialisation differs");
    let mut eng2 = idx2.streaming_lookup();
    let mut eng1 = idx.streaming_lookup();
    for q in queries {
        assert_eq!(eng1.query(q), eng2.query(q));
    }
    (total, positive)
}

fn grid<W: KmerWord>(ks: &[(u32, [u32; 3])], shape_sets: &[Vec<u64>], seed: u64) {
    let seqs = make_text(seed);
    let queries = make_queries(&seqs, seed + 1);
    for &(k, m) in ks {
        for shapes in shape_sets {
            for levels in 1..=3 {
                for &t in &[[64u16, 64, 64], [2, 3, 4]] {
                    for &(use_ht, loc) in &[(true, false), (true, true), (false, false), (false, true)] {
                        let params = Params { k, levels, m, t, threshold: 0, loc, use_ht, shapes: shapes.clone() };
                        if let Ok(g) = params.geometry() {
                            if g.window > W::MAX_BASES
                                || (W::BITS == 128 && g.window <= 32 && shapes.is_empty() && k != 32)
                            {
                                continue;
                            }
                        } else {
                            continue;
                        }
                        let (total, pos) = check_index::<W>(&seqs, &params, &queries);
                        assert!(pos > 0 && pos < total, "degenerate test {params:?}: {pos}/{total}");
                    }
                }
            }
        }
    }
}

#[test]
fn plain_u64() {
    grid::<u64>(&[(15, [7, 9, 11]), (31, [12, 15, 18]), (32, [16, 20, 31]), (9, [3, 4, 5])], &[vec![]], 1);
}

#[test]
fn plain_u128() {
    grid::<u128>(&[(33, [12, 15, 18]), (47, [15, 20, 25]), (63, [20, 25, 31]), (32, [16, 20, 31])], &[vec![]], 2);
}

#[test]
fn shapes_u64() {
    let sets = vec![
        vec![0b1011101101],                          // palindrome, K=3.. (short)
        vec![0x5FFF_FFFD],                           // benchmark palindrome, L=31
        vec![0b11110111011111],                      // non-palindromic
        vec![0b1101111011, 0b110111111, 0b11111101], // three mixed shapes, different runs
        vec![0b111111111111],                        // contiguous shape
    ];
    grid::<u64>(&[(0, [3, 4, 5])], &sets, 3);
}

#[test]
fn shapes_u128() {
    let sets = vec![
        vec![0x7FFF_FFF9],                                         // benchmark shape: L=34
        vec![(1u64 << 40) | 0x7FFF_FFFF],                          // long span, weight 32
        vec![0b1011_1111_1111_1111_1111_1111_1111_1111_1101_0111], // non-palindromic long
        vec![0x5FFF_FFFD, 0x7FFF_FFFB],
    ];
    grid::<u128>(&[(0, [7, 9, 11])], &sets, 4);
}

#[test]
fn threshold_drops_frequent_last_level_kmers() {
    let seqs = make_text(9);
    let text = Text::from_sequences(&seqs);
    let base = Params { k: 15, levels: 1, m: [7, 9, 11], t: [1, 1, 1], ..Params::default() };
    let full = RsHash::<u64>::build(text.clone(), &base).unwrap();
    let filt = RsHash::<u64>::build(text, &Params { threshold: 2, ..base.clone() }).unwrap();
    let a: usize = full.last_level().iter().map(|t| t.len()).sum();
    let b: usize = filt.last_level().iter().map(|t| t.len()).sum();
    assert!(b < a, "threshold should drop repeated last-level k-mers ({b} vs {a})");
    // k-mers occurring once in the text are still found
    let oracle = Oracle::new(&filt);
    for ps in oracle.maps[0].values() {
        if ps.len() == 1 {
            let v: u64 = filt.text().window(ps[0], 15);
            let rc_count = oracle.maps[0].get(&v.rc(15)).map_or(0, |x| x.len());
            if rc_count == 0 {
                assert!(filt.lookup(v));
            }
        }
    }
}

/// Larger text (multi-block rank/select, wide offsets) with default-like
/// parameters, plus planted repeats.
#[test]
fn larger_text() {
    let mut rng = Rng(77);
    let mut seqs: Vec<Vec<u8>> = (0..40)
        .map(|_| {
            let n = 1000 + rng.below(20000) as usize;
            random_seq(&mut rng, n)
        })
        .collect();
    let rep = random_seq(&mut rng, 300);
    for s in seqs.iter_mut().take(20) {
        let at = rng.below((s.len() - 300) as u64) as usize;
        s[at..at + 300].copy_from_slice(&rep);
    }
    let mut queries: Vec<Vec<u8>> = Vec::new();
    for s in seqs.iter().take(10) {
        let a = rng.below((s.len() - 500) as u64) as usize;
        let mut q = s[a..a + 500].to_vec();
        q[250] = b'A';
        queries.push(q.clone());
        queries.push(rc_seq(&q));
    }
    queries.push(rep.clone());
    for _ in 0..10 {
        queries.push(random_seq(&mut rng, 300));
    }
    for (params, _) in [
        (Params { k: 31, levels: 2, m: [16, 19, 23], t: [8, 16, 64], ..Params::default() }, 0),
        (Params { k: 31, levels: 3, m: [13, 15, 17], t: [4, 4, 8], loc: true, ..Params::default() }, 0),
        (
            Params {
                k: 0,
                levels: 2,
                m: [16, 19, 23],
                t: [8, 16, 64],
                loc: true,
                use_ht: false,
                shapes: vec![0x5FFF_FFFD, 0x1FFF_FFFF],
                ..Params::default()
            },
            0,
        ),
    ] {
        let (total, pos) = check_index::<u64>(&seqs, &params, &queries);
        assert!(pos > 0 && pos < total);
    }
    let params = Params { k: 45, levels: 2, m: [19, 23, 23], t: [8, 16, 64], loc: true, ..Params::default() };
    check_index::<u128>(&seqs, &params, &queries);
}
