//! Small deterministic helpers shared by the CLI, tools and tests.

use crate::index::RsHash;
use crate::word::KmerWord;

/// SplitMix64 pseudo-random generator (deterministic, dependency-free).
#[derive(Clone, Debug)]
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[inline]
    pub fn next_u128(&mut self) -> u128 {
        ((self.next_u64() as u128) << 64) | self.next_u64() as u128
    }

    /// Uniform-ish value in `[0, n)`.
    #[inline]
    pub fn below(&mut self, n: u64) -> u64 {
        ((self.next_u64() as u128 * n as u128) >> 64) as u64
    }
}

impl<W: KmerWord> RsHash<W> {
    /// `rand_text_kmers`: `n` windows sampled from the text (inside one
    /// sequence); every other window is reverse-complemented.
    pub fn random_text_windows(&self, n: usize, seed: u64) -> Vec<W> {
        let mut rng = SplitMix64(seed);
        let text = self.text();
        let l = self.window() as u64;
        let mut out = Vec::with_capacity(n);
        if text.len() < l {
            return out;
        }
        let mut tries = 0usize;
        while out.len() < n && tries < 100 * n + 1000 {
            tries += 1;
            let p = crate::text::TEXT_START + rng.below(text.len() - l + 1);
            let (_, _, e) = text.locate(p);
            if p + l > e {
                continue;
            }
            let v: W = text.window(p, l as u32);
            out.push(if out.len() % 2 == 0 { v.rc(l as u32) } else { v });
        }
        out
    }

    /// `rand_kmers`: `n` uniformly random windows.
    pub fn random_windows(&self, n: usize, seed: u64) -> Vec<W> {
        let mut rng = SplitMix64(seed);
        let mask = W::mask(2 * self.window());
        (0..n).map(|_| W::from_u128_trunc(rng.next_u128()) & mask).collect()
    }
}
