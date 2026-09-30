//! Elias–Fano representation of a sorted (non-decreasing) sequence, with the
//! operations RSHash needs from its local sux fork (`EliasFano.hpp`):
//! `contains(x, &rank)`, `rank`, `select`, `select(r, &next)`.
//!
//! Values may be `u64` or `u128` (the last level of an index over windows of
//! more than 32 bases, or of gapped k-mers of weight 32, needs keys wider than
//! 64 bits). Lower-bit widths above 64 are split over two packed vectors.

use crate::bits::{BitVec, CompactVec, RankSelect};
use crate::io::{Reader, Writer, invalid};
use crate::word::KmerWord;
use std::io;
use std::marker::PhantomData;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EliasFano<W: KmerWord> {
    n: usize,
    l: u32,
    lower_lo: CompactVec,
    lower_hi: CompactVec,
    upper: RankSelect,
    _w: PhantomData<W>,
}

/// `floor(log2(x))` for `x > 0`.
fn floor_log2(x: u128) -> u32 {
    127 - x.leading_zeros()
}

impl<W: KmerWord> EliasFano<W> {
    /// Build from sorted values, all `< universe`.
    ///
    /// As in sux, `l = floor(log2(universe / n))` (0 if the quotient is 0).
    pub fn new(values: &[W], universe: u128) -> Self {
        let n = values.len();
        debug_assert!(values.windows(2).all(|w| w[0] <= w[1]), "EF input must be sorted");
        let l = if n == 0 {
            0
        } else {
            let q = universe / n as u128;
            if q == 0 { 0 } else { floor_log2(q).min(W::BITS) }
        };
        let lo_w = l.min(64);
        let hi_w = l.saturating_sub(64);
        let mut lower_lo = CompactVec::new(n, lo_w);
        let mut lower_hi = CompactVec::new(if hi_w > 0 { n } else { 0 }, hi_w);
        let max_high = values.last().map(|&v| Self::high_of(v, l)).unwrap_or(0);
        let mut upper = BitVec::new(n + max_high as usize + 1);
        for (i, &v) in values.iter().enumerate() {
            let low = v.to_u128() & crate::word::mask128(l);
            lower_lo.set(i, low as u64 & crate::word::mask64(lo_w));
            if hi_w > 0 {
                lower_hi.set(i, (low >> 64) as u64);
            }
            upper.set(Self::high_of(v, l) as usize + i);
        }
        Self { n, l, lower_lo, lower_hi, upper: RankSelect::new(upper), _w: PhantomData }
    }

    #[inline(always)]
    fn high_of(v: W, l: u32) -> u64 {
        if l >= 128 { 0 } else { (v.to_u128() >> l) as u64 }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.n
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    #[inline(always)]
    fn lower(&self, i: usize) -> u128 {
        let lo = self.lower_lo.get(i) as u128;
        if self.l > 64 { lo | ((self.lower_hi.get(i) as u128) << 64) } else { lo }
    }

    #[inline(always)]
    fn compose(&self, high: u64, low: u128) -> W {
        let v = if self.l >= 128 { low } else { ((high as u128) << self.l) | low };
        W::from_u128_trunc(v)
    }

    /// The `i`-th value.
    #[inline]
    pub fn select(&self, i: usize) -> W {
        debug_assert!(i < self.n);
        let high = (self.upper.select1(i) - i) as u64;
        self.compose(high, self.lower(i))
    }

    /// `(select(i), select(i+1))`, like sux's `select(r, &next)`.
    #[inline]
    pub fn select_pair(&self, i: usize) -> (W, W) {
        debug_assert!(i + 1 < self.n);
        let p = self.upper.select1(i);
        let a = self.compose((p - i) as u64, self.lower(i));
        // next one after p: usually within the same word
        let w = p / 64;
        let rest = self.upper.bits().words()[w] & !crate::word::mask64((p % 64) as u32 + 1);
        let q = if rest != 0 { w * 64 + rest.trailing_zeros() as usize } else { self.upper.select1(i + 1) };
        let b = self.compose((q - i - 1) as u64, self.lower(i + 1));
        (a, b)
    }

    /// First position in `upper` and first index of the values whose high
    /// part equals `h` (one `select0`, as in the sux fork's `contains`).
    #[inline(always)]
    fn bucket_start(&self, h: u64) -> (usize, usize) {
        if h == 0 {
            (0, 0)
        } else {
            let z = self.upper.select0(h as usize - 1) + 1;
            (z, z - h as usize)
        }
    }

    /// Index of the first value `>= x` if it lies in the high bucket of `x`
    /// (`Ok`), or the index after the bucket (`Err`).
    #[inline(always)]
    fn lower_bound_in_bucket(&self, x: W) -> Result<usize, usize> {
        if self.n == 0 {
            return Err(0);
        }
        let h = Self::high_of(x, self.l);
        let max_high = (self.upper.len() - self.n - 1) as u64;
        if h > max_high {
            return Err(self.n);
        }
        let (mut pos, mut i) = self.bucket_start(h);
        let target = x.to_u128() & crate::word::mask128(self.l);
        let words = self.upper.bits().words();
        loop {
            // bucket ends at the next zero
            if (words[pos / 64] >> (pos % 64)) & 1 == 0 {
                return Err(i);
            }
            if self.lower(i) >= target {
                return Ok(i);
            }
            pos += 1;
            i += 1;
        }
    }

    /// If `x` is present, its (first) index, i.e. `contains(x, rank_out)`.
    #[inline]
    pub fn contains(&self, x: W) -> Option<usize> {
        match self.lower_bound_in_bucket(x) {
            Ok(i) if self.lower(i) == x.to_u128() & crate::word::mask128(self.l) => Some(i),
            _ => None,
        }
    }

    /// Number of values `< x`.
    #[inline]
    pub fn rank(&self, x: W) -> usize {
        match self.lower_bound_in_bucket(x) {
            Ok(i) | Err(i) => i,
        }
    }

    /// Predecessor query: the largest index `i` with `select(i) <= x`, with
    /// `select(i)` and `select(i + 1)` (the latter `None` for the last value).
    /// One `select0` plus local bit scans (used to find the sequence
    /// containing a text position).
    #[inline]
    pub fn pred_pair(&self, x: W) -> Option<(usize, W, Option<W>)> {
        if self.n == 0 {
            return None;
        }
        let h = Self::high_of(x, self.l);
        let max_high = (self.upper.len() - self.n - 1) as u64;
        let target = x.to_u128() & crate::word::mask128(self.l);
        let words = self.upper.bits().words();
        let bit = |p: usize| (words[p / 64] >> (p % 64)) & 1 == 1;
        // (upper position, index) of the predecessor
        let (ppos, pi) = if h > max_high {
            let p = self.upper.prev_one(self.upper.len())?;
            (p, self.n - 1)
        } else {
            let (mut pos, mut i) = self.bucket_start(h);
            let mut best: Option<(usize, usize)> = None;
            while bit(pos) && self.lower(i) <= target {
                best = Some((pos, i));
                pos += 1;
                i += 1;
            }
            match best {
                Some(b) => b,
                None => {
                    // predecessor is the last value of an earlier bucket
                    let (bpos, bi) = self.bucket_start(h);
                    if bi == 0 {
                        return None;
                    }
                    (self.upper.prev_one(bpos)?, bi - 1)
                }
            }
        };
        let v = self.compose((ppos - pi) as u64, self.lower(pi));
        let next = if pi + 1 < self.n {
            let q = self.upper.next_one(ppos + 1);
            Some(self.compose((q - pi - 1) as u64, self.lower(pi + 1)))
        } else {
            None
        };
        Some((pi, v, next))
    }

    /// Size in bits.
    pub fn bit_size(&self) -> usize {
        self.lower_lo.bit_size() + self.lower_hi.bit_size() + self.upper.bit_size()
    }

    pub fn iter(&self) -> impl Iterator<Item = W> + '_ {
        (0..self.n).map(|i| self.select(i))
    }

    pub fn write(&self, w: &mut Writer) -> io::Result<()> {
        w.u32(W::BITS)?;
        w.u64(self.n as u64)?;
        w.u32(self.l)?;
        self.lower_lo.write(w)?;
        self.lower_hi.write(w)?;
        self.upper.write(w)
    }

    pub fn read(r: &mut Reader) -> io::Result<Self> {
        if r.u32()? != W::BITS {
            return Err(invalid("Elias-Fano word size mismatch"));
        }
        let n = r.u64()? as usize;
        let l = r.u32()?;
        let lower_lo = CompactVec::read(r)?;
        let lower_hi = CompactVec::read(r)?;
        let upper = RankSelect::read(r)?;
        if lower_lo.len() != n || upper.count_ones() != n {
            return Err(invalid("Elias-Fano length mismatch"));
        }
        Ok(Self { n, l, lower_lo, lower_hi, upper, _w: PhantomData })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(state: &mut u64) -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *state >> 33
    }

    fn check<W: KmerWord>(vals: &[W], universe: u128, probes: &[W]) {
        let ef = EliasFano::<W>::new(vals, universe);
        assert_eq!(ef.len(), vals.len());
        for (i, &v) in vals.iter().enumerate() {
            assert_eq!(ef.select(i), v);
            if i + 1 < vals.len() {
                assert_eq!(ef.select_pair(i), (v, vals[i + 1]));
            }
        }
        for &x in vals.iter().chain(probes) {
            // predecessor (last index with v <= x)
            let pp = vals.partition_point(|&v| v <= x);
            let expect = if pp == 0 { None } else { Some((pp - 1, vals[pp - 1], vals.get(pp).copied())) };
            assert_eq!(ef.pred_pair(x), expect, "pred_pair({x:?})");
            let rank = vals.partition_point(|&v| v < x);
            assert_eq!(ef.rank(x), rank, "rank({x:?})");
            let expect = if rank < vals.len() && vals[rank] == x { Some(rank) } else { None };
            assert_eq!(ef.contains(x), expect, "contains({x:?})");
        }
    }

    #[test]
    fn ef_u64() {
        let mut st = 1u64;
        for &(n, ubits) in &[(0usize, 10u32), (1, 10), (10, 4), (1000, 36), (1000, 12), (5000, 62), (3000, 20)] {
            let u = 1u128 << ubits;
            let mut vals: Vec<u64> = (0..n).map(|_| (lcg(&mut st) << 31 ^ lcg(&mut st)) % u as u64).collect();
            vals.sort();
            let probes: Vec<u64> = (0..500).map(|_| (lcg(&mut st) << 31 ^ lcg(&mut st)) % u as u64).collect();
            check(&vals, u, &probes);
        }
    }

    #[test]
    fn ef_u64_full_universe() {
        let mut vals = vec![0u64, 5, 5, 5, u64::MAX - 1, u64::MAX];
        vals.sort();
        check(&vals, 1u128 << 64, &[1, 4, 6, u64::MAX - 2]);
    }

    #[test]
    fn ef_u128() {
        let mut st = 9u64;
        for &(n, ubits) in &[(1usize, 126u32), (10, 126), (1000, 126), (1000, 70), (2000, 100), (100, 128)] {
            let u = if ubits == 128 { u128::MAX } else { 1u128 << ubits };
            let rnd = |st: &mut u64| -> u128 {
                let a = (lcg(st) as u128) << 96 | (lcg(st) as u128) << 64 | (lcg(st) as u128) << 32 | lcg(st) as u128;
                a % u
            };
            let mut vals: Vec<u128> = (0..n).map(|_| rnd(&mut st)).collect();
            vals.sort();
            let probes: Vec<u128> = (0..300).map(|_| rnd(&mut st)).collect();
            check(&vals, u, &probes);
        }
    }
}
