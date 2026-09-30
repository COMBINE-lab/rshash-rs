//! Packed 2-bit k-mer words (port of `util.hpp`, the rolling updates of
//! `kmerview`/`longkmerview`, and the BMI2 bit tricks of `shape.hpp`).
//!
//! Encoding (seqan3 `dna4`): A=0, C=1, G=2, T=3, complement = `rank ^ 3`.
//! Base `i` of a window sits at bits `2i..2i+2`, i.e. the first base is in the
//! lowest bits. Windows of up to 32 bases use `u64`, up to 64 bases `u128`.

use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{BitAnd, BitOr, BitXor, Not, Shl, Shr};

/// `compute_mask(n)` from `util.hpp`: the `n` lowest bits set (n may equal 64).
#[inline(always)]
pub const fn mask64(nbits: u32) -> u64 {
    if nbits >= 64 { u64::MAX } else { (1u64 << nbits) - 1 }
}

#[inline(always)]
pub const fn mask128(nbits: u32) -> u128 {
    if nbits >= 128 { u128::MAX } else { (1u128 << nbits) - 1 }
}

/// Software `pext` (parallel bit extract), used when BMI2 is unavailable.
#[inline]
pub fn pext64_soft(x: u64, mut mask: u64) -> u64 {
    let mut res = 0u64;
    let mut bb = 1u64;
    while mask != 0 {
        let low = mask & mask.wrapping_neg();
        if x & low != 0 {
            res |= bb;
        }
        mask ^= low;
        bb <<= 1;
    }
    res
}

/// Software `pdep` (parallel bit deposit).
#[inline]
pub fn pdep64_soft(x: u64, mut mask: u64) -> u64 {
    let mut res = 0u64;
    let mut bb = 1u64;
    while mask != 0 {
        let low = mask & mask.wrapping_neg();
        if x & bb != 0 {
            res |= low;
        }
        mask ^= low;
        bb <<= 1;
    }
    res
}

/// `_pext_u64`.
#[inline(always)]
pub fn pext64(x: u64, mask: u64) -> u64 {
    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
    {
        // SAFETY: guarded by the `bmi2` target feature.
        unsafe { core::arch::x86_64::_pext_u64(x, mask) }
    }
    #[cfg(not(all(target_arch = "x86_64", target_feature = "bmi2")))]
    {
        pext64_soft(x, mask)
    }
}

/// `_pdep_u64`.
#[inline(always)]
pub fn pdep64(x: u64, mask: u64) -> u64 {
    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
    {
        // SAFETY: guarded by the `bmi2` target feature.
        unsafe { core::arch::x86_64::_pdep_u64(x, mask) }
    }
    #[cfg(not(all(target_arch = "x86_64", target_feature = "bmi2")))]
    {
        pdep64_soft(x, mask)
    }
}

/// Reverse complement of a packed k-mer of `k <= 32` bases (`crc` in `util.hpp`).
#[inline(always)]
pub fn rc64(x: u64, k: u32) -> u64 {
    debug_assert!((1..=32).contains(&k));
    let c = !x;
    let mut res = c.swap_bytes();
    const C1: u64 = 0x0f0f_0f0f_0f0f_0f0f;
    const C2: u64 = 0x3333_3333_3333_3333;
    res = ((res & C1) << 4) | ((res & (C1 << 4)) >> 4);
    res = ((res & C2) << 2) | ((res & (C2 << 2)) >> 2);
    res >> (64 - 2 * k)
}

/// Reverse complement of a packed k-mer of `k <= 64` bases.
#[inline(always)]
pub fn rc128(x: u128, k: u32) -> u128 {
    debug_assert!((1..=64).contains(&k));
    let c = !x;
    let mut res = c.swap_bytes();
    const C1: u128 = 0x0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f;
    const C2: u128 = 0x3333_3333_3333_3333_3333_3333_3333_3333;
    res = ((res & C1) << 4) | ((res & (C1 << 4)) >> 4);
    res = ((res & C2) << 2) | ((res & (C2 << 2)) >> 2);
    res >> (128 - 2 * k)
}

/// A packed window of bases: `u64` (<= 32 bases) or `u128` (<= 64 bases).
///
/// Keys derived from windows (canonical k-mers, `pext`-ed gapped k-mers) use
/// the same type.
pub trait KmerWord:
    Copy
    + Clone
    + Eq
    + Ord
    + Hash
    + Debug
    + Default
    + Send
    + Sync
    + 'static
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + BitXor<Output = Self>
    + Not<Output = Self>
    + Shl<u32, Output = Self>
    + Shr<u32, Output = Self>
{
    /// Number of bits of the word.
    const BITS: u32;
    /// Maximum number of bases a window of this type can hold.
    const MAX_BASES: u32 = Self::BITS / 2;
    const ZERO: Self;
    const ONES: Self;

    fn from_u64(x: u64) -> Self;
    fn from_u128_trunc(x: u128) -> Self;
    /// The lowest 64 bits.
    fn low_u64(self) -> u64;
    fn to_u128(self) -> u128;
    /// Low `nbits` bits set; `nbits` may equal `BITS`.
    fn mask(nbits: u32) -> Self;
    /// Reverse complement of a `k`-base window.
    fn rc(self, k: u32) -> Self;
    /// Parallel bit extract.
    fn pext(self, mask: Self) -> Self;
    /// Parallel bit deposit.
    fn pdep(self, mask: Self) -> Self;
    fn count_ones(self) -> u32;
    /// Number of significant bits (`std::bit_width`).
    fn bit_width(self) -> u32 {
        Self::BITS - self.leading_zeros()
    }
    fn leading_zeros(self) -> u32;
}

impl KmerWord for u64 {
    const BITS: u32 = 64;
    const ZERO: Self = 0;
    const ONES: Self = u64::MAX;

    #[inline(always)]
    fn from_u64(x: u64) -> Self {
        x
    }
    #[inline(always)]
    fn from_u128_trunc(x: u128) -> Self {
        x as u64
    }
    #[inline(always)]
    fn low_u64(self) -> u64 {
        self
    }
    #[inline(always)]
    fn to_u128(self) -> u128 {
        self as u128
    }
    #[inline(always)]
    fn mask(nbits: u32) -> Self {
        mask64(nbits)
    }
    #[inline(always)]
    fn rc(self, k: u32) -> Self {
        rc64(self, k)
    }
    #[inline(always)]
    fn pext(self, mask: Self) -> Self {
        pext64(self, mask)
    }
    #[inline(always)]
    fn pdep(self, mask: Self) -> Self {
        pdep64(self, mask)
    }
    #[inline(always)]
    fn count_ones(self) -> u32 {
        u64::count_ones(self)
    }
    #[inline(always)]
    fn leading_zeros(self) -> u32 {
        u64::leading_zeros(self)
    }
}

impl KmerWord for u128 {
    const BITS: u32 = 128;
    const ZERO: Self = 0;
    const ONES: Self = u128::MAX;

    #[inline(always)]
    fn from_u64(x: u64) -> Self {
        x as u128
    }
    #[inline(always)]
    fn from_u128_trunc(x: u128) -> Self {
        x
    }
    #[inline(always)]
    fn low_u64(self) -> u64 {
        self as u64
    }
    #[inline(always)]
    fn to_u128(self) -> u128 {
        self
    }
    #[inline(always)]
    fn mask(nbits: u32) -> Self {
        mask128(nbits)
    }
    #[inline(always)]
    fn rc(self, k: u32) -> Self {
        rc128(self, k)
    }
    /// 128-bit pext as two 64-bit halves joined at `popcount(mask.lo)`
    /// (the `lo_weight` trick of `test/hashtable.cpp`, measured in bits).
    #[inline(always)]
    fn pext(self, mask: Self) -> Self {
        let mlo = mask as u64;
        let mhi = (mask >> 64) as u64;
        let lo = pext64(self as u64, mlo) as u128;
        if mhi == 0 {
            return lo;
        }
        let hi = pext64((self >> 64) as u64, mhi) as u128;
        lo | (hi << mlo.count_ones())
    }
    #[inline(always)]
    fn pdep(self, mask: Self) -> Self {
        let mlo = mask as u64;
        let mhi = (mask >> 64) as u64;
        let lo = pdep64(self as u64, mlo) as u128;
        let hi = pdep64((self >> mlo.count_ones()) as u64, mhi) as u128;
        lo | (hi << 64)
    }
    #[inline(always)]
    fn count_ones(self) -> u32 {
        u128::count_ones(self)
    }
    #[inline(always)]
    fn leading_zeros(self) -> u32 {
        u128::leading_zeros(self)
    }
}

/// Canonical form of a window: `min(fwd, rc)`.
#[inline(always)]
pub fn canonical<W: KmerWord>(fwd: W, rc: W) -> W {
    if fwd < rc { fwd } else { rc }
}

/// Encode an ASCII k-mer (test/debug helper). Non-ACGT characters map to A.
pub fn encode_kmer<W: KmerWord>(s: &[u8]) -> W {
    let mut v = W::ZERO;
    for (i, &c) in s.iter().enumerate() {
        v = v | (W::from_u64(crate::alphabet::rank(c) as u64) << (2 * i as u32));
    }
    v
}

/// Decode a packed window of `k` bases to ASCII (test/debug helper).
pub fn decode_kmer<W: KmerWord>(v: W, k: u32) -> String {
    (0..k).map(|i| b"ACGT"[((v >> (2 * i)).low_u64() & 3) as usize] as char).collect()
}

/// Rolling window over a base (rank) iterator, yielding `(value, value_rev)`
/// for each full window, exactly like `kmerview`/`longkmerview`:
///
/// ```text
/// value     = (value >> 2) | (rank << 2(w-1))
/// value_rev = ((value_rev << 2) | (rank ^ 3)) & mask(2w)
/// ```
pub struct Windows<W: KmerWord, I: Iterator<Item = u8>> {
    ranks: I,
    w: u32,
    shift: u32,
    mask: W,
    value: W,
    value_rev: W,
    filled: u32,
}

impl<W: KmerWord, I: Iterator<Item = u8>> Windows<W, I> {
    pub fn new(ranks: I, w: u32) -> Self {
        assert!(w >= 1 && w <= W::MAX_BASES, "window size {w} out of range");
        Self { ranks, w, shift: 2 * (w - 1), mask: W::mask(2 * w), value: W::ZERO, value_rev: W::ZERO, filled: 0 }
    }
}

impl<W: KmerWord, I: Iterator<Item = u8>> Iterator for Windows<W, I> {
    type Item = (W, W);

    #[inline]
    fn next(&mut self) -> Option<(W, W)> {
        loop {
            let r = self.ranks.next()?;
            self.value = (self.value >> 2) | (W::from_u64(r as u64) << self.shift);
            self.value_rev = ((self.value_rev << 2) | W::from_u64((r ^ 3) as u64)) & self.mask;
            if self.filled + 1 >= self.w {
                return Some((self.value, self.value_rev));
            }
            self.filled += 1;
        }
    }
}

/// Rolling windows over a slice of ranks.
pub fn windows<W: KmerWord>(ranks: &[u8], w: u32) -> Windows<W, std::iter::Copied<std::slice::Iter<'_, u8>>> {
    Windows::new(ranks.iter().copied(), w)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_rc(s: &str) -> String {
        s.bytes()
            .rev()
            .map(|c| match c {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                _ => 'A',
            })
            .collect()
    }

    fn lcg(state: &mut u64) -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *state >> 33
    }

    #[test]
    fn rc_matches_naive() {
        let mut st = 7u64;
        for k in 1..=64u32 {
            for _ in 0..20 {
                let s: String = (0..k).map(|_| b"ACGT"[(lcg(&mut st) & 3) as usize] as char).collect();
                let v128: u128 = encode_kmer(s.as_bytes());
                assert_eq!(decode_kmer(v128.rc(k), k), naive_rc(&s));
                if k <= 32 {
                    let v64: u64 = encode_kmer(s.as_bytes());
                    assert_eq!(decode_kmer(v64.rc(k), k), naive_rc(&s));
                }
            }
        }
    }

    #[test]
    fn pext_pdep_soft_matches_hw() {
        let mut st = 11u64;
        for _ in 0..10000 {
            let x = lcg(&mut st) ^ (lcg(&mut st) << 31);
            let m = lcg(&mut st) ^ (lcg(&mut st) << 29);
            assert_eq!(pext64(x, m), pext64_soft(x, m));
            assert_eq!(pdep64(x, m), pdep64_soft(x, m));
            let x128 = (x as u128) << 64 | (m as u128);
            let m128 = (m as u128) << 64 | (x as u128);
            // pdep(pext(x)) == x & mask
            assert_eq!(x128.pext(m128).pdep(m128), x128 & m128);
        }
    }

    #[test]
    fn rolling_windows() {
        let seq = b"ACGTTGCAAGGCTTAC";
        let ranks: Vec<u8> = seq.iter().map(|&c| crate::alphabet::rank(c)).collect();
        for w in 1..=seq.len() as u32 {
            let got: Vec<(u64, u64)> = windows::<u64>(&ranks, w).collect();
            assert_eq!(got.len(), seq.len() - w as usize + 1);
            for (i, (v, r)) in got.iter().enumerate() {
                let s = &seq[i..i + w as usize];
                assert_eq!(*v, encode_kmer::<u64>(s));
                assert_eq!(*r, v.rc(w));
            }
            let got128: Vec<(u128, u128)> = windows::<u128>(&ranks, w).collect();
            for (a, b) in got.iter().zip(got128.iter()) {
                assert_eq!(a.0 as u128, b.0);
                assert_eq!(a.1 as u128, b.1);
            }
        }
    }
}
