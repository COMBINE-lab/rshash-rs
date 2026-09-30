//! Succinct building blocks: a plain bit vector with rank/select support
//! (replaces sdsl `bit_vector` + sux `SimpleSelect`, and the select structures
//! inside sux's `EliasFano`) and a bit-packed integer vector (port of the
//! pthash-style `bits::compact_vector`).

use crate::io::{Reader, Writer};
use crate::word::{mask64, pdep64};
use std::io;

/// Position of the `r`-th (0-based) set bit of `w`.
#[inline(always)]
pub fn select_in_word(w: u64, r: u32) -> u32 {
    debug_assert!(r < w.count_ones());
    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
    {
        pdep64(1u64 << r, w).trailing_zeros()
    }
    #[cfg(not(all(target_arch = "x86_64", target_feature = "bmi2")))]
    {
        let _ = pdep64;
        let mut w = w;
        for _ in 0..r {
            w &= w - 1;
        }
        w.trailing_zeros()
    }
}

/// Growable bit vector.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BitVec {
    words: Vec<u64>,
    len: usize,
}

impl BitVec {
    pub fn new(len: usize) -> Self {
        Self { words: vec![0; len.div_ceil(64)], len }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> bool {
        debug_assert!(i < self.len);
        (self.words[i / 64] >> (i % 64)) & 1 == 1
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize) {
        debug_assert!(i < self.len);
        self.words[i / 64] |= 1u64 << (i % 64);
    }

    pub fn words(&self) -> &[u64] {
        &self.words
    }

    pub fn count_ones(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    /// Size in bits of the payload.
    pub fn bit_size(&self) -> usize {
        self.words.len() * 64
    }

    pub fn write(&self, w: &mut Writer) -> io::Result<()> {
        w.u64(self.len as u64)?;
        w.vec_u64(&self.words)
    }

    pub fn read(r: &mut Reader) -> io::Result<Self> {
        let len = r.u64()? as usize;
        let words = r.vec_u64()?;
        if words.len() != len.div_ceil(64) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "bit vector length mismatch"));
        }
        Ok(Self { words, len })
    }
}

/// 1 in the low bit of each of seven 9-bit fields (rank9 layout).
const ONES9: u64 = 0x0040_2010_0804_0201;
const MSBS9: u64 = ONES9 << 8;
/// `64 * j` in field `j-1`, for `j = 1..=7`: zeros before word `j` = this - ones.
const WORDBITS9: u64 = 64 | 128 << 9 | 192 << 18 | 256 << 27 | 320 << 36 | 384 << 45 | 448 << 54;

/// Fieldwise `x <= y` on seven 9-bit fields (Vigna's `ULEQ_STEP_9`).
#[inline(always)]
fn uleq9(x: u64, y: u64) -> u64 {
    ((((y | MSBS9).wrapping_sub(x & !MSBS9)) | (x ^ y)) ^ (x & !y)) & MSBS9
}

const BLOCK_BITS: usize = 512;
const BLOCK_WORDS: usize = BLOCK_BITS / 64;
const SAMPLE: usize = 512;

/// Bit vector with constant-time rank and (practically) constant-time
/// select1/select0.
///
/// Layout: cumulative popcounts every 512 bits plus a sample of the block
/// containing every 512-th one (resp. zero). `select` jumps to the sampled
/// block, binary-searches the blocks up to the next sample, then scans at most
/// eight words.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RankSelect {
    bits: BitVec,
    /// Ones before block `b`, for `b in 0..=num_blocks`.
    block_ones: Vec<u64>,
    /// Per block, the ones before word `j` (1..=7) of the block, 9 bits each
    /// at bit `9*(j-1)` (rank9 layout).
    sub: Vec<u64>,
    /// Block containing the `(i*SAMPLE)`-th one.
    sel1: Vec<u32>,
    /// Block containing the `(i*SAMPLE)`-th zero.
    sel0: Vec<u32>,
    /// Exact-position samples of the ones.
    pos1: Samples,
    /// Exact-position samples of the zeros.
    pos0: Samples,
    ones: usize,
}

/// Sampling rate of exact positions for the fast select path.
const POS_SAMPLE: usize = 256;
/// Secondary sampling rate inside long spans.
const SUB_SAMPLE: usize = 32;
/// Max words between two samples for the linear-scan fast path.
const MAX_SCAN_WORDS: u64 = 16;

/// Two-level exact position samples (in the spirit of sux's `SimpleSelect`):
/// the position of every 256th target bit, plus the offsets of every 32nd
/// target bit inside spans longer than [`MAX_SCAN_WORDS`] words.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Samples {
    /// Position of every `POS_SAMPLE`-th target bit, then a sentinel (`len`).
    pos: Vec<u64>,
    /// Per primary sample: start in `sub`, or `u32::MAX` for short spans.
    sub_idx: Vec<u32>,
    /// Offsets from the primary sample of every `SUB_SAMPLE`-th target bit.
    sub: Vec<u32>,
}

impl Samples {
    fn build(bits: &BitVec, ones: bool) -> Self {
        let mut pos = Vec::new();
        let mut chunk: Vec<u64> = Vec::with_capacity(POS_SAMPLE / SUB_SAMPLE);
        let mut chunks: Vec<Vec<u64>> = Vec::new();
        let mut c = 0usize;
        for (wi, &w) in bits.words.iter().enumerate() {
            let valid = if (wi + 1) * 64 <= bits.len() { 64 } else { (bits.len() - wi * 64) as u32 };
            let mut x = if ones { w } else { !w } & mask64(valid);
            while x != 0 {
                let p = (wi * 64) as u64 + x.trailing_zeros() as u64;
                x &= x - 1;
                if c.is_multiple_of(POS_SAMPLE) {
                    if c > 0 {
                        chunks.push(std::mem::take(&mut chunk));
                    }
                    pos.push(p);
                }
                if c.is_multiple_of(SUB_SAMPLE) {
                    chunk.push(p);
                }
                c += 1;
            }
        }
        if c > 0 {
            chunks.push(chunk);
        }
        pos.push(bits.len() as u64);
        let mut sub_idx = Vec::with_capacity(chunks.len());
        let mut sub = Vec::new();
        for (s, ch) in chunks.iter().enumerate() {
            if (pos[s + 1] >> 6) - (pos[s] >> 6) > MAX_SCAN_WORDS && pos[s + 1] - pos[s] < u32::MAX as u64 {
                sub_idx.push(sub.len() as u32);
                // always POS_SAMPLE / SUB_SAMPLE entries (the last chunk is padded
                // with the span end)
                sub.extend(ch.iter().map(|&p| (p - pos[s]) as u32));
                sub.extend(std::iter::repeat_n((pos[s + 1] - pos[s]) as u32, POS_SAMPLE / SUB_SAMPLE - ch.len()));
            } else {
                sub_idx.push(u32::MAX);
            }
        }
        Self { pos, sub_idx, sub }
    }

    /// Start position and remaining count for a linear scan to the `r`-th
    /// target bit, or `None` if the span is too long (use the directory).
    #[inline(always)]
    fn start(&self, r: usize) -> Option<(u64, u32)> {
        let s = r / POS_SAMPLE;
        let p0 = self.pos[s];
        let p1 = self.pos[s + 1];
        if (p1 >> 6) - (p0 >> 6) <= MAX_SCAN_WORDS {
            return Some((p0, (r % POS_SAMPLE) as u32));
        }
        let si = self.sub_idx[s];
        if si == u32::MAX {
            return None;
        }
        let k = (r % POS_SAMPLE) / SUB_SAMPLE;
        let q0 = p0 + self.sub[si as usize + k] as u64;
        let q1 = if k + 1 < POS_SAMPLE / SUB_SAMPLE { p0 + self.sub[si as usize + k + 1] as u64 } else { p1 };
        if (q1 >> 6) - (q0 >> 6) <= MAX_SCAN_WORDS { Some((q0, (r % SUB_SAMPLE) as u32)) } else { None }
    }

    fn bit_size(&self) -> usize {
        64 * self.pos.len() + 32 * (self.sub_idx.len() + self.sub.len())
    }
}

impl RankSelect {
    pub fn new(bits: BitVec) -> Self {
        let nblocks = bits.len().div_ceil(BLOCK_BITS);
        let mut block_ones = Vec::with_capacity(nblocks + 1);
        let mut sub = Vec::with_capacity(nblocks);
        let mut sel1 = Vec::new();
        let mut sel0 = Vec::new();
        let mut ones = 0usize;
        for b in 0..nblocks {
            block_ones.push(ones as u64);
            let zeros_before = b * BLOCK_BITS - ones;
            let lo = b * BLOCK_WORDS;
            let hi = ((b + 1) * BLOCK_WORDS).min(bits.words.len());
            let mut block_count = 0usize;
            let mut packed = 0u64;
            for (j, w) in bits.words[lo..hi].iter().enumerate() {
                if j > 0 {
                    packed |= (block_count as u64) << (9 * (j - 1));
                }
                block_count += w.count_ones() as usize;
            }
            for j in (hi - lo).max(1)..BLOCK_WORDS {
                packed |= (block_count as u64) << (9 * (j - 1));
            }
            sub.push(packed);
            let block_len = ((b + 1) * BLOCK_BITS).min(bits.len()) - b * BLOCK_BITS;
            let block_zeros = block_len - block_count;
            // samples whose target lies in this block
            while sel1.len() * SAMPLE < ones + block_count {
                sel1.push(b as u32);
            }
            while sel0.len() * SAMPLE < zeros_before + block_zeros {
                sel0.push(b as u32);
            }
            ones += block_count;
        }
        block_ones.push(ones as u64);
        sel1.push(nblocks as u32);
        sel0.push(nblocks as u32);
        let pos1 = Samples::build(&bits, true);
        let pos0 = Samples::build(&bits, false);
        Self { bits, block_ones, sub, sel1, sel0, pos1, pos0, ones }
    }

    #[inline(always)]
    pub fn bits(&self) -> &BitVec {
        &self.bits
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.bits.len()
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.bits.is_empty()
    }

    #[inline(always)]
    pub fn count_ones(&self) -> usize {
        self.ones
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> bool {
        self.bits.get(i)
    }

    /// Number of ones in `[0, pos)`.
    #[inline]
    pub fn rank1(&self, pos: usize) -> usize {
        debug_assert!(pos <= self.len());
        let b = pos / BLOCK_BITS;
        let mut r = self.block_ones[b] as usize;
        let wend = pos / 64;
        let j = wend - b * BLOCK_WORDS;
        if j > 0 {
            r += ((self.sub[b] >> (9 * (j - 1))) & 0x1FF) as usize;
        }
        if !pos.is_multiple_of(64) {
            r += (self.bits.words[wend] & mask64((pos % 64) as u32)).count_ones() as usize;
        }
        r
    }

    #[inline]
    pub fn rank0(&self, pos: usize) -> usize {
        pos - self.rank1(pos)
    }

    /// Position of the `r`-th one (0-based).
    #[inline]
    pub fn select1(&self, r: usize) -> usize {
        if let Some((p0, mut rem)) = self.pos1.start(r) {
            let mut wi = (p0 >> 6) as usize;
            let mut w = self.bits.words[wi] & !mask64((p0 & 63) as u32);
            loop {
                let c = w.count_ones();
                if rem < c {
                    return wi * 64 + select_in_word(w, rem) as usize;
                }
                rem -= c;
                wi += 1;
                w = self.bits.words[wi];
            }
        }
        self.select1_blocks(r)
    }

    #[inline(never)]
    fn select1_blocks(&self, r: usize) -> usize {
        debug_assert!(r < self.ones, "select1({r}) out of range ({})", self.ones);
        let s = r / SAMPLE;
        let mut lo = self.sel1[s] as usize;
        let mut hi = self.sel1[s + 1] as usize; // block containing a later one (or end)
        // find the last block b in [lo, hi] with block_ones[b] <= r
        while hi - lo > 8 {
            let mid = (lo + hi) / 2;
            if self.block_ones[mid] as usize <= r {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        while lo + 1 < self.block_ones.len() && self.block_ones[lo + 1] as usize <= r {
            lo += 1;
        }
        let rem = r - self.block_ones[lo] as usize;
        // word within the block: last j with ones-before-word-j <= rem
        let sub = self.sub[lo];
        // number of words j in 1..=7 whose ones-before count is <= rem
        let j = uleq9(sub, rem as u64 * ONES9).count_ones() as usize;
        let before = if j == 0 { 0 } else { ((sub >> (9 * (j - 1))) & 0x1FF) as usize };
        let wi = lo * BLOCK_WORDS + j;
        wi * 64 + select_in_word(self.bits.words[wi], (rem - before) as u32) as usize
    }

    /// Position of the first one at or after `pos` (`len()` if none).
    #[inline(always)]
    pub fn next_one(&self, pos: usize) -> usize {
        let words = self.bits.words();
        let mut wi = pos / 64;
        if wi >= words.len() {
            return self.len();
        }
        let mut w = words[wi] & !mask64((pos % 64) as u32);
        loop {
            if w != 0 {
                return (wi * 64 + w.trailing_zeros() as usize).min(self.len());
            }
            wi += 1;
            if wi >= words.len() {
                return self.len();
            }
            w = words[wi];
        }
    }

    /// Position of the last one strictly before `pos`, if any.
    #[inline(always)]
    pub fn prev_one(&self, pos: usize) -> Option<usize> {
        if pos == 0 {
            return None;
        }
        let words = self.bits.words();
        let mut wi = (pos - 1) / 64;
        let mut w = words[wi] & mask64(((pos - 1) % 64) as u32 + 1);
        loop {
            if w != 0 {
                return Some(wi * 64 + 63 - w.leading_zeros() as usize);
            }
            if wi == 0 {
                return None;
            }
            wi -= 1;
            w = words[wi];
        }
    }

    /// `(select1(r), select1(r + 1))`.
    #[inline(always)]
    pub fn select1_pair(&self, r: usize) -> (usize, usize) {
        let p = self.select1(r);
        (p, self.next_one(p + 1))
    }

    /// Position of the `r`-th zero (0-based).
    #[inline]
    pub fn select0(&self, r: usize) -> usize {
        if let Some((p0, mut rem)) = self.pos0.start(r) {
            let mut wi = (p0 >> 6) as usize;
            let mut w = !self.bits.words[wi] & !mask64((p0 & 63) as u32);
            loop {
                let c = w.count_ones();
                if rem < c {
                    return wi * 64 + select_in_word(w, rem) as usize;
                }
                rem -= c;
                wi += 1;
                w = !self.bits.words[wi];
            }
        }
        self.select0_blocks(r)
    }

    #[inline(never)]
    fn select0_blocks(&self, r: usize) -> usize {
        debug_assert!(r < self.len() - self.ones, "select0({r}) out of range");
        let zeros_before = |b: usize| b * BLOCK_BITS - self.block_ones[b] as usize;
        let s = r / SAMPLE;
        let mut lo = self.sel0[s] as usize;
        let mut hi = self.sel0[s + 1] as usize;
        while hi - lo > 8 {
            let mid = (lo + hi) / 2;
            if zeros_before(mid) <= r {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let nblocks = self.block_ones.len() - 1;
        while lo + 1 < nblocks && zeros_before(lo + 1) <= r {
            lo += 1;
        }
        let rem = r - zeros_before(lo);
        let zsub = WORDBITS9 - self.sub[lo];
        let j = uleq9(zsub, rem as u64 * ONES9).count_ones() as usize;
        let before = if j == 0 { 0 } else { ((zsub >> (9 * (j - 1))) & 0x1FF) as usize };
        let wi = lo * BLOCK_WORDS + j;
        wi * 64 + select_in_word(!self.bits.words[wi], (rem - before) as u32) as usize
    }

    /// Size in bits including the rank/select support.
    pub fn bit_size(&self) -> usize {
        self.bits.bit_size()
            + 64 * (self.block_ones.len() + self.sub.len())
            + self.pos1.bit_size()
            + self.pos0.bit_size()
            + 32 * (self.sel1.len() + self.sel0.len())
    }

    pub fn write(&self, w: &mut Writer) -> io::Result<()> {
        self.bits.write(w)
    }

    pub fn read(r: &mut Reader) -> io::Result<Self> {
        Ok(Self::new(BitVec::read(r)?))
    }
}

/// Fixed-width bit-packed vector of integers of `width <= 64` bits.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CompactVec {
    words: Vec<u64>,
    len: usize,
    width: u32,
    mask: u64,
}

impl CompactVec {
    pub fn new(len: usize, width: u32) -> Self {
        assert!(width <= 64);
        // one extra word so that `get` may always read two words
        let nwords = (len * width as usize).div_ceil(64) + 1;
        Self { words: vec![0; nwords], len, width, mask: mask64(width) }
    }

    pub fn from_slice(values: &[u64], width: u32) -> Self {
        let mut v = Self::new(values.len(), width);
        for (i, &x) in values.iter().enumerate() {
            v.set(i, x);
        }
        v
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline(always)]
    pub fn width(&self) -> u32 {
        self.width
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> u64 {
        debug_assert!(i < self.len);
        if self.width == 0 {
            return 0;
        }
        let pos = i * self.width as usize;
        let wi = pos / 64;
        let off = (pos % 64) as u32;
        let mut v = self.words[wi] >> off;
        if off + self.width > 64 {
            v |= self.words[wi + 1] << (64 - off);
        }
        v & self.mask
    }

    pub fn set(&mut self, i: usize, x: u64) {
        debug_assert!(i < self.len);
        debug_assert!(x & !self.mask == 0, "value {x} does not fit in {} bits", self.width);
        if self.width == 0 {
            return;
        }
        let pos = i * self.width as usize;
        let wi = pos / 64;
        let off = (pos % 64) as u32;
        self.words[wi] = (self.words[wi] & !(self.mask << off)) | (x << off);
        if off + self.width > 64 {
            let spill = 64 - off;
            self.words[wi + 1] = (self.words[wi + 1] & !(self.mask >> spill)) | (x >> spill);
        }
    }

    /// Size in bits of the payload.
    pub fn bit_size(&self) -> usize {
        self.len * self.width as usize
    }

    pub fn write(&self, w: &mut Writer) -> io::Result<()> {
        w.u64(self.len as u64)?;
        w.u32(self.width)?;
        w.vec_u64(&self.words)
    }

    pub fn read(r: &mut Reader) -> io::Result<Self> {
        let len = r.u64()? as usize;
        let width = r.u32()?;
        let words = r.vec_u64()?;
        if width > 64 || words.len() != (len * width as usize).div_ceil(64) + 1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "compact vector length mismatch"));
        }
        Ok(Self { words, len, width, mask: mask64(width) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(state: &mut u64) -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *state >> 33
    }

    #[test]
    fn rank_select_random_densities() {
        let mut st = 3u64;
        for &(len, density) in &[
            (0usize, 50u64),
            (1, 50),
            (100, 50),
            (5000, 1),
            (5000, 99),
            (100_000, 50),
            (200_000, 2),
            (200_000, 98),
            (70_000, 0),
            (70_000, 100),
        ] {
            let mut bv = BitVec::new(len);
            let mut ones = Vec::new();
            let mut zeros = Vec::new();
            for i in 0..len {
                if lcg(&mut st) % 100 < density {
                    bv.set(i);
                    ones.push(i);
                } else {
                    zeros.push(i);
                }
            }
            let rs = RankSelect::new(bv);
            assert_eq!(rs.count_ones(), ones.len());
            for (r, &p) in ones.iter().enumerate() {
                assert_eq!(rs.select1(r), p);
                assert_eq!(rs.rank1(p), r);
            }
            for (r, &p) in zeros.iter().enumerate() {
                assert_eq!(rs.select0(r), p);
                assert_eq!(rs.rank0(p), r);
            }
            assert_eq!(rs.rank1(len), ones.len());
        }
    }

    #[test]
    fn rank_select_skewed() {
        // dense and sparse regions, long runs of zeros and ones
        let mut st = 9u64;
        let len = 300_000;
        let mut bv = BitVec::new(len);
        let mut ones = Vec::new();
        let mut zeros = Vec::new();
        for i in 0..len {
            let density = match i / 30_000 {
                0 => 999,
                1 => 1,
                2 => 0,
                3 => 1000,
                4 => 500,
                5 => 20,
                _ => 990,
            };
            if lcg(&mut st) % 1000 < density {
                bv.set(i);
                ones.push(i);
            } else {
                zeros.push(i);
            }
        }
        let rs = RankSelect::new(bv);
        for (r, &p) in ones.iter().enumerate() {
            assert_eq!(rs.select1(r), p);
        }
        for (r, &p) in zeros.iter().enumerate() {
            assert_eq!(rs.select0(r), p);
        }
    }

    #[test]
    fn compact_vec_roundtrip() {
        let mut st = 5u64;
        for width in 0..=64u32 {
            let vals: Vec<u64> = (0..300).map(|_| (lcg(&mut st) ^ (lcg(&mut st) << 32)) & mask64(width)).collect();
            let cv = CompactVec::from_slice(&vals, width);
            for (i, &v) in vals.iter().enumerate() {
                assert_eq!(cv.get(i), v, "width {width} index {i}");
            }
        }
    }
}

#[cfg(test)]
mod span_tests {
    use super::*;

    #[test]
    fn ef_like_spans() {
        // EF upper bits of 589k random values in 2^36 with l = 16
        let n = 589_021usize;
        let mut st = 1u64;
        let mut vals: Vec<u64> = (0..n)
            .map(|_| {
                st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                (st >> 28) & ((1 << 36) - 1)
            })
            .collect();
        vals.sort();
        let mut bv = BitVec::new(n + (1 << 20) + 1);
        for (i, v) in vals.iter().enumerate() {
            bv.set((v >> 16) as usize + i);
        }
        let rs = RankSelect::new(bv);
        let slow = rs.pos0.pos.windows(2).filter(|w| (w[1] >> 6) - (w[0] >> 6) > MAX_SCAN_WORDS).count();
        eprintln!("pos0 samples {} long {}", rs.pos0.pos.len(), slow);
        for r in 0..rs.len() - rs.count_ones() {
            assert_eq!(rs.select0(r), rs.select0_blocks(r));
        }
    }
}
