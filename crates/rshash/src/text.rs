//! Packed text of all input sequences (port of `pack_dna4_to_uint64`,
//! `mark_sequences`, `get_word64`, `get_base`, `check_overlap`).
//!
//! Layout (as in C++): word 0 is all ones (32 bases of `T` padding), then all
//! sequences are concatenated without separators, 2 bits per base, first base
//! in the lowest bits; the text ends with padding words. The global position
//! of base `j` of the concatenation is `32 + j`. Sequence boundaries are kept
//! in an Elias–Fano sequence (`endpoints`).

use crate::alphabet;
use crate::bits::BitVec;
use crate::ef::EliasFano;
use crate::io::{Reader, Writer, invalid};
use crate::word::KmerWord;
use std::io;

/// Global position of the first base.
pub const TEXT_START: u64 = 32;
/// Trailing padding words (enough for 128-bit reads past the last base).
const TAIL_PADDING: usize = 4;

/// Incrementally packs sequences into a [`Text`].
#[derive(Default)]
pub struct TextBuilder {
    words: Vec<u64>,
    word: u64,
    shift: u32,
    bounds: Vec<u64>,
    len: u64,
    non_acgt: u64,
}

impl TextBuilder {
    pub fn new() -> Self {
        Self { words: vec![u64::MAX], bounds: vec![TEXT_START], ..Default::default() }
    }

    /// Append one sequence (ASCII). Non-ACGT characters become `A` (seqan3 `dna4`).
    pub fn push_ascii(&mut self, seq: &[u8]) {
        for &c in seq {
            if !alphabet::is_acgt(c) {
                self.non_acgt += 1;
            }
            self.push_rank(alphabet::rank(c));
        }
        self.bounds.push(TEXT_START + self.len);
    }

    /// Append one sequence given as ranks.
    pub fn push_ranks(&mut self, ranks: &[u8]) {
        for &r in ranks {
            self.push_rank(r & 3);
        }
        self.bounds.push(TEXT_START + self.len);
    }

    #[inline(always)]
    fn push_rank(&mut self, r: u8) {
        self.word |= (r as u64) << self.shift;
        self.shift += 2;
        self.len += 1;
        if self.shift == 64 {
            self.words.push(self.word);
            self.word = 0;
            self.shift = 0;
        }
    }

    /// Number of non-ACGT characters converted to `A` so far.
    pub fn non_acgt(&self) -> u64 {
        self.non_acgt
    }

    pub fn finish(mut self) -> Text {
        if self.shift != 0 {
            self.words.push(self.word);
        }
        self.words.extend(std::iter::repeat_n(u64::MAX, TAIL_PADDING));
        let universe = *self.bounds.last().unwrap() as u128 + 1;
        let bounds = EliasFano::new(&self.bounds, universe);
        let boundary_blocks = Text::boundary_blocks(&bounds, self.words.len());
        Text { words: self.words, len: self.len, bounds, boundary_blocks }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Text {
    words: Vec<u64>,
    /// Total number of bases (without padding).
    len: u64,
    /// `[32, 32+len0, 32+len0+len1, ..., 32+len]`.
    bounds: EliasFano<u64>,
    /// Derived (not serialised): bit `b` is set if a sequence boundary lies in
    /// positions `[64b, 64b+64)`. Lets most windows skip the Elias–Fano search
    /// in [`Text::within_one_sequence`].
    boundary_blocks: BitVec,
}

impl Text {
    pub fn from_sequences<S: AsRef<[u8]>>(seqs: &[S]) -> Text {
        let mut b = TextBuilder::new();
        for s in seqs {
            b.push_ascii(s.as_ref());
        }
        b.finish()
    }

    /// Total number of bases.
    #[inline(always)]
    pub fn len(&self) -> u64 {
        self.len
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// One past the last global position.
    #[inline(always)]
    pub fn end(&self) -> u64 {
        TEXT_START + self.len
    }

    pub fn num_sequences(&self) -> usize {
        self.bounds.len() - 1
    }

    /// `[start, end)` global positions of sequence `i`.
    #[inline]
    pub fn sequence_bounds(&self, i: usize) -> (u64, u64) {
        self.bounds.select_pair(i)
    }

    pub fn sequence_len(&self, i: usize) -> u64 {
        let (s, e) = self.sequence_bounds(i);
        e - s
    }

    /// Sequence containing global position `pos`: `(id, start, end)`
    /// (`check_overlap` / `endpoints.select(endpoints.rank(pos+1)-1, &end)`).
    ///
    /// Empty sequences never contain a position; the last sequence starting at
    /// or before `pos` is returned.
    #[inline]
    pub fn locate(&self, pos: u64) -> (usize, u64, u64) {
        debug_assert!(pos >= TEXT_START && pos < self.end());
        let (i, s, e) = self.bounds.pred_pair(pos).expect("position before the first sequence");
        (i, s, e.expect("position after the last sequence"))
    }

    fn boundary_blocks(bounds: &EliasFano<u64>, nwords: usize) -> BitVec {
        let mut bv = BitVec::new(nwords.div_ceil(2) + 1);
        for x in bounds.iter() {
            bv.set((x >> 6) as usize);
        }
        bv
    }

    /// Whether `[pos, pos+len)` lies within a single sequence.
    #[inline]
    pub fn within_one_sequence(&self, pos: u64, len: u64) -> bool {
        if pos < TEXT_START || pos + len > self.end() || len == 0 {
            return len == 0 && pos >= TEXT_START && pos <= self.end();
        }
        // a boundary strictly inside the window is some x with pos < x < pos+len
        let (b0, b1) = ((pos + 1) >> 6, (pos + len - 1) >> 6);
        if (b0..=b1).all(|b| !self.boundary_blocks.get(b as usize)) {
            return true;
        }
        let (_, _, e) = self.locate(pos);
        pos + len <= e
    }

    /// 32 bases starting at global position `pos`, first base in the low bits.
    #[inline(always)]
    pub fn word64(&self, pos: u64) -> u64 {
        let block = (pos >> 5) as usize;
        let shift = ((pos & 31) << 1) as u32;
        let lo = self.words[block];
        if shift == 0 {
            return lo;
        }
        let hi = self.words[block + 1];
        (lo >> shift) | (hi << (64 - shift))
    }

    /// The base (rank) at global position `pos`.
    #[inline(always)]
    pub fn base(&self, pos: u64) -> u64 {
        (self.words[(pos >> 5) as usize] >> ((pos & 31) << 1)) & 3
    }

    /// `w` bases starting at `pos` (unmasked above `2w` bits is *not* allowed:
    /// the result is masked).
    #[inline(always)]
    pub fn window<W: KmerWord>(&self, pos: u64, w: u32) -> W {
        if W::BITS == 64 {
            W::from_u64(self.word64(pos)) & W::mask(2 * w)
        } else {
            let lo = self.word64(pos) as u128;
            let v = if w > 32 { lo | ((self.word64(pos + 32) as u128) << 64) } else { lo };
            W::from_u128_trunc(v) & W::mask(2 * w)
        }
    }

    /// Bases `[start, end)` as ranks.
    pub fn ranks(&self, start: u64, end: u64) -> Vec<u8> {
        (start..end).map(|p| self.base(p) as u8).collect()
    }

    /// Size in bits of the packed text.
    pub fn bit_size(&self) -> usize {
        self.words.len() * 64
    }

    pub fn bounds_bit_size(&self) -> usize {
        self.bounds.bit_size()
    }

    pub fn write(&self, w: &mut Writer) -> io::Result<()> {
        w.u64(self.len)?;
        w.vec_u64(&self.words)?;
        self.bounds.write(w)
    }

    pub fn read(r: &mut Reader) -> io::Result<Self> {
        let len = r.u64()?;
        let words = r.vec_u64()?;
        let bounds = EliasFano::read(r)?;
        if (words.len() as u64) < (TEXT_START + len).div_ceil(32) + TAIL_PADDING as u64 || bounds.is_empty() {
            return Err(invalid("text length mismatch"));
        }
        let boundary_blocks = Text::boundary_blocks(&bounds, words.len());
        Ok(Self { words, len, bounds, boundary_blocks })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::word::{decode_kmer, encode_kmer};

    #[test]
    fn layout_and_access() {
        let seqs = ["ACGTACGTTTGACCA", "", "GGGATTACA", "C"];
        let t = Text::from_sequences(&seqs);
        assert_eq!(t.len(), 25);
        assert_eq!(t.num_sequences(), 4);
        assert_eq!(t.sequence_bounds(0), (32, 47));
        assert_eq!(t.sequence_bounds(1), (47, 47));
        assert_eq!(t.sequence_bounds(2), (47, 56));
        assert_eq!(t.sequence_bounds(3), (56, 57));
        assert_eq!(t.locate(47), (2, 47, 56));
        assert_eq!(t.locate(46), (0, 32, 47));
        assert_eq!(t.locate(56), (3, 56, 57));
        // padding before the first base is T
        assert_eq!(t.base(31), 3);
        let all: String = seqs.concat();
        for i in 0..all.len() {
            for w in 1..=(all.len() - i).min(40) as u32 {
                let v: u128 = t.window(32 + i as u64, w);
                assert_eq!(decode_kmer(v, w), &all[i..i + w as usize]);
                if w <= 32 {
                    let v: u64 = t.window(32 + i as u64, w);
                    assert_eq!(v, encode_kmer::<u64>(&all.as_bytes()[i..i + w as usize]));
                }
            }
        }
        assert!(t.within_one_sequence(32, 15));
        assert!(!t.within_one_sequence(33, 15));
        assert!(t.within_one_sequence(47, 9));
        // exhaustive against locate
        for p in 32..t.end() {
            for l in 1..=(t.end() - p) {
                let (_, _, e) = t.locate(p);
                assert_eq!(t.within_one_sequence(p, l), p + l <= e, "p={p} l={l}");
            }
        }
    }
}
