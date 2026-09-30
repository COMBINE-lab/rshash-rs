//! Random (hashed, canonical) minimisers.
//!
//! The score of an m-mer is `min(h(fwd) & mask, h(rc) & mask)` with
//! `h = mixer_64` seeded per level; the minimiser of a k-mer is the minimum
//! score over its `k - m + 1` m-mers. The score itself (not the m-mer) is the
//! key stored in the per-level Elias–Fano set `R`.
//!
//! * Query side: [`MinimizerParams::find`] / [`MinimizerParams::update`] are
//!   exact ports of `find_minimiser` / `update_minimiser` (`rshash.hpp`), which
//!   track the leftmost tied minimum (`left`, offset from the k-mer start) and
//!   the rightmost tied minimum (`right`, offset from the k-mer end).
//! * Build side: [`min_occurrences`] replaces `xor_minimiser_and_positions`
//!   and [`window_minimizers`] replaces `xor_minimiser_and_skmer_positions`.
//!   See DEVIATIONS.md (D1) for the differences.

use crate::hash::Mixer64;
use crate::text::Text;
use crate::word::{KmerWord, mask64};
use std::collections::VecDeque;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MinimizerParams {
    /// k-mer (kernel) length.
    pub k: u32,
    /// Minimiser length.
    pub m: u32,
    pub mask: u64,
    pub hasher: Mixer64,
}

/// State of the query-side rolling minimiser.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MinState {
    pub value: u64,
    /// Start of the leftmost minimal m-mer within the k-mer.
    pub left: u32,
    /// `k - m - start` of the rightmost minimal m-mer.
    pub right: u32,
}

impl MinState {
    /// Whether the leftmost and rightmost minimal m-mers differ.
    #[inline(always)]
    pub fn has_tie(&self, k: u32, m: u32) -> bool {
        self.left != k - m - self.right
    }
}

impl MinimizerParams {
    pub fn new(k: u32, m: u32, seed: u64) -> Self {
        assert!((1..=31).contains(&m) && m <= k, "invalid minimiser length m={m} for k={k}");
        Self { k, m, mask: mask64(2 * m), hasher: Mixer64::new(seed) }
    }

    /// Number of m-mers per k-mer (`span`).
    #[inline(always)]
    pub fn span(&self) -> u32 {
        self.k - self.m + 1
    }

    /// Size of the minimiser universe, `4^m`.
    pub fn universe(&self) -> u128 {
        1u128 << (2 * self.m)
    }

    #[inline(always)]
    pub fn score(&self, mmer: u64, mmer_rc: u64) -> u64 {
        let a = self.hasher.hash(mmer) & self.mask;
        let b = self.hasher.hash(mmer_rc) & self.mask;
        a.min(b)
    }

    /// `find_minimiser`: full scan of the k-mer, from the last m-mer to the first.
    #[inline(never)]
    pub fn find<W: KmerWord>(&self, kmer: W, kmer_rc: W) -> MinState {
        let (k, m) = (self.k, self.m);
        let mut mmer = (kmer >> (2 * (k - m))).low_u64() & self.mask;
        let mut mmer_rc = kmer_rc.low_u64() & self.mask;
        let mut st = MinState { value: self.score(mmer, mmer_rc), left: k - m, right: 0 };
        for i in 1..=(k - m) {
            mmer = (kmer >> (2 * (k - m - i))).low_u64() & self.mask;
            mmer_rc = (kmer_rc >> (2 * i)).low_u64() & self.mask;
            let h = self.score(mmer, mmer_rc);
            if h < st.value {
                st = MinState { value: h, left: k - m - i, right: i };
            } else if h == st.value {
                st.left = k - m - i;
            }
        }
        st
    }

    /// `update_minimiser`: slide by one base (the new m-mer is the last one).
    #[inline(always)]
    pub fn update<W: KmerWord>(&self, kmer: W, kmer_rc: W, st: &mut MinState) {
        if st.left == 0 {
            *st = self.find(kmer, kmer_rc);
            return;
        }
        st.left -= 1;
        let mmer = (kmer >> (2 * (self.k - self.m))).low_u64() & self.mask;
        let mmer_rc = kmer_rc.low_u64() & self.mask;
        let h = self.score(mmer, mmer_rc);
        if h < st.value {
            *st = MinState { value: h, left: self.k - self.m, right: 0 };
        } else if h == st.value {
            st.right = 0;
        } else {
            st.right += 1;
        }
    }
}

/// Sliding-window minimum with leftmost *and* rightmost ties, in amortised
/// O(1) and without branches on the data (the "two stacks" / van Herk–Gil–
/// Werman scheme used by `simd-minimizers`, Groot Koerkamp & Martayan 2025).
///
/// Elements are packed as `(score << 64) | pos` (leftmost tie wins) and
/// `(score << 64) | !pos` (rightmost tie wins), so each push costs a few
/// `u128` minima; every `w` pushes the finished block is turned into suffix
/// minima. Scores are kept exact (no truncation), so the result is identical
/// to [`MinimizerParams::find`].
#[derive(Clone, Debug, Default)]
pub struct SlidingLr {
    w: usize,
    idx: usize,
    ring_l: Vec<u128>,
    ring_r: Vec<u128>,
    pre_l: u128,
    pre_r: u128,
}

impl SlidingLr {
    pub fn new(w: usize) -> Self {
        assert!(w > 0);
        Self { w, idx: 0, ring_l: vec![u128::MAX; w], ring_r: vec![u128::MAX; w], pre_l: u128::MAX, pre_r: u128::MAX }
    }

    #[inline(always)]
    pub fn reset(&mut self) {
        self.idx = 0;
        self.pre_l = u128::MAX;
        self.pre_r = u128::MAX;
    }

    /// Push the score of the m-mer at (absolute) position `pos`; returns the
    /// packed leftmost and rightmost minima of the last `w` pushes (valid once
    /// `w` elements have been pushed since [`Self::reset`]).
    #[inline(always)]
    pub fn push(&mut self, score: u64, pos: u64) -> (u128, u128) {
        let kl = ((score as u128) << 64) | pos as u128;
        let kr = ((score as u128) << 64) | (!pos) as u128;
        let w = self.w;
        // SAFETY: idx < w == ring lengths
        unsafe {
            *self.ring_l.get_unchecked_mut(self.idx) = kl;
            *self.ring_r.get_unchecked_mut(self.idx) = kr;
        }
        self.pre_l = self.pre_l.min(kl);
        self.pre_r = self.pre_r.min(kr);
        self.idx += 1;
        if self.idx == w {
            self.idx = 0;
            for i in (0..w - 1).rev() {
                self.ring_l[i] = self.ring_l[i].min(self.ring_l[i + 1]);
                self.ring_r[i] = self.ring_r[i].min(self.ring_r[i + 1]);
            }
            self.pre_l = u128::MAX;
            self.pre_r = u128::MAX;
        }
        // SAFETY: idx < w
        let (sl, sr) = unsafe { (*self.ring_l.get_unchecked(self.idx), *self.ring_r.get_unchecked(self.idx)) };
        (self.pre_l.min(sl), self.pre_r.min(sr))
    }
}

impl MinimizerParams {
    /// Decode the result of [`SlidingLr::push`] for the k-mer whose first
    /// m-mer is at position `start`.
    #[inline(always)]
    pub fn decode_lr(&self, (l, r): (u128, u128), start: u64) -> MinState {
        let left = (l as u64 - start) as u32;
        let rightmost = (!(r as u64) - start) as u32;
        MinState { value: (l >> 64) as u64, left, right: self.k - self.m - rightmost }
    }

    /// Score of the m-mer at offset `i` of a k-mer (and its rc).
    #[inline(always)]
    pub fn score_at<W: KmerWord>(&self, kmer: W, kmer_rc: W, i: u32) -> u64 {
        let mmer = (kmer >> (2 * i)).low_u64() & self.mask;
        let mmer_rc = (kmer_rc >> (2 * (self.k - self.m - i))).low_u64() & self.mask;
        self.score(mmer, mmer_rc)
    }
}

/// Rolling m-mer scores over text positions `[start, end)`.
struct MmerScores<'a> {
    text: &'a Text,
    p: &'a MinimizerParams,
    pos: u64,
    end: u64,
    fwd: u64,
    rc: u64,
    shift: u32,
}

impl<'a> MmerScores<'a> {
    fn new(text: &'a Text, p: &'a MinimizerParams, start: u64, end: u64) -> Self {
        let mut s = Self { text, p, pos: start, end, fwd: 0, rc: 0, shift: 2 * (p.m - 1) };
        for _ in 0..p.m - 1 {
            s.push_base();
        }
        s
    }

    #[inline(always)]
    fn push_base(&mut self) {
        let b = self.text.base(self.pos);
        self.fwd = (self.fwd >> 2) | (b << self.shift);
        self.rc = ((self.rc << 2) | (b ^ 3)) & self.p.mask;
        self.pos += 1;
    }

    /// Score of the next m-mer.
    #[inline(always)]
    fn next(&mut self) -> Option<u64> {
        if self.pos >= self.end {
            return None;
        }
        self.push_base();
        Some(self.p.score(self.fwd, self.rc))
    }
}

/// Call `f(pos, minimiser)` for every k-mer starting at global position `pos`
/// in `[start, end)` (k-mers must lie fully inside the range).
///
/// Replaces `xor_minimiser_and_skmer_positions`, which is only used to learn
/// the minimiser of each k-mer.
pub fn window_minimizers(text: &Text, p: &MinimizerParams, start: u64, end: u64, mut f: impl FnMut(u64, u64)) {
    if end < start + p.k as u64 {
        return;
    }
    let w = p.span() as u64;
    let mut scores = MmerScores::new(text, p, start, end);
    let mut dq: VecDeque<(u64, u64)> = VecDeque::with_capacity(w as usize + 1);
    let mut i = 0u64;
    while let Some(h) = scores.next() {
        while dq.back().is_some_and(|&(_, v)| v >= h) {
            dq.pop_back();
        }
        dq.push_back((i, h));
        if i + 1 >= w {
            let j = i + 1 - w;
            while dq.front().unwrap().0 < j {
                dq.pop_front();
            }
            f(start + j, dq.front().unwrap().1);
        }
        i += 1;
    }
}

/// Call `f(minimiser, pos)` for every m-mer occurrence (global position `pos`
/// of its first base) that is a minimum — possibly tied — of at least one
/// k-mer lying fully in `[start, end)`. Occurrences are reported once, in
/// increasing position order.
///
/// Replaces `xor_minimiser_and_positions` (D1: tied occurrences carry their
/// true value, and ranges shorter than `k` yield nothing).
pub fn min_occurrences(text: &Text, p: &MinimizerParams, start: u64, end: u64, mut f: impl FnMut(u64, u64)) {
    if end < start + p.k as u64 {
        return;
    }
    let w = p.span() as usize;
    let n_k = (end - start - p.k as u64 + 1) as usize; // number of k-mers
    let ring_len = (w + 1).next_power_of_two();
    let ring_mask = ring_len - 1;
    let mut ring = vec![0u64; ring_len];
    // sliding minimum of scores (per k-mer): indices into ring
    let mut dq_min: VecDeque<(usize, u64)> = VecDeque::with_capacity(w + 1);
    // sliding maximum over the k-mer minima of the k-mers containing an m-mer
    let mut dq_max: VecDeque<(usize, u64)> = VecDeque::with_capacity(w + 1);

    let mut scores = MmerScores::new(text, p, start, end);
    let mut i = 0usize;
    while let Some(h) = scores.next() {
        ring[i & ring_mask] = h;
        while dq_min.back().is_some_and(|&(_, v)| v >= h) {
            dq_min.pop_back();
        }
        dq_min.push_back((i, h));
        if i + 1 >= w {
            let j = i + 1 - w; // k-mer index, also the m-mer whose k-mers are now all known
            while dq_min.front().unwrap().0 < j {
                dq_min.pop_front();
            }
            let mj = dq_min.front().unwrap().1;
            while dq_max.back().is_some_and(|&(_, v)| v <= mj) {
                dq_max.pop_back();
            }
            dq_max.push_back((j, mj));
            // m-mer j is contained in k-mers [j-w+1, j]
            while dq_max.front().unwrap().0 + w <= j {
                dq_max.pop_front();
            }
            if ring[j & ring_mask] == dq_max.front().unwrap().1 {
                f(ring[j & ring_mask], start + j as u64);
            }
        }
        i += 1;
    }
    // trailing m-mers n_k..n_m, contained in k-mers [p-w+1, n_k-1]
    let n_m = i;
    for pm in n_k..n_m {
        while dq_max.front().unwrap().0 + w <= pm {
            dq_max.pop_front();
        }
        if ring[pm & ring_mask] == dq_max.front().unwrap().1 {
            f(ring[pm & ring_mask], start + pm as u64);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::text::Text;

    fn lcg(state: &mut u64) -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *state >> 33
    }

    fn random_seq(st: &mut u64, n: usize, alphabet: usize) -> String {
        (0..n).map(|_| b"ACGT"[(lcg(st) as usize) % alphabet] as char).collect()
    }

    fn brute_scores(text: &Text, p: &MinimizerParams, pos: u64) -> Vec<u64> {
        (0..p.span() as u64)
            .map(|i| {
                let mm: u64 = text.window(pos + i, p.m);
                p.score(mm, mm.rc(p.m))
            })
            .collect()
    }

    #[test]
    fn query_minimiser_matches_brute_force() {
        let mut st = 42u64;
        for &(k, m, alpha) in
            &[(31u32, 18u32, 4usize), (21, 11, 2), (15, 15, 4), (47, 20, 4), (63, 31, 2), (32, 7, 1), (5, 2, 2)]
        {
            let seq = random_seq(&mut st, 400, alpha);
            let text = Text::from_sequences(&[&seq]);
            let p = MinimizerParams::new(k, m, 0x1234 + k as u64);
            let mut rolling: Option<MinState> = None;
            for pos in 32..(32 + 400 - k as u64 + 1) {
                let kmer: u128 = text.window(pos, k);
                let kmer_rc = kmer.rc(k);
                let sc = brute_scores(&text, &p, pos);
                let min = *sc.iter().min().unwrap();
                let leftmost = sc.iter().position(|&x| x == min).unwrap() as u32;
                let rightmost = sc.iter().rposition(|&x| x == min).unwrap() as u32;
                let expect = MinState { value: min, left: leftmost, right: k - m - rightmost };
                assert_eq!(p.find(kmer, kmer_rc), expect);
                if k <= 32 {
                    assert_eq!(p.find(kmer as u64, kmer_rc as u64), expect);
                }
                let r = match rolling {
                    Some(mut r) => {
                        p.update(kmer, kmer_rc, &mut r);
                        r
                    }
                    None => p.find(kmer, kmer_rc),
                };
                assert_eq!(r, expect, "k={k} m={m} pos={pos}");
                rolling = Some(r);
            }
        }
    }

    #[test]
    fn sliding_lr_matches_find() {
        let mut st = 5u64;
        for &(k, m, alpha) in &[(31u32, 18u32, 4usize), (21, 11, 2), (15, 15, 4), (47, 20, 2), (9, 3, 1), (32, 8, 2)] {
            let seq = random_seq(&mut st, 500, alpha);
            let text = Text::from_sequences(&[&seq]);
            let p = MinimizerParams::new(k, m, 3);
            let w = p.span() as usize;
            let mut lr = SlidingLr::new(w);
            for wi in 0..(500 - k as u64 + 1) {
                let kmer: u128 = text.window(32 + wi, k);
                let rc = kmer.rc(k);
                let res = if wi == 0 {
                    let mut r = (0, 0);
                    for i in 0..w as u32 {
                        r = lr.push(p.score_at(kmer, rc, i), i as u64);
                    }
                    r
                } else {
                    lr.push(p.score_at(kmer, rc, k - m), wi + w as u64 - 1)
                };
                assert_eq!(p.decode_lr(res, wi), p.find(kmer, rc), "k={k} m={m} wi={wi}");
            }
        }
    }

    #[test]
    fn build_views_match_brute_force() {
        let mut st = 7u64;
        for &(k, m, alpha, n) in &[
            (31u32, 18u32, 4usize, 500usize),
            (21, 11, 2, 300),
            (15, 15, 4, 100),
            (47, 20, 2, 300),
            (9, 3, 1, 50),
            (10, 4, 4, 9),
            (10, 4, 4, 10),
        ] {
            let seq = random_seq(&mut st, n, alpha);
            let text = Text::from_sequences(&[&seq]);
            let p = MinimizerParams::new(k, m, 99);
            let (s, e) = (32u64, 32 + n as u64);
            let mut expect_occ = std::collections::BTreeSet::new();
            let mut expect_win = Vec::new();
            if n >= k as usize {
                for pos in s..=(e - k as u64) {
                    let sc = brute_scores(&text, &p, pos);
                    let min = *sc.iter().min().unwrap();
                    expect_win.push((pos, min));
                    for (i, &x) in sc.iter().enumerate() {
                        if x == min {
                            expect_occ.insert((pos + i as u64, min));
                        }
                    }
                }
            }
            let mut win = Vec::new();
            window_minimizers(&text, &p, s, e, |pos, v| win.push((pos, v)));
            assert_eq!(win, expect_win);
            let mut occ = Vec::new();
            min_occurrences(&text, &p, s, e, |v, pos| occ.push((pos, v)));
            assert!(occ.windows(2).all(|w| w[0].0 < w[1].0));
            let occ: std::collections::BTreeSet<_> = occ.into_iter().collect();
            assert_eq!(occ, expect_occ, "k={k} m={m} n={n}");
        }
    }
}
