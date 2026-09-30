//! The RSHash index (port of `class RSHash`, `build.cpp` and `io.cpp`).

use crate::bits::{BitVec, CompactVec, RankSelect};
use crate::ef::EliasFano;
use crate::hash::SEEDS;
use crate::io::{Reader, Writer, invalid};
use crate::last_level::LastTable;
use crate::minimizer::{MinimizerParams, min_occurrences, window_minimizers};
use crate::params::{Geometry, ParamError, Params};
use crate::text::Text;
use crate::word::{KmerWord, canonical};
use std::io;

/// One minimiser level: `R` (kept minimisers), `S` (bucket delimiters) and the
/// offsets (global text positions of the minimiser occurrences).
#[derive(Clone, Debug)]
pub struct Level {
    pub mp: MinimizerParams,
    pub threshold: u16,
    pub r: EliasFano<u64>,
    pub s: RankSelect,
    pub offsets: CompactVec,
}

impl Level {
    /// Occurrence index range `[start, end)` of the bucket with rank `rank`.
    #[inline(always)]
    pub fn bucket(&self, rank: usize) -> (usize, usize) {
        self.s.select1_pair(rank)
    }

    /// Number of distinct kept minimisers.
    pub fn num_minimizers(&self) -> usize {
        self.r.len()
    }

    /// Number of stored minimiser occurrences.
    pub fn num_occurrences(&self) -> usize {
        self.offsets.len()
    }

    fn write(&self, w: &mut Writer) -> io::Result<()> {
        self.r.write(w)?;
        self.s.write(w)?;
        self.offsets.write(w)
    }

    fn read(r: &mut Reader, mp: MinimizerParams, threshold: u16) -> io::Result<Self> {
        let rr = EliasFano::read(r)?;
        let s = RankSelect::read(r)?;
        let offsets = CompactVec::read(r)?;
        if s.count_ones() != rr.len() + 1 || s.len() != offsets.len() + 1 {
            return Err(invalid("level size mismatch"));
        }
        Ok(Self { mp, threshold, r: rr, s, offsets })
    }
}

/// A shape placed in the window, with masks converted to the window type.
#[derive(Clone, Copy, Debug)]
pub struct ShapeW<W: KmerWord> {
    pub w_mask: W,
    /// Span `[start, end)` of the shape within the window.
    pub start: u32,
    pub end: u32,
    pub weight: u32,
    /// Palindromic and centred in the window: `pext(rc(X)) = rc(pext(X))`, so
    /// the last level can store canonical keys and probe once.
    pub canonical: bool,
}

/// Build statistics (reported by `print_info`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BuildStats {
    /// Number of kernel k-mers in the text.
    pub text_kmers: u64,
    /// Number of kernel k-mers that ended up in the last level.
    pub last_level_kmers: u64,
    /// Non-ACGT input characters converted to `A`.
    pub non_acgt: u64,
}

/// The RSHash dictionary over windows of type `W` (`u64` for windows of up to
/// 32 bases, `u128` up to 63).
#[derive(Clone, Debug)]
pub struct RsHash<W: KmerWord> {
    pub(crate) params: Params,
    pub(crate) geo: Geometry,
    pub(crate) text: Text,
    pub(crate) levels: Vec<Level>,
    /// One table per shape (one table without shapes).
    pub(crate) last: Vec<LastTable<W>>,
    pub(crate) shapes: Vec<ShapeW<W>>,
    pub(crate) kernel_mask: W,
    pub(crate) window_mask: W,
    pub(crate) stats: BuildStats,
}

#[derive(Debug, thiserror::Error)]
pub enum BuildError {
    #[error(transparent)]
    Params(#[from] ParamError),
    #[error("window of {0} bases does not fit the {1}-bit word type")]
    Word(u32, u32),
}

/// A (minimiser, position) pair of the build-time list, ordered by
/// minimiser then position.
trait Occ: Copy + Ord + Send {
    fn new(v: u64, pos: u64) -> Self;
    fn value(&self) -> u64;
    fn pos(&self) -> u64;
}

/// 12-byte pair for texts whose positions fit in 32 bits (`MinimizerInfo32`).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Occ32 {
    v_hi: u32,
    v_lo: u32,
    pos: u32,
}

impl Occ for Occ32 {
    #[inline(always)]
    fn new(v: u64, pos: u64) -> Self {
        Self { v_hi: (v >> 32) as u32, v_lo: v as u32, pos: pos as u32 }
    }
    #[inline(always)]
    fn value(&self) -> u64 {
        ((self.v_hi as u64) << 32) | self.v_lo as u64
    }
    #[inline(always)]
    fn pos(&self) -> u64 {
        self.pos as u64
    }
}

impl Occ for (u64, u64) {
    #[inline(always)]
    fn new(v: u64, pos: u64) -> Self {
        (v, pos)
    }
    #[inline(always)]
    fn value(&self) -> u64 {
        self.0
    }
    #[inline(always)]
    fn pos(&self) -> u64 {
        self.1
    }
}

/// `get_minimizers` + `filter_freq_minimizers` + `mark_occurences`: all
/// minimiser occurrences of the fragments, sorted, keeping minimisers with at
/// most `threshold` occurrences. Returns the kept minimisers, the bucket
/// delimiters `S` and the offsets.
fn collect_minimizers<O: Occ>(
    text: &Text,
    mp: &MinimizerParams,
    frags: &[Fragment],
    threshold: usize,
    pos_width: u32,
    log: &mut dyn FnMut(&str),
    l: usize,
) -> (Vec<u64>, BitVec, CompactVec) {
    let mut occ: Vec<O> = Vec::new();
    for f in frags {
        min_occurrences(text, mp, f.start, f.end, |v, pos| occ.push(O::new(v, pos)));
    }
    log(&format!("level {}: sorting {} minimizer occurrences...", l + 1, occ.len()));
    occ.sort_unstable();
    // groups of equal minimisers
    let groups = |mut f: Box<dyn FnMut(usize, usize) + '_>| {
        let mut i = 0;
        while i < occ.len() {
            let v = occ[i].value();
            let mut j = i + 1;
            while j < occ.len() && occ[j].value() == v {
                j += 1;
            }
            if j - i <= threshold {
                f(i, j);
            }
            i = j;
        }
    };
    let (mut nkeys, mut npos) = (0usize, 0usize);
    groups(Box::new(|i, j| {
        nkeys += 1;
        npos += j - i;
    }));
    let mut keys = Vec::with_capacity(nkeys);
    let mut s = BitVec::new(npos + 1);
    s.set(0);
    let mut offsets = CompactVec::new(npos, pos_width);
    let mut k = 0usize;
    groups(Box::new(|i, j| {
        keys.push(occ[i].value());
        for o in &occ[i..j] {
            offsets.set(k, o.pos());
            k += 1;
        }
        s.set(k);
    }));
    (keys, s, offsets)
}

/// A maximal run of consecutive k-mers `[start, end)` (global positions, the
/// k-mers are those fully inside the range) handed from one level to the next
/// (`SkmerInfo`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fragment {
    pub start: u64,
    pub end: u64,
}

impl<W: KmerWord> RsHash<W> {
    fn derived(params: &Params) -> Result<(Geometry, Vec<ShapeW<W>>, W, W), BuildError> {
        let geo = params.geometry()?;
        if geo.window > W::MAX_BASES {
            return Err(BuildError::Word(geo.window, W::BITS));
        }
        let shapes = geo
            .shapes
            .as_ref()
            .map(|s| {
                s.shapes
                    .iter()
                    .map(|a| ShapeW {
                        w_mask: W::from_u128_trunc(a.w_mask),
                        start: a.start,
                        end: a.end,
                        weight: a.shape.weight,
                        canonical: a.shape.is_palindrome && a.start == s.length - a.end,
                    })
                    .collect()
            })
            .unwrap_or_default();
        let kernel_mask = geo.shapes.as_ref().map_or(W::mask(2 * geo.kernel), |s| W::from_u128_trunc(s.kernel_mask));
        let window_mask = W::mask(2 * geo.window);
        Ok((geo, shapes, kernel_mask, window_mask))
    }

    fn level_params(params: &Params, geo: &Geometry, l: usize) -> MinimizerParams {
        MinimizerParams::new(geo.kernel, params.m[l], SEEDS[l])
    }

    /// Build the index over `text` (`RSHash::build`).
    pub fn build(text: Text, params: &Params) -> Result<Self, BuildError> {
        Self::build_with_log(text, params, &mut |_| {})
    }

    /// Build, reporting progress messages to `log`.
    pub fn build_with_log(text: Text, params: &Params, log: &mut dyn FnMut(&str)) -> Result<Self, BuildError> {
        let (geo, shapes, kernel_mask, window_mask) = Self::derived(params)?;
        let k = geo.kernel as u64;
        let pos_width = 64 - text.end().leading_zeros();

        let mut stats = BuildStats::default();
        let mut frags: Vec<Fragment> = (0..text.num_sequences())
            .map(|i| text.sequence_bounds(i))
            .filter(|&(s, e)| e - s >= k)
            .map(|(start, end)| Fragment { start, end })
            .collect();
        stats.text_kmers = frags.iter().map(|f| f.end - f.start - k + 1).sum();

        let mut levels = Vec::with_capacity(params.levels as usize);
        for l in 0..params.levels as usize {
            let mp = Self::level_params(params, &geo, l);
            log(&format!("level {}: computing minimizers and positions (m={}, t={})...", l + 1, mp.m, params.t[l]));
            // (minimiser, position) pairs; 12-byte pairs when positions fit in
            // 32 bits (MinimizerInfo32 in C++)
            let (keys, s, offsets) = if text.end() <= u32::MAX as u64 {
                collect_minimizers::<Occ32>(&text, &mp, &frags, params.t[l] as usize, pos_width, log, l)
            } else {
                collect_minimizers::<(u64, u64)>(&text, &mp, &frags, params.t[l] as usize, pos_width, log, l)
            };
            log(&format!("level {}: build R ({} minimizers, {} occurrences)...", l + 1, keys.len(), offsets.len()));
            let r = EliasFano::new(&keys, mp.universe());
            drop(keys);
            let s = RankSelect::new(s);
            let level = Level { mp, threshold: params.t[l], r, s, offsets };

            // get_frequent_skmers: maximal runs of k-mers whose minimiser is not in R
            log(&format!("level {}: collecting k-mers with frequent minimizers...", l + 1));
            let mut next = Vec::new();
            for f in &frags {
                let mut run_start: Option<u64> = None;
                let mut last = u64::MAX;
                let mut last_in = false;
                window_minimizers(&text, &mp, f.start, f.end, |pos, v| {
                    let inr = if v == last { last_in } else { level.r.contains(v).is_some() };
                    last = v;
                    last_in = inr;
                    match (inr, run_start) {
                        (false, None) => run_start = Some(pos),
                        (true, Some(s)) => {
                            next.push(Fragment { start: s, end: pos - 1 + k });
                            run_start = None;
                        }
                        _ => {}
                    }
                });
                if let Some(s) = run_start {
                    next.push(Fragment { start: s, end: f.end });
                }
            }
            frags = next;
            levels.push(level);
        }

        stats.last_level_kmers = frags.iter().map(|f| f.end - f.start - k + 1).sum();
        log(&format!("build last level ({} k-mers)...", stats.last_level_kmers));
        let last = Self::build_last_level(&text, params, &geo, &shapes, &frags, pos_width);

        Ok(Self { params: params.clone(), geo, text, levels, last, shapes, kernel_mask, window_mask, stats })
    }

    fn build_last_level(
        text: &Text,
        params: &Params,
        geo: &Geometry,
        shapes: &[ShapeW<W>],
        frags: &[Fragment],
        pos_width: u32,
    ) -> Vec<LastTable<W>> {
        let k = geo.kernel;
        if shapes.is_empty() {
            let mut pairs: Vec<(W, u64)> = Vec::new();
            for f in frags {
                let ranks_start = f.start;
                let mut v: W = text.window(ranks_start, k);
                let mut rc = v.rc(k);
                let mask = W::mask(2 * k);
                let shift = 2 * (k - 1);
                let mut q = f.start;
                loop {
                    pairs.push((canonical(v, rc), q));
                    if q + k as u64 >= f.end {
                        break;
                    }
                    let b = text.base(q + k as u64);
                    v = (v >> 2) | (W::from_u64(b) << shift);
                    rc = ((rc << 2) | W::from_u64(b ^ 3)) & mask;
                    q += 1;
                }
            }
            vec![LastTable::build(pairs, 2 * k, params.use_ht, params.loc, params.threshold, pos_width)]
        } else {
            let o = geo.overlap as u64;
            let l = geo.window;
            shapes
                .iter()
                .map(|sh| {
                    let mut pairs: Vec<(W, u64)> = Vec::new();
                    for f in frags {
                        let (_, seq_s, seq_e) = text.locate(f.start);
                        // kernel k-mers q in [f.start, f.end - k]; window start p = q - O
                        for q in f.start..=(f.end - k as u64) {
                            let p = q - o;
                            if p + (sh.start as u64) < seq_s || p + (sh.end as u64) > seq_e {
                                continue;
                            }
                            let win: W = text.window(p, l);
                            let key = win.pext(sh.w_mask);
                            let key = if sh.canonical { canonical(key, win.rc(l).pext(sh.w_mask)) } else { key };
                            pairs.push((key, p));
                        }
                    }
                    LastTable::build(pairs, 2 * sh.weight, params.use_ht, params.loc, params.threshold, pos_width)
                })
                .collect()
        }
    }

    pub fn params(&self) -> &Params {
        &self.params
    }

    pub fn geometry(&self) -> &Geometry {
        &self.geo
    }

    pub fn text(&self) -> &Text {
        &self.text
    }

    pub fn levels(&self) -> &[Level] {
        &self.levels
    }

    pub fn last_level(&self) -> &[LastTable<W>] {
        &self.last
    }

    pub fn stats(&self) -> &BuildStats {
        &self.stats
    }

    pub fn set_non_acgt(&mut self, n: u64) {
        self.stats.non_acgt = n;
    }

    /// Window length `L` (= k without shapes): the query unit.
    #[inline(always)]
    pub fn window(&self) -> u32 {
        self.geo.window
    }

    /// Kernel length `K` (= k without shapes).
    #[inline(always)]
    pub fn kernel(&self) -> u32 {
        self.geo.kernel
    }

    pub fn has_locate(&self) -> bool {
        self.params.loc
    }

    pub fn num_sequences(&self) -> usize {
        self.text.num_sequences()
    }

    /// The window of `L` bases starting at global position `pos` (`access`).
    pub fn access(&self, pos: u64) -> W {
        self.text.window(pos, self.geo.window)
    }

    /// Kernel (fwd, rc) of a window (fwd, rc).
    #[inline(always)]
    pub(crate) fn kernel_of(&self, v: W, rc: W) -> (W, W) {
        if self.shapes.is_empty() {
            (v, rc)
        } else {
            let sh = 2 * self.geo.overlap;
            ((v & self.kernel_mask) >> sh, (rc & self.kernel_mask) >> sh)
        }
    }

    // ---------------------------------------------------------------- io

    const MAGIC: &'static [u8; 8] = b"RSHASH01";
    const VERSION: u32 = 1;

    pub fn write(&self, out: &mut dyn io::Write) -> io::Result<()> {
        let mut w = Writer::new(out);
        w.bytes(Self::MAGIC)?;
        w.u32(Self::VERSION)?;
        w.u32(W::BITS)?;
        self.params.write(&mut w)?;
        w.u64(self.stats.text_kmers)?;
        w.u64(self.stats.last_level_kmers)?;
        w.u64(self.stats.non_acgt)?;
        self.text.write(&mut w)?;
        w.u32(self.levels.len() as u32)?;
        for l in &self.levels {
            l.write(&mut w)?;
        }
        w.u32(self.last.len() as u32)?;
        for t in &self.last {
            t.write(&mut w)?;
        }
        Ok(())
    }

    /// Read the header up to (and including) the word size.
    pub fn read_header(r: &mut dyn io::Read) -> io::Result<u32> {
        let mut rd = Reader::new(r);
        let mut magic = [0u8; 8];
        rd.bytes(&mut magic)?;
        if &magic != Self::MAGIC {
            return Err(invalid("not an rshash-rs index (bad magic)"));
        }
        let version = rd.u32()?;
        if version != Self::VERSION {
            return Err(invalid(&format!("unsupported index version {version}")));
        }
        rd.u32()
    }

    /// Read the body after [`Self::read_header`] returned `W::BITS`.
    pub fn read_body(input: &mut dyn io::Read) -> io::Result<Self> {
        let mut r = Reader::new(input);
        let params = Params::read(&mut r)?;
        let (geo, shapes, kernel_mask, window_mask) = Self::derived(&params).map_err(|e| invalid(&e.to_string()))?;
        let stats = BuildStats { text_kmers: r.u64()?, last_level_kmers: r.u64()?, non_acgt: r.u64()? };
        let text = Text::read(&mut r)?;
        let nlevels = r.u32()? as usize;
        if nlevels != params.levels as usize {
            return Err(invalid("level count mismatch"));
        }
        let mut levels = Vec::with_capacity(nlevels);
        for l in 0..nlevels {
            levels.push(Level::read(&mut r, Self::level_params(&params, &geo, l), params.t[l])?);
        }
        let ntables = r.u32()? as usize;
        if ntables != shapes.len().max(1) {
            return Err(invalid("last-level table count mismatch"));
        }
        let mut last = Vec::with_capacity(ntables);
        for _ in 0..ntables {
            last.push(LastTable::read(&mut r)?);
        }
        Ok(Self { params, geo, text, levels, last, shapes, kernel_mask, window_mask, stats })
    }

    pub fn read(input: &mut dyn io::Read) -> io::Result<Self> {
        let bits = Self::read_header(input)?;
        if bits != W::BITS {
            return Err(invalid("index word size mismatch"));
        }
        Self::read_body(input)
    }
}
