//! Queries: point lookup/locate (port of `lookup1..3` / `locate1..3`),
//! streaming lookup (`streaming_lookup1..3`) and streaming locate
//! (`streaming_locate1..3`).
//!
//! A query window `W` of `L` bases is present if some text window `T` lying
//! in one input sequence matches `W` or its reverse complement: exactly for
//! contiguous k-mers, and on the care positions of at least one shape (whose
//! span must lie in one sequence) for gapped k-mers (D8).

use crate::alphabet;
use crate::index::{Level, RsHash};
use crate::minimizer::MinState;
use crate::text::TEXT_START;
use crate::word::{KmerWord, Windows, canonical};

/// Maximum number of shapes per index.
pub const MAX_SHAPES: usize = 8;

/// A located occurrence: the text window starting at global position `pos`
/// matches the query forward (`forward`) or reverse-complemented.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Hit {
    pub pos: u64,
    pub forward: bool,
}

/// A query window with everything derived from it.
#[derive(Clone, Copy, Debug)]
pub struct QueryKeys<W: KmerWord> {
    pub v: W,
    pub rc: W,
    pub kernel: W,
    pub kernel_rc: W,
    /// `pext(v, w_mask_i)` per shape (`shapes_fwd`).
    pub fwd: [W; MAX_SHAPES],
    /// `pext(rc, w_mask_i)` per shape (`shapes_rev`).
    pub rev: [W; MAX_SHAPES],
}

impl<W: KmerWord> QueryKeys<W> {
    pub fn empty() -> Self {
        Self {
            v: W::ZERO,
            rc: W::ZERO,
            kernel: W::ZERO,
            kernel_rc: W::ZERO,
            fwd: [W::ZERO; MAX_SHAPES],
            rev: [W::ZERO; MAX_SHAPES],
        }
    }
}

impl<W: KmerWord> RsHash<W> {
    #[inline(always)]
    pub fn query_keys(&self, v: W, rc: W) -> QueryKeys<W> {
        let mut q = QueryKeys::empty();
        self.fill_query_keys(v, rc, &mut q);
        q
    }

    /// Fill `q` in place (avoids re-initialising the per-shape arrays).
    #[inline(always)]
    pub fn fill_query_keys(&self, v: W, rc: W, q: &mut QueryKeys<W>) {
        let (kernel, kernel_rc) = self.kernel_of(v, rc);
        q.v = v;
        q.rc = rc;
        q.kernel = kernel;
        q.kernel_rc = kernel_rc;
        for (i, s) in self.shapes.iter().enumerate() {
            q.fwd[i] = v.pext(s.w_mask);
            q.rev[i] = rc.pext(s.w_mask);
        }
    }

    /// Check a text window `win` starting at `p` against the query, in the
    /// given orientation. Returns the bounds of the sequence of the match.
    #[inline(always)]
    fn match_window(&self, win: W, p: i64, forward: bool, q: &QueryKeys<W>) -> Option<(u64, u64)> {
        if self.shapes.is_empty() {
            let target = if forward { q.v } else { q.rc };
            if win == target && p >= TEXT_START as i64 {
                let p = p as u64;
                let l = self.geo.window as u64;
                if p + l <= self.text.end() {
                    let (_, s, e) = self.text.locate(p);
                    if p + l <= e {
                        return Some((s, e));
                    }
                }
            }
            None
        } else {
            for (i, sh) in self.shapes.iter().enumerate() {
                let target = if forward { q.fwd[i] } else { q.rev[i] };
                if win.pext(sh.w_mask) == target {
                    let s = p + sh.start as i64;
                    let len = (sh.end - sh.start) as u64;
                    if s >= TEXT_START as i64 && (s as u64) + len <= self.text.end() {
                        let (_, ss, se) = self.text.locate(s as u64);
                        if s as u64 + len <= se {
                            return Some((ss, se));
                        }
                    }
                }
            }
            None
        }
    }

    /// Read the window whose kernel starts at `q` (if readable).
    #[inline(always)]
    fn window_at_kernel(&self, q: i64) -> Option<(W, i64)> {
        if q < TEXT_START as i64 || q as u64 + self.geo.kernel as u64 > self.text.end() {
            return None;
        }
        let p = q - self.geo.overlap as i64;
        Some((self.text.window(p as u64, self.geo.window), p))
    }

    /// Candidate kernel starts `(q, forward)` for the minimiser occurrence at
    /// `o`, in C++ check order: left fwd, left rc, [right fwd, right rc].
    #[inline(always)]
    fn candidates(&self, level: &Level, o: u64, st: &MinState) -> ([(i64, bool); 4], usize) {
        let km = (self.geo.kernel - level.mp.m) as i64;
        let o = o as i64;
        let offset = o - km; // offsets[i] = stored + 1 - span
        let mut c = [(0i64, true); 4];
        c[0] = (offset + km - st.left as i64, true);
        c[1] = (offset + st.left as i64, false);
        if st.has_tie(self.geo.kernel, level.mp.m) {
            c[2] = (offset + st.right as i64, true);
            c[3] = (offset + km - st.right as i64, false);
            (c, 4)
        } else {
            (c, 2)
        }
    }

    /// Last level membership (`lookup_last_level`).
    #[inline]
    pub(crate) fn last_contains(&self, q: &QueryKeys<W>) -> bool {
        if self.shapes.is_empty() {
            self.last[0].contains(canonical(q.v, q.rc))
        } else {
            self.last.iter().zip(&self.shapes).enumerate().any(|(i, (t, sh))| {
                if sh.canonical {
                    t.contains(canonical(q.fwd[i], q.rev[i]))
                } else {
                    t.contains(q.fwd[i]) || t.contains(q.rev[i])
                }
            })
        }
    }

    /// Point lookup of one window (`lookup` for a single k-mer).
    pub fn lookup(&self, v: W) -> bool {
        let v = v & self.window_mask;
        let q = self.query_keys(v, v.rc(self.geo.window));
        self.lookup_keys(&q)
    }

    pub fn lookup_keys(&self, q: &QueryKeys<W>) -> bool {
        for level in &self.levels {
            let st = level.mp.find(q.kernel, q.kernel_rc);
            if let Some(rank) = level.r.contains(st.value) {
                let (s, e) = level.bucket(rank);
                return (s..e).any(|i| self.check_occurrence(level, level.offsets.get(i), &st, q));
            }
        }
        self.last_contains(q)
    }

    /// `check<level>` for one occurrence (with the unitig-boundary check, D4).
    #[inline]
    fn check_occurrence(&self, level: &Level, o: u64, st: &MinState, q: &QueryKeys<W>) -> bool {
        let (cands, n) = self.candidates(level, o, st);
        cands[..n]
            .iter()
            .any(|&(kq, fwd)| self.window_at_kernel(kq).is_some_and(|(win, p)| self.matches(win, p, fwd, q)))
    }

    /// Like [`Self::match_window`] without computing the sequence bounds.
    #[inline(always)]
    fn matches(&self, win: W, p: i64, forward: bool, q: &QueryKeys<W>) -> bool {
        if self.shapes.is_empty() {
            let target = if forward { q.v } else { q.rc };
            win == target && p >= TEXT_START as i64 && self.text.within_one_sequence(p as u64, self.geo.window as u64)
        } else {
            self.shapes.iter().enumerate().any(|(i, sh)| {
                let target = if forward { q.fwd[i] } else { q.rev[i] };
                let s = p + sh.start as i64;
                win.pext(sh.w_mask) == target
                    && s >= TEXT_START as i64
                    && self.text.within_one_sequence(s as u64, (sh.end - sh.start) as u64)
            })
        }
    }

    /// Count the windows of `kmers` present in the index (`lookup(vector)`).
    pub fn lookup_many(&self, kmers: &[W]) -> u64 {
        kmers.iter().map(|&v| self.lookup(v) as u64).sum()
    }

    /// All occurrences of a window (`locate` for a single k-mer).
    pub fn locate(&self, v: W, out: &mut Vec<Hit>) {
        let v = v & self.window_mask;
        let q = self.query_keys(v, v.rc(self.geo.window));
        let start = out.len();
        let mut done = false;
        for level in &self.levels {
            let st = level.mp.find(q.kernel, q.kernel_rc);
            if let Some(rank) = level.r.contains(st.value) {
                self.locate_bucket(level, rank, &st, &q, out);
                done = true;
                break;
            }
        }
        if !done {
            self.locate_last(&q, out);
        }
        out[start..].sort_unstable();
        let mut w = start;
        for i in start..out.len() {
            if i == start || out[i] != out[w - 1] {
                out[w] = out[i];
                w += 1;
            }
        }
        out.truncate(w);
    }

    fn locate_bucket(&self, level: &Level, rank: usize, st: &MinState, q: &QueryKeys<W>, out: &mut Vec<Hit>) {
        let (s, e) = level.bucket(rank);
        for i in s..e {
            let (cands, n) = self.candidates(level, level.offsets.get(i), st);
            for &(kq, fwd) in &cands[..n] {
                if let Some((win, p)) = self.window_at_kernel(kq)
                    && self.matches(win, p, fwd, q)
                {
                    out.push(Hit { pos: p as u64, forward: fwd });
                }
            }
        }
    }

    fn locate_last(&self, q: &QueryKeys<W>, out: &mut Vec<Hit>) {
        let mut buf = Vec::new();
        if self.shapes.is_empty() {
            self.last[0].positions(canonical(q.v, q.rc), &mut buf);
            for &p in &buf {
                let win: W = self.text.window(p, self.geo.window);
                if win == q.v {
                    out.push(Hit { pos: p, forward: true });
                }
                if win == q.rc {
                    out.push(Hit { pos: p, forward: false });
                }
            }
        } else {
            for (i, (t, sh)) in self.last.iter().zip(&self.shapes).enumerate() {
                if sh.canonical {
                    buf.clear();
                    t.positions(canonical(q.fwd[i], q.rev[i]), &mut buf);
                    for &p in &buf {
                        let key = self.text.window::<W>(p, self.geo.window).pext(sh.w_mask);
                        if key == q.fwd[i] {
                            out.push(Hit { pos: p, forward: true });
                        }
                        if key == q.rev[i] {
                            out.push(Hit { pos: p, forward: false });
                        }
                    }
                } else {
                    buf.clear();
                    t.positions(q.fwd[i], &mut buf);
                    out.extend(buf.iter().map(|&p| Hit { pos: p, forward: true }));
                    buf.clear();
                    t.positions(q.rev[i], &mut buf);
                    out.extend(buf.iter().map(|&p| Hit { pos: p, forward: false }));
                }
            }
        }
    }

    /// Map a hit to `(sequence id, offset of the window start within it)`.
    /// With shapes the window may start before its sequence (negative offset)
    /// when only the shape's span lies in the sequence.
    pub fn hit_position(&self, hit: &Hit) -> (usize, i64) {
        let kernel_start = hit.pos + self.geo.overlap as u64;
        let (id, s, _) = self.text.locate(kernel_start);
        (id, hit.pos as i64 - s as i64)
    }

    /// A streaming lookup engine (reusable across query sequences).
    pub fn streaming_lookup(&self) -> StreamingLookup<'_, W> {
        StreamingLookup::new(self)
    }

    /// A streaming locate engine (reusable across query sequences).
    pub fn streaming_locate(&self) -> StreamingLocate<'_, W> {
        StreamingLocate::new(self)
    }
}

const INF: u64 = u64::MAX;

/// Per-level state of the streaming lookup.
struct LevelCache<W: KmerWord> {
    rolling: bool,
    st: MinState,
    current: u64,
    current_neg: u64,
    n: usize,
    /// `offsets[i] = stored + 1 - span`, i.e. the first possible kernel start.
    offsets: Vec<i64>,
    /// `span` windows per occurrence, starting at `offsets[i] - O`.
    buffer: Vec<W>,
    /// Bounds of the sequence containing each occurrence (locate only,
    /// `fill_buffer2`).
    bounds: Vec<(u64, u64)>,
}

impl<W: KmerWord> LevelCache<W> {
    fn new(level: &Level) -> Self {
        let span = level.mp.span() as usize;
        let t = level.threshold as usize;
        Self {
            rolling: false,
            st: MinState::default(),
            current: INF,
            current_neg: INF,
            n: 0,
            offsets: vec![0; t],
            buffer: vec![W::ZERO; t * span],
            bounds: Vec::new(),
        }
    }
}

/// Streaming lookup (Algorithm 1 of the paper): consecutive windows of a
/// query first try to *extend* the previous match by one base in the text;
/// otherwise the rolling minimisers of each level are looked up, caching the
/// decoded bucket of the current minimiser and the last negative minimiser.
pub struct StreamingLookup<'a, W: KmerWord> {
    idx: &'a RsHash<W>,
    levels: Vec<LevelCache<W>>,
    /// Number of windows found by extension (`extensions`).
    pub extensions: u64,
}

/// The text window matched by the previous query window, which the next one
/// may extend (`text_pos`, `forward`, `sequence_begin/end`, `text_kmer(_rc)`).
/// Kept in locals of the query loop so it lives in registers.
#[derive(Clone, Copy, Debug)]
struct Anchor<W: KmerWord> {
    forward: bool,
    /// Last base of the window (forward) or first base (reverse complement).
    text_pos: u64,
    seq_start: u64,
    seq_end: u64,
    /// The matched text window (forward) or its reverse-complement orientation.
    win: W,
}

impl<'a, W: KmerWord> StreamingLookup<'a, W> {
    pub fn new(idx: &'a RsHash<W>) -> Self {
        let levels = idx.levels.iter().map(LevelCache::new).collect();
        Self { idx, levels, extensions: 0 }
    }

    /// Number of windows of `seq` (ASCII) present in the index.
    pub fn query(&mut self, seq: &[u8]) -> u64 {
        self.query_with(seq, |_, _| {})
    }

    /// Like [`Self::query`], calling `f(window index, found)` per window.
    pub fn query_with(&mut self, seq: &[u8], f: impl FnMut(usize, bool)) -> u64 {
        if self.idx.shapes.is_empty() { self.run::<false>(seq, f) } else { self.run::<true>(seq, f) }
    }

    /// The query loop, monomorphised for contiguous k-mers / shapes.
    #[inline(always)]
    fn run<const SHAPES: bool>(&mut self, seq: &[u8], mut f: impl FnMut(usize, bool)) -> u64 {
        let idx = self.idx;
        let l = idx.geo.window as usize;
        if seq.len() < l {
            return 0;
        }
        let levels = &mut self.levels[..];
        for c in levels.iter_mut() {
            c.rolling = false;
        }
        let shift = 2 * (l as u32 - 1);
        let mask = idx.window_mask;
        let mut q = QueryKeys::empty();
        let (mut v, mut rc) = (W::ZERO, W::ZERO);
        for &c in &seq[..l - 1] {
            let r = alphabet::rank(c) as u64;
            v = (v >> 2) | (W::from_u64(r) << shift);
            rc = ((rc << 2) | W::from_u64(r ^ 3)) & mask;
        }
        let mut anchor: Option<Anchor<W>> = None;
        let (mut hits, mut ext) = (0u64, 0u64);
        let mut rolling_valid = true;
        for (wi, &c) in seq[l - 1..].iter().enumerate() {
            let r = alphabet::rank(c) as u64;
            v = (v >> 2) | (W::from_u64(r) << shift);
            rc = ((rc << 2) | W::from_u64(r ^ 3)) & mask;
            if SHAPES {
                // extension compares the per-shape keys (kernel only if needed)
                for (i, sh) in idx.shapes.iter().enumerate() {
                    q.fwd[i] = v.pext(sh.w_mask);
                    q.rev[i] = rc.pext(sh.w_mask);
                }
            }
            let extended = match anchor.as_mut() {
                Some(a) => {
                    if SHAPES {
                        extend::<W, true>(idx, a, &q)
                    } else {
                        extend_plain(idx, a, v, rc)
                    }
                }
                None => false,
            };
            let hit = if extended {
                ext += 1;
                // the rolling minimisers are recomputed after an extension
                rolling_valid = false;
                true
            } else {
                q.v = v;
                q.rc = rc;
                if SHAPES {
                    let (k, krc) = idx.kernel_of(v, rc);
                    q.kernel = k;
                    q.kernel_rc = krc;
                } else {
                    q.kernel = v;
                    q.kernel_rc = rc;
                }
                if !rolling_valid {
                    for c in levels.iter_mut() {
                        c.rolling = false;
                    }
                    rolling_valid = true;
                }
                let (hit, a) = cascade(idx, levels, &q);
                anchor = a;
                hit
            };
            hits += hit as u64;
            f(wi, hit);
        }
        self.extensions += ext;
        hits
    }
}

/// `extend_in_text` for contiguous k-mers: compare the one new base.
#[inline(always)]
fn extend_plain<W: KmerWord>(idx: &RsHash<W>, a: &mut Anchor<W>, v: W, rc: W) -> bool {
    let text = &idx.text;
    if a.forward {
        a.text_pos += 1;
        a.text_pos < a.seq_end && text.base(a.text_pos) == (v >> (2 * (idx.geo.window - 1))).low_u64()
    } else {
        if a.text_pos <= a.seq_start {
            return false;
        }
        a.text_pos -= 1;
        text.base(a.text_pos) == rc.low_u64() & 3
    }
}

/// `extend_in_text`: does the text continue the previous match by one base?
#[inline(always)]
fn extend<W: KmerWord, const SHAPES: bool>(idx: &RsHash<W>, a: &mut Anchor<W>, q: &QueryKeys<W>) -> bool {
    let text = &idx.text;
    let l = idx.geo.window;
    if !SHAPES {
        if a.forward {
            a.text_pos += 1;
            a.text_pos < a.seq_end && text.base(a.text_pos) == (q.v >> (2 * (l - 1))).low_u64()
        } else {
            if a.text_pos <= a.seq_start {
                return false;
            }
            a.text_pos -= 1;
            text.base(a.text_pos) == q.rc.low_u64() & 3
        }
    } else if a.forward {
        a.text_pos += 1;
        let b = text.base(a.text_pos);
        a.win = (a.win >> 2) | (W::from_u64(b) << (2 * (l - 1)));
        let p = a.text_pos + 1 - l as u64;
        idx.shapes.iter().enumerate().any(|(i, sh)| {
            a.win.pext(sh.w_mask) == q.fwd[i] && p + sh.start as u64 >= a.seq_start && p + sh.end as u64 <= a.seq_end
        })
    } else {
        if a.text_pos == 0 {
            return false;
        }
        a.text_pos -= 1;
        let b = text.base(a.text_pos);
        a.win = ((a.win << 2) | W::from_u64(b)) & idx.window_mask;
        let p = a.text_pos;
        idx.shapes.iter().enumerate().any(|(i, sh)| {
            a.win.pext(sh.w_mask) == q.rev[i] && p + sh.start as u64 >= a.seq_start && p + sh.end as u64 <= a.seq_end
        })
    }
}

/// The level cascade for a window that could not be extended. Returns
/// whether the window is present and the match to extend from (none for
/// last-level hits, as in C++). The common path (rolling update, cached
/// minimiser) is inlined into the query loop; rare paths are out of line.
#[inline(always)]
fn cascade<W: KmerWord>(idx: &RsHash<W>, levels: &mut [LevelCache<W>], q: &QueryKeys<W>) -> (bool, Option<Anchor<W>>) {
    for l in 0..levels.len() {
        let level = &idx.levels[l];
        let c = &mut levels[l];
        if c.rolling {
            level.mp.update(q.kernel, q.kernel_rc, &mut c.st);
        } else {
            c.st = level.mp.find(q.kernel, q.kernel_rc);
            c.rolling = true;
        }
        let v = c.st.value;
        let resolved = v == c.current || (v != c.current_neg && fetch_bucket(idx, level, c, v));
        if resolved {
            let anchor = lookup_buffer(idx, level, &levels[l], q);
            for c2 in &mut levels[l + 1..] {
                c2.rolling = false;
            }
            for c2 in &mut levels[..l] {
                c2.current_neg = c2.st.value;
            }
            return (anchor.is_some(), anchor);
        }
    }
    for c in levels.iter_mut() {
        c.current_neg = c.st.value;
    }
    (last_level_contains(idx, q), None)
}

#[inline(never)]
fn last_level_contains<W: KmerWord>(idx: &RsHash<W>, q: &QueryKeys<W>) -> bool {
    idx.last_contains(q)
}

/// If minimiser `v` is kept at this level, decode its bucket into the cache.
#[inline(never)]
fn fetch_bucket<W: KmerWord>(idx: &RsHash<W>, level: &Level, c: &mut LevelCache<W>, v: u64) -> bool {
    match level.r.contains(v) {
        Some(rank) => {
            fill_buffer(idx, level, rank, c, false);
            c.current = v;
            true
        }
        None => false,
    }
}

#[allow(clippy::collapsible_if)]
/// `lookup_buffer` + `check_minimiser_pos(2)`: check the left (and, on
/// ties, right) candidates of every cached occurrence, in C++ order.
#[inline(never)]
fn lookup_buffer<W: KmerWord>(
    idx: &RsHash<W>,
    level: &Level,
    c: &LevelCache<W>,
    q: &QueryKeys<W>,
) -> Option<Anchor<W>> {
    let st = c.st;
    let span = level.mp.span() as usize;
    let km = (idx.geo.kernel - level.mp.m) as usize;
    let tie = st.has_tie(idx.geo.kernel, level.mp.m);
    // candidate buffer indices in C++ order: left fwd, left rc, right fwd, right rc
    let (j0, j1) = (span - 1 - st.left as usize, st.left as usize);
    let (j2, j3) = (st.right as usize, km - st.right as usize);
    let o_ = idx.geo.overlap as i64;
    let n = c.n;
    assert!(c.offsets.len() >= n && c.buffer.len() >= n * span);
    let found = |i: usize, j: usize, fwd: bool| -> Option<Anchor<W>> {
        // SAFETY: i < n, j < span, and the buffer holds n * span windows
        let win = unsafe { *c.buffer.get_unchecked(i * span + j) };
        let p = unsafe { *c.offsets.get_unchecked(i) } + j as i64 - o_;
        let (s, e) = idx.match_window(win, p, fwd, q)?;
        let text_pos = if fwd { (p + idx.geo.window as i64 - 1) as u64 } else { p as u64 };
        Some(Anchor { forward: fwd, text_pos, seq_start: s, seq_end: e, win })
    };
    if idx.shapes.is_empty() {
        // cheap equality first; the bounds check only runs on a match
        let (v, rc) = (q.v, q.rc);
        for i in 0..n {
            let b = i * span;
            // SAFETY: as above
            let w = |j: usize| unsafe { *c.buffer.get_unchecked(b + j) };
            if w(j0) == v
                && let Some(a) = found(i, j0, true)
            {
                return Some(a);
            }
            if w(j1) == rc
                && let Some(a) = found(i, j1, false)
            {
                return Some(a);
            }
            if tie {
                if w(j2) == v
                    && let Some(a) = found(i, j2, true)
                {
                    return Some(a);
                }
                if w(j3) == rc
                    && let Some(a) = found(i, j3, false)
                {
                    return Some(a);
                }
            }
        }
        None
    } else {
        for i in 0..n {
            if let Some(a) = found(i, j0, true).or_else(|| found(i, j1, false)) {
                return Some(a);
            }
            if tie && let Some(a) = found(i, j2, true).or_else(|| found(i, j3, false)) {
                return Some(a);
            }
        }
        None
    }
}

/// `fill_buffer<level>` / `fill_buffer2<level>`: decode all windows around
/// each occurrence of the bucket (and, for locate, the bounds of the sequence
/// containing the occurrence).
#[allow(clippy::explicit_counter_loop)]
fn fill_buffer<W: KmerWord>(idx: &RsHash<W>, level: &Level, rank: usize, c: &mut LevelCache<W>, with_bounds: bool) {
    let (s, e) = level.bucket(rank);
    let n = e - s;
    let span = level.mp.span() as usize;
    if c.offsets.len() < n {
        c.offsets.resize(n, 0);
        c.buffer.resize(n * span, W::ZERO);
    }
    if with_bounds {
        c.bounds.clear();
    }
    let km = (idx.geo.kernel - level.mp.m) as i64;
    let o_ = idx.geo.overlap as i64;
    let l = idx.geo.window;
    let shift = 2 * (l - 1);
    for i in 0..n {
        let o = level.offsets.get(s + i);
        let off = o as i64 - km;
        c.offsets[i] = off;
        if with_bounds {
            let (_, bs, be) = idx.text.locate(o);
            c.bounds.push((bs, be));
        }
        let p0 = off - o_;
        let buf = &mut c.buffer[i * span..(i + 1) * span];
        if p0 >= 0 {
            let p0 = p0 as u64;
            let mut win: W = idx.text.window(p0, l);
            buf[0] = win;
            // the following bases, 32 at a time (as `get_word64` in C++)
            let mut next = p0 + l as u64;
            let mut bits = idx.text.word64(next);
            let mut avail = 32u32;
            for slot in buf.iter_mut().skip(1) {
                if avail == 0 {
                    bits = idx.text.word64(next);
                    avail = 32;
                }
                win = (win >> 2) | (W::from_u64(bits & 3) << shift);
                *slot = win;
                bits >>= 2;
                avail -= 1;
                next += 1;
            }
        } else {
            for (j, slot) in buf.iter_mut().enumerate() {
                let p = (p0 + j as i64).max(0) as u64;
                *slot = idx.text.window(p, l);
            }
        }
    }
    c.n = n;
}

/// Streaming locate: rolling minimisers per level (no extension), reporting
/// all occurrences of every window. As in `streaming_locate`, the decoded
/// bucket of the current minimiser (windows and sequence bounds of each
/// occurrence) is cached across consecutive windows.
pub struct StreamingLocate<'a, W: KmerWord> {
    idx: &'a RsHash<W>,
    levels: Vec<LevelCache<W>>,
    hits: Vec<Hit>,
}

impl<'a, W: KmerWord> StreamingLocate<'a, W> {
    pub fn new(idx: &'a RsHash<W>) -> Self {
        Self { idx, levels: idx.levels.iter().map(LevelCache::new).collect(), hits: Vec::new() }
    }

    /// Locate every window of `seq` (ASCII); `f(window index, hits)` is called
    /// for each window. Returns the number of windows with at least one hit.
    pub fn query(&mut self, seq: &[u8], mut f: impl FnMut(usize, &[Hit])) -> u64 {
        let idx = self.idx;
        for c in &mut self.levels {
            c.rolling = false;
        }
        let mut found = 0u64;
        let windows = Windows::<W, _>::new(seq.iter().map(|&c| alphabet::rank(c)), idx.geo.window);
        let mut q = QueryKeys::empty();
        for (wi, (v, rc)) in windows.enumerate() {
            idx.fill_query_keys(v, rc, &mut q);
            self.hits.clear();
            let mut resolved = false;
            for l in 0..self.levels.len() {
                let level = &idx.levels[l];
                let c = &mut self.levels[l];
                if c.rolling {
                    level.mp.update(q.kernel, q.kernel_rc, &mut c.st);
                } else {
                    c.st = level.mp.find(q.kernel, q.kernel_rc);
                    c.rolling = true;
                }
                let val = c.st.value;
                let hit_level = if val == c.current {
                    true
                } else if val != c.current_neg {
                    match level.r.contains(val) {
                        Some(rank) => {
                            fill_buffer(idx, level, rank, c, true);
                            c.current = val;
                            true
                        }
                        None => false,
                    }
                } else {
                    false
                };
                if hit_level {
                    Self::locate_buffer(idx, level, c, &q, &mut self.hits);
                    for c2 in &mut self.levels[l + 1..] {
                        c2.rolling = false;
                    }
                    for c2 in &mut self.levels[..l] {
                        c2.current_neg = c2.st.value;
                    }
                    resolved = true;
                    break;
                }
            }
            if !resolved {
                for c in &mut self.levels {
                    c.current_neg = c.st.value;
                }
                idx.locate_last(&q, &mut self.hits);
            }
            if self.hits.len() > 1 {
                self.hits.sort_unstable();
                self.hits.dedup();
            }
            found += !self.hits.is_empty() as u64;
            f(wi, &self.hits);
        }
        found
    }

    /// `locate_buffer` + `report_minimiser_pos(2)`: every candidate of every
    /// occurrence, checked against the cached windows and sequence bounds.
    #[inline]
    fn locate_buffer(idx: &RsHash<W>, level: &Level, c: &LevelCache<W>, q: &QueryKeys<W>, out: &mut Vec<Hit>) {
        let st = c.st;
        let span = level.mp.span() as usize;
        let km = (idx.geo.kernel - level.mp.m) as usize;
        let o_ = idx.geo.overlap as i64;
        let l = idx.geo.window as i64;
        let mut js = [(0usize, true); 4];
        js[0] = (span - 1 - st.left as usize, true);
        js[1] = (st.left as usize, false);
        let ncand = if st.has_tie(idx.geo.kernel, level.mp.m) {
            js[2] = (st.right as usize, true);
            js[3] = (km - st.right as usize, false);
            4
        } else {
            2
        };
        for i in 0..c.n {
            let off = c.offsets[i];
            let (bs, be) = c.bounds[i];
            let (bs, be) = (bs as i64, be as i64);
            let buf = &c.buffer[i * span..(i + 1) * span];
            for &(j, fwd) in &js[..ncand] {
                let win = buf[j];
                let p = off + j as i64 - o_;
                let ok = if idx.shapes.is_empty() {
                    win == if fwd { q.v } else { q.rc } && p >= bs && p + l <= be
                } else {
                    idx.shapes.iter().enumerate().any(|(si, sh)| {
                        win.pext(sh.w_mask) == if fwd { q.fwd[si] } else { q.rev[si] }
                            && p + sh.start as i64 >= bs
                            && p + sh.end as i64 <= be
                    })
                };
                if ok {
                    out.push(Hit { pos: p as u64, forward: fwd });
                }
            }
        }
    }
}
