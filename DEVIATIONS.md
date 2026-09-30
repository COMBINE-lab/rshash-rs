# Deviations from the C++ implementation

The port follows the C++ RSHash on branch `multishapes` (HEAD `b273cea`,
"longer kmers").

- **What stays the same:** everything that defines the data structure.
- **What was fixed:** behaviour that was a bug, undefined, or clearly unfinished on that branch now follows the *intended* semantics.

This file lists every such difference, so the fixes can be reported upstream. File and line references are to `rshash/source/`.

## Kept identical

These parts define the structure, so they match C++ exactly:

- **Alphabet and text layout.**
  - Alphabet is seqan3 `dna4`: A=0, C=1, G=2, T=3; complement is `^3`; non-ACGT input becomes A.
  - k-mer bit layout: the first base is in the lowest bits.
  - The text is a packed concatenation of sequences with a leading 32-base `T` padding word. Global position = 32 + offset.
- **Hashing.**
  - Minimiser hash `mixer_64`: `x * 0x517cc1b727220a95 ^ MurmurHash2_64(seed)`.
  - Per-level seeds are `1`, `0x296DBD3332568C64` and `0xE59A385F0376C9F6`. The magics were checked against the C++ build.
  - The canonical m-mer score is `min(h(fwd) & mask, h(rc) & mask)`. The score itself is the key stored in `R`.
- **Query-side minimiser** (`find_minimiser` / `update_minimiser`): leftmost/rightmost tie tracking, scan order, and rolling update are all ported verbatim.
- **Level cascade.**
  - A k-mer belongs to the first level whose bucket for its minimiser was kept (≤ `t_l` occurrences). There is no fallback to deeper levels.
  - Deeper levels only see the maximal runs of k-mers whose minimisers were frequent.
  - `R` is Elias–Fano over universe `4^m`. `S` is a bitvector of bucket delimiters. Offsets are bit-packed global positions of the minimiser occurrences.
- **Candidate positions** (`check` / `check_minimiser_pos(2)`):
  - The left candidate is checked forward and reverse-complement.
  - The right candidate is also checked when the minimum is tied at different positions.
  - Candidates are checked in the same order.
- **Streaming lookup** (Algorithm 1):
  - Extension in the text by one base.
  - Per-level rolling minimisers.
  - Caches for the current bucket (decoded windows) and the last negative minimiser.
  - Last-level hits do not seed extensions.
- **Shape geometry.**
  - Shapes are integer masks: bit *i* = base *i* of the window.
  - The kernel is the shortest of the shapes' longest runs of ones; ties go to the lowest-bit run.
  - Shapes are aligned so that their runs start at the kernel.
  - Minimisers are computed on the kernel.
- **Command line.**
  - Subcommands `build | lookup | locate | bench`, with the same options and defaults: `-k 31 --level 2 --m1 18 --m2 21 --m3 23 --t1/2/3 64 -t 0 --loc --ht --shapes`.
  - `--ht` is inverted as in C++: it selects Elias–Fano for the last level.
  - The same report labels are kept, so `benchmarks/benchmark_rshash.sh` still greps them.

## Deviations

### D1 — Build-time minimiser occurrences (`xor_minimiser_and_positions`)

C++ issues:

- When draining tied minima, the view emits the tie's position with a **stale** minimiser value (`cached.minimiser_value` is not updated). In the first window that value is 0; later it is the previous minimiser.
- These junk (value, position) pairs inflate bucket counts. That shifts minimisers across the thresholds `t_l`.
- A range shorter than the window emits one bogus `{0, 0}` item.
- The view reads one element past the end of its input.

The port emits every m-mer occurrence that is a (possibly tied) minimum of at least one k-mer, each exactly once, with its true value (`minimizer::min_occurrences`). Short ranges emit nothing.

Buckets are fully sorted by (minimiser, position). kxsort's insertion-sort cutoff in C++ compares only the minimiser, so positions within a bucket were only partially ordered.

Consequences:

- Membership answers are unchanged: the C++ query also checks the rightmost tied occurrence, which is emitted correctly.
- Per-level statistics differ slightly. For example, on C. elegans `no minimiser1` is 11,584,506 vs 11,584,496.
- The number of *extensions* differs by roughly 0.001%, because the first matching occurrence in a bucket can differ.

The super-k-mer view (`xor_minimiser_and_skmer_positions`) was only used to learn each k-mer's minimiser. It is replaced by a sliding-window minimum (`minimizer::window_minimizers`), which yields the same minimiser values.

### D2 — Threshold widths

C++ issue: `--t1/--t2/--t3` are parsed as `uint16_t`, but the `RSHash` constructor takes `uint8_t` (`rshash.hpp:185`). So `--t1 256` silently becomes 0.

The port keeps them `u16` end to end.

### D3 — Point lookup uses the level's own `m`

C++ issue: `check<level>` compares against `k - m1 - right` at every level (`lookup.cpp:306`). Levels 2 and 3 therefore skip the right tie candidate, or check the wrong one.

The port uses the level's own `m`.

### D4 — Point lookup respects sequence boundaries and uses all shapes

C++ point lookup (`lookup1..3`) has three problems:

- It does not check that the candidate window lies in a single input sequence. Windows spanning two concatenated sequences are reported as present. Streaming lookup does check this.
- It uses only `shapes[0]`.
- With shapes, it passes uninitialised `shapes_fwd`/`shapes_rev` arrays to the last level.

In the port, point and streaming lookup share the exact same membership semantics.

### D5 — Last level without hash tables (`--ht`)

C++ issues:

- `lookup_last_level` returns `false` when `!use_ht` (`lookup.cpp:98-108`, marked "todo"). So every k-mer stored in the last level is missed.
- `locate_last_level` queries `r5`, which is never built.
- With shapes, the `!use_ht` build keys `r4` on plain canonical windows, while the query probes with gapped keys.

In the port, the Elias–Fano last level works for lookup and locate, with or without shapes (one table per shape). It supports keys wider than 64 bits.

### D6 — Locate

C++ issues:

- `streaming_locate` never returns positions (every `positions.push_back` is commented out).
- The hash-table last level is a no-op, so locate misses every last-level k-mer. That is the default `--loc` build.
- `streaming_locate1/2` compute the shape kernel from an uninitialised `kmer`.
- The shape paths use the never-initialised member `overlap` and the unaligned `shape.mask`.
- `locate1..3` accumulate into an uninitialised counter.

The port's `locate` / `StreamingLocate` report every occurrence as a `Hit { pos, forward }`, deduplicated:

- `pos` is the text position where the window starts.
- `forward` says whether the window matched the query or its reverse complement.

`rshash locate --output` writes `(query, window, sequence id, offset, strand)` rows. `num_positions` counts real hits.

On E. coli (k=21, `--loc`), C++ reports 8,087,576 positive k-mers and Rust reports 8,288,636, which equals `lookup` and the brute-force oracle. The difference is exactly the last-level k-mers.

### D7 — Last-level positions and the `-t` filter

C++ issues:

- Last-level positions are truncated to `uint32_t` (`process_freq_kmers`, `FlatMap`).
- For shapes, the position is taken at the fragment start rather than at the start of the (extended) window.
- `-t` is ignored for shapes.
- `count_kmers` counts only with `shapes[0]` ("todo: multiple shapes").

The port:

- stores window starts at full width (`bit_width(text length)` bits);
- applies `-t` to every table, per shape.

`-t` semantics are unchanged: keys occurring `>= t` times in the last level are dropped.

### D8 — Strand semantics of gapped k-mers

**Intended semantics.** A query window `W` matches under shape `S` if some text window `T`, whose shape span lies in one sequence, satisfies either:

- `pext(T, S) == pext(W, S)` (forward), or
- `pext(T, S) == pext(rc(W), S)` (reverse complement).

**C++ behaviour.**

- Levels 1–3 already verify exactly this (`check_minimiser_pos`).
- The last level instead stores the canonical key `min(pext(T,S), pext(rc(T),S))` over the shape's *own* window, and probes with `min(pext(W,w_mask), pext(rc(W),w_mask))` over the *aligned* window.
- Those two keys agree only when the shape is a palindrome and sits centred in the aligned window.
- The palindrome check enforcing this in commit `56c98bc` is commented out again in `b273cea`.

**Port.** The last level stores forward keys and probes with both `pext(W)` and `pext(rc(W))`. When a shape is palindromic and centred, it stores canonical keys and probes once, which is provably equivalent.

**Parity evidence** (E. coli, 100k reads):

- `--shapes 15855` (palindromic, off-centre kernel): C++ reports 6,983,785 positive windows. Rust and the oracle both report 10,548,210.
- `--shapes 268435445` (non-palindromic): C++ misses 99 windows.

### D9 — Window overlap for multiple shapes

C++ issue: `Shapes32.overlap = max_i max(overlap_left_i, overlap_right_i)`. A shape whose longest run is longer than the kernel then extends past the window length `L = K + 2·O` (its `w_mask` spills beyond `L`).

The port uses `O = max_i max(overlap_right_i, overlap_left_i + run_i − K)`:

- Every aligned shape fits in the window.
- The kernel stays centred, so kernel minimisers stay strand-symmetric.
- When every run equals `K`, this equals the C++ value.

### D10 — Robustness fixes

**Empty level crash.** A level that receives no k-mers indexes `minimizers[0]` of an empty vector (`filter_freq_minimizers`).

- Example: E. coli with default `-k 31 --level 2` segfaults in `get_minimizers<2>` because no level-1 minimiser is frequent.
- The port handles empty levels.

**Short and empty sequences.**

- C++ problems:
  - `size - k + 1` underflows in `mark_sequences` and in the CLI's k-mer counts.
  - The minimiser view emits a bogus item for sequences shorter than `k`.
  - `get_frequent_skmers` reads an uninitialised `cur_freq`.
- The port:
  - skips short sequences;
  - supports empty ones, with duplicate endpoints allowed;
  - counts `max(0, len - L + 1)` windows per query.

**Out-of-bounds reads.**

- `extend_in_text` reads the base before its bound check, and `kmerview` dereferences the end iterator.
- The port checks bounds first.

**Uninitialised state.**

- In plain mode the sentinel `Shape32` leaves `kernel_length`, overlaps and `w_mask` uninitialised. The member `overlap` is serialised but never set.
- The port represents plain mode as "no shapes" (`Option<Shapes>`).

**Report fixes.**

- `print_info` uses `bit_width(text.size()*32)` as the offset width, but the offsets are stored with `bit_width(endpoints.size())`.
- `number_unitigs` counts the two sentinels.
- `text length` includes padding.
- The port reports actual sizes.

**`bench` negatives.** With shapes, negatives are drawn with the window length `L`; C++ uses the kernel length.

### D11 — Generalisations

**Longer windows.** Windows (k or `L`) of up to 63 bases are supported by the index itself, with `u128` words chosen automatically at build time.

- In C++, `Shape64`, `Shapes64` and `longkmerview` are only used by `test/hashtable.cpp`.
- The limits stay the same: shape weight ≤ 32, and m ≤ 31 (C++ allows m = 32, which overflows `1ULL << 2m`).

**More shapes.** Up to 8 shapes are allowed; C++ dispatches at most 3.

**Shape syntax.**

- `--shapes` takes u64 values in decimal, `0x…` or `0b…` notation, repeated or comma-separated. C++ accepts `uint32` decimal only.
- The shape set is printed at build time.

### D12 — Index file format

- The port uses its own little-endian format: magic `RSHASH01`, a version, and parameters, followed by every component.
- Hash-table keys are written sorted, so files are deterministic and re-serialisation is byte-identical.
- Select and hash structures are rebuilt on load.
- C++ cereal files cannot be read. Byte compatibility was not a goal, because the C++ files contain uninitialised fields and hash-table iteration order.

### D13 — Implementation substitutions (no semantic effect)

| C++ | Port |
|---|---|
| sux `EliasFano` fork | `ef::EliasFano`: same `l = floor(log2(u/n))` layout, `contains → rank`, `select`, `select(r, &next)`, predecessor |
| sux `SimpleSelect`/`SimpleSelectZeroHalf`, sdsl `bit_vector` | `bits::RankSelect`: rank9-style directory, plus a sux-style interleaved select inventory (position of every 256th target bit and 16-bit offsets of every 64th) built only for the select direction each structure uses |
| pthash `compact_vector` | `bits::CompactVec` |
| gtl `flat_hash_set`/`FlatMap` | std hash set/map with a murmur3-`fmix64` hasher, behind a derived one-hash bit filter (8–16 bits/key) that screens absent keys; locate buckets as bitvector plus packed positions |
| seqan3 `sequence_file_input` | needletail (FASTA/FASTQ, gz/bz2/xz/zstd) |
| kxsort | `sort_unstable` on 12-byte (minimiser, position) pairs when positions fit in 32 bits (as `MinimizerInfo32`), else 16-byte pairs |
| cereal | own serialisation (D12) |

Two more differences:

- The streaming lookup engine keeps its caches across query records; C++ reallocates them per record. Results are unaffected.
- `rshash lookup` and `locate` load the index *before* the queries. On hosts with transparent huge pages this gives the index arrays better page placement, which measurably speeds up queries.
- Report labels match C++ where their meaning is the same:
  - `no distinct minimiser{i}` prints the number of occurrences; C++ prints `s.size()`, which is occurrences + 1.
  - `no freq kmers` prints the number of distinct last-level keys, as in C++.
  - An extra line, `last level kmer occurrences`, gives the number of k-mer positions routed to the last level.
