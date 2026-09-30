# rshash-rs

A pure Rust implementation of the **RSHash** k-mer dictionary, ported from the upstream C++ implementation at **[jonsmcode/rshash](https://github.com/jonsmcode/rshash)** (branch `multishapes`, commit [`b273cea`](https://github.com/jonsmcode/rshash/commit/b273ceac0d2ef8499ea474ed1d1bb93e02ee43e1)).

RSHash is a compressed, exact k-mer dictionary:
- It replaces the minimal perfect hash functions of SSHash with Elias–Fano bitvectors that support rank and select.
- It handles skewed minimiser distributions with up to three minimiser levels plus a last level for the remaining k-mers.
- It answers membership queries (point and streaming) and, optionally, *locate* queries.

This port includes the gapped k-mer / k-mer pattern (**shapes**) support of the upstream `multishapes` branch. It is also generic over windows of up to 63 bases.
- **Kept identical to C++:** everything that defines the data structure (hash functions and seeds, minimiser rules, level cascade, text layout, shape alignment).
- **Fixed:** bugs and unfinished code paths on that branch now follow their intended semantics.
- **Documented:** every difference is listed in [DEVIATIONS.md](DEVIATIONS.md).

## Citation

If you use RSHash, please cite the original paper:

> Jonas Schulte-Mattler and Knut Reinert. **RSHash: a fast and space-efficient hash table for k-mers.**
> *Bioinformatics* 42(Suppl 2):btag440, 2026. doi:[10.1093/bioinformatics/btag440](https://doi.org/10.1093/bioinformatics/btag440)
> (PMID [42635238](https://pubmed.ncbi.nlm.nih.gov/42635238/), [PMC13501323](https://pmc.ncbi.nlm.nih.gov/articles/PMC13501323/))

```bibtex
@article{schulte-mattler2026rshash,
  author  = {Schulte-Mattler, Jonas and Reinert, Knut},
  title   = {{RSHash}: a fast and space-efficient hash table for k-mers},
  journal = {Bioinformatics},
  volume  = {42},
  number  = {Supplement_2},
  pages   = {btag440},
  year    = {2026},
  doi     = {10.1093/bioinformatics/btag440}
}
```

## Build

```bash
cargo build --release           # binaries: target/release/{rshash, rshash-tools}
cargo test --release            # unit tests + brute-force exactness tests
```

`.cargo/config.toml` builds with `-C target-cpu=native`, like the C++ `-march=native`, so the BMI2 `pext`/`pdep` instructions are used. Portable fallbacks exist.

## Command line

The command line matches the C++ binary:

```bash
# contiguous 31-mers, 2 minimiser levels
rshash build  -i genome.fa -d genome.rsh -k 31 --level 2 --m1 18 --m2 21 --t1 64 --t2 64
rshash lookup -d genome.rsh -q reads.fq            # streaming lookup
rshash lookup -d genome.rsh -q reads.fq --dump per_read.tsv
rshash locate -d genome.rsh -q reads.fq --output hits.tsv   # needs --loc at build
rshash bench  -d genome.rsh                        # random positive/negative point lookups

# gapped k-mers: one or more shapes (repeat --shapes or use commas)
rshash build -i genome.fa -d gapped.rsh --shapes 0x5FFFFFFD --shapes 0x1FFFFFFF --loc
```

| Option | Default | Meaning |
|---|---|---|
| `-k` | 31 | k-mer length (1–63); ignored with `--shapes` |
| `--level` | 2 | number of minimiser levels (1–3) |
| `--m1 --m2 --m3` | 18 21 23 | minimiser length per level (≤ 31, ≤ kernel length) |
| `--t1 --t2 --t3` | 64 | maximum bucket size per level; minimisers occurring more often move to the next level |
| `-t` | 0 | last level: drop k-mers occurring at least `t` times (0 = keep all) |
| `--loc` | off | build locate support |
| `--ht` | off | last level as Elias–Fano instead of hash tables (inverted, as in C++) |
| `--shapes` | – | shape masks: decimal, `0x…` or `0b…` |

### Shapes

Bit *i* of a shape (from the least significant bit) marks base *i* of the window as a care position. Because windows store their first base in the lowest bits, the binary literal reads right-to-left in sequence order.

For example, `0x5FFFFFFD` = `1011111111111111111111111111101` is a palindromic shape:

- window length 31;
- 27-base contiguous kernel;
- weight 29.

With several shapes:

- The kernel `K` is the shortest of their longest runs of ones.
- Every shape is aligned so that its run starts at the kernel.
- The common window has `L = K + 2·O` bases.
- Minimisers are computed on the kernel only.

A window matches if, for at least one shape, its care positions equal those of an indexed text window whose shape span lies in one input sequence. The text window may match the query or the query's reverse complement.

### Tools (`rshash-tools`)

Ports of the upstream `test/` programs:

| Tool | Purpose |
|---|---|
| `hashtable` | Brute-force hash-set oracle. `--aligned` uses the index's shape semantics; without it, the C++ tool's per-shape, forward-only semantics. |
| `datasets` | Plants `-n` random reads of length `-l` into a text, then mutates them with `-e`% errors ("artificial ds for gapped kmers"). |
| `preproc` | Splits sequences at non-ACGT characters, keeping runs of at least k bases. |
| `ministats` | Minimiser frequency / k-mer coverage table. |

## Library

```rust
use rshash::{AnyIndex, Params, input};

let (text, _non_acgt) = input::read_text("genome.fa")?;
let params = Params { k: 31, levels: 2, m: [18, 21, 23], t: [64, 64, 64], loc: true, ..Params::default() };
let index = AnyIndex::build(text, &params)?;     // u64 windows up to 32 bases, u128 up to 63
index.save_file("genome.rsh")?;

rshash::with_index!(&index, idx => {
    let mut engine = idx.streaming_lookup();     // reusable across reads
    let hits = engine.query(b"ACGT...");         // number of windows present
    let mut loc = idx.streaming_locate();
    loc.query(b"ACGT...", |window, hits| { /* Hit { pos, forward } */ });
});
```

`RsHash<W>` is the index over `u64` or `u128` windows. `AnyIndex` picks the word type from the window length.

## Layout

```
crates/rshash/src/        library (one module per C++ component)
  alphabet.rs  word.rs    dna4 encoding, packed windows, rc, pext/pdep (util.hpp, kmerview)
  hash.rs                 MurmurHash2_64 + mixer_64 (minimiser_views.hpp)
  shape.rs                shapes and their alignment (shape.hpp)
  text.rs                 packed text + sequence endpoints (build.cpp, rshash.hpp)
  bits.rs  ef.rs          rank/select bitvector, packed ints, Elias–Fano (EliasFano.hpp, compact_vector.hpp)
  minimizer.rs            query-side and build-side minimisers
  index.rs                RSHash, build, serialisation (rshash.hpp, build.cpp, io.cpp)
  last_level.rs           hash-set / hash-map / Elias–Fano last level (flat_map.hpp)
  query.rs                point/streaming lookup and locate (lookup.cpp, locate.cpp)
  stats.rs                print_info report
crates/rshash/tests/      brute-force exactness tests
crates/rshash-cli/        `rshash` binary (main.cpp)
crates/rshash-tools/      `rshash-tools` binary (test/*.cpp)
parity/                   C++ oracle build, data preparation, parity harness
```

## Testing

- **Unit tests.** They cover:
  - the MurmurHash magics against golden values (verified against the C++ build);
  - reverse complements and rolling windows against naive implementations;
  - `pext`/`pdep` fallbacks against BMI2;
  - rank/select on uniform and skewed densities;
  - Elias–Fano (`contains`, `rank`, `select`, predecessor) including 128-bit universes;
  - shape geometry;
  - query-side and build-side minimisers against brute force.
- **Exactness tests** (`crates/rshash/tests/exactness.rs`). Every answer is compared with a brute-force oracle, per window:
  - point lookup;
  - streaming lookup;
  - point locate;
  - streaming locate.

  The grid covers:
  - levels 1–3;
  - large and tiny thresholds (so deep levels and the last level are populated);
  - hash-table vs Elias–Fano last level, with and without `--loc`;
  - k ∈ {9, 15, 31, 32, 33, 47, 63};
  - single, multiple, palindromic, non-palindromic, off-centre and >32-base shapes.

  The inputs include repeats, low-complexity sequences, both strands, empty and short sequences, and chimeric queries across sequence boundaries. Save/load is also checked: answers must be identical and re-serialisation byte-identical.

### Parity with the C++ implementation

```bash
parity/prepare_data.sh [--large]                 # E. coli (+ C. elegans) and sampled reads
module load gcc/14.2.0 cmake/3.29.2              # C++20 toolchain for the reference
git -C ../rshash submodule update --init --recursive
cmake -S parity/cpp -B parity/work/cpp-build -DCMAKE_BUILD_TYPE=Release -DRSHASH_SRC=$PWD/../rshash
cmake --build parity/work/cpp-build -j
parity/parity.sh                                 # builds both indexes, compares every read
```

`parity/cpp` builds the unmodified upstream sources plus a small driver (`cpp_dump.cpp`) that prints per-read hit counts. `parity.sh` compares those counts with `rshash lookup --dump` and with the `hashtable` oracle.

Results on E. coli K-12 split into 3,036 unitig-like chunks, with 100k 150-bp reads (1% errors, 20% random):

| config | oracle | C++ | Rust | reads differing |
|---|---:|---:|---:|---:|
| `-k 31 --level 1` | 6,873,370 | 6,873,370 | 6,873,370 | 0 |
| `-k 31 --level 2 --t1 4 --t2 8` | 6,873,370 | 6,873,370 | 6,873,370 | 0 |
| `-k 31 --level 3 --t1 2 --t2 2 --t3 4` | 6,873,370 | 6,873,370 | 6,873,370 | 0 |
| `-k 25 --level 2 --t1 3 --t2 3` | 7,696,248 | 7,696,248 | 7,696,248 | 0 |
| `-k 32 --level 1` | 6,743,444 | 6,743,444 | 6,743,444 | 0 |
| `-k 21 --loc` | 8,288,636 | 8,288,636 | 8,288,636 | 0 |
| `--shapes 0x5FFFFFFD` (levels 1, 2) | 7,013,210 | 7,013,210 | 7,013,210 | 0 |
| `--shapes 15855` (off-centre palindrome, D8) | 10,548,210 | 6,983,785 | 10,548,210 | expected |
| `--shapes 268435445` (non-palindromic, D8) | 7,185,032 | 7,184,933 | 7,185,032 | expected |

The same comparison on C. elegans (100 Mbp genome, 500k reads, default parameters) gives 35,164,753 positive k-mers in all three. No read differs.

### Performance

Measured single-threaded on an AMD EPYC 9555, with the upstream C++ built with gcc 14.2 `-O3 -march=native`:

| | C++ | Rust |
|---|---:|---:|
| E. coli, `-k 31 --level 1`, streaming lookup (ns/k-mer) | 9.5 | 10.6 |
| E. coli, 3 levels, streaming lookup (ns/k-mer) | 19.8 | 20.4 |
| E. coli, `--shapes 0x5FFFFFFD`, streaming lookup (ns/k-mer) | 12.6 | 14.4 |
| E. coli, space (bits/k-mer, `-k 31 --level 1`) | 7.79 | 7.79 |
| C. elegans, build time / peak memory | 3.3 s / 431 MB | 4.0 s / 342 MB |
| C. elegans, space (bits/k-mer) | 7.66 | 7.71 |
| C. elegans, streaming lookup (ns/k-mer) | 22.5 | 24.7 |
| E. coli, random point lookups, positive / negative (ns) | 98 / 68 | 110 / 59 |

Locate isn't directly comparable. The C++ version neither returns positions nor reports last-level hits (D6).

## License

rshash-rs is distributed under the BSD 3-Clause license (see [LICENSE](LICENSE)). It is a port of the upstream RSHash, which is © 2026 Jonas Schulte-Mattler and distributed under the MIT license. That license is reproduced in [LICENSE-UPSTREAM](LICENSE-UPSTREAM).
