//! Auxiliary tools (ports of the upstream `test/` programs):
//!
//! * `hashtable` — brute-force hash-set oracle (`test/hashtable.cpp`);
//! * `datasets`  — plant random reads into a text and mutate them
//!   (`test/test_datasets.cpp`, "artificial ds for gapped kmers");
//! * `preproc`   — split sequences at non-ACGT characters (`test/preprocess.cpp`);
//! * `ministats` — minimiser frequency / k-mer coverage table (`test/mini_stats.cpp`).

#![allow(clippy::manual_checked_ops)]

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use rshash::alphabet;
use rshash::hash::SEEDS;
use rshash::input;
use rshash::minimizer::{MinimizerParams, min_occurrences, window_minimizers};
use rshash::shape::Shapes;
use rshash::text::Text;
use rshash::util::SplitMix64;
use rshash::word::{KmerWord, windows};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn parse_shape(s: &str) -> Result<u64, String> {
    let s = s.trim();
    let r = if let Some(h) = s.strip_prefix("0x") {
        u64::from_str_radix(&h.replace('_', ""), 16)
    } else if let Some(b) = s.strip_prefix("0b") {
        u64::from_str_radix(&b.replace('_', ""), 2)
    } else {
        s.parse::<u64>()
    };
    r.map_err(|e| format!("invalid shape '{s}': {e}"))
}

#[derive(Parser)]
#[command(name = "rshash-tools", version, about = "RSHash auxiliary tools")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Hash-set oracle: count query k-mers (or gapped k-mers) present in the input.
    Hashtable {
        #[arg(short = 'i', long = "input")]
        input: PathBuf,
        #[arg(short = 'q', long = "query")]
        query: PathBuf,
        #[arg(short = 'k', long = "k-mer", default_value_t = 31)]
        k: u32,
        #[arg(long = "shapes", value_parser = parse_shape, value_delimiter = ',', num_args = 1..)]
        shapes: Vec<u64>,
        /// Use the index's semantics: aligned shape windows, both strands,
        /// shape spans within one sequence (num_positive_kmers then equals
        /// `rshash lookup`). Default: the C++ tool's semantics (per-shape
        /// windows, forward strand only).
        #[arg(long = "aligned")]
        aligned: bool,
    },
    /// Plant `n` random reads of length `l` into the text, then mutate the reads.
    Datasets {
        #[arg(short = 'i', long = "input")]
        input: PathBuf,
        #[arg(short = 'l', long = "length", default_value_t = 100)]
        length: usize,
        #[arg(short = 'n', long = "number", default_value_t = 10000)]
        number: usize,
        /// error rate in percent
        #[arg(short = 'e', long = "error", default_value_t = 5)]
        error: u64,
        #[arg(long = "seed", default_value_t = 1)]
        seed: u64,
    },
    /// Split sequences at non-ACGT characters, keeping runs of at least k bases.
    Preproc {
        #[arg(short = 'i', long = "input")]
        input: PathBuf,
        #[arg(short = 'o', long = "output")]
        output: PathBuf,
        #[arg(short = 'k', long = "k-mer")]
        k: usize,
    },
    /// Minimiser frequency distribution and k-mer coverage per threshold.
    Ministats {
        #[arg(short = 'i', long = "input")]
        input: PathBuf,
        #[arg(short = 'k', long = "kmer", default_value_t = 31)]
        k: u32,
        #[arg(short = 'm', long = "mini", default_value_t = 18)]
        m: u32,
        #[arg(short = 't', long = "thres", default_value_t = 64)]
        t: usize,
    },
}

fn main() {
    let cli = Cli::parse();
    let r = match cli.cmd {
        Cmd::Hashtable { input, query, k, shapes, aligned } => hashtable(&input, &query, k, &shapes, aligned),
        Cmd::Datasets { input, length, number, error, seed } => datasets(&input, length, number, error, seed),
        Cmd::Preproc { input, output, k } => preproc(&input, &output, k),
        Cmd::Ministats { input, k, m, t } => ministats(&input, k, m, t),
    };
    if let Err(e) = r {
        eprintln!("[ERROR] {e:#}");
        std::process::exit(1);
    }
}

fn read(path: &Path) -> Result<Vec<Vec<u8>>> {
    input::read_sequences(path).with_context(|| format!("reading {}", path.display()))
}

fn ranks(s: &[u8]) -> Vec<u8> {
    s.iter().map(|&c| alphabet::rank(c)).collect()
}

fn write_fasta(path: &Path, seqs: &[Vec<u8>], id: impl Fn(usize) -> String) -> Result<()> {
    let mut w =
        std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?);
    for (i, s) in seqs.iter().enumerate() {
        writeln!(w, ">{}", id(i))?;
        w.write_all(s)?;
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}

// ------------------------------------------------------------------ hashtable

fn hashtable(input: &Path, query: &Path, k: u32, shapes: &[u64], aligned: bool) -> Result<()> {
    println!("loading text...");
    let text = read(input)?;
    println!("building hashtable...");
    let queries_loaded;
    let t0;
    let (mut kmers, mut found_kmers, mut found_queries) = (0u64, 0u64, 0u64);
    if shapes.is_empty() {
        if k > 63 {
            bail!("k must be <= 63");
        }
        let mut ht: HashSet<u128> = HashSet::new();
        for s in &text {
            for (v, r) in windows::<u128>(&ranks(s), k) {
                ht.insert(v.min(r));
            }
        }
        println!("loading queries...");
        let queries = read(query)?;
        queries_loaded = queries.len();
        println!("querying...");
        t0 = Instant::now();
        for q in &queries {
            let mut qf = false;
            for (v, r) in windows::<u128>(&ranks(q), k) {
                let f = ht.contains(&v.min(r));
                qf |= f;
                found_kmers += f as u64;
                kmers += 1;
            }
            found_queries += qf as u64;
        }
    } else if !aligned {
        // C++ semantics: per-shape window of the shape's own length, forward keys
        let sh = Shapes::new(shapes)?;
        let mut hts: Vec<HashSet<u128>> = vec![HashSet::new(); sh.shapes.len()];
        for (i, a) in sh.shapes.iter().enumerate() {
            for s in &text {
                for (v, _) in windows::<u128>(&ranks(s), a.shape.length) {
                    hts[i].insert(v.pext(a.shape.mask));
                }
            }
        }
        println!("loading queries...");
        let queries = read(query)?;
        queries_loaded = queries.len();
        println!("querying...");
        t0 = Instant::now();
        for q in &queries {
            let mut found = false;
            let r = ranks(q);
            for (i, a) in sh.shapes.iter().enumerate() {
                if found {
                    break;
                }
                for (v, _) in windows::<u128>(&r, a.shape.length) {
                    if hts[i].contains(&v.pext(a.shape.mask)) {
                        found_kmers += 1;
                        found = true;
                    }
                    kmers += 1;
                }
            }
            found_queries += found as u64;
        }
    } else {
        // index semantics: aligned windows of L bases, both strands
        let sh = Shapes::new(shapes)?;
        let t = Text::from_sequences(&text);
        let l = sh.length;
        let mut hts: Vec<HashSet<u128>> = vec![HashSet::new(); sh.shapes.len()];
        for (i, a) in sh.shapes.iter().enumerate() {
            for sid in 0..t.num_sequences() {
                let (s, e) = t.sequence_bounds(sid);
                let span = (a.end - a.start) as u64;
                if e - s < span {
                    continue;
                }
                for sp in s..=e - span {
                    let win: u128 = t.window(sp - a.start as u64, l);
                    hts[i].insert(win.pext(a.w_mask));
                }
            }
        }
        println!("loading queries...");
        let queries = read(query)?;
        queries_loaded = queries.len();
        println!("querying...");
        t0 = Instant::now();
        for q in &queries {
            let mut qf = false;
            for (v, r) in windows::<u128>(&ranks(q), l) {
                let f = sh
                    .shapes
                    .iter()
                    .enumerate()
                    .any(|(i, a)| hts[i].contains(&v.pext(a.w_mask)) || hts[i].contains(&r.pext(a.w_mask)));
                qf |= f;
                found_kmers += f as u64;
                kmers += 1;
            }
            found_queries += qf as u64;
        }
    }
    let ns = t0.elapsed().as_nanos() as f64;
    println!("==== query report:");
    println!("num_kmers = {kmers}");
    println!("num_reads = {queries_loaded}");
    println!("num_positive_kmers = {found_kmers} ({}%)", found_kmers as f64 / kmers as f64 * 100.0);
    println!("time_per_kmer = {}", ns / kmers as f64);
    println!("time_per_read = {}", ns / queries_loaded as f64);
    println!("num_found_queries = {found_queries} ({}%)", found_queries as f64 / queries_loaded as f64 * 100.0);
    Ok(())
}

// ------------------------------------------------------------------- datasets

fn datasets(input: &Path, length: usize, n: usize, error: u64, seed: u64) -> Result<()> {
    let mut rng = SplitMix64(seed);
    println!("building random sequences...");
    let mut reads: Vec<Vec<u8>> =
        (0..n).map(|_| (0..length).map(|_| b"ACGT"[rng.below(4) as usize]).collect()).collect();
    println!("loading text...");
    let mut text: Vec<Vec<u8>> = read(input)?.into_iter().map(|s| alphabet::to_ascii(&ranks(&s))).collect();
    let text_length: usize = text.iter().map(|s| s.len()).sum();
    println!("inserting random sequences into text...");
    if n > 0 {
        let distance = text_length / n;
        let (mut seq_no, mut position) = (0usize, 0usize);
        for r in &reads {
            if seq_no >= text.len() {
                bail!("text too short for {n} reads of length {length}");
            }
            if text[seq_no].len() < length {
                bail!("text sequence shorter than read length");
            }
            if position + length > text[seq_no].len() {
                seq_no += 1;
                position = 0;
                if seq_no >= text.len() || text[seq_no].len() < length {
                    bail!("text sequence shorter than read length");
                }
            }
            text[seq_no][position..position + length].copy_from_slice(r);
            position += distance;
        }
    }
    let out_text = PathBuf::from(format!("{}.test.fasta", input.display()));
    write_fasta(&out_text, &text, |i| format!("seq={i}"))?;
    println!("writing errors into random sequences...");
    for r in &mut reads {
        for c in r.iter_mut() {
            if rng.below(100) < error {
                *c = b"ACGT"[(alphabet::rank(*c) ^ 3) as usize];
            }
        }
    }
    let out_reads = PathBuf::from(format!("{}.reads.fasta", input.display()));
    write_fasta(&out_reads, &reads, |i| format!("seq={i}"))?;
    Ok(())
}

// -------------------------------------------------------------------- preproc

fn preproc(input: &Path, output: &Path, k: usize) -> Result<()> {
    println!("loading text...");
    let seqs = read(input)?;
    println!("parsing sequences...");
    let (mut valid_length, mut valid_kmers, mut invalid, mut total) = (0u64, 0u64, 0u64, 0u64);
    let mut out: Vec<Vec<u8>> = Vec::new();
    for rec in &seqs {
        let mut cur: Vec<u8> = Vec::new();
        for &c in rec {
            if alphabet::is_acgt(c) {
                cur.push(b"ACGT"[alphabet::rank(c) as usize]);
                valid_length += 1;
            } else {
                invalid += 1;
                if cur.len() >= k {
                    valid_kmers += (cur.len() - k + 1) as u64;
                    out.push(std::mem::take(&mut cur));
                }
                cur.clear();
            }
        }
        if cur.len() >= k {
            valid_kmers += (cur.len() - k + 1) as u64;
            out.push(cur);
        }
        total += rec.len() as u64;
    }
    println!("no sequences: {}", seqs.len());
    println!("input text length: {total}");
    println!("valid text length: {valid_length}");
    println!("valid kmers: {valid_kmers}");
    println!("invalid kmers: {invalid} {}%", invalid as f64 / valid_length as f64 * 100.0);
    println!("output sequences: {}", out.len());
    println!("saving output...");
    write_fasta(output, &out, |i| format!("seq_{i}"))
}

// ------------------------------------------------------------------ ministats

fn ministats(input: &Path, k: u32, m: u32, t: usize) -> Result<()> {
    let (text, _) = input::read_text(input)?;
    let mp = MinimizerParams::new(k, m, SEEDS[0]);
    let mut freq: HashMap<u64, u32> = HashMap::new();
    let mut kmers = 0u64;
    for i in 0..text.num_sequences() {
        let (s, e) = text.sequence_bounds(i);
        min_occurrences(&text, &mp, s, e, |v, _| *freq.entry(v).or_default() += 1);
        kmers += (e - s + 1).saturating_sub(k as u64);
    }
    // histogram of kept minimisers (frequency <= t)
    let mut counter = vec![0u64; t];
    let mut no_minimizers = 0u64;
    for &c in freq.values() {
        if (c as usize) <= t {
            counter[c as usize - 1] += 1;
            no_minimizers += 1;
        }
    }
    // k-mers whose minimiser has frequency <= i+1
    let mut covered = vec![0u64; t + 1];
    for i in 0..text.num_sequences() {
        let (s, e) = text.sequence_bounds(i);
        window_minimizers(&text, &mp, s, e, |_, v| {
            let c = freq[&v] as usize;
            if c <= t {
                covered[c] += 1;
            }
        });
    }
    println!("frequency, superkmers, kmers covered");
    let mut cum = 0u64;
    for j in 0..t {
        cum += covered[j + 1];
        println!(
            "{},{},{}",
            j + 1,
            counter[j] as f64 / no_minimizers.max(1) as f64 * 100.0,
            cum as f64 / kmers.max(1) as f64 * 100.0
        );
    }
    Ok(())
}
