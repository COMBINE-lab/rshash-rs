//! `rshash` command line (port of `source/main.cpp`).
//!
//! Same positional command, options, defaults and report labels as the C++
//! `rsindex` binary:
//!
//! ```text
//! rshash build  -i text.fa -d index.rsh [-k 31] [--level 2] [--m1 18 --m2 21 --m3 23]
//!               [--t1 64 --t2 64 --t3 64] [-t 0] [--loc] [--ht] [--shapes S ...]
//! rshash lookup -d index.rsh -q reads.fq
//! rshash locate -d index.rsh -q reads.fq [--output hits.tsv]
//! rshash bench  -d index.rsh
//! ```

use anyhow::{Context, Result, bail};
use clap::{Parser, ValueEnum};
use rshash::query::Hit;
use rshash::{AnyIndex, KmerWord, Params, RsHash, input, with_index};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq)]
enum Cmd {
    Build,
    Lookup,
    Locate,
    Bench,
}

fn parse_shape(s: &str) -> Result<u64, String> {
    let s = s.trim();
    let r = if let Some(h) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u64::from_str_radix(&h.replace('_', ""), 16)
    } else if let Some(b) = s.strip_prefix("0b").or_else(|| s.strip_prefix("0B")) {
        u64::from_str_radix(&b.replace('_', ""), 2)
    } else {
        s.parse::<u64>()
    };
    r.map_err(|e| format!("invalid shape '{s}': {e}"))
}

#[derive(Parser, Debug)]
#[command(name = "rshash", version, about = "RSHash: a k-mer dictionary based on rank/select bitvectors (Rust port)")]
struct Args {
    /// command: build, lookup, locate, bench
    #[arg(value_enum)]
    cmd: Cmd,
    /// input file (FASTA/FASTQ, optionally compressed)
    #[arg(short = 'i', long = "input")]
    input: Option<PathBuf>,
    /// query file (FASTA/FASTQ, optionally compressed)
    #[arg(short = 'q', long = "query")]
    query: Option<PathBuf>,
    /// index file
    #[arg(short = 'd', long = "dict")]
    dict: PathBuf,
    /// k-mer length (ignored with --shapes)
    #[arg(short = 'k', long = "k-mer", default_value_t = 31)]
    k: u32,
    /// number of minimiser levels (1..3)
    #[arg(long = "level", default_value_t = 2)]
    level: u32,
    /// minimiser1 length
    #[arg(long = "m1", default_value_t = 18)]
    m1: u32,
    /// minimiser2 length
    #[arg(long = "m2", default_value_t = 21)]
    m2: u32,
    /// minimiser3 length
    #[arg(long = "m3", default_value_t = 23)]
    m3: u32,
    /// threshold1 (max minimiser bucket size of level 1)
    #[arg(long = "t1", default_value_t = 64)]
    t1: u16,
    /// threshold2
    #[arg(long = "t2", default_value_t = 64)]
    t2: u16,
    /// threshold3
    #[arg(long = "t3", default_value_t = 64)]
    t3: u16,
    /// max k-mer/shape frequency threshold on the last level (0 = off)
    #[arg(short = 't', default_value_t = 0)]
    t: u16,
    /// enable locate
    #[arg(long = "loc")]
    loc: bool,
    /// do not use hashtable on last level (use Elias-Fano instead)
    #[arg(long = "ht")]
    ht: bool,
    /// shape values (repeat the option or separate by commas; decimal, 0x.., 0b..);
    /// bit i = base i of the window (first base = least significant bit)
    #[arg(long = "shapes", value_parser = parse_shape, value_delimiter = ',', num_args = 1..)]
    shapes: Vec<u64>,
    /// locate: write every hit as TSV (query, window, sequence, offset, strand)
    #[arg(long = "output")]
    output: Option<PathBuf>,
    /// lookup: write per-record counts as TSV (record, length, windows, hits),
    /// the format of parity/cpp/cpp_dump.cpp
    #[arg(long = "dump")]
    dump: Option<PathBuf>,
    /// bench: number of k-mers per round
    #[arg(long = "bench-kmers", default_value_t = 1_000_000)]
    bench_kmers: usize,
}

fn main() {
    let args = Args::parse();
    let r = match args.cmd {
        Cmd::Build => build(&args),
        Cmd::Lookup => lookup(&args),
        Cmd::Locate => locate(&args),
        Cmd::Bench => bench(&args),
    };
    if let Err(e) = r {
        eprintln!("[ERROR] {e:#}");
        std::process::exit(1);
    }
}

fn build(args: &Args) -> Result<()> {
    let Some(input_path) = &args.input else { bail!("provide input file.") };
    let params = Params {
        k: args.k,
        levels: args.level,
        m: [args.m1, args.m2, args.m3],
        t: [args.t1, args.t2, args.t3],
        threshold: args.t,
        loc: args.loc,
        use_ht: !args.ht,
        shapes: args.shapes.clone(),
    };
    let geo = params.geometry()?;
    if let Some(s) = &geo.shapes {
        println!("{s}");
    }
    println!("loading text...");
    let t0 = Instant::now();
    let (text, non_acgt) = input::read_text(input_path).with_context(|| format!("reading {}", input_path.display()))?;
    println!("text length: {}", text.len());
    println!("no sequences: {}", text.num_sequences());
    if non_acgt > 0 {
        eprintln!("warning: {non_acgt} non-ACGT characters were converted to A");
    }
    println!("building dict...");
    let mut index = AnyIndex::build_with_log(text, &params, &mut |m| println!("{m}"))?;
    index.set_non_acgt(non_acgt);
    println!("{}", index.report());
    index.save_file(&args.dict).with_context(|| format!("writing {}", args.dict.display()))?;
    eprintln!("build time: {:.2}s", t0.elapsed().as_secs_f64());
    Ok(())
}

fn load(args: &Args) -> Result<AnyIndex> {
    AnyIndex::load_file(&args.dict).with_context(|| format!("loading {}", args.dict.display()))
}

fn num_windows(len: usize, l: u32) -> u64 {
    (len + 1).saturating_sub(l as usize) as u64
}

fn lookup(args: &Args) -> Result<()> {
    let Some(qpath) = &args.query else { bail!("provide query file.") };
    // load the index before the queries (see `locate`)
    let index = load(args)?;
    println!("loaded index...");
    println!("loading queries...");
    let queries = input::read_sequences(qpath).with_context(|| format!("reading {}", qpath.display()))?;
    println!("querying...");
    let l = index.window();
    let mut dump: Option<Box<dyn Write>> = match &args.dump {
        Some(p) => Some(Box::new(std::io::BufWriter::new(std::fs::File::create(p)?))),
        None => None,
    };
    let (kmers, found, extensions, ns) = with_index!(&index, idx => {
        let mut eng = idx.streaming_lookup();
        let t0 = Instant::now();
        let mut found = 0u64;
        let mut kmers = 0u64;
        for (i, q) in queries.iter().enumerate() {
            let h = eng.query(q);
            let w = num_windows(q.len(), l);
            if let Some(d) = dump.as_mut() {
                writeln!(d, "{i}\t{}\t{w}\t{h}", q.len())?;
            }
            found += h;
            kmers += w;
        }
        (kmers, found, eng.extensions, t0.elapsed().as_nanos())
    });
    if let Some(d) = dump.as_mut() {
        d.flush()?;
    }
    println!("==== query report:");
    println!("num_kmers = {kmers}");
    println!("num_positive_kmers = {found} ({}%)", found as f64 / kmers as f64 * 100.0);
    println!("time_per_kmer = {}", ns as f64 / kmers as f64);
    println!("extensions = {extensions}");
    Ok(())
}

fn locate(args: &Args) -> Result<()> {
    let Some(qpath) = &args.query else { bail!("provide query file.") };
    // load the index before the queries: its large arrays then get better
    // (transparent huge page) placement, measurably speeding up queries
    println!("loading dict...");
    let index = load(args)?;
    println!("loading queries...");
    let queries = input::read_sequences(qpath).with_context(|| format!("reading {}", qpath.display()))?;
    if !index.has_locate() {
        bail!("index does not support locate");
    }
    let mut out: Option<Box<dyn Write>> = match &args.output {
        Some(p) => Some(Box::new(std::io::BufWriter::new(std::fs::File::create(p)?))),
        None => None,
    };
    println!("locating...");
    let l = index.window();
    let (kmers, found, positions, ns) = with_index!(&index, idx => {
        let mut eng = idx.streaming_locate();
        let t0 = Instant::now();
        let (mut found, mut kmers, mut positions) = (0u64, 0u64, 0u64);
        for (qi, q) in queries.iter().enumerate() {
            found += eng.query(q, |wi, hits: &[Hit]| {
                positions += hits.len() as u64;
                if let Some(o) = out.as_mut() {
                    for h in hits {
                        let (sid, off) = idx.hit_position(h);
                        let _ = writeln!(o, "{qi}\t{wi}\t{sid}\t{off}\t{}", if h.forward { '+' } else { '-' });
                    }
                }
            });
            kmers += num_windows(q.len(), l);
        }
        (kmers, found, positions, t0.elapsed().as_nanos())
    });
    if let Some(o) = out.as_mut() {
        o.flush()?;
    }
    println!("==== query report:");
    println!("num_kmers = {kmers}");
    println!("num_positive_kmers = {found}");
    println!("num_positions = {positions}");
    println!("time_per_kmer = {}", ns as f64 / kmers as f64);
    Ok(())
}

/// Repeat timed rounds like the C++ `bench` until the running error is below
/// 5% of the mean (at least 10, at most 51 rounds). Returns (mean, sd, last count).
fn bench_rounds(label: &str, mut round_fn: impl FnMut(u64) -> (u64, f64)) -> (f64, f64, u64) {
    let mut error = 0.0;
    let mut times = Vec::new();
    let mut sum = 0.0;
    let mut round = 0u64;
    let mut found = 0;
    while (round < 10 || error / round as f64 > 0.05 * (sum / round as f64)) && round <= 50 {
        let (f, ns) = round_fn(round + 1);
        found = f;
        sum += ns;
        round += 1;
        error += ((sum / round as f64) - ns).abs();
        times.push(ns);
        println!(
            "round {round} {label}found {found} time per kmer: {ns}, avg: {}, error: {}",
            sum / round as f64,
            error / round as f64
        );
    }
    let mean = times.iter().sum::<f64>() / times.len() as f64;
    let avg = sum / round as f64;
    let sd = (times.iter().map(|x| (x - avg) * (x - avg)).sum::<f64>() / times.len() as f64).sqrt();
    (mean, sd, found)
}

fn bench_index<W: KmerWord>(idx: &RsHash<W>, n: usize) {
    println!("bench pos lookup...");
    let (mean, sd, found) = bench_rounds("", |r| {
        let kmers = idx.random_text_windows(n, r);
        let t0 = Instant::now();
        let f = idx.lookup_many(&kmers);
        (f, t0.elapsed().as_nanos() as f64 / kmers.len() as f64)
    });
    println!("==== positive lookup:");
    println!("num_kmers = {n}");
    println!("num_positive_kmers = {found} ({}%)", found as f64 / n as f64 * 100.0);
    println!("pos_time_per_kmer = {mean}");
    println!("pos_time_per_kmer_variance = {sd}");

    println!("bench neg lookup...");
    let (mean, sd, found) = bench_rounds("", |r| {
        let kmers = idx.random_windows(n, r);
        let t0 = Instant::now();
        let f = idx.lookup_many(&kmers);
        (f, t0.elapsed().as_nanos() as f64 / kmers.len() as f64)
    });
    println!("==== negative lookup:");
    println!("num_kmers = {n}");
    println!("num_negative_kmers = {found} ({}%)", found as f64 / n as f64 * 100.0);
    println!("neg_time_per_kmer = {mean}");
    println!("neg_time_per_kmer_variance = {sd}");

    if idx.has_locate() {
        println!("bench locate...");
        let mut hits = Vec::new();
        for (label, positive) in [("pos", true), ("neg", false)] {
            println!("bench {label} locate...");
            let (mean, sd, found) = bench_rounds("", |r| {
                let kmers = if positive { idx.random_text_windows(n, r) } else { idx.random_windows(n, r) };
                let t0 = Instant::now();
                let mut f = 0u64;
                for &v in &kmers {
                    hits.clear();
                    idx.locate(v, &mut hits);
                    f += hits.len() as u64;
                }
                (f, t0.elapsed().as_nanos() as f64 / kmers.len() as f64)
            });
            println!("==== locate:");
            println!("num_kmers = {n}");
            let what = if positive { "positive" } else { "negative" };
            println!("num_{what}_positions = {found} ({}%)", found as f64 / n as f64 * 100.0);
            println!("{label}_locate_time_per_kmer = {mean}");
            println!("{label}_locate_time_per_kmer_variance = {sd}");
        }
    }
}

fn bench(args: &Args) -> Result<()> {
    println!("loading dict...");
    let index = load(args)?;
    with_index!(&index, idx => bench_index(idx, args.bench_kmers));
    Ok(())
}
