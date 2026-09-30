//! Rust port of **RSHash** (Schulte-Mattler & Reinert), a compressed exact
//! k-mer dictionary that replaces SSHash's minimal perfect hash functions by
//! Elias–Fano bitvectors with rank/select support, with a layered
//! multi-minimiser scheme for skewed minimiser distributions — including the
//! gapped k-mer (shape) support of the upstream `multishapes` branch.
//!
//! ```
//! use rshash::{AnyIndex, Params, text::Text};
//! let text = Text::from_sequences(&["ACGTACGTTTGACCAGGATTACAGATTACCA"]);
//! let params = Params { k: 15, levels: 1, m: [7, 9, 11], ..Params::default() };
//! let index = AnyIndex::build(text, &params).unwrap();
//! assert_eq!(index.streaming_lookup_count(b"TTGACCAGGATTACA").0, 1);
//! ```

#![allow(clippy::type_complexity, clippy::large_enum_variant, clippy::needless_range_loop)]

pub mod alphabet;
pub mod bits;
pub mod ef;
pub mod hash;
pub mod index;
pub mod input;
pub mod io;
pub mod last_level;
pub mod minimizer;
pub mod params;
pub mod query;
pub mod shape;
pub mod stats;
pub mod text;
pub mod util;
pub mod word;

pub use index::{BuildError, RsHash};
pub use params::Params;
pub use query::{Hit, StreamingLocate, StreamingLookup};
pub use shape::Shapes;
pub use word::KmerWord;

use std::io::{Read, Write};
use text::Text;

/// An index over `u64` windows (up to 32 bases) or `u128` windows (up to 63).
#[derive(Clone, Debug)]
pub enum AnyIndex {
    W64(RsHash<u64>),
    W128(RsHash<u128>),
}

/// Dispatch a generic expression over both index variants.
#[macro_export]
macro_rules! with_index {
    ($any:expr, $idx:ident => $body:expr) => {
        match $any {
            $crate::AnyIndex::W64($idx) => $body,
            $crate::AnyIndex::W128($idx) => $body,
        }
    };
}

impl AnyIndex {
    /// Build, choosing the window type from the window length.
    pub fn build(text: Text, params: &Params) -> Result<Self, BuildError> {
        Self::build_with_log(text, params, &mut |_| {})
    }

    pub fn build_with_log(text: Text, params: &Params, log: &mut dyn FnMut(&str)) -> Result<Self, BuildError> {
        let g = params.geometry()?;
        Ok(if g.window <= 32 {
            AnyIndex::W64(RsHash::build_with_log(text, params, log)?)
        } else {
            AnyIndex::W128(RsHash::build_with_log(text, params, log)?)
        })
    }

    pub fn window(&self) -> u32 {
        with_index!(self, i => i.window())
    }

    pub fn kernel(&self) -> u32 {
        with_index!(self, i => i.kernel())
    }

    pub fn params(&self) -> &Params {
        with_index!(self, i => i.params())
    }

    pub fn has_locate(&self) -> bool {
        with_index!(self, i => i.has_locate())
    }

    pub fn report(&self) -> String {
        with_index!(self, i => i.report())
    }

    pub fn set_non_acgt(&mut self, n: u64) {
        with_index!(self, i => i.set_non_acgt(n))
    }

    /// Streaming lookup of one sequence: `(hits, extensions)`.
    /// (For many queries, keep a [`StreamingLookup`] engine instead.)
    pub fn streaming_lookup_count(&self, seq: &[u8]) -> (u64, u64) {
        with_index!(self, i => {
            let mut e = i.streaming_lookup();
            let h = e.query(seq);
            (h, e.extensions)
        })
    }

    /// Point lookup of an ASCII window of exactly `window()` bases.
    pub fn lookup_ascii(&self, kmer: &[u8]) -> bool {
        assert_eq!(kmer.len(), self.window() as usize);
        with_index!(self, i => i.lookup(word::encode_kmer(kmer)))
    }

    /// Point locate of an ASCII window of exactly `window()` bases.
    pub fn locate_ascii(&self, kmer: &[u8]) -> Vec<Hit> {
        assert_eq!(kmer.len(), self.window() as usize);
        let mut out = Vec::new();
        with_index!(self, i => i.locate(word::encode_kmer(kmer), &mut out));
        out
    }

    pub fn save(&self, out: &mut dyn Write) -> std::io::Result<()> {
        with_index!(self, i => i.write(out))
    }

    pub fn load(input: &mut dyn Read) -> std::io::Result<Self> {
        let bits = RsHash::<u64>::read_header(input)?;
        match bits {
            64 => Ok(AnyIndex::W64(RsHash::read_body(input)?)),
            128 => Ok(AnyIndex::W128(RsHash::read_body(input)?)),
            _ => Err(io::invalid("bad word size")),
        }
    }

    pub fn save_file(&self, path: impl AsRef<std::path::Path>) -> std::io::Result<()> {
        let f = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(f);
        self.save(&mut w)?;
        w.flush()
    }

    pub fn load_file(path: impl AsRef<std::path::Path>) -> std::io::Result<Self> {
        let f = std::fs::File::open(path)?;
        let mut r = std::io::BufReader::new(f);
        Self::load(&mut r)
    }
}
