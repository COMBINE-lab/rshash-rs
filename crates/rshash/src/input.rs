//! FASTA/FASTQ input (replaces seqan3 `sequence_file_input<dna4>`): plain or
//! compressed, multi-line records joined, one text sequence per record,
//! non-ACGT characters converted to `A`.

use crate::text::{Text, TextBuilder};
use std::path::Path;

pub type Error = needletail::errors::ParseError;

/// Call `f(sequence)` for every record of a FASTA/FASTQ file.
pub fn for_each_record(path: impl AsRef<Path>, mut f: impl FnMut(&[u8])) -> Result<(), Error> {
    let mut reader = needletail::parse_fastx_file(path)?;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        let seq = rec.seq();
        f(&seq);
    }
    Ok(())
}

/// Read all records into a [`Text`]; returns the text and the number of
/// non-ACGT characters that were converted to `A`.
pub fn read_text(path: impl AsRef<Path>) -> Result<(Text, u64), Error> {
    let mut b = TextBuilder::new();
    for_each_record(path, |s| b.push_ascii(s))?;
    let n = b.non_acgt();
    Ok((b.finish(), n))
}

/// Read all records into memory.
pub fn read_sequences(path: impl AsRef<Path>) -> Result<Vec<Vec<u8>>, Error> {
    let mut v = Vec::new();
    for_each_record(path, |s| v.push(s.to_vec()))?;
    Ok(v)
}
