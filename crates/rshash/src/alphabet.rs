//! seqan3 `dna4` alphabet: A=0, C=1, G=2, T=3 (case-insensitive).
//!
//! As in seqan3's char conversion, `U`/`u` maps to T and every other
//! character (N, IUPAC codes, ...) maps to A (rank 0).

const fn build_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
}

static RANK: [u8; 256] = build_table();

/// Rank of an ASCII nucleotide.
#[inline(always)]
pub fn rank(c: u8) -> u8 {
    RANK[c as usize]
}

/// Whether `c` is one of `ACGTUacgtu` (i.e. converts without loss).
#[inline(always)]
pub fn is_acgt(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u')
}

/// Convert an ASCII sequence to ranks.
pub fn to_ranks(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&c| rank(c)).collect()
}

/// Convert ranks back to ASCII.
pub fn to_ascii(ranks: &[u8]) -> Vec<u8> {
    ranks.iter().map(|&r| b"ACGT"[(r & 3) as usize]).collect()
}
