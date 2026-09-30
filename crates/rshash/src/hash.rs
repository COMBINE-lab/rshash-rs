//! Minimiser hashing (port of `MurmurHash2_64` and `mixer_64` from
//! `minimiser_views.hpp`).

/// Per-level seeds (`rshash.hpp:16-18`).
pub const SEEDS: [u64; 3] = [1, 0x296D_BD33_3256_8C64, 0xE59A_385F_0376_C9F6];

/// MurmurHash64A (Austin Appleby), as used by `mixer_64::seed`.
pub fn murmur_hash2_64(key: &[u8], seed: u64) -> u64 {
    const M: u64 = 0xc6a4_a793_5bd1_e995;
    const R: u32 = 47;
    let len = key.len();
    let mut h: u64 = seed ^ (len as u64).wrapping_mul(M);

    let nblocks = len / 8;
    for i in 0..nblocks {
        let mut k = u64::from_le_bytes(key[8 * i..8 * i + 8].try_into().unwrap());
        k = k.wrapping_mul(M);
        k ^= k >> R;
        k = k.wrapping_mul(M);
        h ^= k;
        h = h.wrapping_mul(M);
    }

    let tail = &key[8 * nblocks..];
    let rem = len & 7;
    if rem > 0 {
        for i in (0..rem).rev() {
            h ^= (tail[i] as u64) << (8 * i);
        }
        h = h.wrapping_mul(M);
    }

    h ^= h >> R;
    h = h.wrapping_mul(M);
    h ^= h >> R;
    h
}

/// `mixer_64`: `hash(x) = (x * 0x517cc1b727220a95) ^ magic`,
/// with `magic = MurmurHash2_64(&seed, 8, 0)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Mixer64 {
    magic: u64,
}

impl Mixer64 {
    pub fn new(seed: u64) -> Self {
        Self { magic: murmur_hash2_64(&seed.to_le_bytes(), 0) }
    }

    #[inline(always)]
    pub fn hash(&self, x: u64) -> u64 {
        x.wrapping_mul(0x517c_c1b7_2722_0a95) ^ self.magic
    }

    pub fn magic(&self) -> u64 {
        self.magic
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn magics() {
        // Golden values; cross-checked against the C++ build in parity/.
        assert_eq!(Mixer64::new(SEEDS[0]).magic(), 0x8fbb_8d81_5c9e_092e);
        assert_eq!(Mixer64::new(SEEDS[1]).magic(), 0x7171_fede_0774_eda2);
        assert_eq!(Mixer64::new(SEEDS[2]).magic(), 0x0a53_d287_bb29_1669);
    }
}
