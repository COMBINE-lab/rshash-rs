//! The last level: k-mers (or gapped k-mers) whose minimisers were too
//! frequent at every minimiser level (port of `last_level`,
//! `process_freq_kmers`, `count_kmers`, `build_level<4>` and `FlatMap`).
//!
//! One table per shape (one table without shapes). Keys:
//! * contiguous k-mers: canonical `min(fwd, rc)`;
//! * shapes: forward keys `pext(window, w_mask)` (D8), probed with both the
//!   query and its reverse complement.
//!
//! Storage: a hash set (`use_ht && !loc`), a hash map key → bucket
//! (`use_ht && loc`), or an Elias–Fano set whose rank is the bucket
//! (`!use_ht`). With `loc`, buckets delimit the text positions of each key.

use crate::bits::{BitVec, CompactVec, RankSelect};
use crate::ef::EliasFano;
use crate::io::{Reader, Writer, invalid};
use crate::word::KmerWord;
use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};
use std::io;

/// Fast hasher for packed k-mer keys (murmur3 `fmix64` finaliser, which
/// mixes all input bits into both the low bits used for the bucket index and
/// the high bits used for the control byte).
#[derive(Default, Clone, Copy)]
pub struct WordHasher {
    state: u64,
}

const HASH_K: u64 = 0x9E37_79B9_7F4A_7C15;

impl Hasher for WordHasher {
    #[inline(always)]
    fn finish(&self) -> u64 {
        self.state
    }
    #[inline(always)]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.write_u64(b as u64);
        }
    }
    #[inline(always)]
    fn write_u64(&mut self, x: u64) {
        let mut h = self.state.rotate_left(5) ^ x;
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        h ^= h >> 33;
        self.state = h.wrapping_add(HASH_K);
    }
    #[inline(always)]
    fn write_u128(&mut self, x: u128) {
        self.write_u64(x as u64);
        self.write_u64((x >> 64) as u64);
    }
}

pub type WordBuildHasher = BuildHasherDefault<WordHasher>;

#[derive(Clone, Debug)]
pub enum Keys<W: KmerWord> {
    Set(HashSet<W, WordBuildHasher>),
    /// key → bucket index
    Map(HashMap<W, u32, WordBuildHasher>),
    Ef(EliasFano<W>),
}

/// One last-level table.
#[derive(Clone, Debug)]
pub struct LastTable<W: KmerWord> {
    pub keys: Keys<W>,
    /// Bucket delimiters (one bit per position + 1, ones at bucket starts).
    pub buckets: Option<RankSelect>,
    /// Text positions (window starts), grouped by bucket.
    pub positions: Option<CompactVec>,
    /// Width in bits of a key (for the Elias–Fano universe).
    pub key_bits: u32,
    /// Derived (not serialised) one-hash bit filter over the keys (8–16 bits
    /// per key): most probes of absent keys, the common case in streaming
    /// lookups, stop here without touching the table.
    filter: KeyFilter,
}

/// One-hash bit filter (a Bloom filter with a single hash function).
#[derive(Clone, Debug, Default)]
struct KeyFilter {
    bits: Vec<u64>,
    mask: u64,
}

#[inline(always)]
fn fmix64(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    h ^ (h >> 33)
}

impl KeyFilter {
    fn new<W: KmerWord>(keys: impl Iterator<Item = W>, n: usize) -> Self {
        if n == 0 {
            return Self::default();
        }
        let nbits = (8 * n).next_power_of_two().max(64);
        let mut f = Self { bits: vec![0; nbits / 64], mask: nbits as u64 - 1 };
        for k in keys {
            let h = Self::hash(k) & f.mask;
            f.bits[(h / 64) as usize] |= 1 << (h % 64);
        }
        f
    }

    #[inline(always)]
    fn hash<W: KmerWord>(k: W) -> u64 {
        let x = k.to_u128();
        fmix64((x as u64) ^ ((x >> 64) as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
    }

    #[inline(always)]
    fn may_contain<W: KmerWord>(&self, k: W) -> bool {
        if self.bits.is_empty() {
            return false;
        }
        let h = Self::hash(k) & self.mask;
        // SAFETY: h <= mask < 64 * bits.len()
        (unsafe { *self.bits.get_unchecked((h / 64) as usize) } >> (h % 64)) & 1 == 1
    }
}

impl<W: KmerWord> LastTable<W> {
    /// Build from `(key, position)` pairs. Keys occurring `>= threshold` times
    /// are dropped (unless `threshold == 0`).
    pub fn build(
        mut pairs: Vec<(W, u64)>,
        key_bits: u32,
        use_ht: bool,
        loc: bool,
        threshold: u16,
        pos_width: u32,
    ) -> Self {
        pairs.sort_unstable();
        pairs.dedup();
        let mut keys: Vec<W> = Vec::new();
        let mut counts: Vec<u32> = Vec::new();
        let mut positions: Vec<u64> = Vec::new();
        let mut i = 0;
        while i < pairs.len() {
            let key = pairs[i].0;
            let mut j = i;
            while j < pairs.len() && pairs[j].0 == key {
                j += 1;
            }
            let occ = j - i;
            if threshold == 0 || occ < threshold as usize {
                keys.push(key);
                counts.push(occ as u32);
                if loc {
                    positions.extend(pairs[i..j].iter().map(|p| p.1));
                }
            }
            i = j;
        }
        drop(pairs);

        let keys_struct = if use_ht {
            if loc {
                let mut map = HashMap::with_capacity_and_hasher(keys.len(), WordBuildHasher::default());
                for (b, &k) in keys.iter().enumerate() {
                    map.insert(k, b as u32);
                }
                Keys::Map(map)
            } else {
                let mut set = HashSet::with_capacity_and_hasher(keys.len(), WordBuildHasher::default());
                set.extend(keys.iter().copied());
                Keys::Set(set)
            }
        } else {
            Keys::Ef(EliasFano::new(&keys, 1u128.checked_shl(key_bits).unwrap_or(u128::MAX)))
        };

        let (buckets, positions) = if loc {
            let mut bv = BitVec::new(positions.len() + 1);
            bv.set(0);
            let mut j = 0usize;
            for &c in &counts {
                j += c as usize;
                bv.set(j);
            }
            (Some(RankSelect::with_samples(bv, true, false)), Some(CompactVec::from_slice(&positions, pos_width)))
        } else {
            (None, None)
        };

        let filter = KeyFilter::new(keys.iter().copied(), keys.len());
        Self { keys: keys_struct, buckets, positions, key_bits, filter }
    }

    #[inline]
    pub fn contains(&self, key: W) -> bool {
        if !self.filter.may_contain(key) {
            return false;
        }
        match &self.keys {
            Keys::Set(s) => s.contains(&key),
            Keys::Map(m) => m.contains_key(&key),
            Keys::Ef(ef) => ef.contains(key).is_some(),
        }
    }

    /// Text positions (window starts) stored for `key` (requires `loc`).
    pub fn positions(&self, key: W, out: &mut Vec<u64>) {
        let (Some(buckets), Some(pos)) = (&self.buckets, &self.positions) else { return };
        if !self.filter.may_contain(key) {
            return;
        }
        let b = match &self.keys {
            Keys::Set(_) => return,
            Keys::Map(m) => match m.get(&key) {
                Some(&b) => b as usize,
                None => return,
            },
            Keys::Ef(ef) => match ef.contains(key) {
                Some(b) => b,
                None => return,
            },
        };
        let (s, e) = buckets.select1_pair(b);
        out.extend((s..e).map(|i| pos.get(i)));
    }

    pub fn len(&self) -> usize {
        match &self.keys {
            Keys::Set(s) => s.len(),
            Keys::Map(m) => m.len(),
            Keys::Ef(ef) => ef.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn num_positions(&self) -> usize {
        self.positions.as_ref().map_or(0, |p| p.len())
    }

    /// Approximate size in bits.
    pub fn bit_size(&self) -> usize {
        let keys = match &self.keys {
            // one control byte + key per slot (as `print_info` does for gtl)
            Keys::Set(s) => s.capacity() * (W::BITS as usize + 8),
            Keys::Map(m) => m.capacity() * (W::BITS as usize + 32 + 8),
            Keys::Ef(ef) => ef.bit_size(),
        };
        keys + 64 * self.filter.bits.len()
            + self.buckets.as_ref().map_or(0, |b| b.bit_size())
            + self.positions.as_ref().map_or(0, |p| p.bit_size())
    }

    fn sorted_keys(&self) -> Vec<W> {
        let mut keys: Vec<W> = match &self.keys {
            Keys::Set(s) => s.iter().copied().collect(),
            Keys::Map(m) => {
                let mut v: Vec<(u32, W)> = m.iter().map(|(&k, &b)| (b, k)).collect();
                v.sort_unstable();
                return v.into_iter().map(|x| x.1).collect();
            }
            Keys::Ef(ef) => return ef.iter().collect(),
        };
        keys.sort_unstable();
        keys
    }

    pub fn write(&self, w: &mut Writer) -> io::Result<()> {
        w.u32(self.key_bits)?;
        match &self.keys {
            Keys::Set(_) | Keys::Map(_) => {
                w.u8(if matches!(self.keys, Keys::Set(_)) { 0 } else { 1 })?;
                // keys in sorted (= bucket) order: deterministic files
                let keys = self.sorted_keys();
                w.u64(keys.len() as u64)?;
                for k in keys {
                    write_word(w, k)?;
                }
            }
            Keys::Ef(ef) => {
                w.u8(2)?;
                ef.write(w)?;
            }
        }
        w.bool(self.buckets.is_some())?;
        if let (Some(b), Some(p)) = (&self.buckets, &self.positions) {
            b.write(w)?;
            p.write(w)?;
        }
        Ok(())
    }

    pub fn read(r: &mut Reader) -> io::Result<Self> {
        let key_bits = r.u32()?;
        let keys = match r.u8()? {
            tag @ (0 | 1) => {
                let n = r.u64()? as usize;
                let mut keys = Vec::with_capacity(n);
                for _ in 0..n {
                    keys.push(read_word::<W>(r)?);
                }
                if tag == 0 {
                    let mut set = HashSet::with_capacity_and_hasher(n, WordBuildHasher::default());
                    set.extend(keys);
                    Keys::Set(set)
                } else {
                    let mut map = HashMap::with_capacity_and_hasher(n, WordBuildHasher::default());
                    for (b, k) in keys.into_iter().enumerate() {
                        map.insert(k, b as u32);
                    }
                    Keys::Map(map)
                }
            }
            2 => Keys::Ef(EliasFano::read(r)?),
            _ => return Err(invalid("bad last-level tag")),
        };
        let (buckets, positions) = if r.bool()? {
            (Some(RankSelect::read_with_samples(r, true, false)?), Some(CompactVec::read(r)?))
        } else {
            (None, None)
        };
        let mut t = Self { keys, buckets, positions, key_bits, filter: KeyFilter::default() };
        let n = t.len();
        t.filter = KeyFilter::new(t.sorted_keys().into_iter(), n);
        Ok(t)
    }
}

pub fn write_word<W: KmerWord>(w: &mut Writer, x: W) -> io::Result<()> {
    if W::BITS == 64 { w.u64(x.low_u64()) } else { w.u128(x.to_u128()) }
}

pub fn read_word<W: KmerWord>(r: &mut Reader) -> io::Result<W> {
    Ok(if W::BITS == 64 { W::from_u64(r.u64()?) } else { W::from_u128_trunc(r.u128()?) })
}
