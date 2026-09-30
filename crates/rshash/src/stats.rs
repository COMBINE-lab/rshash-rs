//! Space report (port of `RSHash::print_info`; same labels, so the upstream
//! benchmark scripts can grep them).

use crate::index::RsHash;
use crate::last_level::Keys;
use crate::word::KmerWord;
use std::fmt::Write;

impl<W: KmerWord> RsHash<W> {
    /// The `====== report ======` block printed after `build`.
    pub fn report(&self) -> String {
        let mut o = String::new();
        let n = self.text.len();
        let tk = self.stats.text_kmers.max(1) as f64;
        let per = |bits: usize| bits as f64 / tk;

        let mut no_min = [0usize; 3];
        let mut no_occ = [0usize; 3];
        let mut r_bits = [0usize; 3];
        let mut s_bits = [0usize; 3];
        let mut off_bits = [0usize; 3];
        let mut universe = [1f64; 3];
        for (l, lv) in self.levels.iter().enumerate() {
            no_min[l] = lv.num_minimizers();
            no_occ[l] = lv.num_occurrences();
            r_bits[l] = lv.r.bit_size();
            s_bits[l] = lv.s.bit_size();
            off_bits[l] = lv.offsets.bit_size();
            universe[l] = lv.mp.universe() as f64;
        }
        let last_keys: usize = self.last.iter().map(|t| t.len()).sum();
        let last_bits: usize = self.last.iter().map(|t| t.bit_size()).sum();
        let last_pos: usize = self.last.iter().map(|t| t.num_positions()).sum();
        let last_pos_bits: usize = self.last.iter().map(|t| t.positions.as_ref().map_or(0, |p| p.bit_size())).sum();
        let r4_bits: usize = self.last.iter().map(|t| if let Keys::Ef(ef) = &t.keys { ef.bit_size() } else { 0 }).sum();
        let s4_bits: usize = self.last.iter().map(|t| t.buckets.as_ref().map_or(0, |b| b.bit_size())).sum();

        let _ = writeln!(o, "====== report ======");
        let _ = writeln!(o, "text length: {n}");
        let _ = writeln!(o, "textkmers: {}", self.stats.text_kmers);
        for l in 0..3 {
            let i = l + 1;
            let _ = writeln!(o, "no minimiser{i}: {}", no_min[l]);
            let _ = writeln!(o, "no distinct minimiser{i}: {}", no_occ[l]);
            let avg = if no_min[l] > 0 { no_occ[l] as f64 / no_min[l] as f64 } else { 0.0 };
            let _ = writeln!(o, "avg superkmers{i}: {avg}");
        }
        // as in C++: distinct keys stored in the last level
        let _ = writeln!(o, "no freq kmers: {} {}%", last_keys, last_keys as f64 / tk * 100.0);
        let _ = writeln!(o, "last level kmer occurrences: {}", self.stats.last_level_kmers);
        for l in 0..3 {
            let _ = writeln!(o, "density r{}: {}%", l + 1, no_min[l] as f64 / universe[l] * 100.0);
        }
        for l in 0..3 {
            let _ = writeln!(o, "density s{}: {}%", l + 1, no_min[l] as f64 / (no_occ[l] + 1) as f64 * 100.0);
        }
        let _ = writeln!(o, "\nspace per kmer in bit:");
        let _ = writeln!(o, "text: {}", per(2 * n as usize));
        let _ = writeln!(o, "endpoints: {}", per(self.text.bounds_bit_size()));
        for l in 0..3 {
            let _ = writeln!(o, "offsets{}: {}", l + 1, per(off_bits[l]));
        }
        let _ = writeln!(o, "offsets4: {}", per(last_pos_bits));
        let _ = writeln!(o, "Last_level: {}", per(last_bits));
        for l in 0..3 {
            let _ = writeln!(o, "R_{}: {}", l + 1, per(r_bits[l]));
        }
        let _ = writeln!(o, "R_4: {}", per(r4_bits));
        for l in 0..3 {
            let _ = writeln!(o, "S_{}: {}", l + 1, per(s_bits[l]));
        }
        let _ = writeln!(o, "S_4: {}", per(s4_bits));
        let _ = last_pos;
        let total = 2 * n as usize
            + self.text.bounds_bit_size()
            + off_bits.iter().sum::<usize>()
            + r_bits.iter().sum::<usize>()
            + s_bits.iter().sum::<usize>()
            + last_bits;
        let _ = write!(o, "total: {}", per(total));
        o
    }
}
