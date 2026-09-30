// Component micro-benchmark of the C++ streaming lookup, mirroring
// crates/rshash/examples/components.rs. Accesses RSHash internals.
#include <chrono>
#include <iostream>
#define private public
#include "rshash.hpp"
#undef private

struct traits : seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};

int main(int argc, char **argv)
{
    RSHash index;
    index.load(argv[1]);
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> qs;
    seqan3::sequence_file_input<traits> fin{argv[2]};
    for (auto &r : fin) {
        seqan3::bitpacked_sequence<seqan3::dna4> q;
        q.assign(r.sequence().begin(), r.sequence().end());
        qs.push_back(std::move(q));
    }
    using clk = std::chrono::high_resolution_clock;
    if (argc > 3) { // streaming only, repeated
        uint64_t ext = 0, s = 0, n = 0;
        for (auto &q : qs) n += q.size() - index.window_size + 1;
        auto t = clk::now();
        for (int rep = 0; rep < std::atoi(argv[3]); rep++)
            for (auto &q : qs) s += index.streaming_lookup(q, ext);
        std::cout << "streaming only: " << std::chrono::duration<double, std::nano>(clk::now() - t).count() / n / std::atoi(argv[3]) << " ns/kmer (" << s << ")\n";
        return 0;
    }
    uint64_t n = 0, acc = 0;
    auto t = clk::now();
    for (auto &q : qs)
        for (auto &&w : q | rshash::views::kmerview({.window_size = index.window_size})) { acc ^= w.value ^ w.value_rev; n++; }
    double base = std::chrono::duration<double, std::nano>(clk::now() - t).count() / n;
    std::cout << "windows only:        " << base << " ns/kmer (" << acc << ")\n";

    t = clk::now();
    uint64_t s = 0;
    for (auto &q : qs) {
        bool rolling = false; uint64_t m; size_t l, r;
        for (auto &&w : q | rshash::views::kmerview({.window_size = index.window_size})) {
            if (rolling) index.update_minimiser<1>(w.value, w.value_rev, m, l, r);
            else { m = index.find_minimiser<1>(w.value, w.value_rev, l, r); rolling = true; }
            s += m;
        }
    }
    std::cout << "+ rolling minimiser: " << std::chrono::duration<double, std::nano>(clk::now() - t).count() / n << " ns/kmer (" << s << ")\n";

    t = clk::now();
    s = 0;
    for (auto &q : qs) {
        bool rolling = false; uint64_t m, last = ~0ULL, rank; size_t l, r;
        for (auto &&w : q | rshash::views::kmerview({.window_size = index.window_size})) {
            if (rolling) index.update_minimiser<1>(w.value, w.value_rev, m, l, r);
            else { m = index.find_minimiser<1>(w.value, w.value_rev, l, r); rolling = true; }
            if (m != last) { s += index.r1.contains(m, rank); last = m; }
        }
    }
    std::cout << "+ R1.contains:       " << std::chrono::duration<double, std::nano>(clk::now() - t).count() / n << " ns/kmer (" << s << ")\n";

    t = clk::now();
    s = 0;
    for (auto &q : qs) {
        for (auto &&w : q | rshash::views::kmerview({.window_size = index.window_size})) {
            size_t l, r;
            s += index.find_minimiser<1>(w.value, w.value_rev, l, r);
        }
    }
    std::cout << "find every window:   " << std::chrono::duration<double, std::nano>(clk::now() - t).count() / n << " ns/kmer (" << s << ")\n";

    uint64_t ext = 0;
    t = clk::now();
    s = 0;
    for (auto &q : qs) s += index.streaming_lookup(q, ext);
    std::cout << "streaming lookup:    " << std::chrono::duration<double, std::nano>(clk::now() - t).count() / n << " ns/kmer (" << s << ")\n";
}
