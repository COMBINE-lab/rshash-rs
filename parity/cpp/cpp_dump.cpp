// Parity driver: loads a C++ RSHash index and prints, for every query record,
// "record<TAB>length<TAB>windows<TAB>hits" using streaming_lookup. With
// "--magic" it prints the mixer_64 magics of the three level seeds.
#include <iostream>
#include "rshash.hpp"

struct dump_traits : seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};

int main(int argc, char **argv)
{
    if (argc == 2 && std::string(argv[1]) == "--magic") {
        for (uint64_t s : {seed1, seed2, seed3}) {
            mixer_64 m(s);
            std::printf("%016llx %016llx\n", (unsigned long long) s, (unsigned long long) (m.hash(0)));
        }
        return 0;
    }
    if (argc != 3) {
        std::cerr << "usage: rshash_dump <index> <queries> | --magic\n";
        return 1;
    }
    RSHash index;
    index.load(argv[1]);
    const uint64_t L = index.getshapes().shapes[0].value != std::numeric_limits<uint32_t>::max()
                           ? index.getshapes().length : index.getk();
    seqan3::sequence_file_input<dump_traits> fin{argv[2]};
    uint64_t i = 0, total = 0, ext = 0;
    for (auto &rec : fin) {
        seqan3::bitpacked_sequence<seqan3::dna4> q;
        q.assign(rec.sequence().begin(), rec.sequence().end());
        uint64_t windows = q.size() >= L ? q.size() - L + 1 : 0;
        uint64_t hits = windows ? index.streaming_lookup(q, ext) : 0;
        total += hits;
        std::cout << i++ << '\t' << q.size() << '\t' << windows << '\t' << hits << '\n';
    }
    std::cerr << "total hits " << total << " extensions " << ext << '\n';
    return 0;
}
