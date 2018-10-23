#include <chrono>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include <sdsl/int_vector.hpp>

#include "common.h"
#include "algo1.hpp"
#include "algo2.hpp"
#include "algo3.hpp"
#include "algo4.hpp"

using namespace seqan;

template <typename TIV, typename TRng>
void randText(TIV & iv, TRng & rng, unsigned const len)
{
    resize(iv, len);
    for (uint64_t i = 0; i < len; ++i)
        iv[i] = Dna(rng() % 4); // 0 is reserved for sentinel character!
}

int main(int argc, char *argv[])
{
    auto now = std::chrono::system_clock::now();
    auto seed = 1532123185108119664;//std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    std::cout << "Seed: " << seed << '\n';
    std::mt19937_64 rng(seed);

    unsigned const threads = 1;//omp_get_num_threads();
    constexpr unsigned errors = 0;

    while (true)
    {
        typedef StringSet<String<Dna>, Owner<ConcatDirect<> > > TGenome;
        TGenome genome;

        // TODO: bug when k-mer * overlap >= (?) text length

        std::cout << "String lengths: ";
        unsigned stringsetsize = (rng() % 3) + 1;
        for (unsigned ss = 0; ss < stringsetsize; ++ss)
        {
            unsigned textLength = (rng() % 25) + 15;
            DnaString chr;
            randText(chr, rng, textLength);
            appendValue(genome, chr);
        //     std::cout << textLength << ' ';
        }
        // std::cout << '\n';
        // appendValue(genome, DnaString("ACAGTCAAAGCTCGAG"));
        // appendValue(genome, DnaString("TCACTAGTCAGTGCT"));
        // appendValue(genome, DnaString("GTCAGCTCGAG"));
        // appendValue(genome, DnaString("TCACTAGTCA"));

        Index<TGenome, TIndexConfig> index(genome);
        indexCreate(index, FibreSALF());
        auto & text = indexText(index).concat;

        for (unsigned j = 0; j < 40; ++j)
        {
            unsigned length = (rand() % 10) + 4;
            unsigned reads = std::max(1, static_cast<int>(length) - 1 - (rand() % 3));
            unsigned overlap = length - reads + 1;

            std::cout << "Length (K): " << length << ", Overlap: " << overlap << ", Reads: " << reads << '\n';

            std::vector<uint8_t> c1(seqan::length(text) - length + 1);
            std::vector<uint8_t> c2(seqan::length(text) - length + 1);
            std::vector<uint8_t> c3(seqan::length(text) - length + 1);
            std::vector<uint8_t> c4(seqan::length(text) - length + 1);

            // #pragma omp sections
            {
                { runAlgoTrivial<errors>(index, text, length, c1, 0, threads); }
                // #pragma omp section
                { runAlgo2<errors>(index, text, length, c2, overlap, threads); }
                // #pragma omp section
                { runAlgo3<errors>(index, text, length, c3, overlap, threads); }
                // #pragma omp section
                { runAlgo4<errors>(index, text, length, c4, overlap, threads); }
            }

            if (c1 != c2 || c1 != c3 || c1 != c4)
            {
                std::cout << "L : " << length << '\n';
                for (unsigned ss = 0; ss < stringsetsize; ++ss)
                {
                    std::cout << genome[ss] << '\n';
                }

                std::cout << "    ";
                for (unsigned i = 0; i < seqan::length(genome[0]); ++i)
                    std::cout << genome[0][i] << "  ";
                for (unsigned i = 0; i < seqan::length(genome[1]); ++i)
                    std::cout << genome[1][i] << "  ";
                std::cout << "\n    ";
                for (unsigned i = 0; i < seqan::length(text); ++i)
                    std::cout << text[i] << "  ";
                std::cout << "\n   ";
                for (unsigned i = 0; i < c1.size(); ++i)
                    std::cout << ((i < 10) ? " " : "") << i << " ";
                std::cout << "\nc1: ";
                for (unsigned i = 0; i < c1.size(); ++i)
                    std::cout << (unsigned)c1[i] << "  ";
                std::cout << "\nc2: ";
                for (unsigned i = 0; i < c2.size(); ++i)
                    std::cout << (unsigned)c2[i] << "  ";
                std::cout << "\nc3: ";
                for (unsigned i = 0; i < c3.size(); ++i)
                    std::cout << (unsigned)c3[i] << "  ";
                std::cout << "\nc4: ";
                for (unsigned i = 0; i < c4.size(); ++i)
                    std::cout << (unsigned)c4[i] << "  ";
                std::cout << '\n';
                exit(1);
            }
        }
    }
}
