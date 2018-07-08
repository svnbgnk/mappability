#include <chrono>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include <sdsl/int_vector.hpp>

#include "common.h"
#include "algo1.hpp"
#include "algo1_approx.hpp"
#include "algo2.hpp"

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
    auto seed = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    std::cout << "Seed: " << seed << '\n';
    std::mt19937_64 rng(seed);

    unsigned const threads = 1;//omp_get_num_threads();
    constexpr unsigned errors = 0;
    unsigned threshold = 999;

    while (true)
    {
        unsigned textLength = 10000;

        DnaString genome;
        randText(genome, rng, textLength);
        Index<DnaString, TIndexConfig> index(genome);
        indexCreate(index, FibreSALF());
        auto & text = indexText(index);

        for (unsigned j = 2; j < 3; ++j)
        {
            unsigned length = j;

            std::vector<uint8_t> c1(seqan::length(text) - length + 1);
            runAlgoTrivial<errors>(index, genome, length, c1, 0, threads);

            std::vector<uint8_t> c2(seqan::length(text) - length + 1);
            runAlgo1_approx<errors>(index, genome, length, c2, threads, threshold);

            if (c1 != c2)
            {
                std::cout << "L : " << length << '\n';
                std::cout << "T : " << genome << '\n';

                std::cout << "c1: ";
                for (unsigned i = 0; i < c1.size(); ++i)
                    std::cout << (unsigned)c1[i] << ' ';
                std::cout << "\nc2: ";
                for (unsigned i = 0; i < c2.size(); ++i)
                    std::cout << (unsigned)c2[i] << ' ';
                std::cout << '\n';
                exit(1);
            }
            std::cout << "." << std::flush;

            // for (unsigned overlap = 0; overlap <= length - errors - 2; ++overlap) // because there have to be enough characters for the infix using search schemes
            // {
            //     sdsl::int_vector<16> c2(seqan::length(text) - length + 1);
            //     runAlgo2<errors>(index, genome, length, c2, length - overlap, threads);
            //
            //     for (unsigned i = 0; i < c1.size(); ++i)
            //     {
            //         if (c1[i] != c2[i])
            //         {
            //             std::cout << "\nOverlap: " << overlap << ", Length: " << length << '\n';
            //             for (unsigned j = 0; j < seqan::length(genome); ++j)
            //                 std::cout << genome[j] << ' ';
            //             std::cout << '\n';
            //             for (unsigned j = 0; j < c1.size(); ++j)
            //                 std::cout << c1[j] << ' ';
            //             std::cout << '\n';
            //             for (unsigned j = 0; j < c2.size(); ++j)
            //                 std::cout << c2[j] << ' ';
            //             std::cout << '\n';
            //             exit(1);
            //         }
            //     }
            //     std::cout << "." << std::flush;
            // }
        }
    }
}
