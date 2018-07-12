#ifndef COMMON_AUXILLARY_H_
#define COMMON_AUXILLARY_H_

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace seqan;

sdsl::bit_vector create_random_bit_v(int length);

template <typename T> 
void printv(T a);

inline bool file_exists (const std::string& name);

template <size_t nbrBlocks, size_t N>
void print_search_scheme(std::array<OptimalSearch<nbrBlocks>, N> & searchsscheme);

template <typename TText, typename TIndex, typename TIndexSpec>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
              bool const fwd);

template <typename TText, typename TIndex, typename TIndexSpec>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              int const number_of_indeces,
              bool const fwd);

template <typename TText, typename TIndex, typename TIndexSpec>
void print_fullsa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
              bool const fwd);

void printbit(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, Pair<uint8_t, Pair<uint32_t, uint32_t>> brange);

void printPair(pair<uint32_t, uint32_t> p);

namespace seqan{
template <typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle, typename TIndex,
          size_t nbrBlocks,
          typename TDir>
void directSearchDummy(TDelegateD & delegateDirect,
                  Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > olditer,
                  TNeedle const & needle,
                  vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, 
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  Pair<uint8_t, Pair<uint32_t, uint32_t>> brange,
                  TDir const & /**/);
}

#endif
