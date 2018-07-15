#ifndef COMMON_AUXILLARY_H_
#define COMMON_AUXILLARY_H_

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace seqan;

struct majorCaseParameters{
    bool nomappability = true;
    bool directsearch = true;
    bool compmappable = true;
    bool suspectunidirectional = true;
    
    //binaryNumber
    int stepcheck = 4;
    int distancetoblockend = 2;
    
    int directsearch_th = 2;
    float filter_th = 0.5;
    
    float flipdensity = 0.5;
    
    int intervalsize = 3;
    
    void print(){
        cout << "Cases Enabled: " << "\n";
        cout << nomappability << " " << directsearch << " " << compmappable << " " << suspectunidirectional << "\n";
        cout << "Params: " << "\n";
        cout << "stepcheck: " << stepcheck << "\n";
        cout << "distancetoblockend: " << distancetoblockend << "\n";
        cout << "directsearch_th: " << directsearch_th << "\n";
        cout << "filter_th: " << filter_th << "\n";
        cout << "flipdensity: " << flipdensity << "\n";
        cout << "intervalsize: " << intervalsize << "\n";
    }
};


struct myGlobalParameters{
public:
    bool startUnidirectional = false;
    majorCaseParameters normal;
    majorCaseParameters uni;
    
    
    void print(){
        normal.print();
        uni.print();
    }
};

extern int global;
extern myGlobalParameters params;



/*
struct trackCases{
public:
    int returncodes[10];
};*/

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
    
void testglobal();
    
enum class ReturnCode {
	NOMAPPABILITY, DIRECTSEARCH, COMPMAPPABLE, ONEDIRECTION, MAPPABLE, FINISHED, UNIDIRECTIONAL, SUSPECTUNIDIRECTIONAL, FILTER, ERROR
};

enum class BV {
	RIGHT = 0, MIDDLE = 1, LEFT = 2
};


typedef String<Dna, Alloc<>> TString;
typedef StringSet<TString, Owner<ConcatDirect<> > > TText;
typedef Index<TText, TIndexConfig> MyIndex;
    
 /*   
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
                  TDir const & );*/
}

#endif
