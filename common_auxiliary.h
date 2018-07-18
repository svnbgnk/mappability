#ifndef COMMON_AUXILLARY_H_
#define COMMON_AUXILLARY_H_

#include <sdsl/bit_vectors.hpp>

namespace seqan{

struct hit{
    bool rev;
    Pair <unsigned, unsigned> occ;
    uint8_t errors;
    DnaString read;
};

template <size_t nbrBlocks, size_t N>
constexpr inline void _optimalSearchSchemeSetMapParams(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{
    for (OptimalSearch<nbrBlocks> & s : ss){
        int bsize = s.pi.size();
        uint8_t min = s.pi[0];
        uint8_t max = s.pi[0];
        // maybe < N?
        for(int i = 0; i < bsize; ++i){
            if(min > s.pi[i])
                min = s.pi[i];
            if(max < s.pi[i])
                max = s.pi[i];
            s.min[i] = min;
            s.max[i] = max;
        }
        uint8_t lastValue = s.pi[bsize - 1];
        int k = bsize - 2;
        while(k >= 0){
            if(s.pi[k] == lastValue - 1 || s.pi[k] == lastValue + 1)
            {
                lastValue = s.pi[k];
                --k;
            }else{
                s.startUniDir = k + 1;
                break;
            }
        }
        s.chronBL[s.pi[0] - 1]  = s.blocklength[0];
        for(int j = 1; j < bsize; ++j)
            s.chronBL[s.pi[j] - 1] = s.blocklength[j] -  s.blocklength[j - 1];
        for(int j = 1; j < bsize; ++j)
            s.chronBL[j] += s.chronBL[j - 1];
        
        s.revChronBL[s.pi[bsize - 1] - 1]  = s.blocklength[bsize - 1] - s.blocklength[bsize - 2];
        for(int i = bsize - 2; i >= 0; --i){
            s.revChronBL[s.pi[i] - 1] = s.blocklength[i] - ((i > 0) ? s.blocklength[i - 1] : 0);
        }
        for(int i = bsize - 2; i >= 0; --i)
            s.revChronBL[i] += s.revChronBL[i + 1];  
    }
}

    
enum class ReturnCode {
	NOMAPPABILITY, DIRECTSEARCH, COMPMAPPABLE, ONEDIRECTION, MAPPABLE, FINISHED, UNIDIRECTIONAL, SUSPECTUNIDIRECTIONAL, FILTER, ERROR
};

enum class BV {
	RIGHT = 0, MIDDLE = 1, LEFT = 2
};

template <typename TVector, typename TVSupport>
std::vector<std::pair<uint32_t, uint32_t>> getConsOnes(std::vector<std::pair<TVector, TVSupport>> & bitvectors, 
                                             Pair<uint8_t, Pair<uint32_t, uint32_t>> inside_bit_interval,
                                             int const intervalsize);


typedef String<Dna, Alloc<>> TString;
typedef StringSet<TString, Owner<ConcatDirect<> > > TText;
typedef Index<TText, TIndexConfig> MyIndex;
typedef sdsl::bit_vector TBitvector;
typedef sdsl::rank_support_v<> TSupport;
// typedef sdsl::rank_support_v5<> TSupport;
    
}

#endif
