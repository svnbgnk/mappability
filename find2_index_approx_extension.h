#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

// #include <seqan/arg_parse.h>
// #include <seqan/seq_io.h>
#include <seqan/index.h>
#include <sdsl/bit_vectors.hpp>
#include "common.h"


namespace seqan{
/*    
template <size_t N>
struct OptimalSearch
{
    std::array<uint8_t, N> pi; // order of the blocks. permutation of [1..n]
    std::array<uint8_t, N> l; // minimum number of errors at the end of the corresponding block
    std::array<uint8_t, N> u; // maximum number of errors at the end of the corresponding block

    std::array<uint32_t, N> blocklength; // cumulated
    std::array<uint32_t, N> chronblocklength;
    values / prefix sums
    std::array<uint8_t, N> min;
    std::array<uint8_t, N> max;
    
    
    uint32_t startPos; //wrong oder so startPos still get 0
    uint8_t startUniDir;
    
    
};



template <typename TVoidType>
struct OptimalSearchSchemes<0, 2, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 3> VALUE
    {{
        { {{2, 1, 3, 4}}, {{0, 0, 1, 1}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, {{2, 1, 1, 1}}, {{2, 2, 3, 4}}, 2, 0 },
        { {{3, 2, 1, 4}}, {{0, 0, 0, 0}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, {{3, 2, 1, 1}}, {{3, 3, 3, 4}}, 3, 0 },
        { {{4, 3, 2, 1}}, {{0, 0, 0, 2}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, {{4, 3, 2, 1}}, {{4, 4, 4, 4}}, 0, 0 }
    }};
};
    */
enum class ReturnCode {
	NOMAPPABILITY, DIRECTSEARCH, COMPMAPPABLE, ONEDIRECTION, MAPPABLE, FINISHED, UNIDIRECTIONAL, FIRSTTIMEUNIDIRECTIONAL, ERROR
};

template <size_t nbrBlocks, size_t N>
constexpr inline void _optimalSearchSchemeSetMapParams(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{
    for (OptimalSearch<nbrBlocks> & s : ss){
        uint8_t min = s.pi[0];
        uint8_t max = s.pi[0];
        // maybe < N?
        for(int i = 0; i < s.pi.size(); ++i){
            if(min > s.pi[i])
                min = s.pi[i];
            if(max < s.pi[i])
                max = s.pi[i];
            s.min[i] = min;
            s.max[i] = max;
        }
        uint8_t lastValue = s.pi[s.pi.size() - 1];
        int k = s.pi.size() - 2;
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
        s.chronblocklength[s.pi[0] - 1]  = s.blocklength[0];
        for(int j = 1; j < s.blocklength.size(); ++j)
            s.chronblocklength[s.pi[j] - 1] = s.blocklength[j] -  s.blocklength[j - 1];
        for(int j = 1; j < s.blocklength.size(); ++j)
            s.chronblocklength[j] += s.chronblocklength[j - 1];
    }
}

/*    
template< size_t N>
uint8_t mymin(std::array <uint8_t, N> v, uint8_t end)
{
    uint8_t min = v[0];
    for(uint8_t i = 1; i < end; ++i){
        if(v[i] < min)
            min = v[i];
    }
    return min;
}

template< size_t N>
uint8_t mymax(std::array <uint8_t, N> v, uint8_t end){
    uint8_t max = v[0];
    for(uint8_t i = 1; i < end; ++i)
    {
        if(v[i] > max)
            max = v[i];
    }
    return max;
}*/

template <typename TText, typename TIndex, typename TIndexSpec>
void print_fullsa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
              bool const fwd)
{
    int size = seqan::length(iter.fwdIter.index->sa);
    uint32_t number_of_indeces = size - bitvectors[0].first.size();
    int noi = number_of_indeces;
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(!fwd){
        for(uint32_t i = 0; i < size; ++i){
            int seq = iter.revIter.index->sa[i].i1;
            int sa = iter.revIter.index->sa[i].i2;
            cout << i << "(" << seq << ")" << ": " << sa + sequenceLengths[seq] << endl;
        }
    }else{
        for(uint32_t i = 0; i < size; ++i){
            int seq = iter.fwdIter.index->sa[i].i1;
            int sa = iter.fwdIter.index->sa[i].i2;
            cout << i << "(" << seq << ")" << ": " << sa + sequenceLengths[seq] << endl;
        }
    }
}



template <typename TText, typename TIndex, typename TIndexSpec>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
              bool const fwd)
{
    int size = seqan::length(iter.fwdIter.index->sa);
    Pair<uint32_t, uint32_t> dirrange = (fwd) ? range(iter.fwdIter) : range(iter.revIter);
    uint32_t number_of_indeces = size - bitvectors[0].first.size();
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(!fwd){
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.revIter.index->sa[i].i1;
            int sa = iter.revIter.index->sa[i].i2;
            cout << i << ": " << iter.revIter.index->sa[i] << endl;
        }
    }else{
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.fwdIter.index->sa[i].i1;
            int sa = iter.fwdIter.index->sa[i].i2;
            cout << i << ": " << sa + sequenceLengths[seq] << endl;
        }
    }
}

template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks, typename TDir>
Pair<uint8_t, Pair<uint32_t, uint32_t>> get_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                        vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,    
                                        OptimalSearch<nbrBlocks> const & s,
                                        uint8_t const blockIndex,
                                        TDir const & /**/) 
{
//     printv(s.pi);
    
    Pair<uint32_t, uint32_t> dirrange = (std::is_same<TDir, Rev>::value) ? range(iter.fwdIter) : range(iter.revIter);
    uint8_t needed_bitvector;
    if (std::is_same<TDir, Rev>::value)
        needed_bitvector = s.min[blockIndex] - 1;//mymin(s.pi, blockIndex) - 1;
    else
        needed_bitvector = bitvectors.size() - s.max[blockIndex];// + 1 - 1//mymax(s.pi, blockIndex) - 1;
    
    //TODO find out why this does not work my fix this so i dont have to use bitvectors as an input
    //uint32_t number_of_indeces = countSequences(iter.fwdIter.index);
    cout << "maxe: " << (int)s.max[blockIndex] << endl;
    uint32_t number_of_indeces = seqan::length(iter.fwdIter.index->sa) - bitvectors[needed_bitvector].first.size();
    cout << "shift_" << (int)needed_bitvector << "_bitvector" << endl;
    cout << "Sa Range: " << dirrange << endl;
    dirrange.i1 = dirrange.i1 - number_of_indeces;
    dirrange.i2 = dirrange.i2 - number_of_indeces;
    cout << "Bit Range: " << dirrange << endl;
    cout << ((std::is_same<TDir, Rev>::value) ? "reverse case" : "forward case") << endl;
    Pair<uint8_t, Pair<uint32_t, uint32_t>> brange(needed_bitvector, dirrange);
    return brange;
}

void printPair(pair<uint32_t, uint32_t> p){
    cout << "<" << p.first << ", " << p.second << ">";
}

void printb(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, Pair<uint8_t, Pair<uint32_t, uint32_t>> brange){
    sdsl::bit_vector const & rb = bitvectors[brange.i1].first;
    cout << "bitv " << (int)brange.i1 << "brange 1: " << brange.i2.i1 << "  brange 2: " << brange.i2.i2 << endl;
    for(int i = brange.i2.i1; i < brange.i2.i2; ++i)
        cout << i << " Bit: " << rb[i] << endl;
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
bool filterCurrent_interval(TDelegate & delegate,
                     TDelegateD & delegateDirect,
                     Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > & iter,
                     TNeedle const & needle,
                     vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                     uint32_t const needleLeftPos,
                     uint32_t const needleRightPos,
                     uint8_t const errors,
                     OptimalSearch<nbrBlocks> const & s,
                     uint8_t const blockIndex,
                     Pair<uint8_t, Pair<uint32_t, uint32_t>> brange,
                     TDir const & /**/)
{
    uint32_t intervalfilter_size = 3;
    float threshold = 0.5;  //equal or more than half zeroes  
    //TODO may add this as additional input?
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
    sdsl::bit_vector & b = bitvectors[brange.i1].first;
    sdsl::rank_support_v<> & rb = bitvectors[brange.i1].second;
    rb.set_vector(&b);
    
    if(rb(rb.size()) / static_cast<float>(rb.size() - 1) <= threshold){
        uint32_t startPos = brange.i2.i1, endPos = brange.i2.i2;
        for(uint32_t i = startPos; i < endPos; ++i){
            if(b[i] != 0)
                break; 
            ++startPos;
        }
        for(uint32_t i = startPos; i < endPos; ++i){
            if(b[endPos - 1 - i] != 0)
                break;
            --endPos;
        }
        if(startPos > endPos)
            cout << "Error bit vector has only zeroes this should have been checked by check_interval" << endl; 
        cout << "Size: " << endPos - startPos << endl;
        cout << "startPos: " << startPos << " endPos: " << endPos << endl;
        
        vector<pair<uint32_t, uint32_t>> consOnes;
        uint32_t k = startPos;
        uint32_t startOneInterval = startPos;
        while(k < endPos){
            uint32_t interval = 0;
            //TODO delete second condition it should end with 1
            while(b[k + interval] == 0 && (k + interval) < endPos){
                ++interval;
            }
            if(interval >= intervalfilter_size){
                consOnes.push_back(make_pair(startOneInterval, k));
                startOneInterval = k + interval;
            }
            k += interval;
            interval = 0;
            ++k;
        }
        consOnes.push_back(make_pair(startOneInterval, k));
        
        for(int i = 0; i < consOnes.size(); ++i){
            printPair(consOnes[i]);
            if (std::is_same<TDir, Rev>::value){
                iter.revIter.vDesc.range.i1 = consOnes[i].first;
                iter.revIter.vDesc.range.i2 = consOnes[i].second;
                //TODO call function with unidirectional iter
                _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir());
            }
            else
            {
                //TODO does it everything above work for reverse?
                //TODO call function with unidirectional iter
                iter.fwdIter.vDesc.range.i1 = consOnes[i].first;
                iter.fwdIter.vDesc.range.i2 = consOnes[i].second;
                _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir());
            }
        }
        cout << endl;
        return(true);
    }
    else
    {
        return(false);
    }    
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
bool filter_interval(TDelegate & delegate,
                     TDelegateD & delegateDirect,
                     Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > & iter,
                     TNeedle const & needle,
                     vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                     uint32_t const needleLeftPos,
                     uint32_t const needleRightPos,
                     uint8_t const errors,
                     OptimalSearch<nbrBlocks> const & s,
                     uint8_t const blockIndex,
                     Pair<uint8_t, Pair<uint32_t, uint32_t>> brange,
                     TDir const & /**/)
{
    uint32_t intervalfilter_size = 3;
    float threshold = 0.5;  //equal or more than half zeroes     
    sdsl::bit_vector & b = bitvectors[brange.i1].first;
    sdsl::rank_support_v<> & rb = bitvectors[brange.i1].second;
    rb.set_vector(&b);
    
    if(rb(rb.size()) / static_cast<float>(rb.size() - 1) <= threshold){
        uint32_t startPos = brange.i2.i1, endPos = brange.i2.i2;
        for(uint32_t i = startPos; i < endPos; ++i){
            if(b[i] != 0)
                break; 
            ++startPos;
        }
        for(uint32_t i = startPos; i < endPos; ++i){
            if(b[endPos - 1 - i] != 0)
                break;
            --endPos;
        }
        if(startPos > endPos)
            cout << "Error bit vector has only zeroes this should have been checked by check_interval" << endl; 
        cout << "Size: " << endPos - startPos << endl;
        cout << "startPos: " << startPos << " endPos: " << endPos << endl;
        
        vector<pair<uint32_t, uint32_t>> consOnes;
        uint32_t k = startPos;
        uint32_t startOneInterval = startPos;
        while(k < endPos){
            uint32_t interval = 0;
            //TODO delete second condition it should end with 1
            while(b[k + interval] == 0 && (k + interval) < endPos){
                ++interval;
            }
            if(interval >= intervalfilter_size){
                consOnes.push_back(make_pair(startOneInterval, k));
                startOneInterval = k + interval;
            }
            k += interval;
            interval = 0;
            ++k;
        }
        consOnes.push_back(make_pair(startOneInterval, k));
        
        for(int i = 0; i < consOnes.size(); ++i){
            printPair(consOnes[i]);
            if (std::is_same<TDir, Rev>::value){
                iter.revIter.vDesc.range.i1 = consOnes[i].first;
                iter.revIter.vDesc.range.i2 = consOnes[i].second;
                //TODO already checked if we are in the last block so no need for std::min
                //TODO call function with unidirectional iter
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Rev());
            }
            else
            {
                //TODO does it everything above work for reverse?
                iter.fwdIter.vDesc.range.i1 = consOnes[i].first;
                iter.fwdIter.vDesc.range.i2 = consOnes[i].second;
                //TODO call function with unidirectional iter
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Fwd());
            }
        }
        cout << endl;
        return(true);
    }
    else
    {
        return(false);
    }    
}

template<size_t nbrBlocks>
ReturnCode check_interval(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                           Pair<uint8_t, Pair<uint32_t, uint32_t>> brange,
                           OptimalSearch<nbrBlocks> const & s)
{
    sdsl::rank_support_v<> & rb = bitvectors[brange.i1].second;
    rb.set_vector(&bitvectors[brange.i1].first);
    uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
//     sdsl::rank_support_v<> rbi(& bi);
    if(ivalOne == 0)
        return ReturnCode::NOMAPPABILITY;
    
    /*
    if(ivalOne <= 3){ // add additional constrains from chris?
        return ReturnCode::DIRECTSEARCH;
    }
    if(ivalOne == (brange.i2.i2 - brange.i2.i1))
        return ReturnCode::COMPMAPPABLE;
    //Only Case FIRSTTIMEUNIDIRECTIONAL doesnt get covered here
//     if(s.startUniDir >= blockIndex)
//         return ReturnCode::UNIDIRECTIONAL
    */
    return ReturnCode::MAPPABLE;
}

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
                  TDir const & /**/)
{
    cout << "Test Dummy Function: " << endl;
    
//     for(int i = 0; i < seqan::length(iter.index->sa); ++i)
//         cout << iter.index->sa[i] << endl; 
    
//     cout << "Length of Sequences" << endl;
//     for(int i = 0; i < 4; ++i)
//         cout << olditer.fwdIter.index->sa[i] << endl;
    cout << "Ranges" << endl;
    cout <<  iter.vDesc.range << endl;
    cout << "Old Iter ranges" << endl;
    cout << olditer.fwdIter.vDesc.range << endl;
    cout << "Reverse" << endl;
    cout << olditer.revIter.vDesc.range << endl;
        //TODO dont care about different sequences in the index in the moment
    int needleSize = needleRightPos- needleLeftPos - 1;
    cout << "Size of searched needle: " << needleSize << endl;
    vector<int> saveOppositePositions;

    
    int number_of_indeces = seqan::length(olditer.fwdIter.index->sa) - bitvectors[0].first.size();
    vector<int> sequenceLengths(number_of_indeces, 0); 
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[olditer.fwdIter.index->sa[i].i1] = olditer.fwdIter.index->sa[i].i2;
    
    printv(sequenceLengths);
    
    auto fwdr = olditer.fwdIter.vDesc.range;
    auto revr = olditer.revIter.vDesc.range;
    
    if(std::is_same<TDir, Rev>::value){
        cout << "Reverse case therefore need to save fwd sa position" << endl;
        for(int i = fwdr.i1; i < fwdr.i2; ++i){
            cout << olditer.fwdIter.index->sa[i] << endl;
        }
        cout << "other sa" << endl;
        for(int i = revr.i1; i < revr.i2; ++i){
            cout << olditer.revIter.index->sa[i] << endl;
        }
        
        for(int i = revr.i1; i < revr.i2; ++i){
            int seq = olditer.revIter.index->sa[i].i1;
            cout << "Seq Length" << sequenceLengths[seq] << endl;
            int findn = sequenceLengths[seq] - static_cast<int>(olditer.revIter.index->sa[i].i2) - needleSize;
            cout << "findn  " << findn << endl;
            int k = fwdr.i1;
            while(k < fwdr.i2){
                if(findn == olditer.fwdIter.index->sa[k].i2)
                    break;
                ++k;
            }
            saveOppositePositions.push_back(k);
        }     
    }else{
        cout << "Forward case therefore need to save rev sa position" << endl;
        for(int i = revr.i1; i < revr.i2; ++i){
            cout << olditer.revIter.index->sa[i] << endl;
        }
        cout << "other sa" << endl;
        for(int i = fwdr.i1; i < fwdr.i2; ++i){
            cout << olditer.fwdIter.index->sa[i] << endl;
        }
        
        for(int i = fwdr.i1; i < fwdr.i2; ++i){
            int seq = olditer.fwdIter.index->sa[i].i1;
            int findn = sequenceLengths[seq] - olditer.fwdIter.index->sa[i].i2 - needleSize;
            cout << olditer.fwdIter.index->sa[i].i2 << " equivalent to  " << findn << endl;
            int k = revr.i1;
            while(k < revr.i2){
                if(findn == olditer.revIter.index->sa[k].i2)
                    break;
                ++k;
            }
            saveOppositePositions.push_back(k);
        }
    }
        
        
    cout << "Saved Positions" << endl;
    for(int i = 0; i < saveOppositePositions.size(); ++i)
        cout << saveOppositePositions[i] << endl;
    cout << "End" << endl;
    
}

template <typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
void directSearch(TDelegateD & delegateDirect,
                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, 
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  Pair<uint8_t, Pair<uint32_t, uint32_t>> brange,
                  TDir const & /**/)
{
    //TODO  stop using both needle pos ?
    bool reverse = std::is_same<TDir, Rev>::value;
    StringSet<DnaString> const & genome = (reverse) ? indexText(*iter.fwdIter.index) : indexText(*iter.revIter.index);
    vector<Pair<uint16_t, uint32_t>> hitsv;
    vector<uint8_t> errorsv;
    
    for(int i = 0; i < brange.i2.i2 - brange.i2.i1; ++i){
        cout << "Direct Search" << endl;
        cout << "NLP " << needleLeftPos <<  endl;
        cout << "NRP " <<  needleRightPos <<  endl;
        if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
            uint8_t errors2 = errors;
            bool valid = true;
            Pair<uint16_t, uint32_t> sa_info;
            uint32_t startPos;
            if(reverse){
                sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + i];
                startPos = sa_info.i2 - needleLeftPos;
            }
            else
            {
                sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + i];
                startPos = sa_info.i2 + needleRightPos - 1;
            }
            // iter.fwdIter.vDesc.range.i1 is not the same brange.i2.i1 since sentinels are at the beginning!!!
            cout <<  "Sa info" <<  sa_info <<  endl;
            cout << "StartPos " << startPos << endl;
            //search remaining blocks
            for(int j = blockIndex; j < s.pi.size(); ++j){
                int blockStart = (s.pi[j] - 2 == 0) ? 0 : s.chronblocklength[s.pi[j] - 2];
                if(reverse)
                    cout << "searching Parts:" << startPos + blockStart << " - " << startPos + s.chronblocklength[s.pi[j] - 1] << "; ";
                else
                    cout << "searching Parts:" << startPos - blockStart << " - " << startPos - s.chronblocklength[s.pi[j] - 1] << "; ";
                // compare bases to needle
                for(int k = blockStart; k <  s.chronblocklength[s.pi[j] - 1]; ++k){
                    int sign = (reverse) ? 1 : -1;
                    if(needle[k] != genome[sa_info.i1][startPos + sign * k])
                        ++errors2;
                }
                if(errors2 < s.l[j] || errors2 > s.u[j]){
                    cout << "Triggered: " << (int)errors2 << endl;
                    valid = false;
                    break;
                }
            }
            if(valid){
                cout << "Hiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiit" << endl;
                cout << (int)errors2 << endl;
                uint32_t occ = (reverse) ? startPos : seqan::length(genome[sa_info.i1]) - startPos - 1;
                hitsv.push_back(Pair<uint16_t,uint32_t>(sa_info.i1, occ));
                cout << "Hit occ: " << hitsv[hitsv.size() - 1] << endl;
                errorsv.push_back(errors2);
            }
        }
    }
    delegateDirect(hitsv, needle, errorsv);
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
ReturnCode checkCurrentMappability(TDelegate & delegate,
                            TDelegateD & delegateDirect,
                            Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                            TNeedle const & needle,
                            vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,    
                            uint32_t const needleLeftPos,
                            uint32_t const needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            TDir const & /**/)
{
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval = get_bitvector_interval(iter, s, blockIndex, TDir());
    ReturnCode rcode = check_interval(bitvectors, bit_interval, s);
//   cout << "Return code: " << (int)rcode << endl;
        
    if(rcode == ReturnCode::NOMAPPABILITY)
        return ReturnCode::FINISHED;        
        
    if(rcode == ReturnCode::DIRECTSEARCH){
        //search directly in Genome
        //TODO I only need left value for rev and right value for fwd so delte one input?
        directSearch(delegateDirect, iter, needle, bitvectors, needleRightPos - 1 , needleRightPos - 1, errors, s, blockIndex, bit_interval, TDir());
        return ReturnCode::FINISHED;
    }
    
    //TODO check if we are in unidirection case in this case COMPMAPPABLE should no be used
    if(rcode == ReturnCode::COMPMAPPABLE){
        uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
        _optimalSearchSchemeChildren(delegate, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), HammingDistance());
        return ReturnCode::FINISHED;
    }
    //TODO check if we are in one direction case  and apply filter interval
    if(rcode == ReturnCode::UNIDIRECTIONAL){
        bool did_filter = filterCurrent_interval(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s,blockIndex, bit_interval, TDir());
        if(did_filter)
            return ReturnCode::FINISHED;
    }
    if(s.startUniDir == blockIndex)
        return ReturnCode::FIRSTTIMEUNIDIRECTIONAL;
    return ReturnCode::MAPPABLE;
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
ReturnCode checkMappability(TDelegate & delegate,
                            TDelegateD & delegateDirect,
                            Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                            TNeedle const & needle,
                            vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,    
                            uint32_t const needleLeftPos,
                            uint32_t const needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            bool const goToRight2,
                            TDir const & /**/)
{
    uint32_t infixPosLeft, infixPosRight;
    bool finished = false;

    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it reverse
        infixPosLeft = needleRightPos - 1;
        infixPosRight = needleLeftPos + s.blocklength[blockIndex - 1] - 1;
        //TODO maybe this one is the same as in the fwd case?
        finished = infixPosLeft == 0 && infixPosRight + 2 == length(needle) + 1;
    }
    else{
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        infixPosLeft = needleRightPos - s.blocklength[blockIndex - 1] - 1;
        infixPosRight = needleLeftPos - 1;
        finished = infixPosLeft == 0 && needleRightPos == length(needle) + 1;
    }

    //check if we are done with the needle    
    if(!finished){
        Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
        if(goToRight2)
            bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, Rev());
        else
            bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, Fwd());
        
        cout << "Printttt" << endl;
        if(goToRight2)
            print_sa(iter, bitvectors, true);
        else
            print_sa(iter, bitvectors, false);
//         cout << "Print full sa" << endl;
//         print_fullsa(iter, bitvectors, true);
        printb(bitvectors, bit_interval);
        
        
        ReturnCode rcode = check_interval(bitvectors, bit_interval, s);
//      cout << "Return code: " << (int)rcode << endl;
        
        if(rcode == ReturnCode::NOMAPPABILITY){
            cout << "NOMAPPABILITYNOMAPPABILITYNOMAPPABILITYNOMAPPABILITYNOMAPPABILITYNOMAPPABILITY" << endl;
            return ReturnCode::FINISHED;        
        }
        if(rcode == ReturnCode::DIRECTSEARCH){
            //search directly in Genome
            //TODO I only need left value for rev and right value for fwd so delte one input?
            if(goToRight2){
                directSearch(delegateDirect, iter, needle, bitvectors, needleRightPos - 1 , needleRightPos - 1, errors, s, blockIndex, bit_interval, Rev());
//                 directSearchDummy(delegateDirect, iter.fwdIter, iter, needle, bitvectors, needleLeftPos , infixPosRight + 2, errors, s, blockIndex, bit_interval, Rev());
            }else{
                directSearch(delegateDirect, iter, needle, bitvectors, needleRightPos - 1 , needleRightPos - 1, errors, s, blockIndex, bit_interval, Fwd());
//                 directSearchDummy(delegateDirect, iter.fwdIter, iter, needle, bitvectors, infixPosLeft , needleRightPos, errors, s, blockIndex, bit_interval, Fwd());
            }
            return ReturnCode::FINISHED;
        }
        //TODO check if we are in unidirection case in this case COMPMAPPABLE should no be used
        if(rcode == ReturnCode::COMPMAPPABLE){
            if(std::is_same<TDir, Rev>::value){
                if (goToRight2)
                {
                    _optimalSearchScheme(delegate, iter, needle, needleLeftPos, infixPosRight + 2, errors, s, blockIndex, Rev(), HammingDistance());
                }
                else
                {
                    _optimalSearchScheme(delegate, iter, needle, needleLeftPos, infixPosRight + 2, errors, s, blockIndex, Fwd(), HammingDistance());
                }
            }
            else
            { 
                if (goToRight2)
                {
                    _optimalSearchScheme(delegate, iter, needle, infixPosLeft, needleRightPos, errors, s, blockIndex, Rev(), HammingDistance());
                }
                else
                {
                    _optimalSearchScheme(delegate, iter, needle, infixPosLeft, needleRightPos, errors, s, blockIndex, Fwd(), HammingDistance());
                }
            }
            return ReturnCode::FINISHED;
        }
        if(rcode == ReturnCode::MAPPABLE && s.startUniDir >= blockIndex){
            bool did_filter;
            if(goToRight2){
                did_filter = filter_interval(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s,blockIndex, bit_interval, Rev());
            }
            else
            {
                did_filter = filter_interval(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s,blockIndex, bit_interval, Fwd());
            }
            
            //TODO iter changes from iter -> iter.revIter now has the Type: ??? or do it in the function filter_interval
            if(did_filter)
                return ReturnCode::FINISHED;        
        }
        if(s.startUniDir == blockIndex)
            return ReturnCode::FIRSTTIMEUNIDIRECTIONAL;
    }
    return ReturnCode::MAPPABLE;
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeChildren(TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,   
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const & /**/)
{
    bool goToRight = std::is_same<TDir, Rev>::value;
    if (goDown(iter, TDir()))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()),
                                   needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1]);
            if (minErrorsLeftInBlock > 0 && charsLeft + delta < minErrorsLeftInBlock + 1u) // charsLeft - 1 < minErrorsLeftInBlock - delta
            {
                continue;
            }
            int32_t needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t needleRightPos2 = needleRightPos + goToRight;

            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
            
//                 ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, goToRight2, TDir());
//                 //TDir() is in this case Rev() ...
//                 if(rcode == ReturnCode::FINISHED)
//                     return;
//                 if(rcode == ReturnCode::FIRSTTIMEUNIDIRECTIONAL){
//                     //TODO call same with unidirectional iter?? then stop 
//                 }
                
                if (goToRight2)
                {
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev());
                }
                else
                {
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd());
                }
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, TDir());
            }
        } while (goRight(iter, TDir()));
    }
}   

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeExact(TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,    
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const & /**/)
{
    // not in last block and next Block is larger then current block
    bool goToRight2 = (blockIndex < s.pi.size() - 1) && s.pi[blockIndex + 1] > s.pi[blockIndex];
    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it reverse
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

        cout << "SearchScheme: " << endl;
        printv(s.pi);
        cout << "Needleeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee" << endl;
         for(int i = infixPosLeft; i < infixPosRight + 1; ++i)
            cout << needle[i];
        cout << endl;
        
        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
            return;

        //TODO Check if we are Done here?
        //TODO does something go wrong in checkMappability if blockIndex is not increased by one? then do Todo before
        cout << "checkMappability Rev Index" << endl;
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), goToRight2, TDir());
        //TDir() is in this case Rev() ...
        if(rcode == ReturnCode::FINISHED)
            return;
        if(rcode == ReturnCode::FIRSTTIMEUNIDIRECTIONAL){
            //TODO call same with unidirectional iter?? what if iter is already unidirectional then stop!!
        }
            
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd());
        }
    }
    else
    {
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        int32_t infixPosLeft = needleRightPos - s.blocklength[blockIndex] - 1;
        int32_t infixPosRight = needleLeftPos - 1;

        cout << "SearchScheme: " << endl;
        printv(s.pi);
        cout << "Needleeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee(Fwd Idx)" << endl;
        for(int i = infixPosLeft; i < infixPosRight + 1; ++i)
            cout << needle[i];
        cout << endl;
        
        while (infixPosRight >= infixPosLeft)
        {
            if (!goDown(iter, needle[infixPosRight], TDir()))
                return;
            --infixPosRight;
        }

        //TODO Check if we are Done here?
        cout << "checkMappability Fwd Index" << endl;
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), goToRight2, TDir());
        //TDir() is in this case Fwd() ...
        if(rcode == ReturnCode::FINISHED)
            return;
        
        if(rcode == ReturnCode::FIRSTTIMEUNIDIRECTIONAL){
            //TODO call same with unidirectional iter?? what if iter is already unidirectional then stop!!
        }
        
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd());
        }
    }
}
    

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,   
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,                             
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 TDir const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

//     cout << "Step: " << needleRightPos - needleLeftPos - 1 << "    ss: "; printv(s.pi); cout << endl;
    // Done. (Last step)
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        delegate(iter, needle, errors);
    }
    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
    }
    // Approximate search in current block.
    else
    {
   /*   
        //TODO implement faster mod 
        if((needleRightPos - needleLeftPos - 1) % 4 == 0){
            //TODO Check if we are Done here?
            //TODO stop doing on loop to much (enter checkCurrentMappability a second time)
            bool goToRight2 = std::is_same<TDir, Rev>::value;
            ReturnCode rcode = checkCurrentMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, goToRight2, TDir());
            // blockIndex is different here!!! goToRight2 is missing
            if(rcode == ReturnCode::FINISHED)
                return;
        }*/
        
        _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir());
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                                 OptimalSearch<nbrBlocks> const & s)
{
    bool initialDirection = s.pi[1] > s.pi[0];
    if(initialDirection)
        _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Rev());
    else
        _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Fwd());
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks, size_t N>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss)
{
    for (auto & s : ss)
        _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s);
}  

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors)
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeSetMapParams(scheme);
    print_search_scheme(scheme);
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, scheme);
}
  
  

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TParallelTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
     TParallelTag const & /**/)
{
    typedef typename Iterator<StringSet<TNeedle, TStringSetSpec> const, Rooted>::Type TNeedleIt;
    typedef typename Reference<TNeedleIt>::Type                                       TNeedleRef;
    iterate(needles, [&](TNeedleIt const & needleIt)
    {
        TNeedleRef needle = value(needleIt);
        find<minErrors, maxErrors>(delegate, delegateDirect, index, needle, bitvectors);
    },
    Rooted(), TParallelTag());
}    
    
template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors)
{
    find<minErrors, maxErrors>(delegate, delegateDirect, index, needles, bitvectors, Serial());
}

}

#endif
