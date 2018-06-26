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

    std::array<uint32_t, N> blocklength; // cumulated values / prefix sums
    std::array<uint8_t, N> min;
    std::array<uint8_t, N> max;
    uint8_t startUniDir;
    uint32_t startPos;
};



template <typename TVoidType>
struct OptimalSearchSchemes<0, 2, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 3> VALUE
    {{
        { {{2, 1, 3, 4}}, {{0, 0, 1, 1}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, {{2, 1, 1, 1}}, {{2, 2, 3, 4}}, 0, 2 },
        { {{3, 2, 1, 4}}, {{0, 0, 0, 0}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, {{3, 2, 1, 1}}, {{3, 3, 3, 4}}, 0, 3 },
        { {{4, 3, 2, 1}}, {{0, 0, 0, 2}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, {{4, 3, 2, 1}}, {{4, 4, 4, 4}}, 0, 0 }
    }};
};
    */
enum class ReturnCode {
	NOMAPPABILITY, DIRECTSEARCH, COMPMAPPABLE, ONEDIRECTION, MAPPABLE, FINISHED, ERROR
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
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              Pair<uint32_t, uint32_t>  const & range)
{
    int number_of_indeces = countSequences(iter.fwdIter.index);
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    for(uint32_t i = range.i1; i < range.i2; ++i){
        int seq = iter.fwdIter.index->sa[i].i1;
        int sa = iter.fwdIter.index->sa[i].i2;
        cout << i << ": " << sa + sequenceLengths[seq] << endl;
    }
}

template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks>
Pair<uint8_t, Pair<uint32_t, uint32_t>> get_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                        OptimalSearch<nbrBlocks> const & s,
                                        uint8_t const blockIndex,
                                        bool const goToRight2) 
{
    Pair<uint32_t, uint32_t> dirrange = (goToRight2) ? range(iter.fwdIter) : range(iter.revIter);
    uint8_t needed_bitvector;
//     sdsl::bit_vector b(dirrange.i2 - dirrange.i1, 0);
    if (goToRight2)
        needed_bitvector = s.min[blockIndex];//mymin(s.pi, blockIndex) - 1;
    else
        needed_bitvector = s.max[blockIndex];//mymax(s.pi, blockIndex);
    uint32_t number_of_indeces = countSequences(iter.fwdIter.index);
//     cout << "Min Element: " << (int)min << endl;
    cout << "shift_" << (int)needed_bitvector << "_bitvector" << endl;
    
//     for(int i = dirrange.i1; i < dirrange.i2; ++i)
//         b[i - dirrange.i1] = bitvectors.at(needed_bitvector)[i + number_of_indeces];
    dirrange.i1 = dirrange.i1 + number_of_indeces;
    dirrange.i2 = dirrange.i2 + number_of_indeces;
    Pair<uint8_t, Pair<uint32_t, uint32_t>> brange(needed_bitvector, dirrange); 
    return brange;
}

void printPair(pair<uint32_t, uint32_t> p){
    cout << "<" << p.first << ", " << p.second << ">";
}

void printb(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, Pair<uint8_t, Pair<uint32_t, uint32_t>> brange){
    sdsl::bit_vector const & rb = bitvectors[brange.i1].first;
    cout << "Ranksupport vector size " << rb.size() << endl;
    cout << "bitv " << (int)brange.i1 << "brange 1 " << brange.i2.i1 << "brange 2 " << brange.i2.i2 << endl;
    for(int i = brange.i2.i1; i < brange.i2.i2; ++i)
        cout << i << " Bit: " << rb[i] << endl;
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks>
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
                     bool const goToRight)
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
        
        vector<pair<uint32_t, uint32_t>> consZeroes;
        uint32_t k = startPos;
        uint32_t startOneInterval = startPos;
        while(k < endPos){
            uint32_t interval = 0;
            //TODO delete second condition it should end with 1
            while(b[k + interval] == 0 && (k + interval) < endPos){
                ++interval;
            }
            if(interval >= intervalfilter_size){
                consZeroes.push_back(make_pair(startOneInterval, k));
                startOneInterval = k + interval;
            }
            k += interval;
            interval = 0;
            ++k;
        }
        consZeroes.push_back(make_pair(startOneInterval, k));
        
        for(int i = 0; i < consZeroes.size(); ++i){
            printPair(consZeroes[i]);
            if (goToRight){
                iter.revIter.vDesc.range.i1 = consZeroes[i].first;
                iter.revIter.vDesc.range.i2 = consZeroes[i].second;
                //TODO already checked if we are in the last block so no need for std::min
                //TODO call function with unidirectional iter
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev());
            }
            else
            {
                //TODO does it everything above work for reverse?
                iter.fwdIter.vDesc.range.i1 = consZeroes[i].first;
                iter.fwdIter.vDesc.range.i2 = consZeroes[i].second;
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd());
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

ReturnCode check_interval(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                           Pair<uint8_t, Pair<uint32_t, uint32_t>> brange)
{
    sdsl::rank_support_v<> & rb = bitvectors[brange.i1].second;
    rb.set_vector(&bitvectors[brange.i1].first);
    uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
//     sdsl::rank_support_v<> rbi(& bi);
    if(ivalOne == 0)
        return ReturnCode::NOMAPPABILITY;
    if(ivalOne <= 3){ //TODO add additional constrains from chris
        return ReturnCode::DIRECTSEARCH;
    }
    if(ivalOne == (brange.i2.i2 - brange.i2.i1))
        return ReturnCode::COMPMAPPABLE; 
    return ReturnCode::MAPPABLE;
}

template <typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks>
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
                  bool const goToRight2)
{
    //TODO  stop using both needle pos ?
    bool reverse = goToRight2;
    StringSet<DnaString> const & genome = (reverse) ? indexText(*iter.fwdIter.index) : indexText(*iter.revIter.index);
    vector<Pair<uint16_t, uint32_t>> hitsv;
    vector<uint8_t> errorsv;
    
    //TODO export to _optimalSearchSchemeSetMapParams and/or set OptimalSearch Struct for values??
    std::vector<int> bl (s.pi.size() + 1,0);
    cout << "Print blocklengths" << endl;
    bl[s.pi[0]]  = s.blocklength[0];
    for(int j = 1; j < s.blocklength.size(); ++j)
        bl[s.pi[j]] = s.blocklength[j] -  s.blocklength[j - 1];
    for(int j = 1; j < bl.size(); ++j)
        bl[j] += bl[j - 1];
    for(int j = 0; j < bl.size(); ++j)
        cout << bl[j] << " ";
    cout << endl;
    
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
                if(reverse)
                    cout << "searching Parts:" << startPos + bl[s.pi[j] - 1] << " - " << startPos + bl[s.pi[j]] << "; ";
                else
                    cout << "searching Parts:" << startPos - bl[s.pi[j] - 1] << " - " << startPos - bl[s.pi[j]] << "; ";
                // compare bases to needle
                for(int k = bl[s.pi[j] - 1]; k <  bl[s.pi[j]]; ++k){
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
    //Are sure we can not be finished
    //TODO transform goToRight2 to TDir in functions
    bool goToRight2 = (std::is_same<TDir, Rev>::value);
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval = get_bitvector_interval(iter, s, blockIndex, goToRight2);
    ReturnCode rcode = check_interval(bitvectors, bit_interval);
//   cout << "Return code: " << (int)rcode << endl;
        
    if(rcode == ReturnCode::NOMAPPABILITY)
        return ReturnCode::FINISHED;        
        
    if(rcode == ReturnCode::DIRECTSEARCH){
        //search directly in Genome
        //TODO can remove std::min since i already checked if we were finished
        //TODO I only need left value for rev and right value for fwd so delte one input?
        directSearch(delegateDirect, iter, needle, bitvectors, needleRightPos - 1 , needleRightPos - 1, errors, s, blockIndex, bit_interval, goToRight2);
        return ReturnCode::FINISHED;
    }
    
    //TODO check if we are in unidirection case in this case COMPMAPPABLE should no be used
    if(rcode == ReturnCode::COMPMAPPABLE){
        uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
        _optimalSearchSchemeChildren(delegate, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), HammingDistance());
        return ReturnCode::FINISHED;
    }
    //TODO check if we are in one direction case  and apply filter interval
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
                            bool const newBlock,
                            TDir const & /**/)
{
    uint32_t infixPosLeft, infixPosRight;
    bool finished = false;
    //TODO Does something go wrong if blockIndex == 0 (in case we come from _optimalSearchScheme doing approxiamte search -> check inf blockIndex == 0 or newBlock skip first condition (only need these values for newBlock condition))
//     if(newBlock){
    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it reverse
        infixPosLeft = needleRightPos - 1;
        infixPosRight = needleLeftPos + s.blocklength[blockIndex - 1] - 1;
        //TODO maybe this one is the same as in the fwd case?
        finished = infixPosLeft != 0 || infixPosRight + 2 != length(needle) + 1;
    }
    else{
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        infixPosLeft = needleRightPos - s.blocklength[blockIndex - 1] - 1;
        infixPosRight = needleLeftPos - 1;
        finished = infixPosLeft != 0 || needleRightPos != length(needle) + 1;
    }
//     }
    //check if we are done with the needle
    // new check mappability of the next block
    if(!finished){
        Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval = get_bitvector_interval(iter, s, blockIndex, goToRight2);
        ReturnCode rcode = check_interval(bitvectors, bit_interval);
//      cout << "Return code: " << (int)rcode << endl;
        
        if(rcode == ReturnCode::NOMAPPABILITY)
            return ReturnCode::FINISHED;        
        
        if(rcode == ReturnCode::DIRECTSEARCH){
            //search directly in Genome
            //TODO can remove std::min since i already checked if we were finished
            //TODO I only need left value for rev and right value for fwd so delte one input?
            directSearch(delegateDirect, iter, needle, bitvectors, needleRightPos - 1 , needleRightPos - 1, errors, s, blockIndex, bit_interval, goToRight2);
            return ReturnCode::FINISHED;
        }
        if(newBlock){
            //TODO this can also happen if we are not in a newBlock
            //TODO no const on blockIndex ?
            //TODO remove newBlock condition
            if(rcode == ReturnCode::COMPMAPPABLE){
                if(std::is_same<TDir, Rev>::value){
                    if (goToRight2)
                    {
                        _optimalSearchScheme(delegate, iter, needle, needleLeftPos, infixPosRight + 2, errors, s, std::min(blockIndex + 0, static_cast<uint8_t>(s.u.size()) - 1), Rev(), HammingDistance());
                    }
                    else
                    {
                        _optimalSearchScheme(delegate, iter, needle, needleLeftPos, infixPosRight + 2, errors, s, std::min(blockIndex + 0, static_cast<uint8_t>(s.u.size()) - 1), Fwd(), HammingDistance());
                    }
                }
                else
                { 
                    if (goToRight2)
                    {
                        _optimalSearchScheme(delegate, iter, needle, infixPosLeft, needleRightPos, errors, s, std::min(blockIndex + 0, static_cast<uint8_t>(s.u.size()) - 1), Rev(), HammingDistance());
                    }
                    else
                    {
                        _optimalSearchScheme(delegate, iter, needle, infixPosLeft, needleRightPos, errors, s, std::min(blockIndex + 0, static_cast<uint8_t>(s.u.size()) - 1), Fwd(), HammingDistance());
                    }
                }
                return ReturnCode::FINISHED;
            }
            if(rcode == ReturnCode::MAPPABLE && s.startUniDir >= blockIndex){ 
                bool did_filter = filter_interval(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s,blockIndex, bit_interval, goToRight2);
                //TODO iter changes from iter -> iter.revIter now has the Type: ??? or do it in the function filter_interval
                if(did_filter)
                    return ReturnCode::FINISHED;
            }
        }
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
/*                
                //TODO Check if we are Done here?
                //TODO does something go wrong in checkMappability if blockIndex is not increased by one? then do Todo before
                //TODO check all Input parameters again!!
                ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, goToRight2, true, TDir());
                //TDir() is in this case Rev() ...
                if(rcode == ReturnCode::FINISHED)
                    return;
                
                */
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

        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
            return;

        //TODO Check if we are Done here?
        //TODO does something go wrong in checkMappability if blockIndex is not increased by one? then do Todo before
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), goToRight2, true, TDir());
        //TDir() is in this case Rev() ...
        if(rcode == ReturnCode::FINISHED)
            return;
        
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

        while (infixPosRight >= infixPosLeft)
        {
            if (!goDown(iter, needle[infixPosRight], TDir()))
                return;
            --infixPosRight;
        }
        
        //TODO Check if we are Done here?
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), goToRight2, true, TDir());
        //TDir() is in this case Fwd() ...
        if(rcode == ReturnCode::FINISHED)
            return;
        
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
        if((needleRightPos - needleLeftPos - 1) % 4 == 0){
            //TODO Check if we are Done here?
            //TODO stop doing on loop to much (enter checkMappability a second time)
            bool goToRight2 = std::is_same<TDir, Rev>::value;
            ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, goToRight2, false, TDir());
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
    _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Rev());
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
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
//     auto & testIndex = ;
//     auto & genome2 = indexText(*it.fwdIter.index);
// //     auto & testIndex = index(it)
//     for(int i = 0; i < seqan::length(genome2[1]); ++i){
//         cout << genome2[1][i];
//     }
//     cout << endl;

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
