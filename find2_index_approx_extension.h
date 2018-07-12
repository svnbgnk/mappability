#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

#include <seqan/index.h>
#include <sdsl/bit_vectors.hpp>
#include "common.h"
#include "common_auxiliary.h"
#include "find2_index_approx_unidirectional.h"
#include "find2_index_approx_compmappable.h"

template <typename TIter>
struct isBidirectionalIter
{
     static constexpr bool VALUE = false;
};

template <typename TText, typename TIndex, typename TIndexSpec>
struct isBidirectionalIter<Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > >
{
     static constexpr bool VALUE = true;
};

namespace seqan{
/*    
template <size_t N>
struct OptimalSearch
{
    std::array<uint8_t, N> pi; // order of the blocks. permutation of [1..n]
    std::array<uint8_t, N> l; // minimum number of errors at the end of the corresponding block
    std::array<uint8_t, N> u; // maximum number of errors at the end of the corresponding block

    std::array<uint32_t, N> blocklength; // cumulated values / prefix sums
    std::array<uint32_t, N> chronBL;
    std::array<uint32_t, N> revChronBL;
    std::array<uint8_t, N> min;
    std::array<uint8_t, N> max;
    uint32_t startPos; //wrong position so i still get 0 from initialization
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
	NOMAPPABILITY, DIRECTSEARCH, COMPMAPPABLE, ONEDIRECTION, MAPPABLE, FINISHED, UNIDIRECTIONAL, SUSPECTUNIDIRECTIONAL, ERROR
};

enum class BV {
	RIGHT = 0, MIDDLE = 1, LEFT = 2
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

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
void filter_interval(TDelegate & delegate,
                     TDelegateD & delegateDirect,
                     Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                     TNeedle const & needle,
                     vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                     uint32_t const needleLeftPos,
                     uint32_t const needleRightPos,
                     uint8_t const errors,
                     OptimalSearch<nbrBlocks> const & s,
                     uint8_t const blockIndex,
                     Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                     TDir const & /**/)
{
    uint32_t intervalfilter_size = 3;      
    sdsl::bit_vector & b = bitvectors[inside_bit_interval.i1].first;
/*    
 * //TODO need bitvector 
    int n;
    if(std::is_same<TDir, Rev>::value)
        n = 4 + s.max[blockIndex];
    else
        n = s.min[blockIndex];
    uint32_t number_of_indeces = seqan::length(iter.fwdIter.index->sa) - bitvectors[needed_bitvector].first.size();
    cout << "selected bitvector: " << (int)needed_bitvector << endl;
//     cout << "Sa Range: " << dirrange << endl;  
    dirrange.i1 = dirrange.i1 - number_of_indeces;
    dirrange.i2 = dirrange.i2 - number_of_indeces;
    //TODO create correct range 
    */

    cout << "In filterInterval" << endl;
    printbit(bitvectors, inside_bit_interval);
    
    if(std::is_same<TDir, Rev>::value)
        print_sa(iter, bitvectors, true);
    else
        print_sa(iter, bitvectors, false);
    vector<pair<uint32_t, uint32_t>> consOnes;
    
    uint32_t k = inside_bit_interval.i2.i1;
    uint32_t startOneInterval = inside_bit_interval.i2.i1;
    while(k < inside_bit_interval.i2.i2){
        uint32_t interval = 0;
        //TODO delete second condition it should end with 1
        while(b[k + interval] == 0 && (k + interval) < inside_bit_interval.i2.i2){
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
//     consOnes.push_back(make_pair(inside_bit_interval.i2.i1, inside_bit_interval.i2.i2));
    uint32_t noi = seqan::length(iter.fwdIter.index->sa) - bitvectors[0].first.size(); // number_of_indeces
    
    for(int i = 0; i < consOnes.size(); ++i){
         cout << "Print Bit Range";
         printPair(consOnes[i]);
         cout << endl;
        if (std::is_same<TDir, Rev>::value){
            //TODO call DirectSearch here if the interval is to small also use block Index ....?
            iter.revIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.revIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(delegate, delegateDirect, iter.revIter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Rev());
            cout << "reverse case inside" << endl;
            cout << "Print SA: " << endl;
            print_sa(iter, bitvectors, false);
        }
        else
        {
            //TODO does it everything above work for reverse?
            iter.fwdIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.fwdIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(delegate, delegateDirect, iter.fwdIter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Fwd());
            cout << "forward case inside" << endl;
            cout << "Print SA: " << endl;
            print_sa(iter, bitvectors, true);
        }
    } 
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
                  Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                  TDir const & /**/)
{
    cout << "directSearchdefault " <<  endl;
    cout << "NLP: " <<  needleLeftPos <<  endl;
    cout << "NRP: " <<  needleRightPos <<  endl;
    cout << "errors:  " <<  (int)errors <<  endl;
    vector<Pair<uint16_t, uint32_t>> hitsv;
    vector<uint8_t> errorsv;
    auto const & genome = indexText(*iter.fwdIter.index);
    for(int i = 0; i < brange.i2.i2 - brange.i2.i1; ++i){
        cout << "I: " << i << endl;
        if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
            cout << "blockIndex: " << (int)blockIndex << endl;
            uint8_t errors2 = errors;
            bool valid = true;
            Pair<uint16_t, uint32_t> sa_info;
            uint32_t startPos;
            // mappability information is in reverse index order if we use the forward index
            if(std::is_same<TDir, Rev>::value){
                cout << "rev" << endl;
                sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + i];
                startPos = sa_info.i2 - needleLeftPos;
            }
            else{
                cout << "fwd" << endl;
                sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + i];
                //calculate correct starting position of the needle  on the forward index
                startPos = seqan::length(genome[sa_info.i1]) - sa_info.i2 - needleRightPos + 1;
            }
            cout << "StartPos " << startPos << endl;
            //search remaining blocks
            for(int j = blockIndex; j < s.pi.size(); ++j){
                int blockStart = (s.pi[j] - 1 == 0) ? 0 : s.chronBL[s.pi[j] - 2];
                int blockEnd = s.chronBL[s.pi[j] - 1];
                cout << "searching Parts:" << blockStart << " - " << blockEnd << "; ";
                cout << endl;
                // compare bases to needle
                if(std::is_same<TDir, Rev>::value){
                    if(needleRightPos - 1 > blockStart && needleRightPos - 1 < blockEnd){
                        cout << "changing Blockstart, should only happen (once) when we come from checkcurrentmappability!!" << endl;
                        blockStart = needleRightPos - 1;
                        cout << "searching Parts:" << blockStart << " - " << blockEnd << "; ";
                    }
                }else{
                    if(needleLeftPos > blockStart && needleLeftPos < blockEnd){
                        cout << "changing Blockend fwd" << endl;
                        blockEnd = needleLeftPos;
                        cout << "searching Parts:" << blockStart << " - " << blockEnd << "; ";
                    }
                }
                for(int k = blockStart; k <  blockEnd; ++k){
                    //                     if(needle[k] != it)
                    if(needle[k] != genome[sa_info.i1][startPos + k])
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
                uint32_t occ = startPos;
                hitsv.push_back(Pair<uint16_t,uint32_t>(sa_info.i1, occ));
                cout << "Hit occ: " << hitsv[hitsv.size() - 1] << endl;
                errorsv.push_back(errors2);
            }
        }
    }
    delegateDirect(hitsv, needle, errorsv);
}

//TODO load bitvectors inside a struct to make accessing the correct bitvector easier
template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks>
Pair<uint8_t, Pair<uint32_t, uint32_t>> get_bitvector_interval_inside(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                        vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,    
                                        OptimalSearch<nbrBlocks> const & s,
                                        uint8_t const blockIndex,
                                        bool const goToRight2) 
{
    cout << "Get Inside bitvector Interval" << endl;
    cout << "bitvector_interval blockIndex: " << (int)blockIndex << endl;
    cout << ((goToRight2) ? "reverse Index" : "forward index") << endl;
    Pair<uint32_t, uint32_t> dirrange = (goToRight2) ? range(iter.revIter) : range(iter.fwdIter);
    uint8_t needed_bitvector;
    uint8_t size = s.pi.size();
    uint8_t bitvsize = bitvectors.size();
    cout << "max: " << (int)s.max[blockIndex - 1] << endl;
    cout << "bitvsize: " << (int)bitvsize << "end bitvector size" << endl;
    if (goToRight2)
        needed_bitvector = bitvsize - s.max[blockIndex - 1];      
    else
        needed_bitvector = s.min[blockIndex - 1] - 1;


    uint32_t number_of_indeces = seqan::length(iter.fwdIter.index->sa) - bitvectors[needed_bitvector].first.size();
    cout << "selected bitvector inside: " << (int)needed_bitvector << endl;
//     cout << "Sa Range: " << dirrange << endl;  
    dirrange.i1 = dirrange.i1 - number_of_indeces;
    dirrange.i2 = dirrange.i2 - number_of_indeces;
//     cout << "Bit Range: " << dirrange << endl;  
    
    Pair<uint8_t, Pair<uint32_t, uint32_t>> brange(needed_bitvector, dirrange);
    return brange;
}

template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks, typename TDir>
Pair<uint8_t, Pair<uint32_t, uint32_t>> get_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                        vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,    
                                        OptimalSearch<nbrBlocks> const & s,
                                        uint8_t const blockIndex,
                                        TDir const & /**/) 
{
//     printv(s.pi);
    cout << "bitvector_interval blockIndex: " << (int)blockIndex << endl;
    Pair<uint32_t, uint32_t> dirrange = (std::is_same<TDir, Rev>::value) ? range(iter.fwdIter) : range(iter.revIter);
    uint8_t needed_bitvector;
    uint8_t bitvsize = bitvectors.size();
    
    uint8_t size = s.pi.size();
    uint8_t firstE = s.pi[0];
    if(bitvsize == 3){
        if(firstE == size)
            needed_bitvector = 2; //BV::LEFT;
        else if(firstE == 1)
            needed_bitvector = 0; //BV::RIGHT;
        else 
            needed_bitvector = 1; //BV::MIDDLE;
    }else{
        if (std::is_same<TDir, Rev>::value)
                needed_bitvector = s.min[blockIndex] - 1;//mymin(s.pi, blockIndex) - 1;
            else
                needed_bitvector = bitvsize - s.max[blockIndex];// + 1 - 1//mymax(s.pi, blockIndex) - 1; 
    }
    //TODO find out why this does not work
    //uint32_t number_of_indeces = countSequences(iter.fwdIter.index);
//     cout << "maxe: " << (int)s.max[blockIndex] << endl;
    uint32_t number_of_indeces = seqan::length(iter.fwdIter.index->sa) - bitvectors[needed_bitvector].first.size();
    cout << "selected bitvector: " << (int)needed_bitvector << endl;
//     cout << "Sa Range: " << dirrange << endl;  
    dirrange.i1 = dirrange.i1 - number_of_indeces;
    dirrange.i2 = dirrange.i2 - number_of_indeces;
//     cout << "Bit Range: " << dirrange << endl;  
//     cout << ((std::is_same<TDir, Rev>::value) ? "reverse case" : "forward case") << endl;
    Pair<uint8_t, Pair<uint32_t, uint32_t>> brange(needed_bitvector, dirrange);
    return brange;
}

template<typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks>
ReturnCode testUnidirectional(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                          vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex,
                          bool const goToRight2)
{
    
    // allowd flips per intervalSize
    float flipDensity = 1/static_cast<float>(2);
    
    // need bitinterval from inside the pattern to filter according to the mappability form
    //therefore i also need to acces the block before because of that block i got mappability of both sides
    auto bit_interval = get_bitvector_interval_inside(iter, bitvectors, s, blockIndex, goToRight2);
    sdsl::bit_vector & b2 = bitvectors[bit_interval.i1].first;
    
    //squas interval
    uint32_t startPos = bit_interval.i2.i1, endPos = bit_interval.i2.i2;
    
    for(uint32_t i = startPos; i < endPos; ++i){
        if(b2[i] != 0)
            break; 
        ++startPos;
    }
    for(uint32_t i = endPos - 1; i >= startPos; --i){
        if(b2[i] != 0)
            break;
        --endPos;
    }
    
    if(startPos > endPos){
        cout << "Error bit vector has only zeroes this should have been checked by checkinterval" << endl;
        cout << "Size: " << endPos - startPos << endl;
        exit(0);
    }
    // order of bits
    bool last = b2[startPos];
    uint32_t pos = startPos;
    uint32_t count = 0; 
    while(pos < endPos){
        if(b2[pos] != last){
            ++count;
            last = !last;
        }
        ++pos;
    }  
    float ivalSize = brange.i2.i2 - brange.i2.i1;
    // if next condition is true then brange will be modified!!!!
    // it will contain mappability of the bitvector from the other side of already searched needle
    
    // only interested in changes inside the supinterval (startPos - endPos)
    // if 0 got Cutoff that means more changes but is at the same time also good
    // therefore ignore them
    cout << "Test flipdensitey " << endl;
    cout << ivalSize * flipDensity - 1 << endl;
    cout << "count: " << count << endl;

    if(ivalSize * flipDensity - 1 > static_cast<float>(count)){
        cout << "Continue UNIDIRECTIONAL" << endl;
        brange.i1 = bit_interval.i1;
        cout << "New selected bitvector by inside function: " << (int)bit_interval.i1 << endl;
        brange.i2.i1 = startPos;
        brange.i2.i2 = endPos;
        return ReturnCode::UNIDIRECTIONAL;
    }
    return ReturnCode::MAPPABLE;
}

template<size_t nbrBlocks>
ReturnCode checkInterval(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex)
{
    float filter_threshold = 0.5; 
    int directSearch_Threshold = 2;
    sdsl::bit_vector & b = bitvectors[brange.i1].first;
    sdsl::rank_support_v<> & rb = bitvectors[brange.i1].second; 
    rb.set_vector(&b);
    
    uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
    if(ivalOne == 0)
        return ReturnCode::NOMAPPABILITY;

    if(ivalOne < (s.pi.size() - blockIndex - 1) * directSearch_Threshold){ //<4 TODO add additional constrains from chris? //incorparated blockIndex 
        return ReturnCode::DIRECTSEARCH;
    }
    
    if(ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
        return ReturnCode::COMPMAPPABLE;
    //equal or more than half zeroes     
    float ivalSize = brange.i2.i2 - brange.i2.i1;
    if(s.startUniDir <= blockIndex && ivalOne/ ivalSize <= filter_threshold){
        return ReturnCode::SUSPECTUNIDIRECTIONAL;
    }
    return ReturnCode::MAPPABLE;
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
    cout << "checkCurrentMappability:" << endl;
    cout << "NLP: " << (int)needleLeftPos << endl;
    cout << "NRP: " << (int)needleRightPos << endl;
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, TDir());
    ReturnCode rcode = checkInterval(bitvectors, bit_interval, s, blockIndex);
//   cout << "Return code: " << (int)rcode << endl;
    
    cout << "cPrintttt" << endl;
    if(std::is_same<TDir, Rev>::value)
        print_sa(iter, bitvectors, true);
    else
        print_sa(iter, bitvectors, false);
    printbit(bitvectors, bit_interval);
    cout << "cPrintend" << endl;
    
    if(rcode == ReturnCode::NOMAPPABILITY){
        cout << "checkCurrentMappability Stopp" << endl;
        return ReturnCode::FINISHED;        
    }
    if(rcode == ReturnCode::DIRECTSEARCH){
        cout << "checkCurrentMappability DirectSearch" << endl;
        //search directly in Genome
        directSearch(delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir());
        return ReturnCode::FINISHED;
    }
 /*   

    if(rcode == ReturnCode::COMPMAPPABLE){
        uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
        _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), HammingDistance());
        return ReturnCode::FINISHED;
    }
 */
 
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
                            uint32_t const current_needleLeftPos,
                            uint32_t const current_needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            bool const goToRight2,
                            TDir const & /**/)
{
    //TODO probably dont need needleLeftPos and needleRightPos sind there important border wasnt shifted
    bool finished = current_needleLeftPos == 0 && current_needleRightPos == length(needle) + 1;
    //check if we are done with the needle    
    if(!finished){
        Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
        if(goToRight2)
            bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, Rev());
        else
            bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, Fwd());
        
        ReturnCode rcode = checkInterval(bitvectors, bit_interval, s, blockIndex);
//      cout << "Return code: " << (int)rcode << endl;
        
        
        cout << "Printttt" << endl;
        if(goToRight2)
            print_sa(iter, bitvectors, true);
        else
            print_sa(iter, bitvectors, false);
        printbit(bitvectors, bit_interval);
        cout << "Printend" << endl;
        
        if(rcode == ReturnCode::NOMAPPABILITY){
            cout << "NOMAPPABILITYNOMAPPABILITYNOMAPPABILITYNOMAPPABILITYNOMAPPABILITYNOMAPPABILITY" << endl;
            return ReturnCode::FINISHED;        
        }
        if(rcode == ReturnCode::DIRECTSEARCH){
            cout << "DIRECTSEARCHDIRECTSEARCHDIRECTSEARCHDIRECTSEARCHDIRECTSEARCHDIRECTSEARCH" << endl;
            cout << "NLP: " << current_needleLeftPos << endl;
            cout << "NRP: " << current_needleRightPos << endl;
            if(std::is_same<TDir, Rev>::value)
                cout << "Iter start pos fwd" << iter.fwdIter.vDesc.range.i1 << endl;
            else
                cout << "Iter start pos rev" << iter.revIter.vDesc.range.i1 << endl;
            cout << "errors: " << (int)errors << endl;
            //search directly in Genome
            //TODO I only need left value for rev and right value for fwd so delte one input?
            // search in the next blocks only therefore need current error count
            if(goToRight2){
                directSearch(delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Rev());
//                 directSearchDummy(delegateDirect, iter.fwdIter, iter, needle, bitvectors, needleLeftPos , infixPosRight + 2, errors, s, blockIndex, bit_interval, Rev());
            }else{
                directSearch(delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Fwd());
//                 directSearchDummy(delegateDirect, iter.fwdIter, iter, needle, bitvectors, infixPosLeft , needleRightPos, errors, s, blockIndex, bit_interval, Fwd());
            }
            return ReturnCode::FINISHED;
        }
        
        //TODO check if we are in unidirection case in this case COMPMAPPABLE should no be used
        if(rcode == ReturnCode::COMPMAPPABLE){
            cout << "COMPMAPPABLECOMPMAPPABLECOMPMAPPABLECOMPMAPPABLECOMPMAPPABLE" << endl;
            cout << "NLP: " << current_needleLeftPos << endl;
            cout << "NRP: " << current_needleRightPos << endl;
            cout << "errors: " << (int)errors << endl;
            if (goToRight2)
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, Rev(), HammingDistance());
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, Fwd(), HammingDistance());
            }
            return ReturnCode::FINISHED;
        }
        
        if(rcode == ReturnCode::SUSPECTUNIDIRECTIONAL)
            rcode == testUnidirectional(iter, bitvectors, bit_interval, s, blockIndex, goToRight2);
        
        if(rcode == ReturnCode::UNIDIRECTIONAL){
            //range on iter was changed in function before
            if(goToRight2){
                filter_interval(delegate, delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Rev());
            }
            else
            {
                filter_interval(delegate, delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s,blockIndex, bit_interval, Fwd());
            }
            return ReturnCode::FINISHED;        
        }
    }else{
        cout << "Already finished" << endl;
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
    cout << "optimalSearchSchemeChildren" << endl;
    bool goToRight = std::is_same<TDir, Rev>::value;                                 
    if (goDown(iter, TDir()))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            cout << "TestLabelEdge:" << endl;
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()),
                                   needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1]);
            if (minErrorsLeftInBlock > 0 && charsLeft + delta < minErrorsLeftInBlock + 1u) // charsLeft - 1 < minErrorsLeftInBlock - delta
            {
                continue;
            }
            int32_t needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t needleRightPos2 = needleRightPos + goToRight;
            
            cout << "Succesful" << endl;
            cout << "NLP: " << needleLeftPos2 << endl;
            cout << "NRP: " << needleRightPos2 << endl;
            cout << "errors +delta: " << errors + delta << endl;
            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];                
                cout << "checkMappability Call from Children" << endl;
                ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, goToRight2, TDir());
                if(rcode == ReturnCode::FINISHED)
                    continue;
                
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
    bool goToRight2 = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] : s.pi[blockIndex] > s.pi[blockIndex - 1];
    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
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
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, goToRight2, TDir());
        //TDir() is in this case Rev() ...
        if(rcode == ReturnCode::FINISHED)
            return;            
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Rev());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Fwd());
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
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, goToRight2, TDir());
        //TDir() is in this case Fwd() ...
        if(rcode == ReturnCode::FINISHED)
            return;
        
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd());
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
    cout << "optimalSearchScheme" << endl;
    cout << "blockIndex: " << (int)blockIndex << endl;
    cout << "NLP: " << (int)needleLeftPos << endl;
    cout << "NRP: " << (int)needleRightPos << endl;
    cout << "Error: " << (int)errors << endl;
    cout << "SS: ";
    printv(s.pi);
    if(needleRightPos - needleLeftPos > 1)
        cout << "RangeSize: " << iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 << endl;
    
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

//     cout << "Step: " << needleRightPos - needleLeftPos - 1 << "    ss: "; printv(s.pi); cout << endl;
    // Done. (Last step)
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        for (auto occ : getOccurrences(iter)){
        cout << "BiSearch Hits: "<< (Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ)) << endl;
        }
        cout << "Finished bidirectional Search" << endl;
        //last input only matters for unidirectional search
        delegate(iter, needle, errors, false);
    }
    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
    }
    // Approximate search in current block.
    else
    {
      
    //TODO implement faster mod   ((temp3 ) == 0)
    //TODO check if we are not starting in a new block or near it!!
    cout << "test checkCurrentMappability" << endl;
    int step = (needleRightPos - needleLeftPos - 1);
    cout << "Step: " << step << endl;
    
    //Parameters
    int dfbe = 2; //distanceFromBlockEnd
    //0b11 == 4 parameter
    int pblocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
    if(((step & 0b11) == 0) && needleRightPos - needleLeftPos - 1 + dfbe < s.blocklength[blockIndex] && needleRightPos - needleLeftPos - 1 - dfbe > pblocklength){
        cout << "Checking" << endl;
        //TODO stop doing on loop to much (enter checkCurrentMappability a second time)
        ReturnCode rcode = checkCurrentMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
        if(rcode == ReturnCode::FINISHED)
            return;
    }
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
    cout << "New Needle:" << endl; //TODO revert this
    for(int i = 0; i < length(needle); ++i)
        cout << needle[i];
    cout << endl;
    _optimalSearchSchemeSetMapParams(scheme);
    print_search_scheme(scheme); //TODO revert this
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
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec>
void find(const int minErrors,
     const int maxErrors,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors)
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(delegate, delegateDirect, index, needles, bitvectors);
                break;
        case 2: find<0, 2>(delegate, delegateDirect, index, needles, bitvectors);
                break;
        case 3: find<0, 3>(delegate, delegateDirect, index, needles, bitvectors);
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

template <typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
void find(const int minErrors,
     const int maxErrors,
     TDelegate & delegate,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & /**/)
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(delegate, index, needles, TDistanceTag());
                break;
        case 2: find<0, 2>(delegate, index, needles, TDistanceTag());
                break;
        case 3: find<0, 3>(delegate, index, needles, TDistanceTag());
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}


#endif
