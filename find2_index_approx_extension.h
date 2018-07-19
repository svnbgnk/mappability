#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include "global.h"
#include "auxiliary.h"
#include "common_auxiliary.h"
#include "find2_index_approx_unidirectional.h"
#include "find2_index_approx_compmappable.h"
#include "find2_index_approx_start_unidirectional.h"

void testglobal(){
    std::cout << "Lets test global" << std::endl;
    std::cout << global << std::endl;
    params.print();
}


//TODO add //squash interval as function

namespace seqan{
    
template <typename TVector, typename TVSupport>
vector<pair<uint32_t, uint32_t>> getConsOnes(vector<pair<TVector, TVSupport>> & bitvectors, // TODO: const
                                             Pair<uint8_t, Pair<uint32_t, uint32_t>> inside_bit_interval, // TODO: const &
                                             int const intervalsize)
{
    TVector & b = bitvectors[inside_bit_interval.i1].first;
    vector<pair<uint32_t, uint32_t>> consOnes;
    uint32_t k = inside_bit_interval.i2.i1;
    uint32_t startOneInterval = inside_bit_interval.i2.i1;
    while(k < inside_bit_interval.i2.i2){
        uint32_t interval = 0;
        //TODO delete second condition it should end with 1
        while(b[k + interval] == 0 && (k + interval) < inside_bit_interval.i2.i2){
            ++interval;
        }
        if(interval >= intervalsize){
            consOnes.push_back(make_pair(startOneInterval, k));
            startOneInterval = k + interval;
        }
        k += interval;
        interval = 0;
        ++k;
    }
    consOnes.push_back(make_pair(startOneInterval, k));
    return(consOnes);
}

    
// template <size_t N>
// struct OptimalSearch
// {
//     std::array<uint8_t, N> pi; // order of the blocks. permutation of [1..n]
//     std::array<uint8_t, N> l; // minimum number of errors at the end of the corresponding block
//     std::array<uint8_t, N> u; // maximum number of errors at the end of the corresponding block
// 
//     std::array<uint32_t, N> blocklength; // cumulated values / prefix sums
//     std::array<uint32_t, N> chronBL;
//     std::array<uint32_t, N> revChronBL;
//     std::array<uint8_t, N> min;
//     std::array<uint8_t, N> max;
//     uint32_t startPos; //wrong position so i still get 0 from initialization
//     uint8_t startUniDir; 
// };  
    

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
void filter_interval(TDelegate & delegate,
                     TDelegateD & delegateDirect,
                     Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                     TNeedle const & needle,
                     vector<pair<TVector, TVSupport>> & bitvectors,
                     uint32_t const needleLeftPos,
                     uint32_t const needleRightPos,
                     uint8_t const errors,
                     OptimalSearch<nbrBlocks> const & s,
                     uint8_t const blockIndex,
                     Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                     TDir const & )
{  
    vector<pair<uint32_t, uint32_t>> consOnes = getConsOnes(bitvectors, inside_bit_interval, params.normal.intervalsize);
//     consOnes.push_back(make_pair(inside_bit_interval.i2.i1, inside_bit_interval.i2.i2));
    uint32_t noi = seqan::length(iter.fwdIter.index->sa) - bitvectors[0].first.size(); // number_of_indeces
    
    for(int i = 0; i < consOnes.size(); ++i){
        if (std::is_same<TDir, Rev>::value){
            //TODO call DirectSearch here if the interval is to small also use block Index ....?
            iter.revIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.revIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(delegate, delegateDirect, iter.revIter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Rev());
        }
        else
        {
            iter.fwdIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.fwdIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(delegate, delegateDirect, iter.fwdIter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Fwd());
        }
    } 
}

// remove hitsvOutput return vector errors
template <typename TDelegateD,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
void genomeSearch(TDelegateD & delegateDirect,
                  TNeedle const & needle,
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  TDir const & ,
                  auto const & genome,
                  Pair<uint16_t, uint32_t> const & sa_info)
{
    bool valid = true;
    for(int j = blockIndex; j < s.pi.size(); ++j){
        int blockStart = (s.pi[j] - 1 == 0) ? 0 : s.chronBL[s.pi[j] - 2];
        int blockEnd = s.chronBL[s.pi[j] - 1];
        // compare bases to needle
        if(std::is_same<TDir, Rev>::value){
            if(needleRightPos - 1 > blockStart && needleRightPos - 1 < blockEnd)
                blockStart = needleRightPos - 1;
        }
        else
        {
            if(needleLeftPos > blockStart && needleLeftPos < blockEnd)
                blockEnd = needleLeftPos;
        }
        for(int k = blockStart; k <  blockEnd; ++k){
            if(needle[k] != genome[sa_info.i1][sa_info.i2 + k])
                ++errors;
        }
        if(errors < s.l[j] || errors > s.u[j]){
            valid = false;
            break;
        }
    }
    if(valid){
        delegateDirect(sa_info, needle, errors);
    }
}

template <typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
void directSearch(TDelegateD & delegateDirect,
                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  vector<pair<TVector, TVSupport>> & bitvectors, 
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                  TDir const & )
{
//     hitsv.reserve(countOccurrences(iter.fwdIter)); //TODO fix this
    auto const & genome = indexText(*iter.fwdIter.index);
    for(int i = 0; i < brange.i2.i2 - brange.i2.i1; ++i){
        if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
            Pair<uint16_t, uint32_t> sa_info;
            uint32_t startPos;
            // mappability information is in reverse index order if we use the forward index
            if(std::is_same<TDir, Rev>::value){
                sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + i];
                sa_info.i2 = sa_info.i2 - needleLeftPos;
            }
            else
            {
                sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + i];
                //calculate correct starting position of the needle  on the forward index
                sa_info.i2 = seqan::length(genome[sa_info.i1]) - sa_info.i2 - needleRightPos + 1;
            }
            //search remaining blocks
            genomeSearch(delegateDirect, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), genome, sa_info);
        }
    }
//     delegateDirect(hitsv, needle, errorsv);
}

//TODO load bitvectors inside a struct to make accessing the correct bitvector easier
template <typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks>
Pair<uint8_t, Pair<uint32_t, uint32_t>> get_bitvector_interval_inside(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                        vector<pair<TVector, TVSupport>> & bitvectors,    
                                        OptimalSearch<nbrBlocks> const & s,
                                        uint8_t const blockIndex,
                                        bool const goToRight2) 
{
    Pair<uint32_t, uint32_t> dirrange = (goToRight2) ? range(iter.revIter) : range(iter.fwdIter);
    uint8_t needed_bitvector;
    uint8_t size = s.pi.size();
    uint8_t bitvsize = bitvectors.size();
    if (goToRight2)
        needed_bitvector = bitvsize - s.max[blockIndex - 1];      
    else
        needed_bitvector = s.min[blockIndex - 1] - 1;

    //TODO use countSequences if it works
    uint32_t number_of_indeces = seqan::length(iter.fwdIter.index->sa) - bitvectors[needed_bitvector].first.size();
    dirrange.i1 = dirrange.i1 - number_of_indeces;
    dirrange.i2 = dirrange.i2 - number_of_indeces;
    
    Pair<uint8_t, Pair<uint32_t, uint32_t>> brange(needed_bitvector, dirrange);
    return brange;
}

template <typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
Pair<uint8_t, Pair<uint32_t, uint32_t>> get_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                                               vector<pair<TVector, TVSupport>> & bitvectors,    
                                                               OptimalSearch<nbrBlocks> const & s,
                                                               uint8_t const blockIndex,
                                                               TDir const & ) 
{
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
    uint32_t number_of_indeces = seqan::length(iter.fwdIter.index->sa) - bitvectors[needed_bitvector].first.size(); 
    dirrange.i1 = dirrange.i1 - number_of_indeces;
    dirrange.i2 = dirrange.i2 - number_of_indeces;
    Pair<uint8_t, Pair<uint32_t, uint32_t>> brange(needed_bitvector, dirrange);
    return brange;
}

template<typename TText, typename TIndex, typename TIndexSpec,
         typename TVector, typename TVSupport,
         size_t nbrBlocks>
bool testUnidirectionalFilter(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                          vector<pair<TVector, TVSupport>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex,
                          bool const goToRight2)
{
    // need bitinterval from inside the pattern to filter according to the mappability form
    //therefore i also need to acces the block before because of that block i got mappability of both sides
    auto bit_interval = get_bitvector_interval_inside(iter, bitvectors, s, blockIndex, goToRight2);
    TVector & b2 = bitvectors[bit_interval.i1].first;
    
    //squash interval
    uint32_t startPos = bit_interval.i2.i1, endPos = bit_interval.i2.i2;
    
    for(uint32_t i = startPos; i < endPos; ++i){
        if(b2[i] != 0)
            break; 
        ++startPos;
    }
    // TODO: schleifen schöner machen. ohne i
    for(uint32_t i = endPos - 1; i >= startPos; --i){
        if(b2[i] != 0)
            break;
        --endPos;
    }
    
    if(startPos > endPos){
        cout << "Error bit vector has only zeroes this should have been checked by checkinterval" << endl;
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
    // if 0 got Cutoff that means more changes but is also good at the same time
    // therefore ignore them
    // allowd flips per intervalSize
    if(ivalSize * params.normal.invflipdensity - 1 > static_cast<float>(count)){
        brange.i1 = bit_interval.i1;
        brange.i2.i1 = startPos;
        brange.i2.i2 = endPos;
        return true;
    }
    return false;
}

template<typename TVector, typename TVSupport,
         size_t nbrBlocks>
ReturnCode checkInterval(vector<pair<TVector, TVSupport>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex)
{
    TVector & b = bitvectors[brange.i1].first;
    TVSupport & rb = bitvectors[brange.i1].second; 
    rb.set_vector(&b);
    
    uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
    float ivalSize = brange.i2.i2 - brange.i2.i1;
    
    if(params.normal.nomappability && ivalOne == 0)
        return ReturnCode::NOMAPPABILITY;
    // TODO: else raus
    else if(params.normal.directsearch && ivalOne < (s.pi.size() - blockIndex - 1 + params.normal.directsearchblockoffset) * params.normal.directsearch_th)
        return ReturnCode::DIRECTSEARCH;
    
    else if(params.normal.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
        return ReturnCode::COMPMAPPABLE;
    
    //equal or more than half zeroes
    else if(params.normal.suspectunidirectional && s.startUniDir <= blockIndex && ivalOne/ ivalSize <= params.normal.filter_th)
        return ReturnCode::SUSPECTUNIDIRECTIONAL;
        
    else
        return ReturnCode::MAPPABLE;
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
ReturnCode checkCurrentMappability(TDelegate & delegate,
                            TDelegateD & delegateDirect,
                            Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                            TNeedle const & needle,
                            vector<pair<TVector, TVSupport>> & bitvectors,    
                            uint32_t const needleLeftPos,
                            uint32_t const needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            uint8_t const minErrorsLeftInBlock,
                            TDir const & )
{
    //TODO add else ifs
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, TDir());
    ReturnCode rcode = checkInterval(bitvectors, bit_interval, s, blockIndex);
    
    // TODO:
    // switch(rcode) {
    //   case ReturnCode::NOMAPPABILITY
    // }
    
    if(rcode == ReturnCode::NOMAPPABILITY)
        return ReturnCode::FINISHED;        
    //search directly in Genome
    if(rcode == ReturnCode::DIRECTSEARCH){
        directSearch(delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir());
        return ReturnCode::FINISHED;
    }
    if(rcode == ReturnCode::COMPMAPPABLE){
        _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), HammingDistance());
        return ReturnCode::FINISHED;
    }
    return ReturnCode::MAPPABLE;
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
ReturnCode checkMappability(TDelegate & delegate,
                            TDelegateD & delegateDirect,
                            Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                            TNeedle const & needle,
                            vector<pair<TVector, TVSupport>> & bitvectors,    
                            uint32_t const current_needleLeftPos,
                            uint32_t const current_needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            bool const goToRight2,
                            TDir const & )
{
    bool finished = current_needleLeftPos == 0 && current_needleRightPos == length(needle) + 1;
    //check if we are done with the needle    
    if(!finished){
        Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
        if(goToRight2)
            bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, Rev());
        else
            bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, Fwd());
        
        ReturnCode rcode = checkInterval(bitvectors, bit_interval, s, blockIndex);
        
        if(rcode == ReturnCode::NOMAPPABILITY)
            return ReturnCode::FINISHED;   
        if(rcode == ReturnCode::DIRECTSEARCH)
        {
            //search directly in Genome
            // search in the next blocks only therefore need current error count
            if(goToRight2)
            {
                directSearch(delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Rev());
            }
            else
            {
                directSearch(delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Fwd());
            }
            return ReturnCode::FINISHED;
        }
        if(rcode == ReturnCode::COMPMAPPABLE)
        {
            if (goToRight2)
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, Rev(), HammingDistance());
            else
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, Fwd(), HammingDistance());
            return ReturnCode::FINISHED;
        }
        
        if(rcode == ReturnCode::SUSPECTUNIDIRECTIONAL)
        {
            //test unidirectional changes iter range if true
            if(testUnidirectionalFilter(iter, bitvectors, bit_interval, s, blockIndex, goToRight2)){
            //range on iter was changed in function before
                if(goToRight2)
                    filter_interval(delegate, delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Rev());
                else
                    filter_interval(delegate, delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s,blockIndex,    bit_interval, Fwd());
                return ReturnCode::FINISHED;        
            }
        }
    }
    return ReturnCode::MAPPABLE;
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeChildren(TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         vector<pair<TVector, TVSupport>> & bitvectors,   
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const & )                                 
{
    bool goToRight = std::is_same<TDir, Rev>::value;                                 
    if (goDown(iter, TDir()))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()), needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1]);
            if (minErrorsLeftInBlock > 0 && charsLeft + delta < minErrorsLeftInBlock + 1u) 
            {
                continue;
            }
            int32_t needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t needleRightPos2 = needleRightPos + goToRight;
            //finished Block
            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];                
                
                ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, goToRight2, TDir());
                if(rcode == ReturnCode::FINISHED)
                    continue;
                
                if (goToRight2)
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev());
                else
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd());
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
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeExact(TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      vector<pair<TVector, TVSupport>> & bitvectors,    
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const & )
{
    // not in last block and next Block is larger then current block
    bool goToRight2 = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] : s.pi[blockIndex] > s.pi[blockIndex - 1];
    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it reverse
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;
        
        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
            return;
        
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, goToRight2, TDir());
        //TDir() is in this case Rev() ...
        if(rcode == ReturnCode::FINISHED)
            return;            
        if (goToRight2)
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Rev());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Fwd());
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

        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, goToRight2, TDir());
        if(rcode == ReturnCode::FINISHED)
            return;
        
        if (goToRight2)
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd());
    }
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 vector<pair<TVector, TVSupport>> & bitvectors, // TODO: const   
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,                             
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 TDir const & )
{    
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    // Done. (Last step)
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        //last input only matters for unidirectional searches (has to be false in this case)
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
    //int dfbe = 2; //distanceFromBlockEnd
    int pblocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
    int step = (needleRightPos - needleLeftPos - 1);
    if(((step & params.normal.step) == 0) && 
        needleRightPos - needleLeftPos - 1 + params.normal.distancetoblockend < s.blocklength[blockIndex] && static_cast<int>(needleRightPos - needleLeftPos - 1) - params.normal.distancetoblockend > pblocklength)
    {
        ReturnCode rcode = checkCurrentMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir());
        if(rcode == ReturnCode::FINISHED)
            return;
    }
    _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir());
    }
}  

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 vector<pair<TVector, TVSupport>> & bitvectors,
                                 OptimalSearch<nbrBlocks> const & s)
{
    bool initialDirection = s.pi[1] > s.pi[0];
    if(!params.startUnidirectional || s.startUniDir > 0){
        if(initialDirection)
            _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Rev());
        else
            _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Fwd());
    }else{
        if(initialDirection)
            _uniOptimalSearchScheme(delegate, delegateDirect, it.revIter, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Rev());
        else
            _uniOptimalSearchScheme(delegate, delegateDirect, it.fwdIter, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Fwd());
    }
}

//TODO use TDistanceTag in fucntions before
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 OptimalSearch<nbrBlocks> const & s)
{
    bool initialDirection = s.pi[1] > s.pi[0];
    if(initialDirection)
        _optimalSearchScheme(delegate, delegateDirect, it, needle, s.startPos, s.startPos + 1, 0, s, 0, Rev(), HammingDistance());
    else
        _optimalSearchScheme(delegate, delegateDirect, it, needle, s.startPos, s.startPos + 1, 0, s, 0, Fwd(), HammingDistance());

}

  
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks, size_t N>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 vector<pair<TVector, TVSupport>> & bitvectors,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss)
{
    for (auto & s : ss)
        _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s);
}  

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks, size_t N>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss)
{
    for (auto & s : ss)
        _optimalSearchScheme(delegate, delegateDirect, it, needle, s);
}


template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec,
          typename TVector, typename TVSupport>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     vector<pair<TVector, TVSupport>> & bitvectors)
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeSetMapParams(scheme);
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, scheme);
}

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle)
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeSetMapParams(scheme);
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(delegate, delegateDirect, it, needle, scheme);
}


template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TVector, typename TVSupport,
          typename TParallelTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<pair<TVector, TVSupport>> & bitvectors,
     TParallelTag const & )
{
    typedef typename Iterator<StringSet<TNeedle, TStringSetSpec> const, Rooted>::Type TNeedleIt;
    typedef typename Reference<TNeedleIt>::Type                                       TNeedleRef;
    checkTime time;
    iterate(needles, [&](TNeedleIt const & needleIt)
    {
        /*
        if(params.clocking && time.stopnow(params.terminateDuration)){
            cout << "Trying to break/return" << endl;
            for(int i = 0; i < length(needles); ++i){
                TNeedleRef needle = value(needleIt);
//                 cout << "s";
            }
            //return does not work;
        }*/ 
        
        TNeedleRef needle = value(needleIt);
        if(!(params.clocking && time.stopnow(params.terminateDuration))){
            find<minErrors, maxErrors>(delegate, delegateDirect, index, needle, bitvectors);
        }else{
            params.wasStopped = true;
        }
        
    },
    Rooted(), TParallelTag());
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
     TParallelTag const & )
{
    typedef typename Iterator<StringSet<TNeedle, TStringSetSpec> const, Rooted>::Type TNeedleIt;
    typedef typename Reference<TNeedleIt>::Type                                       TNeedleRef;
    iterate(needles, [&](TNeedleIt const & needleIt)
    {
        TNeedleRef needle = value(needleIt); 
        find<minErrors, maxErrors>(delegate, delegateDirect, index, needle);
    },
    Rooted(), TParallelTag());
}  
   

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TVector, typename TVSupport>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<pair<TVector, TVSupport>> & bitvectors)
{
    find<minErrors, maxErrors>(delegate, delegateDirect, index, needles, bitvectors, Serial());
}


template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles)
{
    find<minErrors, maxErrors>(delegate, delegateDirect, index, needles, Serial());
}



template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TVector, typename TVSupport>
void find(const int minErrors,
     const int maxErrors,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index, // TODO: const?
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<pair<TVector, TVSupport>> & bitvectors) // TODO: const?
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
     TDistanceTag const & )
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

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec>
void find(const int minErrors,
     const int maxErrors,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles)
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(delegate, delegateDirect, index, needles);
                break;
        case 2: find<0, 2>(delegate, delegateDirect, index, needles);
                break;
        case 3: find<0, 3>(delegate, delegateDirect, index, needles);
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

}



#endif
