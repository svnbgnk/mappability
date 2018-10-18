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

namespace seqan{

template <typename TVector, typename TVSupport>
inline void getConsOnes(std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                uint32_t const intervalsize,
                std::vector<std::pair<uint32_t, uint32_t>> & consOnesOutput)
{
    TVector & b = bitvectors[inside_bit_interval.i1].first;
    uint32_t k = inside_bit_interval.i2.i1;
    uint32_t startOneInterval = inside_bit_interval.i2.i1;
    while(k < inside_bit_interval.i2.i2){
        uint32_t interval = 0;
        //TODO delete second condition it should end with 1
        while(b[k + interval] == 0 && (k + interval) < inside_bit_interval.i2.i2){
            ++interval;
        }
        if(interval >= intervalsize){
            consOnesOutput.push_back(make_pair(startOneInterval, k));
            startOneInterval = k + interval;
        }
        k += interval;
        interval = 0;
        ++k;
    }
    consOnesOutput.push_back(make_pair(startOneInterval, k));
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void filter_interval(TDelegate & delegate,
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
                     TDir const & ,
                     TDistanceTag const &)
{
    vector<pair<uint32_t, uint32_t>> consOnes;
    getConsOnes(bitvectors, inside_bit_interval, params.normal.intervalsize, consOnes);
    uint32_t noi = countSequences(*iter.fwdIter.index);

    for(uint32_t i = 0; i < consOnes.size(); ++i){
        if (std::is_same<TDir, Rev>::value){
            iter.revIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.revIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(delegate, delegateDirect, iter.revIter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Rev(), TDistanceTag());
        }
        else
        {
            iter.fwdIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.fwdIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(delegate, delegateDirect, iter.fwdIter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Fwd(), TDistanceTag());
        }
    }
}


//TODO remove blockindex ??
template <typename TDelegateD,
          typename TNeedle,
          size_t nbrBlocks>
inline void genomeSearch(TDelegateD & delegateDirect,
                  TNeedle const & needle,
                  uint8_t errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  auto const & genome,
                  Pair<uint16_t, uint32_t> const & sa_info,
                  std::array<uint32_t, nbrBlocks> & blockStarts,
                  std::array<uint32_t, nbrBlocks> & blockEnds)
{
    for(uint32_t j = 0; j < nbrBlocks - blockIndex; ++j){
        // compare bases to needle
        for(uint32_t k = blockStarts[j]; k <  blockEnds[j]; ++k){
            if(needle[k] != genome[sa_info.i1][sa_info.i2 + k]){
                ++errors;
            }
        }
        if(errors < s.l[blockIndex + j] || errors > s.u[blockIndex + j]){
            return;
        }
    }
    delegateDirect(sa_info, needle, errors);

}

template <typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void directSearch(TDelegateD & delegateDirect,
                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  vector<pair<TVector, TVSupport>> & bitvectors,
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                  TDir const & ,
                  TDistanceTag const &)
{
    //TODO implement myers bitvector
    auto const & genome = indexText(*iter.fwdIter.index);
    std::array<uint32_t, nbrBlocks> blockStarts;
    std::array<uint32_t, nbrBlocks> blockEnds;
    std::copy(std::begin(s.blockStarts) + blockIndex, std::end(s.blockStarts), std::begin(blockStarts));
    std::copy(std::begin(s.blockEnds) + blockIndex, std::end(s.blockEnds), std::begin(blockEnds));

    if(std::is_same<TDir, Rev>::value){
        //modify blockstart in case we are still inside a block
        if(needleRightPos - 1 > blockStarts[0] && needleRightPos - 1 < blockEnds[0])
            blockStarts[0] = needleRightPos - 1;

        for(uint32_t i = 0; i < brange.i2.i2 - brange.i2.i1; ++i){
            if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                // mappability information is in reverse index order if we use the forward index
                Pair<uint16_t, uint32_t> sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + i];
                uint32_t const chromlength = length(genome[sa_info.i1]);
                //Info make sure we dont DS search something going over the chromosom edge
                //check left chromosom boundry && check right chromosom boundry
                if(!(needleLeftPos <= sa_info.i2 && chromlength - 1 >= sa_info.i2 - needleLeftPos + length(needle) - 1))
                    continue;

                sa_info.i2 = sa_info.i2 - needleLeftPos;

                //search remaining blocks
                genomeSearch(delegateDirect, needle, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
            }
        }
    }
    else
    {
        //modify blockend in case we are still inside a block
        if(needleLeftPos > blockStarts[0] && needleLeftPos < blockEnds[0])
            blockEnds[0] = needleLeftPos;

        for(uint32_t i = 0; i < brange.i2.i2 - brange.i2.i1; ++i){
            if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                Pair<uint16_t, uint32_t> sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + i];
                uint32_t const chromlength = length(genome[sa_info.i1]);
                //check left chromosom boundry && check right chromosom boundry
                if(!(chromlength - 1 >= sa_info.i2 + needleRightPos - 1 && sa_info.i2 + needleRightPos - 1 >= length(needle) + 1))
                    continue;
                //calculate correct starting position of the needle  on the forward index
                sa_info.i2 = chromlength - sa_info.i2 - needleRightPos + 1;

                //search remaining blocks
                genomeSearch(delegateDirect, needle, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
            }
        }
    }

}

//TODO load bitvectors inside a struct to make accessing the correct bitvector easier
template <typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks>
inline void get_bitvector_interval_inside(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                              vector<pair<TVector, TVSupport>> & bitvectors,
                              OptimalSearch<nbrBlocks> const & s,
                              uint8_t const blockIndex,
                              Pair<uint8_t, Pair<uint32_t, uint32_t>> & brangeOutput,
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

    uint32_t nseq = countSequences(*iter.fwdIter.index);
    dirrange.i1 = dirrange.i1 - nseq;
    dirrange.i2 = dirrange.i2 - nseq;

    brangeOutput.i1 = needed_bitvector;
    brangeOutput.i2 = dirrange;
}


template <typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void get_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                       vector<pair<TVector, TVSupport>> & bitvectors,
                       OptimalSearch<nbrBlocks> const & s,
                       uint8_t const blockIndex,
                       Pair<uint8_t, Pair<uint32_t, uint32_t>> & brangeOutput,
                       TDir const & )
{
    Pair<uint32_t, uint32_t> dirrange = (std::is_same<TDir, Rev>::value) ? range(iter.fwdIter) : range(iter.revIter);
    uint8_t needed_bitvector;
    if (std::is_same<TDir, Rev>::value)
        needed_bitvector = s.min[blockIndex] - 1;
    else
        needed_bitvector = bitvectors.size() - s.max[blockIndex];// + 1 - 1

    uint32_t nseq = countSequences(*iter.fwdIter.index);
    dirrange.i1 = dirrange.i1 - nseq;
    dirrange.i2 = dirrange.i2 - nseq;

    brangeOutput.i1 = needed_bitvector;
    brangeOutput.i2 = dirrange;
}

template<typename TText, typename TIndex, typename TIndexSpec,
         typename TVector, typename TVSupport,
         size_t nbrBlocks>
inline bool testUnidirectionalFilter(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                          vector<pair<TVector, TVSupport>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex,
                          bool const goToRight2)
{
    // need bitinterval from inside the pattern to filter according to the mappability form
    //therefore i also need to acces the block before because of that block i got mappability of both sides
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
    get_bitvector_interval_inside(iter, bitvectors, s, blockIndex, bit_interval, goToRight2);
    TVector & b2 = bitvectors[bit_interval.i1].first;

    //squash interval
    uint32_t startPos = bit_interval.i2.i1, endPos = bit_interval.i2.i2;
    uint32_t startPos2 = startPos;
    uint32_t endPos2 = endPos;

    while(b2[startPos] == 0 && startPos < endPos)
        ++startPos;

    while(b2[endPos - 1] == 0 && endPos > startPos)
        --endPos;

    if(startPos > endPos){
        cout << "Error bit vector has only zeroes this should have been checked by checkinterval" << endl;
        exit(0);
    }

    float ivalSize = brange.i2.i2 - brange.i2.i1;
    uint32_t count = 0;

    if(params.normal.testflipdensity){
        // order of bits
        bool last = b2[startPos];
        uint32_t pos = startPos;
        while(pos < endPos){
            if(b2[pos] != last){
                ++count;
                last = !last;
            }
            ++pos;
        }
    }
    // if next condition is true then brange will be modified!!!!
    // it will contain mappability of the bitvector from the other side of already searched needle

    // only interested in changes inside the supinterval (startPos - endPos)
    // if 0 got Cutoff that means more changes but is also good at the same time
    // therefore ignore them
    // allowd flips per intervalSize
    if(!params.normal.testflipdensity || ivalSize * params.normal.invflipdensity - 1 > static_cast<float>(count)){
        brange.i1 = bit_interval.i1;
        brange.i2.i1 = startPos;
        brange.i2.i2 = endPos;
        return true;
    }
    return false;
}

template<typename TVector, typename TVSupport,
         size_t nbrBlocks>
inline ReturnCode checkInterval(vector<pair<TVector, TVSupport>> & bitvectors,
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

    if(params.normal.directsearch && ivalOne < (s.pi.size() - blockIndex - 1 + params.normal.directsearchblockoffset) * params.normal.directsearch_th)
        return ReturnCode::DIRECTSEARCH;

    if(params.normal.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
        return ReturnCode::COMPMAPPABLE;

    //equal or more than half zeroes
    if(params.normal.suspectunidirectional && s.startUniDir <= blockIndex && ivalOne/ ivalSize <= params.normal.filter_th)
        return ReturnCode::SUSPECTUNIDIRECTIONAL;
    else
        return ReturnCode::MAPPABLE;
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode checkCurrentMappability(TDelegate & delegate,
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
                            TDir const & ,
                            TDistanceTag const &)
{
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
    get_bitvector_interval(iter, bitvectors, s, blockIndex, bit_interval, TDir());

    ReturnCode rcode = checkInterval(bitvectors, bit_interval, s, blockIndex);

    switch(rcode){
        case ReturnCode::NOMAPPABILITY:
            return ReturnCode::FINISHED;

        case ReturnCode::DIRECTSEARCH:
        {
            directSearch(delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::COMPMAPPABLE:
        {
            _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), HammingDistance());
            return ReturnCode::FINISHED;
        }

        default:
            return ReturnCode::MAPPABLE;
    }
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode checkMappability(TDelegate & delegate,
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
                            TDir const & ,
                            TDistanceTag const &)
{
    bool finished = current_needleLeftPos == 0 && current_needleRightPos == length(needle) + 1;
    //check if we are done with the needle
    if(!finished){
        Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
        if(goToRight2)
            get_bitvector_interval(iter, bitvectors, s, blockIndex, bit_interval, Rev());
        else
            get_bitvector_interval(iter, bitvectors, s, blockIndex, bit_interval, Fwd());

        ReturnCode rcode = checkInterval(bitvectors, bit_interval, s, blockIndex);
        //TODO switch
         switch(rcode)
         {
            case ReturnCode::NOMAPPABILITY:
                return ReturnCode::FINISHED;

            case ReturnCode::DIRECTSEARCH:
            {
                //search directly in Genome
                // search in the next blocks only therefore need current error count
                //TODO direction plays no role in this anymore?
                if(goToRight2)
                {
                    directSearch(delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Rev(), TDistanceTag());
                }
                else
                {
                    directSearch(delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Fwd(), TDistanceTag());
                }
                return ReturnCode::FINISHED;
            }
            case ReturnCode::COMPMAPPABLE:
            {
                if (goToRight2)
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, Rev(), HammingDistance());
                else
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, Fwd(), HammingDistance());
                return ReturnCode::FINISHED;
            }

            case ReturnCode::SUSPECTUNIDIRECTIONAL:
            {
                //test unidirectional changes iter range if true
                if(testUnidirectionalFilter(iter, bitvectors, bit_interval, s, blockIndex, goToRight2)){
                    //range on iter was changed in function before
                    if(goToRight2)
                        filter_interval(delegate, delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, Rev(), TDistanceTag());
                    else
                        filter_interval(delegate, delegateDirect, iter, needle, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s,blockIndex, bit_interval, Fwd(), TDistanceTag());
                    return ReturnCode::FINISHED;
                }
            }
         }
    }
    return ReturnCode::MAPPABLE;
}

template <typename TDelegate,
          typename TDelegateDirect,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletion(TDelegate & delegate,
                                         TDelegateDirect & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         vector<pair<TVector, TVSupport>> & bitvectors,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         TDir const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    if (minErrorsLeftInBlock == 0)
    {
        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, Rev(),
                                 EditDistance());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, Fwd(),
                                 EditDistance());
        }
    }

    if (maxErrorsLeftInBlock > 0 && goDown(iter, TDir()))
    {
        do
        {
            _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors + 1, s,
                                         blockIndex, TDir());
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
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
                                         TDir const & ,
                                         TDistanceTag const &)
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
                if (std::is_same<TDistanceTag, EditDistance>::value)
                {
                    //TODO it does not make sense to check mappability for deletions, right?
                    _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2,
                                                 errors + delta, s, blockIndex, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

                    ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, goToRight2, TDir(), TDistanceTag());
                    if(rcode == ReturnCode::FINISHED)
                        continue;

                    if (goToRight2)
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev(), TDistanceTag());
                    else
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd(), TDistanceTag());
                }
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, TDir(), TDistanceTag());
            }

            // Deletion
            if (std::is_same<TDistanceTag, EditDistance>::value)
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex,
                                     TDir(), TDistanceTag());
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
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
                                      TDir const &,
                                      TDistanceTag const &)
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

        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, goToRight2, TDir(), TDistanceTag());
        //TDir() is in this case Rev() ...
        if(rcode == ReturnCode::FINISHED)
            return;
        if (goToRight2)
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Fwd(), TDistanceTag());
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

        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, goToRight2, TDir(), TDistanceTag());
        if(rcode == ReturnCode::FINISHED)
            return;

        if (goToRight2)
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd(), TDistanceTag());
    }
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 vector<pair<TVector, TVSupport>> & bitvectors,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 TDir const & ,
                                 TDistanceTag const &)
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
        _optimalSearchSchemeExact(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
    }
    // Approximate search in current block.
    else
    {
        //TODO put this into _optimalSearchSchemeChildren right at the end
        // Insertion
        if (std::is_same<TDistanceTag, EditDistance>::value)
        {
            bool const goToRight = std::is_same<TDir, Rev>::value;
            int32_t const needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t const needleRightPos2 = needleRightPos + goToRight;

            //if we are at the end of block we need to add possible deletions because _optimalSearchScheme does not check it
            // though checking 1 deletion makes no sense since that would be equivalent to a match
            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, TDir());
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, TDir(), TDistanceTag());
            }
        }

        uint32_t pblocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
        uint32_t step = (needleRightPos - needleLeftPos - 1);
        if(((step & (params.normal.step - 1)) == 0) &&
            needleRightPos - needleLeftPos - 1 + params.normal.distancetoblockend < s.blocklength[blockIndex] && static_cast<int>(needleRightPos - needleLeftPos - 1) - params.normal.distancetoblockend > pblocklength)
        {
            ReturnCode rcode = checkCurrentMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
            if(rcode == ReturnCode::FINISHED)
                return;
        }
        _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 vector<pair<TVector, TVSupport>> & bitvectors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 TDistanceTag const &)
{
    bool initialDirection = s.pi[1] > s.pi[0];
    if(!params.startUnidirectional || s.startUniDir > 0){
        if(initialDirection)
            _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Fwd(), TDistanceTag());
    }/*
    else
    {
        if(initialDirection)
            _uniOptimalSearchScheme(delegate, delegateDirect, it.revIter, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Rev(), HammingDistance());
        else
            _uniOptimalSearchScheme(delegate, delegateDirect, it.fwdIter, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Fwd(), HammingDistance());
    }*/
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 OptimalSearch<nbrBlocks> const & s,
                                 TDistanceTag const &)
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
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 vector<pair<TVector, TVSupport>> & bitvectors,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TDistanceTag const &)
{
    for (auto & s : ss)
        _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, s, TDistanceTag());
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TDistanceTag const &)
{
    for (auto & s : ss)
        _optimalSearchScheme(delegate, delegateDirect, it, needle, s, TDistanceTag());
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     vector<pair<TVector, TVSupport>> & bitvectors,
     std::array<OptimalSearch<nbrBlocks>, N> const & ss,
     TDistanceTag const &)
{
    //TODO there has to be a better way
    auto scheme = ss;/*OptimalSearchSchemes<minErrors, maxErrors>::VALUE;*/
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(delegate, delegateDirect, it, needle, bitvectors, scheme, TDistanceTag());
}

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     TDistanceTag const & )
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(delegate, delegateDirect, it, needle, scheme, TDistanceTag());
}


template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TVector, typename TVSupport,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<pair<TVector, TVSupport>> & bitvectors,
     TDistanceTag const & )
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    calcConstParameters(scheme);
    checkTime time;
    int k = 0;
    uint32_t lastcount = 0;
    while(k < length(needles))
    {
        if(!(params.clocking && time.stopnow(params.terminateDuration))){
            find(delegate, delegateDirect, index, needles[k], bitvectors, scheme, TDistanceTag());
        }else{
            params.wasStopped = true;
        }
        ++k;

        uint32_t currentcount = hits.size() + dhits.size() - lastcount;
        readOccCount.push_back(currentcount);
        lastcount += currentcount;
    }
}

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & )
{
    int k = 0;
    uint32_t lastcount = 0;
    while(k < length(needles))
    {
        find<minErrors, maxErrors>(delegate, delegateDirect, index, needles[k], TDistanceTag());
        ++k;
        uint32_t currentcount = hitsDe.size() + dhitsDe.size() - lastcount;
        readOccCountDeT.push_back(currentcount);
        lastcount += currentcount;

    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TVector, typename TVSupport,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 vector<pair<TVector, TVSupport>> & bitvectors, // cant be const since TVSupport.set_vector(&TVector)
                 TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(delegate, delegateDirect, index, needles, bitvectors, TDistanceTag());
                break;
        case 2: find<0, 2>(delegate, delegateDirect, index, needles, bitvectors, TDistanceTag());
                break;
        case 3: find<0, 3>(delegate, delegateDirect, index, needles, bitvectors, TDistanceTag());
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}
/*
template <typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
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
        case 4: find<0, 4>(delegate, index, needles, TDistanceTag());
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
     const int maxErrors,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(delegate, delegateDirect, index, needles, HammingDistance());
                break;
        case 2: find<0, 2>(delegate, delegateDirect, index, needles, HammingDistance());
                break;
        case 3: find<0, 3>(delegate, delegateDirect, index, needles, HammingDistance());
                break;
        case 4: find<0, 4>(delegate, delegateDirect, index, needles, HammingDistance());
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}*/

}



#endif
