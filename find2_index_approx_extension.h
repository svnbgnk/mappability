#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include "common.h"
#include "find2_index_approx_unidirectional.h"

// class OSSContext

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
            consOnesOutput.push_back(std::make_pair(startOneInterval, k));
            startOneInterval = k + interval;
        }
        k += interval;
        interval = 0;
        ++k;
    }
    consOnesOutput.push_back(std::make_pair(startOneInterval, k));
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void filter_interval(TContex & ossContext,
                            TDelegate & delegate,
                            TDelegateD & delegateDirect,
                            Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                            TNeedle const & needle,
                            uint32_t needleId,
                            std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                            uint32_t const needleLeftPos,
                            uint32_t const needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                            TDir const & ,
                            TDistanceTag const &)
{
    std::vector<std::pair<uint32_t, uint32_t>> consOnes;
    getConsOnes(bitvectors, inside_bit_interval, ossContext.normal.intervalsize, consOnes);
    uint32_t noi = countSequences(*iter.fwdIter.index);

    //TODO shorten this
    for(uint32_t i = 0; i < consOnes.size(); ++i){
        if (std::is_same<TDir, Rev>::value){
            iter.revIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.revIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter.revIter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, false, Rev(), TDistanceTag());
        }
        else
        {
            iter.fwdIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.fwdIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter.fwdIter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, false, Fwd(), TDistanceTag());
        }
    }
}

template <typename TContex,
          typename TDelegateD,
          typename TNeedle,
          size_t nbrBlocks>
inline void genomeSearch(TContex & ossContext,
                         TDelegateD & delegateDirect,
                         TNeedle const & needle,
                         uint32_t needleId,
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
    delegateDirect(ossContext, sa_info, posAdd(sa_info, length(needle)), needle, needleId, errors);

}

template<typename TVector, typename TVSupport>
inline bool checkSinglePos(std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                           Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                           uint32_t offset)
{
    if(bitvectors.size() == 0)
        return true;
    else
        return (bitvectors[brange.i1].first[brange.i2.i1 + offset] == 1);
}

inline void saPosOnFwd(Pair<uint16_t, uint32_t> & sa_info,
                       uint32_t const genomelength,
                       uint32_t const occLength)
{
    sa_info.i2 = genomelength - sa_info.i2 - occLength;
}

template <typename TContex,
          typename TDelegateD,
          typename TString,
          typename TNeedle>
inline void alignmentMyersBitvector(TContex & ossContext,
                                    TDelegateD & delegateDirect,
                                    TNeedle const & needle,
                                    uint32_t needleId,
                                    TString const & n_infix,
                                    TString const & ex_infix,
                                    uint32_t const genomelength,
                                    Pair<uint16_t, uint32_t> const & sa_info,
                                    uint8_t max_e,
                                    uint8_t overlap_l,
                                    uint8_t overlap_r,
                                    uint8_t intDel,
                                    bool usingReverseText)
{
//     bool usingReverseText = true;
    //TODO insert return after each delegate call for only best alignment
    uint16_t needleL = length(needle);
    uint16_t ex_infixL = needleL + overlap_l + overlap_r;


    int initialScore = globalAlignmentScore(ex_infix, needle, MyersBitVector());

 //assume more Insertions (in the read) than deletions
    int ins_initialScore = globalAlignmentScore(n_infix, needle, MyersBitVector());

    if(ins_initialScore >= 0 - 2 * max_e || initialScore >= 0 - overlap_l - overlap_r - max_e + intDel) //MM creates one error D creates one error since now it also align to overlap
    {
        auto sa_info_tmp = sa_info;
        //No Insertions or Deletions
//         cout << "E: " << (int)0 << endl;
        TString const & tmp0 = infix(ex_infix, overlap_l, ex_infixL - overlap_r);
        int errors2 = 0 - globalAlignmentScore(tmp0, needle, MyersBitVector()); //
        if(errors2 <= max_e /*&& compareStartAndEnd(needle, tmp0, errors2)*/){
            if(usingReverseText){
                saPosOnFwd(sa_info_tmp, genomelength, needleL);
            }
            delegateDirect(ossContext, sa_info_tmp, posAdd(sa_info_tmp, length(needle)) , needle, needleId, errors2);
        }

        for(uint8_t e = 1; e <= max_e; ++e){
//             cout << "E: " << (int)e << endl;
            for(uint8_t del = 0; del <= e; ++del){
                //del is number of deletions
                uint8_t ins = e - del; //number of insertions
                sa_info_tmp = sa_info;

                if(del > 1 && ins == 0 || ins > 1 && del == 0){
                //only insertion or deletions
                    int16_t pos = (ins > del) ? 1 : (-1);
                    int16_t m = std::max(del,ins);
                    for(int16_t k = 0; k <= m; ++k)
                    {
                        if(!(0 <= overlap_l + (pos * k) && overlap_r >= 0 - (pos * (m - k))))
                            continue;
                        sa_info_tmp = sa_info;
                        sa_info_tmp.i2 = sa_info_tmp.i2 + (pos * k);
                        TString const & tmp2 = infix(ex_infix, overlap_l + (pos * k), ex_infixL - overlap_r - (pos * (m - k)));
                        errors2 = 0 - globalAlignmentScore(tmp2, needle, MyersBitVector());
                        if(errors2 <= max_e /*&& compareStartAndEnd(needle, tmp2, errors2)*/){
                            uint32_t occLength = length(needle) - (pos * m);
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
                            delegateDirect(ossContext, sa_info_tmp, posAdd(sa_info_tmp, occLength) , needle, needleId, errors2);
                        }
                    }
                }
                else
                {
//                     uint32_t occLength = length(needle) - ins + del;
                    //insertions left and deletion right
                    if(overlap_l >= del){
                        TString const & tmp = infix(ex_infix, overlap_l - del, ex_infixL - overlap_r - ins);
                        sa_info_tmp = sa_info;
                        sa_info_tmp.i2 = sa_info_tmp.i2 - del;
                        errors2 = 0 - globalAlignmentScore(tmp, needle, MyersBitVector());
                        if(errors2 <= max_e /*&& compareStartAndEnd(needle, tmp, errors2)*/){
                            uint32_t occLength = length(needle) - ins + del;
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
                            delegateDirect(ossContext, sa_info_tmp, posAdd(sa_info_tmp, occLength) , needle, needleId, errors2);
                        }
                    }

                    //insertions right and deletion left
                    if(overlap_r >= del){
                        sa_info_tmp = sa_info; //just include del from before into the calculation and delete this
                        TString const & tmp1 = infix(ex_infix, overlap_l + ins, ex_infixL - overlap_r + del);
                        errors2 = 0 - globalAlignmentScore(tmp1, needle, MyersBitVector());
                        sa_info_tmp.i2 = sa_info_tmp.i2 + ins;
                        if(errors2 <= max_e/* && compareStartAndEnd(needle, tmp1, errors2)*/){
                            uint32_t occLength = length(needle) - ins + del;
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
                            delegateDirect(ossContext, sa_info_tmp, posAdd(sa_info_tmp, occLength) , needle, needleId, errors2);
                        }
                    }
                }
            }
        }
    }
}


template <typename TContex,
          typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void directSearch(TContex & ossContext,
                         TDelegateD & delegateDirect,
                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                         TNeedle const & needle,
                         uint32_t needleId,
                         std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                         uint32_t const needleLeftPos,
                         uint32_t const needleRightPos,
                         uint8_t const errors,
                         OptimalSearch<nbrBlocks> const & s,
                         uint8_t const blockIndex,
                         Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                         TDir const & ,
                         TDistanceTag const &)
{
    auto const & genome = indexText(*iter.fwdIter.index);

    if (std::is_same<TDistanceTag, EditDistance>::value){

        //TODO if we are only interested in the best hit call return after delegate calls
        uint16_t needleL = length(needle);
        uint8_t max_e = s.u[s.u.size() - 1];
        uint8_t intIns = 0;
        uint8_t intDel = 0;
        //calculate net sum of internal Insertions - Deletions

        if(repLength(iter) < needleRightPos - needleLeftPos - 1)
            intIns = needleRightPos - needleLeftPos - 1 - repLength(iter);
        else
            intDel = repLength(iter) - (needleRightPos - needleLeftPos - 1);
        uint8_t overlap_l = max_e;
        uint8_t overlap_r = max_e;
        //TODO use intIns to modify overlap_r and overlap_l

//         if(needleLeftPos == 0)
//             overlap_l = intIns;
//         if(needleRightPos == needleL + 1)
//             overlap_r = intIns;

        uint16_t ex_infixL = needleL + overlap_l + overlap_r;
        for(uint32_t r = 0; r < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++r)
        {
    //         if(bitvectors[brange.i1].first[brange.i2.i1 + r] == 1){
            if(checkSinglePos(bitvectors, brange, r)){
                Pair<uint16_t, uint32_t> sa_info;
                uint32_t chromlength;
                if(std::is_same<TDir, Rev>::value){
                    sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + r];
                    chromlength = length(genome[sa_info.i1]);
                    if(!(needleLeftPos + overlap_l <= sa_info.i2  && chromlength - 1 >= sa_info.i2 - needleLeftPos + needleL - 1 + overlap_r))
                        continue;
                    sa_info.i2 = sa_info.i2 - needleLeftPos;
                }
                else
                {
                    sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + r];
                    chromlength = length(genome[sa_info.i1]);
                    if(!(chromlength - 1 >= sa_info.i2 + needleRightPos - 1 + overlap_r && sa_info.i2 + needleRightPos - 1 - overlap_l >= length(needle) + 1))
                        continue;
                    sa_info.i2 = chromlength - sa_info.i2 - needleRightPos + 1;
                }

                DnaString const & ex_infix = infix(genome[sa_info.i1], sa_info.i2 - overlap_l, sa_info.i2 + needleL + overlap_r);
                DnaString const & n_infix = infix(genome[sa_info.i1], sa_info.i2, sa_info.i2 + needleL);

                alignmentMyersBitvector(ossContext, delegateDirect, needle, needleId, n_infix, ex_infix, chromlength, sa_info, max_e, overlap_l, overlap_r, intDel, false);
            }
        }
    }
    else
    {
        std::array<uint32_t, nbrBlocks> blockStarts;
        std::array<uint32_t, nbrBlocks> blockEnds;
        std::copy(std::begin(s.blockStarts) + blockIndex, std::end(s.blockStarts), std::begin(blockStarts));
        std::copy(std::begin(s.blockEnds) + blockIndex, std::end(s.blockEnds), std::begin(blockEnds));

        if(std::is_same<TDir, Rev>::value){
            //modify blockstart in case we are still inside a block
            if(needleRightPos - 1 > blockStarts[0] && needleRightPos - 1 < blockEnds[0])
                blockStarts[0] = needleRightPos - 1;

            for(uint32_t i = 0; i < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++i){
//                 if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                if(checkSinglePos(bitvectors, brange, i)){
                    // mappability information is in reverse index order if we use the forward index
                    Pair<uint16_t, uint32_t> sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + i];
                    uint32_t const chromlength = length(genome[sa_info.i1]);
                    //Info make sure we dont DS search something going over the chromosom edge
                    //check left chromosom boundry && check right chromosom boundry
                    if(!(needleLeftPos <= sa_info.i2 && chromlength - 1 >= sa_info.i2 - needleLeftPos + length(needle) - 1))
                        continue;

                    sa_info.i2 = sa_info.i2 - needleLeftPos;

                    //search remaining blocks
                    genomeSearch(ossContext, delegateDirect, needle, needleId, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
                }
            }
        }
        else
        {
            //modify blockend in case we are still inside a block
            if(needleLeftPos > blockStarts[0] && needleLeftPos < blockEnds[0])
                blockEnds[0] = needleLeftPos;

            for(uint32_t i = 0; i < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++i){
//                 if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                if(checkSinglePos(bitvectors, brange, i)){
                    Pair<uint16_t, uint32_t> sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + i];
                    uint32_t const chromlength = length(genome[sa_info.i1]);
                    //check left chromosom boundry && check right chromosom boundry
                    if(!(chromlength - 1 >= sa_info.i2 + needleRightPos - 1 && sa_info.i2 + needleRightPos - 1 >= length(needle) + 1))
                        continue;
                    //calculate correct starting position of the needle  on the forward index
                    sa_info.i2 = chromlength - sa_info.i2 - needleRightPos + 1;

                    //search remaining blocks
                    genomeSearch(ossContext, delegateDirect, needle, needleId , errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
                }
            }
        }
    }
}

template <typename TText, typename TIndex, typename TIndexSpec,
          typename TDir>
inline void request_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                       uint8_t needed_bitvector,
                                       Pair<uint8_t, Pair<uint32_t, uint32_t>> & brangeOutput,
                                       TDir const & )
{
    Pair<uint32_t, uint32_t> dirrange = (std::is_same<TDir, Rev>::value) ? range(iter.fwdIter) : range(iter.revIter);
    uint32_t nseq = countSequences(*iter.fwdIter.index);
    dirrange.i1 = dirrange.i1 - nseq;
    dirrange.i2 = dirrange.i2 - nseq;

    brangeOutput.i1 = needed_bitvector;
    brangeOutput.i2 = dirrange;
}

//TODO load bitvectors inside a struct to make accessing the correct bitvector easier
template <typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks>
inline void get_bitvector_interval_inside(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                          std::vector<std::pair<TVector, TVSupport>> & bitvectors,
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
                       std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                       OptimalSearch<nbrBlocks> const & s,
                       uint8_t const blockIndex,
                       Pair<uint8_t, Pair<uint32_t, uint32_t>> & brangeOutput,
                       TDir const & )
{
    uint8_t needed_bitvector;
    if (std::is_same<TDir, Rev>::value)
        needed_bitvector = s.min[blockIndex] - 1;
    else
        needed_bitvector = bitvectors.size() - s.max[blockIndex];// + 1 - 1

    request_bitvector_interval(iter, needed_bitvector, brangeOutput, TDir());
}

template<typename TContex,
         typename TText, typename TIndex, typename TIndexSpec,
         typename TVector, typename TVSupport,
         size_t nbrBlocks>
inline bool testUnidirectionalFilter(TContex & ossContext,
                                     Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                     std::vector<std::pair<TVector, TVSupport>> & bitvectors,
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
        std::cerr << "Error bit vector has only zeroes this should have been checked by checkinterval" << "\n";
        exit(0);
    }

    float ivalSize = brange.i2.i2 - brange.i2.i1;
    uint32_t count = 0;

    if(ossContext.normal.testflipdensity){
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
    // it will contain mappability of the bitvector anchored at the other side of already searched needle

    // only interested in changes inside the supinterval (startPos - endPos)
    // allowed flips per intervalSize
    if(!ossContext.normal.testflipdensity || ivalSize * ossContext.normal.invflipdensity - 1 > static_cast<float>(count)){
        brange.i1 = bit_interval.i1;
        brange.i2.i1 = startPos;
        brange.i2.i2 = endPos;
        return true;
    }
    return false;
}

template<typename TContex,
         typename TVector, typename TVSupport,
         size_t nbrBlocks>
inline ReturnCode checkInterval(TContex & ossContext,
                                std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                                OptimalSearch<nbrBlocks> const & s,
                                uint8_t const blockIndex)
{
    TVector & b = bitvectors[brange.i1].first;
    TVSupport & rb = bitvectors[brange.i1].second;
    rb.set_vector(&b);

    uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
    uint32_t ivalSize = brange.i2.i2 - brange.i2.i1;

    if(ossContext.normal.nomappability && ivalOne == 0)
        return ReturnCode::NOMAPPABILITY;

//     ivalOne < (s.pi.size() - blockIndex - 1 + ossContext.normal.directsearchblockoffset) * ossContext.normal.directsearch_th
    if(ossContext.normal.directsearch && ossContext.itvCondition(s, blockIndex, ivalOne))
        return ReturnCode::DIRECTSEARCH;

    if(ossContext.normal.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
        return ReturnCode::COMPMAPPABLE;

    //equal or more than half zeroes
    if(ossContext.normal.suspectunidirectional && s.startUniDir <= blockIndex && static_cast<float>(ivalOne) / static_cast<float>(ivalSize) <= ossContext.normal.filter_th)
        return ReturnCode::SUSPECTUNIDIRECTIONAL;

    return ReturnCode::MAPPABLE;
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode checkCurrentMappability(TContex & ossContext,
                                          TDelegate & delegate,
                                          TDelegateD & delegateDirect,
                                          Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                          TNeedle const & needle,
                                          uint32_t needleId,
                                          std::vector<std::pair<TVector, TVSupport>> & bitvectors,
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

    ReturnCode rcode = checkInterval(ossContext, bitvectors, bit_interval, s, blockIndex);

    switch(rcode){
        case ReturnCode::NOMAPPABILITY:
            return ReturnCode::FINISHED;

        case ReturnCode::DIRECTSEARCH:
        {
            directSearch(ossContext, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::COMPMAPPABLE:
        {
            std::vector<std::pair<TVector, TVSupport>> emtpy_bitvectors;
            _optimalSearchSchemeChildren(ossContext, delegate, delegateDirect, iter, needle, needleId, emtpy_bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        default:
            return ReturnCode::MAPPABLE;
    }
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode checkMappability(TContex & ossContext,
                                   TDelegate & delegate,
                                   TDelegateD & delegateDirect,
                                   Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                   TNeedle const & needle,
                                   uint32_t needleId,
                                   std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                   uint32_t const current_needleLeftPos,
                                   uint32_t const current_needleRightPos,
                                   uint8_t const errors,
                                   OptimalSearch<nbrBlocks> const & s,
                                   uint8_t const blockIndex,
                                   bool const lastEdit,
                                   TDir const & ,
                                   TDistanceTag const &)
{
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
    get_bitvector_interval(iter, bitvectors, s, blockIndex, bit_interval, TDir());

    ReturnCode rcode = checkInterval(ossContext, bitvectors, bit_interval, s, blockIndex);
    switch(rcode)
    {
        case ReturnCode::NOMAPPABILITY:
            return ReturnCode::FINISHED;

        case ReturnCode::DIRECTSEARCH:
        {
            //search directly in Genome
            directSearch(ossContext, delegateDirect, iter, needle, needleId, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::COMPMAPPABLE:
        {
            std::vector<std::pair<TVector, TVSupport>> emtpy_bitvectors;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, emtpy_bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, lastEdit, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::SUSPECTUNIDIRECTIONAL:
        {
            //test unidirectional changes iter range if true
            //TODO modfy functions for TDIR
            bool goToRight2 = std::is_same<TDir, Rev>::value;
            if(testUnidirectionalFilter(ossContext, iter, bitvectors, bit_interval, s, blockIndex, goToRight2)){
                //range on iter was changed in function before
                filter_interval(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
                return ReturnCode::FINISHED;
            }
        }
        default:
            return ReturnCode::MAPPABLE;
    }
}

template <typename TContex,
          typename TDelegate,
          typename TDelegateDirect,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletion(TContex & ossContext,
                                         TDelegate & delegate,
                                         TDelegateDirect & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         uint32_t needleId,
                                         std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         bool const lastEdit,
                                         TDir const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    if (minErrorsLeftInBlock == 0)
    {

        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, lastEdit, Rev(), EditDistance());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, lastEdit, Fwd(), EditDistance());
    }

//     bool not_at_end = std::is_same<TDir, Rev>::value && needleRightPos != length(needle) + 1 || !std::is_same<TDir, Rev>::value && needleLeftPos != 0/* || true*/;

    if (/*not_at_end && */maxErrorsLeftInBlock > 0 && goDown(iter, TDir()))
    {
        do
        {
            _optimalSearchSchemeDeletion(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, true, TDir());
        } while (goRight(iter, TDir()));
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeChildren(TContex & ossContext,
                                         TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         uint32_t needleId,
                                         std::vector<std::pair<TVector, TVSupport>> & bitvectors,
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
            if (!std::is_same<TDistanceTag, EditDistance>::value && minErrorsLeftInBlock > 0 && charsLeft + delta < minErrorsLeftInBlock + 1u)
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
                    //use delta instead of false if no mismatches are allowed
                    _optimalSearchSchemeDeletion(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, false, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, false, Rev(), TDistanceTag());
                    else
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, false, Fwd(), TDistanceTag());
                }
            }
            else
            {
                //if want to disable mismatches at the start and end (!delta || not_at_end) && use delta instead of false
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, false, TDir(), TDistanceTag());
            }

            //Deletion
//             bool not_at_end = std::is_same<TDir, Rev>::value && needleRightPos2 != length(needle) + 1 || !std::is_same<TDir, Rev>::value && needleLeftPos2 != 0/* || true*/;
            if (std::is_same<TDistanceTag, EditDistance>::value/* && not_at_end*/)
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, true, TDir(), TDistanceTag());
        } while (goRight(iter, TDir()));
    }
}

template <typename T>
void printv(T a){
    for(int i = 0; i < a.size(); ++i){
        std::cout << static_cast<int> (a.at(i)) << ", ";
    }
    std::cout << "\n";
}

template <size_t nbrBlocks>
void print_search(OptimalSearch<nbrBlocks> const & s){
        std::cout << "Search sscheme: " << (int)s.id << "\n";
        std::cout << "Permutation: " << "\n";
        printv(s.pi);
        std::cout << "Lower bound: " << "\n";
        printv(s.l);
        std::cout << "Upper bound: " << "\n";
        printv(s.u);
        std::cout << "blockLengths: " << "\n";
        printv(s.blocklength);
        std::cout << "chronblockLengths: " << "\n";
        printv(s.chronBL);
        std::cout << "revchronblockLengths: " << "\n";
        printv(s.revChronBL);
        std::cout << "start Pos: " << "\n";
        std::cout << s.startPos << "\n";
        std::cout << "minMax: " << "\n";
        printv(s.min);
        printv(s.max);
        std::cout << "OneDirection" << "\n" << (int)s.startUniDir << "\n";
        std::cout << "\n";

}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeExact(TContex & ossContext,
                                      TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      uint32_t needleId,
                                      std::vector<std::pair<TVector, TVSupport>> & bitvectors,
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

        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, false, Fwd(), TDistanceTag());
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

        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, false, Fwd(), TDistanceTag());
    }
}

template <typename TContex,
          typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport>
inline void filteredDelegate(TContex & ossContext,
                             TDelegate & delegate,
                             Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                             TNeedle const & needle,
                             uint32_t needleId,
                             std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                             uint8_t const errors)
{
    Pair<uint8_t, Pair<uint32_t, uint32_t>> left_bit_interval;
    request_bitvector_interval(iter, 0, left_bit_interval, Rev());

    uint32_t rangeStart = iter.fwdIter.vDesc.range.i1;
    uint32_t rangeEnd = iter.fwdIter.vDesc.range.i2;
    uint32_t lastStart = 0;
    for(uint32_t i = 0; i < rangeEnd - rangeStart; ++i)
    {
        if(bitvectors[left_bit_interval.i1].first[left_bit_interval.i2.i1 + i] == 0 )
        {
            if(i != lastStart){
                iter.fwdIter.vDesc.range.i1 = rangeStart + lastStart;
                iter.fwdIter.vDesc.range.i2 = rangeStart + i - 1;
                delegate(ossContext, iter, needle, needleId, errors, false);
            }
            lastStart = i + 1;
        }
    }
    if(lastStart < rangeEnd - rangeStart){
        iter.fwdIter.vDesc.range.i1 = rangeStart + lastStart;
        iter.fwdIter.vDesc.range.i2 = rangeStart + rangeEnd - rangeStart;
        delegate(ossContext, iter, needle, needleId, errors, false);
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 bool const lastEdit,
                                 TDir const & ,
                                 TDistanceTag const &)
{
    //add strata tag
    if(ossContext.oneSSBestXMapper){
        bool save = false;
        uint32_t readId = getReadId(needleId, ossContext.readCount);
        if(errors > getCurrentErrors(ossContext.ctx, readId))
            save = true;
        if(isMapped(ossContext.ctx, readId)){
            if(errors > getMinErrors(ossContext.ctx, readId) + ossContext.strata)
                return;
        }
        if(save){
//             std::cout << "saving state: " << (int)errors << "\n";
//             std::cout << iter.fwdIter.vDesc.range << "\t" << needleLeftPos << "\t" << needleRightPos << "\t" << (int)s.id << "\t" << (int)blockIndex << "\n";
            bool right = std::is_same<TDir, Rev>::value;

//             std::cout << (int)needleId << "\t" << (int)errors << "\t" << (int)s.pi[0] << "\t" << right << "\n";

            ossContext.saveState(iter, needleLeftPos, needleRightPos, s.id, blockIndex, right, errors);
            return;
        }
    }

    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
    bool const done = minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1;
    bool const atBlockEnd = (blockIndex > 0) ? needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex - 1] : false;        //is not true if we finished needle
    bool const checkMappa = bitvectors.size() != 0;

    // Done. (Last step)
    if (done)
    {
        //last input only matters for unidirectional searches (has to be false in this case)
        if(/*!lastEdit*/true){
            if(checkMappa){
                filteredDelegate(ossContext, delegate, iter, needle, needleId, bitvectors, errors);
            }
            else
            {
                delegate(ossContext, iter, needle, needleId, errors, false);
            }
        }
        return;
    }

    //TODO add minErrorsLeftInBlock == 0 so only we check once after an insertion
    if(atBlockEnd && checkMappa){
        ReturnCode rcode = checkMappability(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, lastEdit, TDir(), TDistanceTag());
        if(rcode == ReturnCode::FINISHED)
            return;
    }

    // Exact search in current block.
    if (maxErrorsLeftInBlock == 0)
    {
        _optimalSearchSchemeExact(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
    }
    else if(!checkMappa && ossContext.itvConditionComp(iter, needleLeftPos, needleRightPos, errors, s, blockIndex))
    {
        //give emtpy bitvector and bitvector range sine we will not check mappability
        Pair<uint8_t, Pair<uint32_t, uint32_t>> dummy_bit_interval;
         directSearch(ossContext, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, dummy_bit_interval, TDir(), TDistanceTag());
    }

    // Approximate search in current block.
    else
    {
//         bool not_at_end = std::is_same<TDir, Rev>::value && needleRightPos != length(needle) || !std::is_same<TDir, Rev>::value && needleLeftPos != 1/* || true*/;
        // Insertion
        if (std::is_same<TDistanceTag, EditDistance>::value/* && not_at_end*/)
        {
            bool const goToRight = std::is_same<TDir, Rev>::value;
            int32_t const needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t const needleRightPos2 = needleRightPos + goToRight;

            //if we are at the end of block we need to add possible deletions because _optimalSearchScheme does not check it
            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                uint8_t const minErrorsLeftInBlock2 = (s.l[blockIndex] > (errors + 1)) ? (s.l[blockIndex] - (errors + 1)) : 0;
                if (minErrorsLeftInBlock2 == 0)
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

                    if (goToRight2)
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, true, Rev(), TDistanceTag());
                    else
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, true, Fwd(), TDistanceTag());
                }
            }
            else
            {
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, true, TDir(), TDistanceTag());
            }
        }
        //checkCurrentMappability
        uint32_t pblocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
        uint32_t step = (needleRightPos - needleLeftPos - 1);
        if(!atBlockEnd && checkMappa && ossContext.inBlockCheckMappabilityCondition(needleLeftPos, needleRightPos, s, blockIndex))
        {
            ReturnCode rcode = checkCurrentMappability(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
            if(rcode == ReturnCode::FINISHED)
                return;
        }
        _optimalSearchSchemeChildren(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TNeedle,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 TDistanceTag const &)
{
    bool initialDirection = s.pi[1] > s.pi[0];
        if(initialDirection)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Fwd(), TDistanceTag());
}

template<typename TIter, typename TSparseIter>
inline void loadIter(TIter & it, TSparseIter & itsparse){

    it.fwdIter.vDesc.range = itsparse.fwdRange;
    it.revIter.vDesc.range.i1 = itsparse.revRangeStart;
    it.revIter.vDesc.range.i2 = it.revIter.vDesc.range.i1 + it.fwdIter.vDesc.range.i2 - it.fwdIter.vDesc.range.i1;
    it.fwdIter.vDesc.repLen = itsparse.repLen;
    it.revIter.vDesc.repLen = itsparse.repLen;

/*
    it.fwdIter.vDesc.range = itsparse.fwd.range;
    it.fwdIter.vDesc.smaller = itsparse.fwd.smaller;
    it.fwdIter.vDesc.repLen = itsparse.fwd.repLen;
    it.fwdIter.vDesc.lastChar = itsparse.fwd.lastChar;

    it.revIter.vDesc.range = itsparse.rev.range;
    it.revIter.vDesc.smaller = itsparse.rev.smaller;
    it.revIter.vDesc.repLen = itsparse.rev.repLen;
    it.revIter.vDesc.lastChar = itsparse.rev.lastChar;

    it.fwdIter._parentDesc.range = itsparse.fwdp.range;
    it.fwdIter._parentDesc.smaller = itsparse.fwdp.smaller;
    it.fwdIter._parentDesc.repLen = itsparse.fwdp.repLen;
    it.fwdIter._parentDesc.lastChar = itsparse.fwdp.lastChar;

    it.revIter._parentDesc.range = itsparse.revp.range;
    it.revIter._parentDesc.smaller = itsparse.revp.smaller;
    it.revIter._parentDesc.repLen = itsparse.revp.repLen;
    it.revIter._parentDesc.lastChar = itsparse.revp.lastChar;*/
}

template<typename TIter>
inline void loadIter(TIter & it, TIter & stateIter){
    it = stateIter;
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks, size_t N,
          typename TNeedle,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 TDistanceTag const &)
{


    for (auto & s : ss)
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, bitvectors, s, needle, needleId, TDistanceTag());

    if(ossContext.oneSSBestXMapper){
        uint32_t readId = getReadId(needleId, ossContext.readCount);
        setCurrentErrors(ossContext.ctx, readId, 0); //TODO this is not needed

        for(uint8_t e = 1; e < ossContext.states.size() && e <= getMinErrors(ossContext.ctx, readId) + ossContext.strata; ++e){
            setCurrentErrors(ossContext.ctx, readId, e);
//             std::cout << "Read: " << readId << "\n";
            for(int j = 0; j < ossContext.states[e].size(); ++j){
                State<MySparseIter> & state = ossContext.states[e][j];
                MyIter tmp = it;
                loadIter(tmp, state.it);
//                 MyIter tmp = state.it;

/*
 * test creating a sparse Iter
                tmp.fwdIter.vDesc.range = state.it.fwdIter.vDesc.range;
                tmp.revIter.vDesc.range.i1 = state.it.revIter.vDesc.range.i1;
                tmp.revIter.vDesc.range.i2 = tmp.revIter.vDesc.range.i1 + tmp.fwdIter.vDesc.range.i2 - tmp.fwdIter.vDesc.range.i1;
                tmp.fwdIter.vDesc.repLen = state.it.fwdIter.vDesc.repLen;
                tmp.revIter.vDesc.repLen = tmp.fwdIter.vDesc.repLen;*/

                if(state.fwdDirection){
    //                 std::cout << "searching forward" << "\n";
                    _optimalSearchScheme(ossContext, delegate, delegateDirect, tmp/*state.it*/, needle, needleId, bitvectors,  state.nlp,  state.nrp, e, ss[state.sId],  state.blockIndex, false, Rev(), TDistanceTag());
                }else{
    //                 std::cout << "searching backwards" << "\n";
                    _optimalSearchScheme(ossContext, delegate, delegateDirect, tmp/*state.it*/, needle, needleId, bitvectors,  state.nlp,  state.nrp, e, ss[state.sId],  state.blockIndex, false, Fwd(), TDistanceTag());
                }
    //             std::cout << "Finished OSS: " << j << "\n";
            }

/*
            for(int i = e; i < ossContext.states.size(); ++i){
                std::cout << "Errors: " << i << "\t times \t" << ossContext.states[i].size();
                if(i == e)
                    std::cout << "\tsearching";
                std::cout << "\n";
            }

            if((getMinErrors(ossContext.ctx, readId) == e) || (e == 1 && getMinErrors(ossContext.ctx, readId) == 0)){
                std::cout << "Found with read with " << (int)getMinErrors(ossContext.ctx, readId) << " errors" << "\n";
            }
            std::cout << "\n";*/

            ossContext.states[e].clear();
        }

//         std::cout << "\n\n\n";

        for(int i = 0; i < ossContext.states.size(); ++i){
            ossContext.states[i].clear();
        }

    }
}

template <size_t minErrors, size_t maxErrors,
          typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks, size_t N,
          typename TChar, typename TStringSpec,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     std::vector<std::pair<TVector, TVSupport>> & bitvectors,
     std::array<OptimalSearch<nbrBlocks>, N> const & ss,
     String<TChar, TStringSpec> const & needle,
     uint32_t needleId,
     TDistanceTag const & )
{
//     auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    auto scheme = ss;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeComputeChronBlocklength(scheme);
//     Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
     MyIter it(index);
    _optimalSearchScheme(ossContext, delegate, delegateDirect, it, bitvectors, scheme, needle, needleId, TDistanceTag());
}

//TODO create Sortkey ?
bool id_smaller(const hit & x, const hit & y)
{
    return x.readId < y.readId;
}
/*
//TODO make it faster takes 8.4 out of 20 seconds .. canceling the speed up of 1.6 sec
template<typename TContex, typename THit>
void removeBadHits(TContex & ossContext, std::vector<THit> & hits)
{
    std::cout << "Start removing Bad Hits: " << "\n";
    uint32_t correct_errors;
    auto condition = [& correct_errors](const THit & v) { return v.errors > correct_errors;};

    auto start = std::chrono::high_resolution_clock::now();
    uint32_t count = hits.size();
    sort(hits.begin(), hits.end(), id_smaller);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Finished sorting after ID: " << elapsed.count() << "s" << "\n";

    uint32_t lastId = 0;
    for(auto it = hits.begin(); it != hits.end();){
        auto startI = it;
        lastId = (*it).readId;
        correct_errors = getMinErrors(ossContext.ctx, lastId) + ossContext.strata;
        while(lastId == (*it).readId && it != hits.end()){
            ++it;
        }
        it = hits.erase(std::remove_if(startI, it, condition), it);
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "Deleted In Text Verification Hits (which were not fulfilling the strata condition): " << count - hits.size() << "\t" << elapsed.count() << "s" << "\n";
}*/

template<typename TContex, typename hit>
void removeBadHits(TContex & ossContext, std::vector<hit> & hits)
{
    std::cout << "Start removing Bad Hits: " << "\n";
    uint32_t correct_errors;
    std::vector<bool> bad(hits.size());
//     auto condition = [& correct_errors](const THit & v) { return v.errors > correct_errors;};

    auto start = std::chrono::high_resolution_clock::now();
    uint32_t count = hits.size();
    sort(hits.begin(), hits.end(), id_smaller);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Finished sorting after ID: " << elapsed.count() << "s" << "\n";

    uint32_t lastId = 0;
    auto itb = bad.begin();
    for(auto it = hits.begin(); it != hits.end();){
        lastId = (*it).readId;
        correct_errors = getMinErrors(ossContext.ctx, lastId) + ossContext.strata;
        while(lastId == (*it).readId && it != hits.end()){
            *itb = (*it).errors > correct_errors;
            ++itb;
            ++it;
        }
    }
    hits.erase(std::remove_if(hits.begin(), hits.end(),
      [&bad, &hits](auto const &i) { return bad.at(&i - hits.data()); }), hits.end());

    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "Deleted In Text Verification Hits (which were not fulfilling the strata condition): " << count - hits.size() << "\t" << elapsed.count() << "s" << "\n";
}

template <size_t minErrors, size_t maxErrors,
          typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TVector, typename TVSupport,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     std::vector<std::pair<TVector, TVSupport>> & bitvectors, // cant be const since TVSupport.set_vector(&TVector)
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & )
{

    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    calcConstParameters(scheme);
    uint32_t needleId = 0;
    uint32_t lastcount = 0;
    while(needleId < length(needles))
    {
//         find<minErrors, maxErrors>(ossContext, delegate, delegateDirect, index, bitvectors, scheme, needles[needleId], needleId, TDistanceTag());
//         ++needleId;

        bool skip = false;
        if(ossContext.bestXMapper){
            uint32_t readId = getReadId(needleId, ossContext.readCount);
            if(isMapped(
                ossContext.ctx, readId)){
                if(getMinErrors(ossContext.ctx, readId) + ossContext.strata < minErrors){
                    skip = true;
//                     std::cout << "Skip: " << needleId << "\n";
                }
            }
        }
        if(!skip){
            find<minErrors, maxErrors>(ossContext, delegate, delegateDirect, index, bitvectors, scheme, needles[needleId], needleId, TDistanceTag());
        }
        needleId++;

        /*
        if(ossContext.trackReadCount){
            uint32_t currentcount = ossContext.hits.size() + ossContext.dhits.size() - lastcount;
            ossContext.readOccCount.push_back(currentcount);
            lastcount += currentcount;
        }*/

    }
//     takes to much time for editDistance
    // wrong hits found with itv (we dont know before hand how many errors will be reported) //TODO make it optional
    if(ossContext.itv && ossContext.oneSSBestXMapper){
        removeBadHits(ossContext, ossContext.dhits);
    }
}
//



/*
template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    std::vector<std::pair<TBitvector, TSupport>> emtpy_bitvectors;
    switch (maxErrors)
    {
        case 1: find<0, 1>(ossContext, delegate, delegateDirect, index, emtpy_bitvectors, needles, TDistanceTag());
                break;
        case 2: find<0, 2>(ossContext, delegate, delegateDirect, index, emtpy_bitvectors, needles, TDistanceTag());
                break;
        case 3: find<0, 3>(ossContext, delegate, delegateDirect, index, emtpy_bitvectors, needles, TDistanceTag());
                break;
        case 4: find<0, 4>(ossContext, delegate, delegateDirect, index, emtpy_bitvectors, needles, TDistanceTag());
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}*/

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TVector, typename TVSupport,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 std::vector<std::pair<TVector, TVSupport>> & bitvectors, // cant be const since TVSupport.set_vector(&TVector)
                 StringSet<TNeedle, TStringSetSpec> const & needles,

                 TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 2: find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 3: find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 4: find<0, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

//no strata (needed for one Scheme Best X)
template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    std::vector<std::pair<TBitvector, TSupport>> emtpy_bitvectors;
    find(minErrors, maxErrors, ossContext, delegate, delegateDirect, index, needles, emtpy_bitvectors, TDistanceTag());

}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TVector, typename TVSupport,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 int strata,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 std::vector<std::pair<TVector, TVSupport>> & bitvectors, // cant be const since TVSupport.set_vector(&TVector
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    if(strata == 99 || ossContext.oneSSBestXMapper)
        strata = maxErrors;

    switch (maxErrors)
    {
        case 1:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 2:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 3:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 3 :
                {
                    find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 4:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 3 :
                {
                    find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 4 :
                {
                    find<0, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;
        }
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                 exit(1);
    }
    //TODO remove since not used
    if(ossContext.itv && ossContext.bestXMapper){
        removeBadHits(ossContext, ossContext.dhits);
    }
}

//no bitvectors so substitute them
template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 const int strata,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    std::vector<std::pair<TBitvector, TSupport>> emtpy_bitvectors;
    find(minErrors, maxErrors, strata, ossContext, delegate, delegateDirect, index, emtpy_bitvectors, needles, TDistanceTag());
}

// for find2_index_approx.h find function
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
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

}



#endif
