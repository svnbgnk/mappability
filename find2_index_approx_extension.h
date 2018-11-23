#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

#include <iostream>
#include <sdsl/bit_vectors.hpp>
// #include "global.h"
#include "auxiliary.h"
#include "common_auxiliary.h"
#include "find2_index_approx_unidirectional.h"
// #include "find2_index_approx_compmappable.h"
// #include "find2_index_approx_start_unidirectional.h"

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
            consOnesOutput.push_back(make_pair(startOneInterval, k));
            startOneInterval = k + interval;
        }
        k += interval;
        interval = 0;
        ++k;
    }
    consOnesOutput.push_back(make_pair(startOneInterval, k));
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
    delegateDirect(sa_info, posAdd(sa_info, length(needle)), needle, needleId, errors);

}

//find first match without insertion and errors at the beginning and end
template<typename TNeedle,
         typename TString>
inline bool compareStartAndEnd(TNeedle const & needle,
                              TString const & infix,
                              uint8_t const alignment_errors)
{
    // It is impossible without additional alignments to check for Insertion and deletions at the beginning and end therefore they will be allowed
    // In the future just the first best hit will be reported
//     std::cout << "compare start and end" << "\n";
    uint8_t errors = 0;
    auto lastBase = needle[0];
    for(int i = 0; i <= alignment_errors; ++i){
//         std::cout << (char)infix[i] << "\t" << (char)needle[i] << "\n";
        if(infix[i] == needle[i]){
            break;
        }else{
            ++errors;
        }
    }
    if(errors > alignment_errors){
//         std::cout << "Denied" << "\n";
        return false;
    }
//     std::cout << "Errors: " << (int)errors << "\n";
    for(int i = 0; i <= alignment_errors - errors; ++i){
//         std::cout << (char)infix[length(infix) - 1 - i] << "\t" << (char)needle[length(needle) - 1 - i] << "\n";
        if(infix[length(infix) - 1 - i] == needle[length(needle) - 1 - i]){
//             std::cout << "i: " << i << "\n";
//             std::cout << "Accepted" << "\n";
            return true;
        }
    }
//     std::cout << "Denied" << "\n";
    return false;

    //allow also no mismatches at the start and end
//     return(infix[0] == needle[0] && infix[length(infix) - 1] == needle[length(needle) - 1]);
}

template<typename TVector, typename TVSupport>
inline bool checkSinglePos(vector<pair<TVector, TVSupport>> & bitvectors,
                           Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                           uint32_t offset)
{
    if(bitvectors.size() == 0)
        return true;
    else
        return (bitvectors[brange.i1].first[brange.i2.i1 + offset] == 1);
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
//     bool usingReverseText = false;

    //TODO insert return after each delegate call for only best alignment
    uint16_t needleL = length(needle);
    uint16_t ex_infixL = needleL + overlap_l + overlap_r;


    int initialScore = globalAlignmentScore(ex_infix, needle, MyersBitVector());

 //assume more Insertions (in the read) than deletions
    int ins_initialScore = globalAlignmentScore(n_infix, needle, MyersBitVector());

    if(ins_initialScore >= 0 - 2 * max_e || initialScore >= 0 - overlap_l - overlap_r - max_e + intDel) //MM creates one error D creates one error since now it also align to overlap
    {
        auto sa_info_tmp = sa_info;
//         cout << ex_infix << "        ex_infix " << (int)overlap_l << "  " << (int)overlap_r << "\n";
//         cout << needle << "        needle" << "\n";
        //No Insertions or Deletions
//         cout << "E: " << (int)0 << endl;
        TString const & tmp0 = infix(ex_infix, overlap_l, ex_infixL - overlap_r);
        int errors2 = 0 - globalAlignmentScore(tmp0, needle, MyersBitVector()); //
        if(errors2 <= max_e /*&& compareStartAndEnd(needle, tmp0, errors2)*/){
//             std::cout << "c1 " << sa_info << "  " << (int) errors2 << "\n";
//             std::cout << tmp0 << "\n";
            if(usingReverseText){
                saPosOnFwd(sa_info_tmp, genomelength, needleL);
            }
            delegateDirect(sa_info, posAdd(sa_info, length(needle)) , needle, needleId, errors2);
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
//                         std::cout << (int)k << ":" << (int)m-k << "\t" << (int)pos << "\n";
//                         std::cout << (int)overlap_l << ":" << (int)(pos * k) << " :: " << (int)overlap_r << ":" << (int)(pos * (m - k))  << endl;
                        if(!(0 <= overlap_l + (pos * k) && overlap_r >= 0 - (pos * (m - k))))
                            continue;
                        sa_info_tmp = sa_info;
                        sa_info_tmp.i2 = sa_info_tmp.i2 + (pos * k);
                        TString const & tmp2 = infix(ex_infix, overlap_l + (pos * k), ex_infixL - overlap_r - (pos * (m - k)));
                        errors2 = 0 - globalAlignmentScore(tmp2, needle, MyersBitVector());
                        if(errors2 <= max_e /*&& compareStartAndEnd(needle, tmp2, errors2)*/){
//                             std::cout << "c2 " << sa_info_tmp << "  " << (int) errors2 << "\n";
//                             std::cout << tmp2 << "\n";
                            uint32_t occLength = length(needle) - (pos * m);
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
                            delegateDirect(sa_info_tmp, posAdd(sa_info_tmp, occLength) , needle, needleId, errors2);
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
//                             std::cout << "c3 " << sa_info_tmp << "  " << (int) errors2 << "\n";
//                             std::cout << tmp << "\n";
                            uint32_t occLength = length(needle) - ins + del;
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
                            delegateDirect(sa_info_tmp, posAdd(sa_info_tmp, occLength) , needle, needleId, errors2);
                        }
                    }

                    //insertions right and deletion left
                    if(overlap_r >= del){
                        sa_info_tmp = sa_info; //just include del from before into the calculation and delete this
                        TString const & tmp1 = infix(ex_infix, overlap_l + ins, ex_infixL - overlap_r + del);
                        errors2 = 0 - globalAlignmentScore(tmp1, needle, MyersBitVector());
                        sa_info_tmp.i2 = sa_info_tmp.i2 + ins;
                        if(errors2 <= max_e/* && compareStartAndEnd(needle, tmp1, errors2)*/){
//                             std::cout << "c4 " << sa_info_tmp << "  " << (int) errors2 << "\n";
//                             std::cout << tmp1 << "\n";
                            uint32_t occLength = length(needle) - ins + del;
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
                            delegateDirect(sa_info_tmp, posAdd(sa_info_tmp, occLength) , needle, needleId, errors2);

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
    auto const & genome = indexText(*iter.fwdIter.index);

    if (std::is_same<TDistanceTag, EditDistance>::value){
        //TODO put this into a function
        //TODO if we are only interested in the best hit call return after delegate calls
        uint16_t needleL = length(needle);
        uint8_t max_e = s.u[s.u.size() - 1];
        int intIns = 0;
        int intDel = 0;
        //calculate net sum of internal Insertions - Deletions

        if(repLength(iter) < needleRightPos - needleLeftPos - 1)
            intIns = needleRightPos - needleLeftPos - 1 - repLength(iter);
        else
            intDel = repLength(iter) - (needleRightPos - needleLeftPos - 1);
        uint8_t overlap_l = max_e;
        uint8_t overlap_r = max_e;
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

                TString const & ex_infix = infix(genome[sa_info.i1], sa_info.i2 - overlap_l, sa_info.i2 + needleL + overlap_r);
                TString const & n_infix = infix(genome[sa_info.i1], sa_info.i2, sa_info.i2 + needleL);

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
                                vector<pair<TVector, TVSupport>> & bitvectors,
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
            vector<pair<TVector, TVSupport>> emtpy_bitvectors;
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
                                   vector<pair<TVector, TVSupport>> & bitvectors,
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
            vector<pair<TVector, TVSupport>> emtpy_bitvectors;
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
                                         vector<pair<TVector, TVSupport>> & bitvectors,
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

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport>
inline void filteredDelegate(TDelegate & delegate,
                             Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                             TNeedle const & needle,
                             uint32_t needleId,
                             vector<pair<TVector, TVSupport>> & bitvectors,
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
                delegate(iter, needle, needleId, errors, false);
            }
            lastStart = i + 1;
        }
    }
    if(lastStart < rangeEnd - rangeStart){
        iter.fwdIter.vDesc.range.i1 = rangeStart + lastStart;
        iter.fwdIter.vDesc.range.i2 = rangeStart + rangeEnd - rangeStart;
        delegate(iter, needle, needleId, errors, false);
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
                                 vector<pair<TVector, TVSupport>> & bitvectors,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 bool const lastEdit,
                                 TDir const & ,
                                 TDistanceTag const &)
{
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
                filteredDelegate(delegate, iter, needle, needleId, bitvectors, errors);
            }
            else
            {
                delegate(iter, needle, needleId, errors, false);
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
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 vector<pair<TVector, TVSupport>> & bitvectors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 TDistanceTag const &)
{
    bool initialDirection = s.pi[1] > s.pi[0];
//     if(!params.startUnidirectional || s.startUniDir > 0){
        if(initialDirection)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Fwd(), TDistanceTag());
//     }else{
//         if(initialDirection)
//             _uniOptimalSearchScheme(delegate, delegateDirect, it.revIter, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Rev(), TDistanceTag());
//         else
//             _uniOptimalSearchScheme(delegate, delegateDirect, it.fwdIter, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Fwd(), TDistanceTag());
//     }
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 OptimalSearch<nbrBlocks> const & s,
                                 TDistanceTag const &)
{
    vector<pair<TBitvector, TSupport>> emtpy_bitvectors;
    bool initialDirection = s.pi[1] > s.pi[0];
    if(initialDirection)
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, emtpy_bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Rev(), TDistanceTag());
    else
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, emtpy_bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Fwd(), TDistanceTag());
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 vector<pair<TVector, TVSupport>> & bitvectors,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TDistanceTag const &)
{
    for (auto & s : ss)
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, bitvectors, s, TDistanceTag());
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TDistanceTag const &)
{
    for (auto & s : ss)
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, s, TDistanceTag());
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec,
          typename TVector, typename TVSupport,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     uint32_t needleId,
     vector<pair<TVector, TVSupport>> & bitvectors,
     std::array<OptimalSearch<nbrBlocks>, N> const & ss,
     TDistanceTag const &)
{
    //TODO there has to be a better way
    auto scheme = ss;/*OptimalSearchSchemes<minErrors, maxErrors>::VALUE;*/
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, bitvectors, scheme, TDistanceTag());
}

template <size_t minErrors, size_t maxErrors,
          typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     uint32_t needleId,
     TDistanceTag const & )
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, scheme, TDistanceTag());
}


template <size_t minErrors, size_t maxErrors,
          typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TVector, typename TVSupport,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<pair<TVector, TVSupport>> & bitvectors,
     TDistanceTag const & )
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    calcConstParameters(scheme);
//     checkTime time;
    int k = 0;
    uint32_t lastcount = 0;
    while(k < length(needles))
    {
        find(ossContext, delegate, delegateDirect, index, needles[k], k, bitvectors, scheme, TDistanceTag());
        /*
        if(!(params.clocking && time.stopnow(params.terminateDuration))){
            find(ossContext, delegate, delegateDirect, index, needles[k], k, bitvectors, scheme, TDistanceTag());
        }else{
            params.wasStopped = true;
        }*/
        ++k;

        if(ossContext.trackReadCount){
            uint32_t currentcount = ossContext.hits.size() + ossContext.dhits.size() - lastcount;
            ossContext.readOccCount.push_back(currentcount);
            lastcount += currentcount;
        }
    }
}

template <size_t minErrors, size_t maxErrors,
          typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & )
{
    int k = 0;
    uint32_t lastcount = 0;
    while(k < length(needles))
    {
        find<minErrors, maxErrors>(ossContext, delegate, delegateDirect, index, needles[k], k, TDistanceTag());
        ++k;
        if(ossContext.trackReadCount){
            uint32_t currentcount = ossContext.hits.size() + ossContext.dhits.size() - lastcount;
            ossContext.readOccCount.push_back(currentcount);
            lastcount += currentcount;
        }

    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TVector, typename TVSupport,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 vector<pair<TVector, TVSupport>> & bitvectors, // cant be const since TVSupport.set_vector(&TVector)
                 TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(ossContext, delegate, delegateDirect, index, needles, bitvectors, TDistanceTag());
                break;
        case 2: find<0, 2>(ossContext, delegate, delegateDirect, index, needles, bitvectors, TDistanceTag());
                break;
        case 3: find<0, 3>(ossContext, delegate, delegateDirect, index, needles, bitvectors, TDistanceTag());
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

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
    switch (maxErrors)
    {
        case 1: find<0, 1>(ossContext, delegate, delegateDirect, index, needles, TDistanceTag());
                break;
        case 2: find<0, 2>(ossContext, delegate, delegateDirect, index, needles, TDistanceTag());
                break;
        case 3: find<0, 3>(ossContext, delegate, delegateDirect, index, needles, TDistanceTag());
                break;
        case 4: find<0, 4>(ossContext, delegate, delegateDirect, index, needles, TDistanceTag());
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

}



#endif
