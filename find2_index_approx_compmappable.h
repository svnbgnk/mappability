#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_COMPMAPPLE_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_COMPMAPPLE_H_

#include <seqan/align.h>

using namespace std;

namespace seqan{




template <typename TDelegateD,
          typename TString, typename TConfig, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void directSearch(TDelegateD & delegateDirect,
                  Iter<Index<StringSet<TString, TConfig>, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  TDir const & /**/,
                  TDistanceTag const & )
{
    auto const & genome = indexText(*iter.fwdIter.index);

    if (std::is_same<TDistanceTag, EditDistance>::value){
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
        if(needleLeftPos == 0)
            overlap_l = intIns;
        if(needleRightPos == needleL + 1)
            overlap_r = intIns;
        uint16_t ex_infixL = needleL + overlap_l + overlap_r;

        for(uint32_t r = iter.fwdIter.vDesc.range.i1; r < iter.fwdIter.vDesc.range.i2; ++r)
        {
            Pair<uint16_t, uint32_t> sa_info = iter.fwdIter.index->sa[r];
            uint32_t const chromlength = length(genome[sa_info.i1]);
            if(!(needleLeftPos + overlap_l <= sa_info.i2  && chromlength - 1 >= sa_info.i2 - needleLeftPos + needleL - 1 + overlap_r))
                continue;
            sa_info.i2 = sa_info.i2 - needleLeftPos;

            TString const & ex_infix = infix(genome[sa_info.i1], sa_info.i2 - overlap_l, sa_info.i2 + needleL + overlap_r);
            TString const & n_infix = infix(genome[sa_info.i1], sa_info.i2, sa_info.i2 + needleL);

            //Unfortunately, since for example (E = 2) 2 Insertions lead to a score of -4 we cannot conclude the number of errors the read will have at the end. (O Errors would also have a score -4 2 Deletion in the beginning to Insertions at the end)
            //for E = 2 we allow a score -6 since 2 MM + 2 D + 2 I.
            //To Insertion on the otherhand lead to 2D + 2I + 4D -> 8 errors

            int initialScore = globalAlignmentScore(ex_infix, needle, MyersBitVector());

            //assume more Insertions (in the read) than deletions
            int ins_initialScore = globalAlignmentScore(n_infix, needle, MyersBitVector());

            if(ins_initialScore >= 0 - 2 * max_e || initialScore >= 0 - overlap_l - overlap_r - max_e + intDel) //MM creates one error D creates one error since now it also align to overlap
            {
                cout << ex_infix << "        ex_infix " << (int)overlap_l << "  " << (int)overlap_r << "\n";
                cout << needle << "        needle" << "\n";
                //No Insertions or Deletions
                TString const & tmp0 = infix(ex_infix, overlap_l, ex_infixL - overlap_r);
                int errors2 = 0 - globalAlignmentScore(tmp0, needle, MyersBitVector()); //
                if(errors2 <= max_e && tmp0[0] == needle[0] && tmp0[length(tmp0) - 1] == needle[needleL - 1]){
                    std::cout << "c1 " << sa_info << "  " << (int) errors2 << "\n";
                    std::cout << tmp0 << "\n";
                    delegateDirect(sa_info , needle, errors2);
                }

                for(uint8_t e = 1; e <= max_e /*overlap*/; ++e){
//                     cout << "E: " << (int)e << endl;
                    for(uint8_t del = 0; del <= e; ++del){
                        //del is number of deletions
                        uint8_t ins = e - del; //number of insertions
                        auto sa_info_tmp = sa_info;

                        if(del > 1 && ins == 0 || ins > 1 && del == 0){
                        //only insertion or deletions
                            int16_t pos = (ins > del) ? 1 : (-1);
                            int16_t m = std::max(del,ins);
                            for(int16_t k = 0; k <= m; ++k)
                            {
//                                 cout << (int)k << ":" << (int)m-k << "\t" << (int)pos << "\n";
//                                 cout << (int)overlap_l << ":" << (int)(pos * k) << " :: " << (int)overlap_r << ":" << (int)(pos * (m - k))  << endl;
                                if(!(0 < overlap_l + (pos * k) && overlap_r > 0 - (pos * (m - k))))
                                    continue;

                                sa_info_tmp = sa_info;
                                sa_info_tmp.i2 = sa_info_tmp.i2 + (pos * k);
                                TString const & tmp2 = infix(ex_infix, overlap_l + (pos * k), ex_infixL - overlap_r - (pos * (m - k)));
                                errors2 = 0 - globalAlignmentScore(tmp2, needle, MyersBitVector());
                                if(errors2 <= max_e && tmp2[0] == needle[0] && tmp2[length(tmp2) - 1] == needle[needleL - 1]){
                                    std::cout << "c2 " << sa_info_tmp << "  " << (int) errors2 << "\n";
                                    std::cout << tmp2 << "\n";
                                    delegateDirect(sa_info_tmp , needle, errors2);
                                }
                            }
                        }
                        else
                        {
                            //insertions left and deletion right
                            if(overlap_l >= del){
                                TString const & tmp = infix(ex_infix, overlap_l - del, ex_infixL - overlap_r - ins);
                                sa_info_tmp.i2 = sa_info_tmp.i2 - del;
                                errors2 = 0 - globalAlignmentScore(tmp, needle, MyersBitVector());
                                if(errors2 <= max_e && tmp[0] == needle[0] && tmp[length(tmp) - 1] == needle[needleL - 1]){
                                    std::cout << "c3 " << sa_info_tmp << "  " << (int) errors2 << "\n";
                                    std::cout << tmp << "\n";
                                    delegateDirect(sa_info_tmp , needle, errors2);
                                }
                            }

                            //insertions right and deletion left
                            if(overlap_r >= del){
                                sa_info_tmp = sa_info; //just include del from before into the calculation and delete this
                                TString const & tmp1 = infix(ex_infix, overlap_l + ins, ex_infixL - overlap_r + del);
                                errors2 = 0 - globalAlignmentScore(tmp1, needle, MyersBitVector());
                                sa_info_tmp.i2 = sa_info_tmp.i2 + ins;
                                if(errors2 <= max_e && tmp1[0] == needle[0] && tmp1[length(tmp1) - 1] == needle[needleL - 1]){
                                    std::cout << "c4 " << sa_info_tmp << "  " << (int) errors2 << "\n";
                                    std::cout << tmp1 << "\n";
                                    delegateDirect(sa_info_tmp , needle, errors2);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        //cut of blockStarts and Ends that where already checked by the search
        std::array<uint32_t, nbrBlocks> blockStarts;
        std::array<uint32_t, nbrBlocks> blockEnds;
        std::copy(std::begin(s.blockStarts) + blockIndex, std::end(s.blockStarts), std::begin(blockStarts));
        std::copy(std::begin(s.blockEnds) + blockIndex, std::end(s.blockEnds), std::begin(blockEnds));

        if(std::is_same<TDir, Rev>::value){
            if(needleRightPos - 1 > blockStarts[0] && needleRightPos - 1 < blockEnds[0])
                blockStarts[0] = needleRightPos - 1;
        }
        else
        {
            if(needleLeftPos > blockStarts[0] && needleLeftPos < blockEnds[0])
                blockEnds[0] = needleLeftPos;
        }

        //iterate over each block according to search scheme
        for(uint32_t i = iter.fwdIter.vDesc.range.i1; i < iter.fwdIter.vDesc.range.i2; ++i)
        {
            bool valid = true;
            Pair<uint16_t, uint32_t> sa_info = iter.fwdIter.index->sa[i];
            //dont need look at the reverse index in this case since i dont use mappability
            uint32_t const chromlength = length(genome[sa_info.i1]);
            if(!(needleLeftPos <= sa_info.i2 && chromlength - 1 >= sa_info.i2 - needleLeftPos + length(needle) - 1))
                continue;

            sa_info.i2 = sa_info.i2 - needleLeftPos;
            uint8_t errors2 = errors;
            for(uint32_t j = 0; j < nbrBlocks - blockIndex; ++j){
                // compare bases to needle
                for(uint32_t k = blockStarts[j]; k <  blockEnds[j]; ++k){
                    if(needle[k] != genome[sa_info.i1][sa_info.i2 + k])
                        ++errors2;
                }
                if(errors2 < s.l[blockIndex + j] || errors2 > s.u[blockIndex + j]){
                    valid = false;
                    break;
                }
            }
            if(valid)
                delegateDirect(sa_info, needle, errors2);
        }
    }
}

template <typename TDelegate,
          typename TDelegateDirect,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletion(TDelegate & delegate,
                                         TDelegateDirect & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
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
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex2, lastEdit, Rev(), EditDistance());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex2, lastEdit, Fwd(), EditDistance());
    }

    bool not_surface =  std::is_same<TDir, Rev>::value && needleRightPos != length(needle) + 1 || !std::is_same<TDir, Rev>::value && needleLeftPos != 0/* || true*/;

    if (not_surface && maxErrorsLeftInBlock > 0 && goDown(iter, TDir()))
    {
        do
        {
//             std::cout << "Deletionb " << needleLeftPos << " : " << needleRightPos << "(" << (int) errors << "\n";
            _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors + 1, s,
                                         blockIndex, true, TDir());
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeChildren(TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const &,
                                         TDistanceTag const & )
{
    bool goToRight = std::is_same<TDir, Rev>::value;
    if (goDown(iter, TDir()))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()),
                                   needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1]);

            // NOTE (cpockrandt): this might not be optimal yet! we have more edges than in the theoretical model,
            // since we go down an edge before we check whether it can even work out!
            if (!std::is_same<TDistanceTag, EditDistance>::value && minErrorsLeftInBlock > 0 && charsLeft + delta < minErrorsLeftInBlock + 1u) // charsLeft - 1 < minErrorsLeftInBlock - delta
            {
                continue;
            }

            int32_t needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t needleRightPos2 = needleRightPos + goToRight;

            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                if (std::is_same<TDistanceTag, EditDistance>::value)
                {
                    _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, delta, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                    {
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, delta /*false*/, Rev(), TDistanceTag());
                    }
                    else
                    {
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, delta /*false*/, Fwd(), TDistanceTag());
                    }
                }
            }
            else
            {
//                 bool not_surface = repLength(iter) > 1 && s.pi[blockIndex] != 1 || s.pi[blockIndex] != s.pi.size()/* || true*/;
                    bool not_surface =  std::is_same<TDir, Rev>::value && needleRightPos2 != length(needle) + 1 || !std::is_same<TDir, Rev>::value && needleLeftPos2 != 0/* || true*/;
                if(!delta || not_surface) //TODO maybe reverse MM
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,
                                     blockIndex, delta/*true*/, TDir(), TDistanceTag()); //TODO maybe reverse MM
            }

            // Deletion
            bool not_surface =  std::is_same<TDir, Rev>::value && needleRightPos2 != length(needle) + 1 || !std::is_same<TDir, Rev>::value && needleLeftPos2 != 0/* || true*/;
//             std::cout << "SS: " << (int) s.pi[0] << "\n";
//             std::cout << "length: " << repLength(iter) << "  NLP: " << needleLeftPos << "  NRP: " << needleRightPos;

            if (std::is_same<TDistanceTag, EditDistance>::value && not_surface){
//                 std::cout << "Deletion. " << needleLeftPos << " : " << needleRightPos << "(" << (int) errors << "\n";
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, true, TDir(), TDistanceTag());
            }
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeExact(TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const & ,
                                      TDistanceTag const &)
{
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
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                 blockIndex2, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                 blockIndex2, false, Fwd(), TDistanceTag());
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
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                 blockIndex2, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                 blockIndex2, false, Fwd(), TDistanceTag());
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 bool const lastEdit,
                                 TDir const & ,
                                 TDistanceTag const & )
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    // Done.
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        if(!lastEdit /*|| true*/){
/*
            std::cout << "Pos: " << "\n";
            for(uint32_t i = iter.fwdIter.vDesc.range.i1; i < iter.fwdIter.vDesc.range.i2; ++i)
                std::cout << iter.fwdIter.index->sa[i] << std::endl;
            std::cout << "SS: " << (int) s.pi[0] << "\n";
            std::cout << "Piece: " << (int) s.pi[blockIndex] << "\n";
            std::cout << "Error: " << (int) errors << "\n";*/

            delegate(iter, needle, errors, false);

        }
    }
    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(),
                                  TDistanceTag());
    }

    else
    {
        // Insertion
//         bool not_surface = repLength(iter) > 0 && s.pi[blockIndex] != 1 || s.pi[blockIndex] != s.pi.size()/* || true*/;
        bool not_surface =  std::is_same<TDir, Rev>::value && needleRightPos != length(needle) || !std::is_same<TDir, Rev>::value && needleLeftPos != 1/* || true*/;

        if (std::is_same<TDistanceTag, EditDistance>::value && not_surface)
        {
            bool const goToRight = std::is_same<TDir, Rev>::value;
            int32_t const needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t const needleRightPos2 = needleRightPos + goToRight;
//             std::cout << "Insertion. " << needleLeftPos2 << " : " << needleRightPos2 << "(" << (int) errors << "\n";

            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                uint8_t const minErrorsLeftInBlock2 = (s.l[blockIndex] > (errors + 1)) ? (s.l[blockIndex] - (errors + 1)) : 0;
                if (minErrorsLeftInBlock2 == 0)
                {
                    uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, true, Rev(), EditDistance());
                    else
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, true, Fwd(), EditDistance());
                }
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, true, TDir(), TDistanceTag());
            }
        }

        if(params.comp.directsearch && iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1 + params.uni.directsearchblockoffset) * params.comp.directsearch_th)
        {
            directSearch(delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
            return;
        }
        _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex,
                                     minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 TDir const & ,
                                 TDistanceTag const & )
{
    _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, false, TDir(), TDistanceTag());
}

}
#endif
