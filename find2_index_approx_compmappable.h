#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_COMPMAPPLE_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_COMPMAPPLE_H_

#include <seqan/align.h>

using namespace std;

namespace seqan{

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
                                         TDir const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    if (minErrorsLeftInBlock == 0)
    {
        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
        if (goToRight2)
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex2, Rev(), EditDistance());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex2, Fwd(), EditDistance());
    }

    if (maxErrorsLeftInBlock > 0 && goDown(iter, TDir()))
    {
        do
        {
            _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors + 1, s,
                                         blockIndex, TDir());
        } while (goRight(iter, TDir()));
    }
}

//
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
    uint16_t needleL = length(needle);
    uint8_t max_e = s.u[s.u.size() - 1];

    uint8_t overlap_l = max_e;
    uint8_t overlap_r = max_e;
    if(needleLeftPos == 0)
        overlap_l = errors;
    if(needleRightPos == length(needle) + 1)
        overlap_r = errors;

    uint8_t rest_errors = max_e - errors;
    if (std::is_same<TDistanceTag, EditDistance>::value){
        for(uint32_t r = iter.fwdIter.vDesc.range.i1; r < iter.fwdIter.vDesc.range.i2; ++r)
        {
            Pair<uint16_t, uint32_t> sa_info = iter.fwdIter.index->sa[r];
            uint32_t const chromlength = length(genome[sa_info.i1]);
            //TODO use new limits everywhere for EditDistance
            if(!(needleLeftPos + overlap_l <= sa_info.i2  && chromlength - 1 >= sa_info.i2 - needleLeftPos + needleL - 1 + overlap_r))
                continue;
            sa_info.i2 = sa_info.i2 - needleLeftPos;

            TString const & ex_needle = infix(genome[sa_info.i1], sa_info.i2 - overlap_l, sa_info.i2 + needleL + overlap_r);
            TString const & n_infix = infix(genome[sa_info.i1], sa_info.i2, sa_info.i2 + needleL);

            //Unfortunately, since for example (E = 2) 2 Insertions lead to a score of -4 we cannot conclude the number of errors the read will have at the end. (O Errors would also have a score -4 2 Deletion in the beginning to Insertions at the end)
            //for E = 2 we allow a score -6 since 2 MM + 2 D + 2 I.
            //To Insertion on the otherhand lead to 2D + 2I + 4D -> 8 errors

            int initialScore = globalAlignmentScore(ex_needle, needle, MyersBitVector());

            //assume more Insertions (in the read) than deletions
            int ins_initialScore = globalAlignmentScore(n_infix, needle, MyersBitVector());

            //TODO switch Insertion and deltions in the following code!!
            if(ins_initialScore >= 0 - 2 * max_e || initialScore >= 0 - overlap_l - overlap_r + max_e) //MM creates one error D creates one error since now it also align to overlap
            {
                //No Insertions or Deletions
                //TODO use length(ex_needle) to make it easier
                TString const & tmp0 = infix(ex_needle, overlap_l, needleL + overlap_l);
                int errors2 = 0 - globalAlignmentScore(tmp0, needle, MyersBitVector()); //
                if(errors2 <= max_e)
                    delegateDirect(sa_info , needle, errors2);

                for(uint8_t e = 1; e <= max_e /*overlap*/; ++e){
//                     cout << "E: " << (int)e << endl;
                    for(uint32_t i = 0; i <= e; ++i){
                        //i is number of insertions
                        uint32_t d = e - i; //number of deletions
                        auto sa_info_tmp = sa_info;

                        if(i > 1 && d == 0 || d > 1 && i == 0){
                        //only insertion or deletions
                            int pos = (d > i) ? 1 : (-1);
                            int m = std::max(i,d);
                            for(int k = 0; k <= m; ++k)
                            {
                                if(overlap_l < (pos * k) || 0 - (pos * (m - k)) > overlap_r)
                                    continue;
                                sa_info_tmp = sa_info;
                                sa_info_tmp.i2 = sa_info_tmp.i2 + (pos * k);
                                TString const & tmp2 = infix(ex_needle, overlap_l + (pos * k), needleL + overlap_l - (pos * (m - k)));
                                errors2 = 0 - globalAlignmentScore(tmp2, needle, MyersBitVector());
                                if(errors2 <= max_e)
                                    delegateDirect(sa_info_tmp , needle, errors2);
                            }
                        }
                        else
                        {
                            if(overlap_l >= i){

                            //insertions left and deletion right
                            TString const & tmp = infix(ex_needle, overlap_l - i, needleL + overlap_l - d);
                            sa_info_tmp.i2 = sa_info_tmp.i2 - i;
                            errors2 = 0 - globalAlignmentScore(tmp, needle, MyersBitVector());
                            if(errors2 <= max_e)
                                delegateDirect(sa_info_tmp , needle, errors2);
                            }
                            if(overlap_r >= i){
                            //insertions right and deletion left
                            sa_info_tmp = sa_info; //TODO just include i from before into the calculation
                            TString const & tmp1 = infix(ex_needle, overlap_l + d, needleL + overlap_l + i);
                            errors2 = 0 - globalAlignmentScore(tmp1, needle, MyersBitVector());
                            sa_info_tmp.i2 = sa_info_tmp.i2 + d;
                            if(errors2 <= max_e)
                                delegateDirect(sa_info_tmp , needle, errors2);
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
    std::vector<uint32_t> blockStarts(s.pi.size() - blockIndex);
    std::vector<uint32_t> blockEnds(s.pi.size() - blockIndex);
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
                    _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                    {
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev(), TDistanceTag());
                    }
                    else
                    {
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd(), TDistanceTag());
                    }
                }
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,
                                     blockIndex, TDir(), TDistanceTag());
            }

            // Deletion
            if (std::is_same<TDistanceTag, EditDistance>::value)
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, TDir(), TDistanceTag());
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
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                 blockIndex2, Rev(), TDistanceTag());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                 blockIndex2, Fwd(), TDistanceTag());
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
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                 blockIndex2, Rev(), TDistanceTag());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                 blockIndex2, Fwd(), TDistanceTag());
        }
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
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    // Done.
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        delegate(iter, needle, errors, false);
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
        if (std::is_same<TDistanceTag, EditDistance>::value)
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
                    uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, Rev(), EditDistance());
                    else
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, Fwd(), EditDistance());
                }
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, TDir(), TDistanceTag());
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


}
#endif
