#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_COMPMAPPLE_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_COMPMAPPLE_H_
using namespace std;

namespace seqan{

template <typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void directSearch(TDelegateD & delegateDirect,
                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  TDir const & /**/)
{
    auto const & genome = indexText(*iter.fwdIter.index);
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
                if (goToRight2)
                {
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,
                                         blockIndex2, Rev(), TDistanceTag());
                }
                else
                {
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,
                                         blockIndex2, Fwd(), TDistanceTag());
                }
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,
                                     blockIndex, TDir(), TDistanceTag());
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
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev(), TDistanceTag());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd(), TDistanceTag());
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
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev(), TDistanceTag());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd(), TDistanceTag());
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
                                 // TDelegateTV & doInTextVerification,
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
        // if (doInTextVerification(s, blockIndex, iter))
        if(params.comp.directsearch && iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1 + params.uni.directsearchblockoffset) * params.comp.directsearch_th)
        {
            directSearch(delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
            return;
        }
        _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex,
                                     minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}


}
#endif
