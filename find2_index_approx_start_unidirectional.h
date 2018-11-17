#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_START_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_START_UNIDIRECTIONAL_EXTENSION_H_

#include "common_auxiliary.h"

namespace seqan{


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void filter_interval(TDelegate & delegate,
                     TDelegateD & delegateDirect,
                     Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                     TNeedle const & needle,
                     vector<pair<TVector, TVSupport>> & bitvectors,
                     uint32_t const needleLeftPos,
                     uint32_t const needleRightPos,
                     uint8_t const errors,
                     OptimalSearch<nbrBlocks> const & s,
                     uint8_t const blockIndex,
                     Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                     TDir const & /**/,
                     TDistanceTag const &,
                     bool & successfulOutput)
{
    vector<pair<uint32_t, uint32_t>> consOnes;
    getConsOnes(bitvectors, inside_bit_interval, params.startuni.intervalsize, consOnes);
    if(consOnes.size() == 1){
        successfulOutput = false;
        return;
    }
    uint32_t nseq = countSequences(*iter.index);
    for(uint32_t i = 0; i < consOnes.size(); ++i){
        iter.vDesc.range.i1 = consOnes[i].first + nseq;
        iter.vDesc.range.i2 = consOnes[i].second + nseq;
        if (std::is_same<TDir, Rev>::value)
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Rev(), TDistanceTag());
        else
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Fwd(), TDistanceTag());
    }
    successfulOutput = true;
}


//TODO merge with old funtion (test unidirectional) diff: iter and paramenter filpDensity
template<typename TText, typename TConfig, typename TIndexSpec,
         typename TVector, typename TVSupport,
         size_t nbrBlocks,
         typename TDir>
inline bool testFilter(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                          vector<pair<TVector, TVSupport>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex,
                          TDir const & )
{
    // need bitinterval from inside the pattern to filter according to the mappability form
    //therefore i also need to acces the block before because of that block i got mappability of both sides
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
    get_bitvector_interval_inside(iter, bitvectors, s, blockIndex, bit_interval, TDir());
    TVector & b2 = bitvectors[bit_interval.i1].first;

    //squash interval
    uint32_t startPos = bit_interval.i2.i1, endPos = bit_interval.i2.i2;

    while(b2[startPos] == 0 && startPos < endPos)
        ++startPos;

    while(b2[endPos - 1] == 0 && endPos > startPos)
        --endPos;

    if(startPos > endPos){
        std::cout << "Error bit vector has only zeroes this should have been checked by checkinterval" << "\n";
        std::cout << "Size: " << endPos - startPos << endl;
        exit(0);
    }


    float ivalSize = brange.i2.i2 - brange.i2.i1;
    uint32_t count = 0;

    if(params.startuni.testflipdensity){
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
    if(ivalSize * params.startuni.invflipdensity - 1 > static_cast<float>(count)){
        brange.i1 = bit_interval.i1;
        brange.i2.i1 = startPos;
        brange.i2.i2 = endPos;
        return true;
    }
    return false;
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _uniOptimalSearchSchemeChildren(TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         vector<pair<TVector, TVSupport>> & bitvectors,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const & /**/,
                                         TDistanceTag const &)
{
    bool goToRight = std::is_same<TDir, Rev>::value;
    if (goDown(iter))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter),
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

                //TODO remove goToRight2 and input should be TDIR
                if (goToRight2)
                {
                    _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev(), TDistanceTag());
                }
                else
                {
                    _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd(), TDistanceTag());
                }
            }
            else
            {
                _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, TDir(), TDistanceTag());
            }
        } while (goRight(iter));
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _uniOptimalSearchSchemeExact(TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      vector<pair<TVector, TVSupport>> & bitvectors,
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const & /**/,
                                      TDistanceTag const &)
{
    bool goToRight2 = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] : s.pi[blockIndex] > s.pi[blockIndex - 1];
    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it forward
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1)))
            return;

        if (goToRight2)
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Rev(), TDistanceTag());
        else
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Fwd(), TDistanceTag());
    }
    else
    {
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        int32_t infixPosLeft = needleRightPos - s.blocklength[blockIndex] - 1;
        int32_t infixPosRight = needleLeftPos - 1;

        while (infixPosRight >= infixPosLeft)
        {
            if (!goDown(iter, needle[infixPosRight]))
                return;
            --infixPosRight;
        }

        if (goToRight2)
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev(), TDistanceTag());
        else
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd(), TDistanceTag());
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _uniOptimalSearchScheme(TDelegate & delegate,
                                    TDelegateD & delegateDirect,
                                    Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,                                 TNeedle const & needle,
                                    vector<pair<TVector, TVSupport>> & bitvectors,
                                    uint32_t const needleLeftPos,
                                    uint32_t const needleRightPos,
                                    uint8_t const errors,
                                    OptimalSearch<nbrBlocks> const & s,
                                    uint8_t const blockIndex,
                                    TDir const & /**/,
                                    TDistanceTag const &)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    bool done = minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1;
    if(blockIndex > 0 && done || needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex - 1]){
        ReturnCode rcode = uniCheckMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, done, false, TDir(), TDistanceTag());
        if(rcode == ReturnCode::FINISHED)
            return;
    }

    // Exact search in current block.
    if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _uniOptimalSearchSchemeExact(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
    }
    // Approximate search in current block.
    else
    {
    _uniOptimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
    }

}

}

#endif