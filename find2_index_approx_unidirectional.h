#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_


using namespace std;

namespace seqan{


//TODO load bitvectors inside a struct to make accessing the correct bitvector easier
template <typename TText, typename TConfig, typename TIndexSpec,
          typename TVector, typename TVSupport,
          typename TDir,
          size_t nbrBlocks>
inline void get_bitvector_interval_inside(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                          vector<pair<TVector, TVSupport>> & bitvectors,
                                          OptimalSearch<nbrBlocks> const & s,
                                          uint8_t const blockIndex,
                                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brangeOutput,
                                          TDir const & )
{
    Pair<uint32_t, uint32_t> dirrange = range(iter);
    uint8_t needed_bitvector;
    uint8_t size = s.pi.size();
    uint8_t bitvsize = bitvectors.size();
    if (std::is_same<TDir, Rev>::value)
        needed_bitvector = bitvsize - s.max[blockIndex - 1];
    else
        needed_bitvector = s.min[blockIndex - 1] - 1;

    uint32_t nseq = countSequences(*iter.index);
    dirrange.i1 = dirrange.i1 - nseq;
    dirrange.i2 = dirrange.i2 - nseq;
    brangeOutput.i1 = needed_bitvector;
    brangeOutput.i2 = dirrange;
}

//search on unidirectional reverse genome
//sa_info is not const
//TODO calculate the EndPos of the needle maybe than it is easier to merge with the other function
template <typename TDelegateD,
          typename TNeedle,
          size_t nbrBlocks>
inline void genomeSearch(TDelegateD & delegateDirect,
                  TNeedle const & needle,
                  uint8_t errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  auto const & rgenome,
                  Pair<uint16_t, uint32_t> & sa_info,
//                   vector<uint32_t> const & blockStarts,
//                   vector<uint32_t> const & blockEnds,
                  std::array<uint32_t, nbrBlocks> & blockStarts,
                  std::array<uint32_t, nbrBlocks> & blockEnds,
                  bool const unidirectionalOnReverseIndex)
{
    uint32_t needleL = length(needle);
    for(uint8_t j = 0; j < nbrBlocks - blockIndex; ++j){
        for(uint32_t k = blockStarts[j]; k < blockEnds[j]; ++k){
            if(needle[needleL - k - 1] != rgenome[sa_info.i1][sa_info.i2 + k])
                ++errors;
        }
        if(errors < s.l[blockIndex + j] || errors > s.u[blockIndex + j]){
            return;
        }
    }
    sa_info.i2 = seqan::length(rgenome[sa_info.i1]) - sa_info.i2 - needleL;
    delegateDirect(sa_info, needle, errors);
}


template <typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void uniDirectSearch(TDelegateD & delegateDirect,
                  Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  vector<pair<TVector, TVSupport>> & bitvectors,
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                  TDir const & /**/)
{

    //this can also be the reverse genome
    auto const & genome = indexText(*iter.index);
    uint32_t needleL = length(needle);
    uint32_t blocks = s.pi.size();

    if(std::is_same<TDir, Rev>::value){

        std::array<uint32_t, nbrBlocks> blockStarts;
        std::array<uint32_t, nbrBlocks> blockEnds;
        std::copy(std::begin(s.revblockStarts) + blockIndex, std::end(s.revblockStarts), std::begin(blockStarts));
        std::copy(std::begin(s.revblockEnds) + blockIndex, std::end(s.revblockEnds), std::begin(blockEnds));

        for(uint32_t i = 0; i < brange.i2.i2 - brange.i2.i1; ++i){
            //this time i use the mappability from "inside" the needle since i can garantue i am at a blockend
            if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                Pair<uint16_t, uint32_t> sa_info = iter.index->sa[iter.vDesc.range.i1 + i];
                uint32_t chromlength = length(genome[sa_info.i1]);
                // mappability information is this time in reverse index order even if we use reverse index (we get_bitvector_interval_inside)
                //check left chromosom boundry && check right chromosom boundry
                if(!(sa_info.i2 >= needleL - needleRightPos + 1 && chromlength - 1 >= sa_info.i2 + needleRightPos - 2))
                    continue;
                sa_info.i2 = sa_info.i2 - needleL + needleRightPos - 1;
                //use modified genomeSearch in case of reverse index
                genomeSearch(delegateDirect, needle, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds, true);
            }
        }
    }
    else
    {
        std::array<uint32_t, nbrBlocks> blockStarts;
        std::array<uint32_t, nbrBlocks> blockEnds;
        std::copy(std::begin(s.blockStarts) + blockIndex, std::end(s.blockStarts), std::begin(blockStarts));
        std::copy(std::begin(s.blockEnds) + blockIndex, std::end(s.blockEnds), std::begin(blockEnds));

        //this time i use the mappability from "inside" the needle since i can garantue i am at a blockend
        for(uint32_t i = 0; i < brange.i2.i2 - brange.i2.i1; ++i){
            if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                Pair<uint16_t, uint32_t> sa_info = iter.index->sa[iter.vDesc.range.i1 + i];
                uint32_t chromlength = length(genome[sa_info.i1]);
                // mappability information is this time in reverse index order even if we use reverse index (we get_bitvector_interval_inside)
                //check left chromosom boundry && check right chromosom boundry
                if(!(needleLeftPos <= sa_info.i2 && chromlength - 1 >= sa_info.i2 - needleLeftPos + needleL - 1))
                    continue;
                //calculate correct starting position of the needle  on the forward index
                sa_info.i2 = sa_info.i2 - needleLeftPos;
                //use modified genomeSearch in case of forward index
                genomeSearch(delegateDirect, needle, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
            }
        }
    }
}

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeChildren(TDelegate & delegate,
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
                                         TDir const & /**/)
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
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev());
                else
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd());
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, TDir());
            }
        } while (goRight(iter));
    }
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeExact(TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      vector<pair<TVector, TVSupport>> & bitvectors,
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const & /**/)
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
            if (!goDown(iter, needle[infixPosRight]))
                return;
            --infixPosRight;
        }

        if (goToRight2)
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd());
    }
}

/*
inline ReturnCode checkInterval(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          uint8_t const blockSize,
                          bool const done,
                          bool const true,
                          uint8_t const blockIndex)
{
    cout << "unicheckInterval no start uni" << endl;
        TVector & b = bitvectors[brange.i1].first;
        TVSupport & rb = bitvectors[brange.i1].second;
        rb.set_vector(&b);

        uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
        if(params.uni.nomappability && ivalOne == 0)
            return ReturnCode::NOMAPPABILITY;

        if(!done){
            if(params.uni.directsearch && ivalOne < (blockSize - blockIndex - 1 + params.uni.directsearchblockoffset) * params.uni.directsearch_th){ //<4
                return ReturnCode::DIRECTSEARCH;
            }
            if(params.uni.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
                return ReturnCode::COMPMAPPABLE;
        }
        return ReturnCode::MAPPABLE;
}
*/

template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
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

    bool done = minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1;
    if(blockIndex > 0 && done || needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex - 1]){
        ReturnCode rcode = uniCheckMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, done, true, TDir());
        if(rcode == ReturnCode::FINISHED)
            return;
    }

    // Exact search in current block.
    if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
    }
    // Approximate search in current block.
    else
    {
    _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir());
    }
}

}
#endif
