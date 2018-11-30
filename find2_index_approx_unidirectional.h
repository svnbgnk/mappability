#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_


using namespace std;

namespace seqan{

template <typename TText, typename TConfig, typename TIndexSpec,
          typename TVector, typename TVSupport,
          typename TDir,
          size_t nbrBlocks>
inline void get_bitvector_interval_inside(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                          std::vector<std::pair<TVector, TVSupport>> & bitvectors,
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
template <typename TContex,
          typename TDelegateD,
          typename TNeedle,
          size_t nbrBlocks>
inline void genomeSearchRev(TContex & ossContext,
                            TDelegateD & delegateDirect,
                            TNeedle const & needle,
                            uint32_t needleId,
                            uint8_t errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            auto const & rgenome,
                            Pair<uint16_t, uint32_t> & sa_info,
                            std::array<uint32_t, nbrBlocks> & blockStarts,
                            std::array<uint32_t, nbrBlocks> & blockEnds)
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
    delegateDirect(sa_info, posAdd(sa_info, length(needle)), needle, needleId, errors);
}

template <typename TContex,
          typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void uniDirectSearch(TContex & ossContext,
                            TDelegateD & delegateDirect,
                            Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                            TNeedle const & needle,
                            uint32_t needleId,
                            std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                            uint32_t const needleLeftPos,
                            uint32_t const needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                            TDir const & /**/,
                            TDistanceTag const & )
{

    //this can also be the reverse genome
    auto const & genome = indexText(*iter.index);
    uint32_t needleL = length(needle);

    if(std::is_same<TDistanceTag, EditDistance>::value){
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

        if(std::is_same<TDir, Rev>::value){
            TNeedle needle2 = needle;
            DnaStringReverse needleRev(needle2);
            TNeedle needleRevCopy = needleRev; //why is this neccesary MyersBitVector() otherwise return error

            for(uint32_t r = 0; r < iter.vDesc.range.i2 - iter.vDesc.range.i1; ++r){
                if(checkSinglePos(bitvectors, brange, r))
                {
                    Pair<uint16_t, uint32_t> sa_info = iter.index->sa[iter.vDesc.range.i1 + r];
                    uint32_t const chromlength = length(genome[sa_info.i1]);
                    if(!(sa_info.i2 >= needleL - needleRightPos + 1 + overlap_l && chromlength - 1 >= sa_info.i2 + needleRightPos - 2 + overlap_l))
                        continue;
                    sa_info.i2 = sa_info.i2 - needleL + needleRightPos - 1;
                    //since we reverse the needle the overlap in infices must also me switched
                    DnaString ex_infix = infix(genome[sa_info.i1], sa_info.i2 - overlap_r, sa_info.i2 + needleL + overlap_l);
                    DnaString n_infix = infix(genome[sa_info.i1], sa_info.i2, sa_info.i2 + needleL);

                    alignmentMyersBitvector(ossContext, delegateDirect, needleRevCopy, needleId, n_infix, ex_infix, chromlength, sa_info, max_e, overlap_l, overlap_r, intDel, true);
                }
            }
        }
        else
        {
            for(uint32_t r = 0; r < iter.vDesc.range.i2 - iter.vDesc.range.i1; ++r){
                if(checkSinglePos(bitvectors, brange, r)){
                    Pair<uint16_t, uint32_t> sa_info = iter.index->sa[iter.vDesc.range.i1 + r];
                    uint32_t const chromlength = length(genome[sa_info.i1]);
                    if(!(needleLeftPos <= sa_info.i2 - overlap_l  && chromlength - 1 >= sa_info.i2 - needleLeftPos + needleL - 1 + overlap_r))
                        continue;
                    //calculate correct starting position of the needle  on the forward index
                    sa_info.i2 = sa_info.i2 - needleLeftPos;

                    DnaString const & ex_infix = infix(genome[sa_info.i1], sa_info.i2 - overlap_l, sa_info.i2 + needleL + overlap_r);
                    DnaString const & n_infix = infix(genome[sa_info.i1], sa_info.i2, sa_info.i2 + needleL);

                    alignmentMyersBitvector(ossContext, delegateDirect, needle, needleId, n_infix, ex_infix, chromlength, sa_info, max_e, overlap_l, overlap_r, intDel, false);
                }
            }
        }
    }
    else
    {
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
                    genomeSearchRev(ossContext, delegateDirect, needle, needleId, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
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
                    genomeSearch(ossContext, delegateDirect, needle, needleId, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
                }
            }
        }
    }
}

template <typename TContex,
          typename TVector, typename TVSupport>
inline ReturnCode uniCheckInterval(TContex & ossContext,
                                   std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                   Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                                   uint8_t const blockSize,
                                   bool const done,
                                   bool const nofilter,
                                   uint8_t const blockIndex)
{
    TVector & b = bitvectors[brange.i1].first;
    TVSupport & rb = bitvectors[brange.i1].second;
    rb.set_vector(&b);

    uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
    if(ossContext.uni.nomappability && ivalOne == 0)
        return ReturnCode::NOMAPPABILITY;

    if(!done){
        if(ossContext.uni.directsearch && ivalOne < (blockSize - blockIndex - 1 + ossContext.uni.directsearchblockoffset) * ossContext.uni.directsearch_th)
            return ReturnCode::DIRECTSEARCH;
    }
    //this is in the moment used to determine if we have to filter the delegate Call
    if(/*ossContext.uni.compmappable && */ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
        return ReturnCode::COMPMAPPABLE;

    return ReturnCode::MAPPABLE;
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode uniCheckMappability(TContex & ossContext,
                                      TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      uint32_t needleId,
                                      std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      bool const done,
                                      bool const nofilter,
                                      TDir const & ,
                                      TDistanceTag const &)
{
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
    get_bitvector_interval_inside(iter, bitvectors, s, blockIndex + done, bit_interval, TDir());
    ReturnCode rcode = uniCheckInterval(ossContext, bitvectors, bit_interval, s.pi.size(), done, nofilter, blockIndex);

    if(rcode == ReturnCode::NOMAPPABILITY)
        return ReturnCode::FINISHED;
    // Done. (Last step)
    if (done)
    {
        bool rev = std::is_same<TDir, Rev>::value;
        //NOTE we only need to check mappability here to not allow occurences with a frequency higher than allowed
        //if the only concern is speed than reporting all occurrences would be faster
        if(rcode != ReturnCode::COMPMAPPABLE /* add option for speed up here*/)
        {
            uint32_t rangeStart = iter.vDesc.range.i1;
            uint32_t rangeEnd = iter.vDesc.range.i2;
            uint32_t lastStart = 0;
            for(uint32_t i = 0; i < rangeEnd - rangeStart; ++i)
            {
                if(bitvectors[bit_interval.i1].first[bit_interval.i2.i1 + i] == 0 )
                {
                    if(i != lastStart){
                        iter.vDesc.range.i1 = rangeStart + lastStart;
                        iter.vDesc.range.i2 = rangeStart + i - 1;
//                         cout << iter.vDesc.range.i1 << " - " << iter.vDesc.range.i2;
                        delegate(iter, needle, needleId, errors, rev);
                    }
                    lastStart = i + 1;
                }
            }
            if(lastStart < rangeEnd - rangeStart){
                iter.vDesc.range.i1 = rangeStart + lastStart;
                iter.vDesc.range.i2 = rangeStart + rangeEnd - rangeStart;
                delegate(iter, needle, needleId, errors, rev);
            }
        }
        else
        {
            delegate(iter, needle, needleId, errors, rev);
        }
        return ReturnCode::FINISHED;
    }

    if(rcode == ReturnCode::DIRECTSEARCH){
        uniDirectSearch(ossContext, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
        return ReturnCode::FINISHED;
    }
    return ReturnCode::MAPPABLE;
}

template <typename TContex,
          typename TDelegate,
          typename TDelegateDirect,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletion(TContex & ossContext,
                                         TDelegate & delegate,
                                         TDelegateDirect & delegateDirect,
                                         Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
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

    if (/*not_at_end && */maxErrorsLeftInBlock > 0 && goDown(iter))
    {
        do
        {
            _optimalSearchSchemeDeletion(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, true, TDir());
        } while (goRight(iter));
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeChildren(TContex & ossContext,
                                         TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         uint32_t needleId,
                                         std::vector<std::pair<TVector, TVSupport>> & bitvectors,
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
            bool delta = !ordEqual(parentEdgeLabel(iter), needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1]);
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
        } while (goRight(iter));
    }
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeExact(TContex & ossContext,
                                      TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      uint32_t needleId,
                                      std::vector<std::pair<TVector, TVSupport>> & bitvectors,
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
        //search take rest of the block and search it reverse
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1)))
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
            if (!goDown(iter, needle[infixPosRight]))
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
          typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 bool const lastEdit,
                                 TDir const & /**/,
                                 TDistanceTag const &)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
    bool done = minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1;
    bool const atBlockEnd = (blockIndex > 0) ? needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex - 1] : false;

    bool rev = std::is_same<TDir, Rev>::value;
    if(done){
        if(true /*!lastEdit*/){
            delegate(iter, needle, needleId, errors, rev);
        }
        return;
    }

    if(done || atBlockEnd){
        ReturnCode rcode = uniCheckMappability(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, done, true, TDir(), TDistanceTag());
        if(rcode == ReturnCode::FINISHED)
            return;
    }

    // Exact search in current block.
    if (maxErrorsLeftInBlock == 0)
    {
        _optimalSearchSchemeExact(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
    }
    //TODO do ITV search here if we are COMPMAPPABLE (check mappability only very few times so ...)
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

        _optimalSearchSchemeChildren(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}

}
#endif





