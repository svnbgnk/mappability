#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_

#include <seqan/modifier.h>

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
          typename TString,
          typename TNeedle>
inline void alignmentMyersBitvector(TDelegateD & delegateDirect,
                                    TNeedle const & needle,
                                    TString const & n_infix,
                                    TString const & ex_infix,
                                    auto const & genome,
                                    Pair<uint16_t, uint32_t> const & sa_info,
                                    uint8_t max_e,
                                    uint8_t overlap_l,
                                    uint8_t overlap_r)
{
    uint16_t needleL = length(needle);
    uint16_t ex_infixL = needleL + overlap_l + overlap_r;

    int initialScore = globalAlignmentScore(ex_infix, needle, MyersBitVector());
    //assume more Insertions (in the read) than deletions
    int ins_initialScore = globalAlignmentScore(n_infix, needle, MyersBitVector());

    if(ins_initialScore >= 0 - 2 * max_e || initialScore >= 0 - overlap_l - overlap_r + max_e) //MM creates one error D creates one error since now it also align to overlap
    {
        //No Insertions or Deletions
        TString const & tmp0 = infix(ex_infix, overlap_l, ex_infixL - overlap_r);
        int errors2 = 0 - globalAlignmentScore(tmp0, needle, MyersBitVector()); //
        if(errors2 <= max_e)
            delegateDirect(sa_info , needle, errors2);

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
                        if(overlap_l < (pos * k) || 0 - (pos * (m - k)) > overlap_r)
                            continue;
                        sa_info_tmp = sa_info;
                        sa_info_tmp.i2 = sa_info_tmp.i2 + (pos * k);
                        TString const & tmp2 = infix(ex_infix, overlap_l + (pos * k), ex_infixL - overlap_r - (pos * (m - k)));
                        errors2 = 0 - globalAlignmentScore(tmp2, needle, MyersBitVector());
                        if(errors2 <= max_e)
                            delegateDirect(sa_info_tmp , needle, errors2);
                    }
                }
                else
                {
                    //insertions left and deletion right
                    if(overlap_l >= del){
                        TString const & tmp = infix(ex_infix, overlap_l - del, ex_infixL - overlap_r - ins);
                        sa_info_tmp.i2 = sa_info_tmp.i2 - del;
                        errors2 = 0 - globalAlignmentScore(tmp, needle, MyersBitVector());
                        if(errors2 <= max_e)
                            delegateDirect(sa_info_tmp , needle, errors2);
                    }

                    //insertions right and deletion left
                    if(overlap_r >= del){
                        sa_info_tmp = sa_info; //just include del from before into the calculation
                        TString const & tmp1 = infix(ex_infix, overlap_l + ins, ex_infixL - overlap_r + del);
                        errors2 = 0 - globalAlignmentScore(tmp1, needle, MyersBitVector());
                        sa_info_tmp.i2 = sa_info_tmp.i2 + ins;
                        if(errors2 <= max_e)
                            delegateDirect(sa_info_tmp , needle, errors2);
                    }
                }
            }
        }
    }
}

template <typename TDelegateD,
          typename TString, typename TConf, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void uniDirectSearch(TDelegateD & delegateDirect,
                  Iter<Index<StringSet<TString, TConf>, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  vector<pair<TVector, TVSupport>> & bitvectors,
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
        uint16_t needleL = length(needle);
        uint8_t max_e = s.u[s.u.size() - 1];
        uint8_t overlap_l = max_e;
        uint8_t overlap_r = max_e;
        if(needleLeftPos == 0)
            overlap_l = errors;
        if(needleRightPos == needleL + 1)
            overlap_r = errors;
        uint16_t ex_infixL = needleL + overlap_l + overlap_r;


        if(std::is_same<TDir, Rev>::value){
            for(uint32_t r = 0; r < brange.i2.i2 - brange.i2.i1; ++r){
                if(bitvectors[brange.i1].first[brange.i2.i1 + r] == 1){
                    Pair<uint16_t, uint32_t> sa_info = iter.index->sa[r];
                    uint32_t const chromlength = length(genome[sa_info.i1]);
                    if(!(sa_info.i2 >= needleL - needleRightPos + 1 && chromlength - 1 >= sa_info.i2 + needleRightPos - 2))
                        continue;
                    sa_info.i2 = sa_info.i2 - needleL + needleRightPos - 1;
                    //TODO fix this mess
                    TString /*const &*/ n_infix = infix(genome[sa_info.i1], sa_info.i2, sa_info.i2 + needleL);
                    TString /*const &*/ ex_infix = infix(genome[sa_info.i1], sa_info.i2 - overlap_l, sa_info.i2 + needleL + overlap_r);
                    DnaStringReverse n_infix_rev(n_infix);
                    DnaStringReverse ex_infix_rev(ex_infix);
                    TString n_infix_rev_copy = n_infix_rev;
                    TString ex_infix_rev_copy = ex_infix_rev;

                    alignmentMyersBitvector(delegateDirect, needle, n_infix_rev_copy, ex_infix_rev_copy, genome, sa_info, max_e, overlap_l, overlap_r);
                }
            }
        }else{
            for(uint32_t r = 0; r < brange.i2.i2 - brange.i2.i1; ++r){

                if(bitvectors[brange.i1].first[brange.i2.i1 + r] == 1){
                    Pair<uint16_t, uint32_t> sa_info = iter.index->sa[r];
                    uint32_t const chromlength = length(genome[sa_info.i1]);
                    if(!(needleLeftPos + overlap_l <= sa_info.i2  && chromlength - 1 >= sa_info.i2 - needleLeftPos + needleL - 1 + overlap_r))
                        continue;
                     sa_info.i2 = sa_info.i2 - needleLeftPos;
                    TString const & ex_infix = infix(genome[sa_info.i1], sa_info.i2 - overlap_l, sa_info.i2 + needleL + overlap_r);
                    TString const & n_infix = infix(genome[sa_info.i1], sa_info.i2, sa_info.i2 + needleL);
                    alignmentMyersBitvector(delegateDirect, needle, n_infix, ex_infix, genome, sa_info, max_e, overlap_l, overlap_r);

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
}

template <typename TDelegate,
          typename TDelegateDirect,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletion(TDelegate & delegate,
                                         TDelegateDirect & delegateDirect,
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

    if (minErrorsLeftInBlock == 0)
    {

        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

        if (goToRight2)
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, Rev(), EditDistance());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, Fwd(), EditDistance());

    }

    if (maxErrorsLeftInBlock > 0 && goDown(iter))
    {
        do
        {
            _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, TDir());
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
                    _optimalSearchSchemeDeletion(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

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
                _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, TDir(), TDistanceTag());
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
            if (!goDown(iter, needle[infixPosRight]))
                return;
            --infixPosRight;
        }

        if (goToRight2)
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd(), TDistanceTag());
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
          typename TDir,
          typename TDistanceTag>
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
                                 TDir const & /**/,
                                 TDistanceTag const &)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    bool done = minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1;
    //TODO is blockIndex > 0 necessary
    if(blockIndex > 0 && done || needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex - 1]){
        ReturnCode rcode = uniCheckMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, done, true, TDir(), TDistanceTag());
        if(rcode == ReturnCode::FINISHED)
            return;
    }
    //cant use else if here since checking mappability at the blockend can lead to mappable then the code below should executed

    // Exact search in current block.
    if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
    }
    // Approximate search in current block.
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
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

                    if (goToRight2)
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, Rev(), TDistanceTag());
                    else
                        _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, Fwd(), TDistanceTag());
                }
            }
            else
            {
                _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
            }
        }
    }
}

}
#endif
