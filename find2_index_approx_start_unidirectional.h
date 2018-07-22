#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_START_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_START_UNIDIRECTIONAL_EXTENSION_H_

#include "common_auxiliary.h"

namespace seqan{

    
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
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
                     Pair<uint8_t, Pair<uint16_t, uint32_t>> & inside_bit_interval,
                     TDir const & /**/,
                     bool & successfulOutput)
{
    vector<std::pair<uint32_t, uint32_t>> consOnes;
    getConsOnes(bitvectors, inside_bit_interval, params.startuni.intervalsize, consOnes);
    if(consOnes.size() == 1){
        successfulOutput = false;
        return;
    }
    //TODO replace with countSequences when it works
    uint32_t nseq = countSequences(*iter.index);
    for(int i = 0; i < consOnes.size(); ++i){
        iter.vDesc.range.i1 = consOnes[i].first + nseq;
        iter.vDesc.range.i2 = consOnes[i].second + nseq;
        //TODO Maybe Link to old Funtion to not filter multiple times?
        if (std::is_same<TDir, Rev>::value)
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Rev());
        else
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Fwd());
    }
    successfulOutput = true;
}

    
//TODO merge with old funtion (test unidirectional) diff: iter and pramenter filpDensity
template<typename TText, typename TConfig, typename TIndexSpec,
         typename TVector, typename TVSupport,
         size_t nbrBlocks,
         typename TDir>
inline bool testFilter(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                          vector<pair<TVector, TVSupport>> & bitvectors,
                          Pair<uint8_t, Pair<uint16_t, uint32_t>> & brange,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex,
                          TDir const & )
{
    // need bitinterval from inside the pattern to filter according to the mappability form
    //therefore i also need to acces the block before because of that block i got mappability of both sides
    Pair<uint8_t, Pair<uint16_t, uint32_t>> bit_interval;
    get_bitvector_interval_inside(iter, bitvectors, s, blockIndex, bit_interval, TDir());
    TVector & b2 = bitvectors[bit_interval.i1].first;
    
    //squash interval
    uint32_t startPos = bit_interval.i2.i1, endPos = bit_interval.i2.i2;
    
    while(b2[startPos] == 0 && startPos < endPos)
        ++startPos;

    while(b2[endPos - 1] == 0 && endPos > startPos)
        --endPos;
    
    if(startPos > endPos){
        cout << "Error bit vector has only zeroes this should have been checked by checkinterval" << endl;
        cout << "Size: " << endPos - startPos << endl;
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

template <typename TVector, typename TVSupport>
inline ReturnCode checkInterval(vector<pair<TVector, TVSupport>> & bitvectors,
                          Pair<uint8_t, Pair<uint16_t, uint32_t>> & brange,
                          uint8_t const blockSize,
                          bool const done,
                          bool const nofilter,
                          uint8_t const blockIndex)
{
    if(!nofilter){
        TVector & b = bitvectors[brange.i1].first;
        TVSupport & rb = bitvectors[brange.i1].second; 
        rb.set_vector(&b);
    
        uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
        if(params.startuni.nomappability && ivalOne == 0)
            return ReturnCode::NOMAPPABILITY;

        if(!done){
            if(params.startuni.directsearch && ivalOne < (blockSize - blockIndex - 1 + params.startuni.directsearchblockoffset) * params.startuni.directsearch_th){ //<4 
                return ReturnCode::DIRECTSEARCH;
            }    
            if(params.startuni.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
                return ReturnCode::COMPMAPPABLE;
        
            //equal or more than half zeroes     
            float ivalSize = brange.i2.i2 - brange.i2.i1;
            if(params.startuni.suspectunidirectional && ivalOne/ ivalSize <= params.startuni.filter_th){
                return ReturnCode::FILTER;
            }
        }
        return ReturnCode::MAPPABLE;
    }
    else
    {
        TVector & b = bitvectors[brange.i1].first;
        TVSupport & rb = bitvectors[brange.i1].second; 
        rb.set_vector(&b);
    
        uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
        if(params.uni.nomappability && ivalOne == 0)
            return ReturnCode::NOMAPPABILITY;

        if(!done){
            if(params.uni.directsearch && ivalOne < (blockSize - blockIndex - 1 + params.uni.directsearchblockoffset) * params.uni.directsearch_th){
                return ReturnCode::DIRECTSEARCH;
            }    
            if(params.uni.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
                return ReturnCode::COMPMAPPABLE;
        }
        return ReturnCode::MAPPABLE;
    }
    
    
    std::cerr << "Something went wrong!!!" << "\n";
    exit(0);
    return ReturnCode::MAPPABLE;
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir> 
inline ReturnCode uniCheckMappability(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 vector<pair<TVector, TVSupport>> & bitvectors,   
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,                             
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 bool const done,
                                 bool const nofilter,
                                 TDir const & )
{
    Pair<uint8_t, Pair<uint16_t, uint32_t>> bit_interval;
    get_bitvector_interval_inside(iter, bitvectors, s, blockIndex + done, bit_interval, TDir());
    ReturnCode rcode = checkInterval(bitvectors, bit_interval, s.pi.size(), done, nofilter, blockIndex);
    
    if(rcode == ReturnCode::NOMAPPABILITY)
        return ReturnCode::FINISHED;
    // Done. (Last step)
    if (done)
    {
        bool rev = std::is_same<TDir, Rev>::value;
        //TODO it is possible to allowed cheap repeats if always false 
        if(rcode == ReturnCode::MAPPABLE) 
        {
            uint32_t rangeStart = iter.vDesc.range.i1;
            uint32_t rangeEnd = iter.vDesc.range.i2;
            int lastStart = 0;
            for(int i = 0; i < rangeEnd - rangeStart; ++i)
            {
                if(bitvectors[bit_interval.i1].first[bit_interval.i2.i1 + i] == 0 )
                {
                    if(i != lastStart){
                        iter.vDesc.range.i1 = rangeStart + lastStart;
                        iter.vDesc.range.i2 = rangeStart + i - 1;
                        cout << iter.vDesc.range.i1 << " - " << iter.vDesc.range.i2;
                        delegate(iter, needle, errors, rev);
                    }
                    lastStart = i + 1;
                }
            }
            iter.vDesc.range.i1 = rangeStart + lastStart;
            iter.vDesc.range.i2 = rangeStart + rangeEnd - rangeStart;
            delegate(iter, needle, errors, rev);
        }
        else
        {
            delegate(iter, needle, errors, rev);
        }
        return ReturnCode::FINISHED; 
    }
    
    if(rcode == ReturnCode::DIRECTSEARCH){
        uniDirectSearch(delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir());
        return ReturnCode::FINISHED;
    }
    
    if(rcode == ReturnCode::FILTER){
        //test filter also modfied iter range if true;
        if(testFilter(iter, bitvectors, bit_interval, s, blockIndex, TDir())){
            bool successful;
            filter_interval(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir(), successful);
            if(successful)
                return(ReturnCode::FINISHED);
        }
    }
    return ReturnCode::MAPPABLE;
}
  
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
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
                {
                    _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev());
                }
                else
                {
                    _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd());
                }
            }
            else
            {
                _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, TDir());
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
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Rev());
        else
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Fwd());
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
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev());
        else
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd());
    }
}
    
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          typename TVector, typename TVSupport,
          size_t nbrBlocks,
          typename TDir>
inline void _uniOptimalSearchScheme(TDelegate & delegate,
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
        ReturnCode rcode = uniCheckMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, done, false, TDir());
        if(rcode == ReturnCode::FINISHED)
            return;
    }
    
    // Exact search in current block.
    if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _uniOptimalSearchSchemeExact(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
    }
    // Approximate search in current block.
    else
    {
    _uniOptimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir());
    }
    
}   
    
}  

#endif
