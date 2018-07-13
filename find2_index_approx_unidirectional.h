#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_

#include <sdsl/bit_vectors.hpp>
#include "common.h"
#include "common_auxiliary.h"

namespace seqan{

//TODO load bitvectors inside a struct to make accessing the correct bitvector easier
template <typename TText, typename TConfig, typename TIndexSpec,
          typename TDir,
          size_t nbrBlocks>
Pair<uint8_t, Pair<uint32_t, uint32_t>> get_bitvector_interval_inside(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                        vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,    
                                        OptimalSearch<nbrBlocks> const & s,
                                        uint8_t const blockIndex,
                                        TDir const & /**/) 
{
    cout << "Get Inside bitvector Interval unidirectional" << endl;
    cout << "bitvector_interval blockIndex: " << (int)blockIndex << endl;
    cout << ((std::is_same<TDir, Rev>::value) ? "reverse Index" : "forward index") << endl;
    Pair<uint32_t, uint32_t> dirrange = range(iter);
    uint8_t needed_bitvector;
    uint8_t size = s.pi.size();
    uint8_t bitvsize = bitvectors.size();
    cout << "max: " << (int)s.max[blockIndex - 1] << endl;
    cout << "bitvsize: " << (int)bitvsize << "end bitvector size" << endl;
    if (std::is_same<TDir, Rev>::value)
        needed_bitvector = bitvsize - s.max[blockIndex - 1];      
    else
        needed_bitvector = s.min[blockIndex - 1] - 1;


    uint32_t number_of_indeces = seqan::length(iter.index->sa) - bitvectors[needed_bitvector].first.size();
    cout << "selected bitvector inside: " << (int)needed_bitvector << endl;
//     cout << "Sa Range: " << dirrange << endl;  
    dirrange.i1 = dirrange.i1 - number_of_indeces;
    dirrange.i2 = dirrange.i2 - number_of_indeces;
//     cout << "Bit Range: " << dirrange << endl;  
    
    Pair<uint8_t, Pair<uint32_t, uint32_t>> brange(needed_bitvector, dirrange);
    return brange;
}
    
    
//TODO not used remove?
template <typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
void directSearch(TDelegateD & delegateDirect,
                  Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, 
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  Pair<uint8_t, Pair<uint32_t, uint32_t>> const & brange,
                  TDir const & /**/)
{
    cout << "directSearch unidirectional " <<  endl;
    cout << "NLP: " <<  needleLeftPos <<  endl;
    cout << "NRP: " <<  needleRightPos <<  endl;
    cout << "errors:  " <<  (int)errors <<  endl;
    vector<Pair<uint16_t, uint32_t>> hitsv;
    vector<uint8_t> errorsv;
    auto const & genome = indexText(*iter.index);
    for(int i = 0; i < brange.i2.i2 - brange.i2.i1; ++i){
        cout << "I: " << i << endl;
        if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
            cout << "blockIndex: " << (int)blockIndex << endl;
            uint8_t errors2 = errors;
            bool valid = true;
            Pair<uint16_t, uint32_t> sa_info;
            uint32_t startPos;
            // mappability information is in reverse index order if we use the forward index
            if(std::is_same<TDir, Rev>::value){
                cout << "rev" << endl;
                sa_info = iter.index->sa[iter.vDesc.range.i1 + i];
                startPos = sa_info.i2 - needleLeftPos;
            }
            else{
                cout << "fwd" << endl;
                sa_info = iter.index->sa[iter.vDesc.range.i1 + i];
                //calculate correct starting position of the needle  on the forward index
                startPos = seqan::length(genome[sa_info.i1]) - sa_info.i2 - needleRightPos + 1;
            }
            cout << "StartPos " << startPos << endl;
            //search remaining blocks
            for(int j = blockIndex; j < s.pi.size(); ++j){
                int blockStart = (s.pi[j] - 1 == 0) ? 0 : s.chronBL[s.pi[j] - 2];
                int blockEnd = s.chronBL[s.pi[j] - 1];
                cout << "searching Parts:" << blockStart << " - " << blockEnd << "; ";
                cout << endl;
                // compare bases to needle
                if(std::is_same<TDir, Rev>::value){
                    if(needleRightPos - 1 > blockStart && needleRightPos - 1 < blockEnd){
                        cout << "changing Blockstart, should only happen (once) when we come from checkcurrentmappability!!" << endl;
                        blockStart = needleRightPos - 1;
                        cout << "searching Parts:" << blockStart << " - " << blockEnd << "; ";
                    }
                }else{
                    if(needleLeftPos > blockStart && needleLeftPos < blockEnd){
                        cout << "changing Blockend fwd" << endl;
                        blockEnd = needleLeftPos;
                        cout << "searching Parts:" << blockStart << " - " << blockEnd << "; ";
                    }
                }
                for(int k = blockStart; k <  blockEnd; ++k){
                    //                     if(needle[k] != it)
                    if(needle[k] != genome[sa_info.i1][startPos + k])
                        ++errors2;
                }
                if(errors2 < s.l[j] || errors2 > s.u[j]){
                    cout << "Triggered: " << (int)errors2 << endl;
                    valid = false;
                    break;
                }
            }
            if(valid){
                cout << "Hiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiit" << endl;
                cout << (int)errors2 << endl;
                uint32_t occ = startPos;
                hitsv.push_back(Pair<uint16_t,uint32_t>(sa_info.i1, occ));
                cout << "Hit occ: " << hitsv[hitsv.size() - 1] << endl;
                errorsv.push_back(errors2);
            }
        }
    }
    delegateDirect(hitsv, needle, errorsv);
}
    
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeChildren(TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,   
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const & /**/)                                 
{
    cout << "_optimalSearchSchemeChildren" << endl;
    cout << "Error: " << (int)errors << endl;
    cout << "blockIndex: " << (int)blockIndex << endl;
    
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
                
                    if(!(std::is_same<TDir, Rev>::value ^ !goToRight2)){
                        cout << "switching direction in unidirectional function" << endl;
                        cout << "blockIndex: " << blockIndex << endl;
                        printv(s.pi);
                        exit(0);
                    }

                //TODO remove goToRight2 and input should be TDIR 
                if (goToRight2)
                {
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev());
                }
                else
                {
                    _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd());
                }
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
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeExact(TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,    
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const & /**/)
{
    cout << "_optimalSearchSchemeExact" << endl;
    cout << "Error: " << (int)errors << endl;
    cout << "blockIndex: " << (int)blockIndex << endl;

    bool goToRight2 = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] : s.pi[blockIndex] > s.pi[blockIndex - 1];
    
        //sanity check
    if((blockIndex < s.pi.size() - 1) && !(std::is_same<TDir, Rev>::value ^ !goToRight2)){
        cout << "switching direction in unidirectional function" << endl;
        cout << "blockIndex: " << blockIndex << endl;
        printv(s.pi);
        exit(0);
    }
    
    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it forward
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

        cout << "SearchScheme: " << endl;
        printv(s.pi);
        cout << "Needleeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee" << endl;
         for(int i = infixPosLeft; i < infixPosRight + 1; ++i)
            cout << needle[i];
        cout << endl;
        
        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1)))
            return;
        
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Rev());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Fwd());
        }
    }
    else
    {
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        int32_t infixPosLeft = needleRightPos - s.blocklength[blockIndex] - 1;
        int32_t infixPosRight = needleLeftPos - 1;

        cout << "SearchScheme: " << endl;
        printv(s.pi);
        cout << "Needleeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee(Fwd Idx)" << endl;
        for(int i = infixPosLeft; i < infixPosRight + 1; ++i)
            cout << needle[i];
        cout << endl;
        
        while (infixPosRight >= infixPosLeft)
        {
            if (!goDown(iter, needle[infixPosRight]))
                return;
            --infixPosRight;
        }
        
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd());
        }
    }
}
    
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,   
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,                             
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 TDir const & /**/)
{
    cout << "_optimalSearchScheme" << endl;
    cout << "Error: " << (int)errors << endl;
    cout << "blockIndex: " << (int)blockIndex << endl;
    cout << "NLP: " << (int)needleLeftPos << endl;
    cout << "NRP: " << (int)needleRightPos << endl;
    cout << "SS: ";
    printv(s.pi);
    cout << "RangeSize: " << iter.vDesc.range.i2 - iter.vDesc.range.i1 << endl;
    
    
    if (std::is_same<TDir, Rev>::value){
        cout << "Reverse Index" << endl;
    }else{
        cout << "Forward Index" << endl;
    }
//     cout << "SaRange: " << iter.vDesc.range << endl;
//     cout << "Reached End" << endl;


//     int directSearchThreshold = 2;


    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
    
    
    if(needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex]){
        Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval = get_bitvector_interval_inside(iter, bitvectors, s, blockIndex + 1, TDir());
        sdsl::bit_vector & b = bitvectors[bit_interval.i1].first;
        sdsl::rank_support_v<> & rb = bitvectors[bit_interval.i1].second;
        rb.set_vector(&b);
        
        cout << "at blockend with blockIndex: " << (int)blockIndex << " " << s.blocklength[blockIndex]   << endl;
        cout << "PrintUNI" << endl;
        printbit(bitvectors, bit_interval);
        cout << "PrintUend" << endl;
    
        if(rb(bit_interval.i2.i2) - rb(bit_interval.i2.i1) == 0){
            cout << "This is very unlikly depending on parameters" << endl;
            return;
        }
        cout << "Iter Range: " << iter.vDesc.range.i2 - iter.vDesc.range.i1 << endl;

        // Done. (Last step)
        if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
        {
            cout << "All Unidirectional occs: " << endl;
            for (auto occ : getOccurrences(iter))
            {                
                cout << occ.i2 << endl;
                if(std::is_same<TDir, Rev>::value)
                    cout << "still need to calc for fwd Index" << endl;
            }
            //TODO disable true statement to allowed cheap repeats 
            if(rb(bit_interval.i2.i2) - rb(bit_interval.i2.i1) < bit_interval.i2.i2 - bit_interval.i2.i1)
            {
                uint32_t rangeStart = iter.vDesc.range.i1;
                uint32_t rangeEnd = iter.vDesc.range.i2;
                int lastStart = 0;
                for(int i = 0; i < rangeEnd - rangeStart; ++i)
                {
                    if(bitvectors[bit_interval.i1].first[bit_interval.i2.i1 + i] == 0)
                    {
                        if(i != lastStart){
                            cout << "Called delegate on the following range:" << endl;
                            iter.vDesc.range.i1 = rangeStart + lastStart;
                            iter.vDesc.range.i2 = rangeStart + i - 1;
                            cout << iter.vDesc.range.i1 << " - " << iter.vDesc.range.i2;
                            delegate(iter, needle, errors, std::is_same<TDir, Rev>::value);
                        }
                        lastStart = i + 1;
                    }
                }
            }
            else
            {
                cout << "CompMappableInterval" << endl;
                delegate(iter, needle, errors, std::is_same<TDir, Rev>::value);
            }
            cout << "Finished unidirectional Search" << endl;
            return; 
        }
    }
    
  /*  
     //very unlikly
    else if(needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex] && iter.vDesc.range.i2 - iter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1) * directSearchThreshold){
        cout << "start direct Search Unidirectional" << endl;
//         bit_interval = get_bitvector_interval_inside(iter, bitvectors, s, blockIndex, TDir());
        directSearch(delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir());
        return;
    }
    
    */
    
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
