#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_

#include <sdsl/bit_vectors.hpp>
#include "common.h"
#include "common_auxiliary.h"

namespace seqan{

    
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
  /*              
                cout << "checkMappability Call from Children" << endl;
                ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, goToRight2, TDir());
                if(rcode == ReturnCode::FINISHED)
                    return;*/
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
    
    // not in last block and next Block is larger then current block
    bool goToRight2 = (blockIndex < s.pi.size() - 1) && s.pi[blockIndex + 1] > s.pi[blockIndex];
        //sanity check
    if(!(std::is_same<TDir, Rev>::value ^ !goToRight2)){
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
 /*       
        //TODO does something go wrong in checkMappability if blockIndex is not increased by one? then do Todo before
        cout << "checkMappability Rev Index" << endl;
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, goToRight2, TDir());
        //TDir() is in this case Rev() ...
        if(rcode == ReturnCode::FINISHED)
            return;     */ 
        //TODO remove goToRight2 and input should be TDIR 
        if (goToRight2)
        {
            cout << "Start Next Function with Block Index (Rev): " << blockIndex2 << endl;
            _optimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Rev());
        }
        else
        {
            cout << "Start Next Function with Block Index (Fwd): " << blockIndex2 << endl;
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
/*
        cout << "checkMappability Fwd Index" << endl;
        ReturnCode rcode = checkMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, infixPosLeft, needleRightPos, errors, s, blockIndex2, goToRight2, TDir());
        //TDir() is in this case Fwd() ...
        if(rcode == ReturnCode::FINISHED)
            return;
        */
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

    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

//     cout << "Step: " << needleRightPos - needleLeftPos - 1 << "    ss: "; printv(s.pi); cout << endl;
    // Done. (Last step)
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        cout << "Finished unidirectional Search" << endl;
        delegate(iter, needle, errors);
    }
    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
    }
    // Approximate search in current block.
    else
    {
   /*   
        //TODO implement faster mod 
        if((needleRightPos - needleLeftPos - 1) % 4 == 0){
            //TODO Check if we are Done here?
            //TODO stop doing on loop to much (enter checkCurrentMappability a second time)
            bool goToRight2 = std::is_same<TDir, Rev>::value;
            ReturnCode rcode = checkCurrentMappability(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, goToRight2, TDir());
            // blockIndex is different here!!! goToRight2 is missing
            if(rcode == ReturnCode::FINISHED)
                return;
        }*/
        
    _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir());
    }
}    
    
    
    
template <typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle, typename TIndex,
          size_t nbrBlocks,
          typename TDir>
void directSearchDummy(TDelegateD & delegateDirect,
                  Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > olditer,
                  TNeedle const & needle,
                  vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, 
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  Pair<uint8_t, Pair<uint32_t, uint32_t>> brange,
                  TDir const & /**/)
{
    cout << "Test Dummy Function: " << endl;
    
//     for(int i = 0; i < seqan::length(iter.index->sa); ++i)
//         cout << iter.index->sa[i] << endl; 
    
//     cout << "Length of Sequences" << endl;
//     for(int i = 0; i < 4; ++i)
//         cout << olditer.fwdIter.index->sa[i] << endl;
    cout << "Ranges" << endl;
    cout <<  iter.vDesc.range << endl;
    cout << "Old Iter ranges" << endl;
    cout << olditer.fwdIter.vDesc.range << endl;
    cout << "Reverse" << endl;
    cout << olditer.revIter.vDesc.range << endl;
        //TODO dont care about different sequences in the index in the moment
    int needleSize = needleRightPos- needleLeftPos - 1;
    cout << "Size of searched needle: " << needleSize << endl;
    vector<int> saveOppositePositions;

    
    int number_of_indeces = seqan::length(olditer.fwdIter.index->sa) - bitvectors[0].first.size();
    vector<int> sequenceLengths(number_of_indeces, 0); 
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[olditer.fwdIter.index->sa[i].i1] = olditer.fwdIter.index->sa[i].i2;
    
    printv(sequenceLengths);
    
    auto fwdr = olditer.fwdIter.vDesc.range;
    auto revr = olditer.revIter.vDesc.range;
    
    if(std::is_same<TDir, Rev>::value){
        cout << "Reverse case therefore need to save fwd sa position" << endl;
        for(int i = fwdr.i1; i < fwdr.i2; ++i){
            cout << olditer.fwdIter.index->sa[i] << endl;
        }
        cout << "other sa" << endl;
        for(int i = revr.i1; i < revr.i2; ++i){
            cout << olditer.revIter.index->sa[i] << endl;
        }
        
        for(int i = revr.i1; i < revr.i2; ++i){
            int seq = olditer.revIter.index->sa[i].i1;
            cout << "Seq Length" << sequenceLengths[seq] << endl;
            int findn = sequenceLengths[seq] - static_cast<int>(olditer.revIter.index->sa[i].i2) - needleSize;
            cout << "findn  " << findn << endl;
            int k = fwdr.i1;
            while(k < fwdr.i2){
                if(findn == olditer.fwdIter.index->sa[k].i2)
                    break;
                ++k;
            }
            saveOppositePositions.push_back(k);
        }     
    }else{
        cout << "Forward case therefore need to save rev sa position" << endl;
        for(int i = revr.i1; i < revr.i2; ++i){
            cout << olditer.revIter.index->sa[i] << endl;
        }
        cout << "other sa" << endl;
        for(int i = fwdr.i1; i < fwdr.i2; ++i){
            cout << olditer.fwdIter.index->sa[i] << endl;
        }
        
        for(int i = fwdr.i1; i < fwdr.i2; ++i){
            int seq = olditer.fwdIter.index->sa[i].i1;
            int findn = sequenceLengths[seq] - olditer.fwdIter.index->sa[i].i2 - needleSize;
            cout << olditer.fwdIter.index->sa[i].i2 << " equivalent to  " << findn << endl;
            int k = revr.i1;
            while(k < revr.i2){
                if(findn == olditer.revIter.index->sa[k].i2)
                    break;
                ++k;
            }
            saveOppositePositions.push_back(k);
        }
    }
        
        
    cout << "Saved Positions" << endl;
    for(int i = 0; i < saveOppositePositions.size(); ++i)
        cout << saveOppositePositions[i] << endl;
    cout << "End" << endl;
    
}




}
#endif
