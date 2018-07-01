#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_UNIDIRECTIONAL_H_

#include <sdsl/bit_vectors.hpp>
#include "common.h"
#include "common_auxiliary.h"

namespace seqan{

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
