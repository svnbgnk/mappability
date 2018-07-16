#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_START_UNIDIRECTIONAL_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_START_UNIDIRECTIONAL_EXTENSION_H_

#include "auxiliary.h"

namespace seqan{

    
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
void filter_interval(TDelegate & delegate,
                     TDelegateD & delegateDirect,
                     Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                     TNeedle const & needle,
                     vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                     uint32_t const needleLeftPos,
                     uint32_t const needleRightPos,
                     uint8_t const errors,
                     OptimalSearch<nbrBlocks> const & s,
                     uint8_t const blockIndex,
                     Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                     TDir const & /**/)
{ 


    cout << "In start Unidirectional filterInterval" << endl;
    printbit(bitvectors, inside_bit_interval);
    vector<pair<uint32_t, uint32_t>> consOnes = getConsOnes(bitvectors, inside_bit_interval, params.startuni.intervalsize);
    
    uint32_t noi = seqan::length(iter.index->sa) - bitvectors[0].first.size(); // number_of_indeces
    //TODO replace with countSequences when it works
    
    for(int i = 0; i < consOnes.size(); ++i){
        cout << "Print Bit Range";
        printPair(consOnes[i]);
        cout << endl;
        iter.vDesc.range.i1 = consOnes[i].first + noi;
        iter.vDesc.range.i2 = consOnes[i].second + noi;
        if (std::is_same<TDir, Rev>::value){
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Rev());
            //TODO Maybe Link to old Funtion to not filter multiple times?
        }
        else
        {
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, Fwd());
        }
    } 
}

    
    //TODO merge with old funtion (test unidirectional) diff: iter and pramenter filpDensity return value
template<typename TText, typename TConfig, typename TIndexSpec,
         size_t nbrBlocks,
         typename TDir>
bool testFilter(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                          vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex,
                          TDir const & )
{
    // need bitinterval from inside the pattern to filter according to the mappability form
    //therefore i also need to acces the block before because of that block i got mappability of both sides
    auto bit_interval = get_bitvector_interval_inside(iter, bitvectors, s, blockIndex, TDir());
    sdsl::bit_vector & b2 = bitvectors[bit_interval.i1].first;
    
    //squas interval
    uint32_t startPos = bit_interval.i2.i1, endPos = bit_interval.i2.i2;
    
    for(uint32_t i = startPos; i < endPos; ++i){
        if(b2[i] != 0)
            break; 
        ++startPos;
    }
    for(uint32_t i = endPos - 1; i >= startPos; --i){
        if(b2[i] != 0)
            break;
        --endPos;
    }
    
    if(startPos > endPos){
        cout << "Error bit vector has only zeroes this should have been checked by checkinterval" << endl;
        cout << "Size: " << endPos - startPos << endl;
        exit(0);
    }
    // order of bits
    bool last = b2[startPos];
    uint32_t pos = startPos;
    uint32_t count = 0; 
    while(pos < endPos){
        if(b2[pos] != last){
            ++count;
            last = !last;
        }
        ++pos;
    }  
    float ivalSize = brange.i2.i2 - brange.i2.i1;
    // if next condition is true then brange will be modified!!!!
    // it will contain mappability of the bitvector from the other side of already searched needle
    
    // only interested in changes inside the supinterval (startPos - endPos)
    // if 0 got Cutoff that means more changes but is at the same time also good
    // therefore ignore them
    cout << "Test flipdensitey " << endl;
    cout << ivalSize * params.startuni.invflipdensity - 1 << endl;
    cout << "count: " << count << endl;

    if(ivalSize * params.startuni.invflipdensity - 1 > static_cast<float>(count)){
        cout << "FILTER!!!!!!!!!" << endl;
        brange.i1 = bit_interval.i1;
        cout << "New selected bitvector by inside function: " << (int)bit_interval.i1 << endl;
        brange.i2.i1 = startPos;
        brange.i2.i2 = endPos;
        return true;
    }
    cout << "Still continue UNIDIRECTIONAL" << endl;
    return false;
}

ReturnCode checkInterval(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                          uint8_t const blockSize,
                          bool const done,
                          bool const nofilter,
                          uint8_t const blockIndex)
{
    if(!nofilter){
        cout << "unicheckInterval" << endl;
        sdsl::bit_vector & b = bitvectors[brange.i1].first;
        sdsl::rank_support_v<> & rb = bitvectors[brange.i1].second; 
        rb.set_vector(&b);
    
        uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
        if(params.startuni.nomappability && ivalOne == 0)
            return ReturnCode::NOMAPPABILITY;

        if(!done){
            if(params.startuni.directsearch && ivalOne < (blockSize - blockIndex - 1) * params.startuni.directsearch_th){ //<4 
                cout << "UNIDIRECTSEARCHDIRECTSEARCHDIRECTSEARCH" << endl;
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
        cout << "MAPPABLEMAPPABLEMAPPABLEMAPPABLEMAPPABLE" << endl;
        return ReturnCode::MAPPABLE;
    }
    else
    {
        cout << "unicheckInterval no start uni" << endl;
        sdsl::bit_vector & b = bitvectors[brange.i1].first;
        sdsl::rank_support_v<> & rb = bitvectors[brange.i1].second; 
        rb.set_vector(&b);
    
        uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
        if(params.uni.nomappability && ivalOne == 0)
            return ReturnCode::NOMAPPABILITY;

        if(!done){
            if(params.uni.directsearch && ivalOne < (blockSize - blockIndex - 1) * params.uni.directsearch_th){ //<4 
                cout << "UNIDIRECTSEARCHDIRECTSEARCHDIRECTSEARCH" << endl;
                return ReturnCode::DIRECTSEARCH;
            }    
            if(params.uni.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
                return ReturnCode::COMPMAPPABLE;
        }
        cout << "MAPPABLEMAPPABLEMAPPABLEMAPPABLEMAPPABLE" << endl;
        return ReturnCode::MAPPABLE;
    }
    
    
    std::cerr << "Something went wrong!!!" << "\n";
    exit(0);
    return ReturnCode::MAPPABLE;
}


template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir> 
ReturnCode uniCheckMappability(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,   
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,                             
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 bool const done,
                                 bool const nofilter,
                                 TDir const & )
{
    cout << "uniCheckMappability" << endl;
    
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval = get_bitvector_interval_inside(iter, bitvectors, s, blockIndex + done, TDir());
    ReturnCode rcode = checkInterval(bitvectors, bit_interval, s.pi.size(), done, nofilter, blockIndex);
    
    cout << "at blockend with blockIndex: " << (int)blockIndex << endl;
    cout << "PrintUNI" << endl;
      
    printbit(bitvectors, bit_interval);
    cout << "PrintUend" << endl;
    
    //TODO checkUnidirectionalMappability here
    if(rcode == ReturnCode::NOMAPPABILITY){
        cout << "This is very unlikly depending on parameters" << endl;
        return ReturnCode::FINISHED;
    }
    cout << "Iter Range: " << iter.vDesc.range.i2 - iter.vDesc.range.i1 << endl;

    // Done. (Last step)
    if (done)
    {
        bool rev = std::is_same<TDir, Rev>::value;
        cout << "All Unidirectional occs: " << endl;
        cout << "Direction during delegate call" << (rev);
        cout << endl;
        for (auto occ : getOccurrences(iter))
        {                
            cout << occ.i2 << endl;
            if(rev)
                cout << "still need to calc for fwd Index" << endl;
        }
        //TODO disable true statement to allowed cheap repeats 
        //
        if(false) 
        {
            uint32_t rangeStart = iter.vDesc.range.i1;
            uint32_t rangeEnd = iter.vDesc.range.i2;
            int lastStart = 0;
            for(int i = 0; i < rangeEnd - rangeStart; ++i)
            {
                if(bitvectors[bit_interval.i1].first[bit_interval.i2.i1 + i] == 0 )
                {
                    if(i != lastStart){
                        cout << "Called delegate on the following range:" << endl;
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
            cout << "CompMappableInterval" << endl;
            delegate(iter, needle, errors, rev);
        }
        cout << "Finished unidirectional Search" << endl;
        return ReturnCode::FINISHED; 
    }else if(rcode == ReturnCode::DIRECTSEARCH){
        cout << "start direct Search Unidirectional" << endl;
        uniDirectSearch(delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir());
        return ReturnCode::FINISHED;
    }else if(rcode == ReturnCode::FILTER){
        //test filter also modfied iter range if true;
        if(testFilter(iter, bitvectors, bit_interval, s, blockIndex, TDir())){
            cout << "filterstartUniDirectionalInterval" << endl;
            filter_interval(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir());
            return(ReturnCode::FINISHED);
        }
    }
    return ReturnCode::MAPPABLE;
}
  
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _uniOptimalSearchSchemeChildren(TDelegate & delegate,
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
    cout << "_uniOptimalSearchSchemeChildren" << endl;
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
          size_t nbrBlocks,
          typename TDir>
inline void _uniOptimalSearchSchemeExact(TDelegate & delegate,
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
    cout << "_uniOptimalSearchSchemeExact" << endl;
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
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Rev());
        }
        else
        {
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, Fwd());
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
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Rev());
        }
        else
        {
            _uniOptimalSearchScheme(delegate, delegateDirect, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, Fwd());
        }
    }
}
    
template <typename TDelegate, typename TDelegateD,
          typename TText, typename TConfig, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _uniOptimalSearchScheme(TDelegate & delegate,
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
    cout << "_uniOptimalSearchScheme" << endl;
    cout << "Error: " << (int)errors << endl;
    cout << "blockIndex: " << (int)blockIndex << endl;
    cout << "NLP: " << (int)needleLeftPos << endl;
    cout << "NRP: " << (int)needleRightPos << endl;
    cout << "SS: ";
    printv(s.pi);
    if(needleRightPos - needleLeftPos > 1) 
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
    cout << "checkUniMappa?" << endl;
    cout << "current blocklength   " << needleRightPos - needleLeftPos - 1 << endl;
    cout << "calc   " << s.blocklength[blockIndex - 1] << endl;
    bool done = minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1;

    if(blockIndex > 0 && done || needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex - 1]){
        cout << "start checkUnidirectionalMappability" << endl;
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