

#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_COMPMAPPLE_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_COMPMAPPLE_H_

namespace seqan{
    
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
//         delegate(iter, needle, errors);
        cout << "Delegate Call comp:" << endl;
        delegate(iter, needle, errors, false); //TODO put iterator test in the original code?
    }
    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(),
                                  TDistanceTag());
    }
 
    else
    {
        //TODO test for directSearch
        cout << "Test for cDirectSearch:" << endl;
        cout << "PrintSA:" << endl;
        print_sa(iter, 1 , !(std::is_same<TDir, Rev>::value));
        
        int directSearchTrigger = 2;
        if(iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1) * directSearchTrigger){
            cout << "Current Range: " << (int)(iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1) << endl;
            directSearch(delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
            return;
        }
        _optimalSearchSchemeChildren(delegate, delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex,
                                     minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}

template <typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
void directSearch(TDelegateD & delegateDirect,
                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  TDir const & /**/)
{
    cout << "compdirectSearch " <<  endl;
    cout << "NLP: " <<  needleLeftPos <<  endl;
    cout << "NRP: " <<  needleRightPos <<  endl;
    cout << "errors:  " <<  (int)errors <<  endl;
    vector<Pair<uint16_t, uint32_t>> hitsv;
    vector<uint8_t> errorsv;
    if(std::is_same<TDir, Rev>::value){
        auto const & genome = indexText(*iter.fwdIter.index);
        for(int i = iter.fwdIter.vDesc.range.i1; i < iter.fwdIter.vDesc.range.i2; ++i){
            cout << "blockIndex: " << (int)blockIndex << endl;
            uint8_t errors2 = errors;
            bool valid = true;
            Pair<uint16_t, uint32_t> sa_info= iter.fwdIter.index->sa[i];
            uint32_t startPos = sa_info.i2 - needleLeftPos;
//            cout <<  "Sa info" <<  sa_info <<  endl; //TODO redo this
            cout << "StartPos " << startPos << endl;
            //search remaining blocks (also finished the current block if neede)
            for(int j = blockIndex; j < s.pi.size(); ++j){
                int blockStart = (s.pi[j] - 1 == 0) ? 0 : s.chronBL[s.pi[j] - 2];
                int blockEnd = s.chronBL[s.pi[j] - 1];
                cout << "searching Parts:" << blockStart << " - " << blockEnd << "; ";
                cout << endl;
                // compare bases to needle                
                if(std::is_same<TDir, Rev>::value){
                        if(needleRightPos - 1 > blockStart && needleRightPos - 1 < blockEnd){
                            cout << "changing Blockstart rev" << endl;
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
    else
    {
        auto const & rgenome = indexText(*iter.revIter.index);
        for(int i = iter.revIter.vDesc.range.i1; i < iter.revIter.vDesc.range.i2; ++i){
            cout << "Direct Search FWD" << endl;
            cout << "blockIndex: " << (int)blockIndex << endl;
            uint8_t errors2 = errors;
            bool valid = true;
            Pair<uint16_t, uint32_t> sa_info = iter.revIter.index->sa[i];
            uint32_t startPos = sa_info.i2 - (length(needle) - needleRightPos + 1);
            // iter.fwdIter.vDesc.range.i1 is not the same brange.i2.i1 since sentinels are at the beginning!!!
//            cout <<  "Sa info" <<  sa_info <<  endl; //TODO redo this
            cout << "StartPos " << startPos << endl;
            //search remaining blocks             
            uint8_t blocks = s.pi.size();
            for(int j = blockIndex; j < s.pi.size(); ++j){
                int blockStart = (s.pi[j] == s.pi.size()) ? 0 : s.revChronBL[s.pi[j]];
                int blockEnd = s.revChronBL[s.pi[j] - 1];
                cout << "searching Parts:" << length(needle) - blockStart << " - " << length(needle) - blockEnd << "; ";
                cout << endl;
                if(needleLeftPos < length(needle) - blockStart && needleLeftPos > length(needle) - blockEnd){
                    cout << "changing Blockstart, should only happen (once) when we come from checkcurrentmappability!!" << endl;
                    blockStart = length(needle) - needleLeftPos; //- 1 + 1
                    cout << "searching Parts:" << blockStart << " - " << blockEnd << "; " << endl;
                }
                for(int k = blockStart; k < blockEnd; ++k){
                    if(needle[length(needle) - k - 1] != rgenome[sa_info.i1][startPos + k])
                        ++errors2;
                }
                cout << "Errors until now: " << (int)errors2 << endl;
                if(errors2 < s.l[j] || errors2 > s.u[j]){
                    cout << "Triggered: " << (int)errors2 << endl;
                    valid = false;
                    break;
                }
            }
            if(valid){
                cout << "Hiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiit" << endl;
                cout << (int)errors2 << endl;
                uint32_t occ = seqan::length(rgenome[sa_info.i1]) - startPos - length(needle);
                hitsv.push_back(Pair<uint16_t,uint32_t>(sa_info.i1, occ));
                cout << "Hit occ FWD Index: " << hitsv[hitsv.size() - 1] << endl;
                errorsv.push_back(errors2);
            }
        }
    }
    delegateDirect(hitsv, needle, errorsv);
    cout << "compdirectSearchEnd " <<  endl;
    
}

}
#endif

