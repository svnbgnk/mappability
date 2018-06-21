#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

// #include <seqan/arg_parse.h>
// #include <seqan/seq_io.h>
#include <seqan/index.h>
#include <sdsl/bit_vectors.hpp>
#include "common.h"


namespace seqan{
/*    
template <size_t N>
struct OptimalSearch
{
    std::array<uint8_t, N> pi; // order of the blocks. permutation of [1..n]
    std::array<uint8_t, N> l; // minimum number of errors at the end of the corresponding block
    std::array<uint8_t, N> u; // maximum number of errors at the end of the corresponding block

    std::array<uint32_t, N> blocklength; // cumulated values / prefix sums
    std::array<uint32_t, N> min;
    std::array<uint32_t, N> max;
    uint8_t startUniDir;
    uint32_t startPos;
};
    */
enum class ReturnCodes {
	BLACK, RED, BLUE, WHITE
};

template <size_t nbrBlocks, size_t N>
constexpr inline void _optimalSearchSchemeSetMapParams(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{
    for (OptimalSearch<nbrBlocks> & s : ss){
        uint8_t min = s.pi[0];
        uint8_t max = s.pi[0];
        // maybe < N?
        for(int i = 0; i < s.pi.size(); ++i){
            if(min > s.pi[i])
                min = s.pi[i];
            if(max < s.pi[i])
                max = s.pi[i];
            s.min[i] = min;
            s.max[i] = max;
        }
        uint8_t lastValue = s.pi[s.pi.size() - 1];
        int k = s.pi.size() - 2;
        while(k >= 0){
            if(s.pi[k] == lastValue - 1 || s.pi[k] == lastValue + 1)
            {
                lastValue = s.pi[k];
                --k;
            }else{
                s.startUniDir = k + 1;
                break;
            }
        }
    }
}

    
template< size_t N>
uint8_t mymin(std::array <uint8_t, N> v, uint8_t end)
{
    uint8_t min = v[0];
    for(uint8_t i = 1; i < end; ++i){
        if(v[i] < min)
            min = v[i];
    }
    return min;
}

template< size_t N>
uint8_t mymax(std::array <uint8_t, N> v, uint8_t end){
    uint8_t max = v[0];
    for(uint8_t i = 1; i < end; ++i)
    {
        if(v[i] > max)
            max = v[i];
    }
    return max;
}

template <typename TText, typename TIndex, typename TIndexSpec>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              Pair<uint32_t, uint32_t>  const & range)
{
    int number_of_indeces = countSequences(iter.fwdIter.index);
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    for(uint32_t i = range.i1; i < range.i2; ++i){
        int seq = iter.fwdIter.index->sa[i].i1;
        int sa = iter.fwdIter.index->sa[i].i2;
        cout << i << ": " << sa + sequenceLengths[seq] << endl;
    }
}

template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks, typename TDir>
sdsl::bit_vector get_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                        vector<sdsl::bit_vector> const & bitvectors,
                                        OptimalSearch<nbrBlocks> const & s,
                                        uint8_t const blockIndex,
                                        Pair<uint32_t, uint32_t>  const & range,
                                        TDir const & /**/) 
{
    uint8_t needed_bitvector;
    sdsl::bit_vector b(range.i2 - range.i1, 0);
    if (std::is_same<TDir, Rev>::value)
        needed_bitvector = mymin(s.pi, blockIndex + 1) - 1;
    else
        needed_bitvector = mymax(s.pi, blockIndex + 1);
    int number_of_indeces = countSequences(iter.fwdIter.index);
//     cout << "Min Element: " << (int)min << endl;
    cout << "shift_" << (int)needed_bitvector << "_bitvector" << endl;
    for(int i = range.i1; i < range.i2; ++i)
        b[i - range.i1] = bitvectors.at(needed_bitvector)[i + number_of_indeces];
    return b;
}

void printb(sdsl::bit_vector b){
    for(int i = 0; i < b.size(); ++i)
        cout << i << " Bit: " << b[i] << endl;
}

template<typename TText, typename TIndex, typename TIndexSpec>
uint8_t squash_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > & iter,
                     sdsl::bit_vector const & bi)
{
    cout << "Range after squash in Function: " << range(iter.fwdIter) << endl;
    //TODO later directly a access correct range in rank_support_v
    sdsl::rank_support_v<> rbi(& bi);
    if(rbi(bi.size() == 0))
        return 0;
    if(rbi(bi.size()) == bi.size())
        return 1; 
    if(rbi(bi.size()) <= 3){
        return 2;
    }
    return 99;
    /*
    uint32_t startPos = 0, endPos = bi.size();
    for(uint32_t i = 0; i < bi.size(); ++i){
        if(bi[i] == 0)
            ++startPos;
        if(bi[bi.size() - 1 - i] == 0)
            --endPos;
    }
    iter.fwdIter.vDesc.range.i1 = iter.fwdIter.vDesc.range.i1 + startPos; 
    iter.fwdIter.vDesc.range.i2 = iter.fwdIter.vDesc.range.i1 + endPos;
    return 2;*/
}

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeChildren(TDelegate & delegate,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         vector<sdsl::bit_vector> const & bitvectors,  
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const & /**/)
{
    bool goToRight = std::is_same<TDir, Rev>::value;
    if (goDown(iter, TDir()))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()),
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
                if (goToRight2)
                {
                    _optimalSearchScheme(delegate, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s,
                                            blockIndex2, Rev());
                }
                else
                {
                    _optimalSearchScheme(delegate, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s,
                                            blockIndex2, Fwd());
                }
            }
            else
            {
                _optimalSearchScheme(delegate, iter, needle, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s,
                                     blockIndex, TDir());
            }
        } while (goRight(iter, TDir()));
    }
}   

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeExact(TDelegate & delegate,
                                      Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      vector<sdsl::bit_vector> const & bitvectors,  
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const & /**/)
{
    // not in last block and next Block is larger then current block
    bool goToRight2 = (blockIndex < s.pi.size() - 1) && s.pi[blockIndex + 1] > s.pi[blockIndex];
    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it reverse
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
            return;

        cout << "Lets print genome" << endl;
//         auto & genome = indexText(index(iter.fwdIter));
//         cout << "First pos " << genome[0][0] << endl;
//         cout << "Second pos " << genome[0][0] << endl;
        cout << "Range: " << range(iter.fwdIter) << endl;
        cout << "Range2: " << range(iter.revIter) << endl;
        seqan::Pair<uint16_t, uint32_t> r = range(iter.fwdIter);
        cout << "Suffices at from that interval" << endl;
        print_sa(iter, r);
        sdsl::bit_vector bit_interval = get_bitvector_interval(iter, bitvectors, s, blockIndex, r, TDir());
        cout << "Mappability of this suffix interval: " << endl;
        printb(bit_interval);
        uint8_t rcode = squash_interval(iter, bit_interval);
        cout << "Return code: " << (int)rcode << endl;
        cout << "Range after squash: " << iter.fwdIter.vDesc.range << endl;
        cout << "Range2 after squash: " << iter.revIter.vDesc.range << endl;
        if(rcode == 2){
            //search directly in Genome
            for(int i = 0; i < bit_interval.size(); ++i){
                if(bit_interval[i] == 1){
                    // in this case (fwd index) 
                    uint32_t pos = iter.fwdIter.vDesc.range.i1 + i;
                    uint16_t seq = iter.fwdIter.vDesc.range.i2 + i;
                    uint32_t sa = iter.fwdIter.index->sa[pos].i2;                    
                    //get correct sa pos muliple sequences
                    // compare full thing again later split into left and right part
                    //genome(seq,sa -needlePosLeft, sa + length(needle))
//                     iter.fwdIter.text;
                }
            }
            
            return;
        }
        if(rcode == 0)
            return;
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev());
        }
        else
        {
            _optimalSearchScheme(delegate, iter, needle, bitvectors, needleLeftPos, infixPosRight + 2, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd());
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
/*       
        cout << "Range: " << range(iter.revIter) << endl;
        auto r2 = range(iter.revIter);
        uint8_t max = mymax(s.pi, blockIndex + 1);
        cout << "Max Element: " << (int)max << endl;
        
        for(int i = r2.i1; i < r2.i2; ++i){
            cout << bitvectors.at(max)[i] << endl;
        }
        
        cout << "shift" << (int)max << "bitvector" << endl;
        for(int i = 0; i < bitvectors.at(1).size(); ++i){
            cout << i << ": " << bitvectors.at(1)[i] << endl;
        }*/
        
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev());
        }
        else
        {
            _optimalSearchScheme(delegate, iter, needle, bitvectors, infixPosLeft, needleRightPos, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd());
        }
    }
}
    

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 vector<sdsl::bit_vector> const & bitvectors,    
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,                             
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 TDir const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

//     cout << "Step: " << needleRightPos - needleLeftPos - 1 << "    ss: "; printv(s.pi); cout << endl;
    // Done.
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        delegate(iter, needle, errors);
    }
    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
    }
    // Approximate search in current block.
    else
    {
        _optimalSearchSchemeChildren(delegate, iter, needle, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex,
                                     minErrorsLeftInBlock, TDir());
    }
}
  
template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 vector<sdsl::bit_vector> const & bitvectors,
                                 OptimalSearch<nbrBlocks> const & s)
{
     _optimalSearchScheme(delegate, it, needle, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, Rev());
}

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks, size_t N>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 vector<sdsl::bit_vector> const & bitvectors,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss)
{
    for (auto & s : ss)
        _optimalSearchScheme(delegate, it, needle, bitvectors, s);
}  

template <size_t minErrors, size_t maxErrors,
          typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec>
inline void
find(TDelegate & delegate,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     vector<sdsl::bit_vector> const & bitvectors)
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
//     print_search_scheme(scheme);
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    cout << "Print Genome: " << endl;
    auto & genome = indexText(index);

//     auto & testIndex = it.index;
//     auto & testIndex = index(it)

//     auto & genome2 = indexText(index(it.fwdIter));
//     auto & genome2 = indexText(it.fwdIter.index);
    cout << genome[0][0] << endl;

    _optimalSearchScheme(delegate, it, needle, bitvectors, scheme);
}
  
  
  

template <size_t minErrors, size_t maxErrors,
          typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TParallelTag>
inline void
find(TDelegate & delegate,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<sdsl::bit_vector> const & bitvectors,
     TParallelTag const & /**/)
{
    typedef typename Iterator<StringSet<TNeedle, TStringSetSpec> const, Rooted>::Type TNeedleIt;
    typedef typename Reference<TNeedleIt>::Type                                       TNeedleRef;
    iterate(needles, [&](TNeedleIt const & needleIt)
    {
        TNeedleRef needle = value(needleIt);
        find<minErrors, maxErrors>(delegate, index, needle, bitvectors);
    },
    Rooted(), TParallelTag());
}    
    
template <size_t minErrors, size_t maxErrors,
          typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec>
inline void
find(TDelegate & delegate,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     vector<sdsl::bit_vector> const & bitvectors)
{
    find<minErrors, maxErrors>(delegate, index, needles, bitvectors, Serial());
}

}

#endif
