using namespace seqan;

template <unsigned errors, typename TIndex, typename TText, typename TContainer>
inline void runAlgoTrivial(TIndex & index, TText const & text, TContainer & c, SearchParams const & searchParams)
{
    typedef typename TContainer::value_type value_type;

    constexpr uint64_t max_val = std::numeric_limits<typename TContainer::value_type>::max();
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, searchParams.length);
    _optimalSearchSchemeComputeChronBlocklength(scheme);

    uint64_t textLength = seqan::length(text);

    #pragma omp parallel for schedule(dynamic, 1000000) num_threads(searchParams.threads)
    for (uint64_t i = 0; i < textLength - searchParams.length + 1; ++i)
    {
        value_type hits = 0;
        auto delegate = [&hits](auto const &it, auto const & /*read*/, unsigned const /*errors*/) {
//             for(uint32_t i = it.fwdIter.vDesc.range.i1; i < it.fwdIter.vDesc.range.i2; ++i)
//                 std::cout << it.fwdIter.index->sa[i] << std::endl;
            if ((uint64_t) hits + countOccurrences(it) <= max_val)
                hits += countOccurrences(it);
            else
                hits = max_val;
        };

        auto delegateDirect = [&hits](auto const & pos, DnaString const & /*needle*/, uint8_t const /*errors*/)
        {
            if((uint64_t) hits + 1 <= max_val)
                ++hits;
        };

        auto const & needle = infix(text, i, i + searchParams.length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        _optimalSearchScheme(delegate, delegateDirect, it, needle, scheme, EditDistance());
        c[i] = hits;
//         std::cout << "Needle: " <<  needle << "\n";
//         std::cout << (int)hits << "\n";
    }

    resetLimits(indexText(index), c, searchParams.length);
}
