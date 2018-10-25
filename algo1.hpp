using namespace seqan;

bool test_smaller(int x, int y){
    return(x < y);
}


bool test_same(int x, int y){
    return(x == y);
}


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

        std::vector<hit> myhits;

/*
        auto delegate = [&hits](auto const &it, auto const & , unsigned const ) {
            if ((uint64_t) hits + countOccurrences(it) <= max_val)
                hits += countOccurrences(it);
            else
                hits = max_val;
        };


        auto delegateDirect = [&hits](auto const & pos, DnaString const & , uint8_t const )
        {
            if((uint64_t) hits + 1 <= max_val)
                ++hits;
        };
        */
        auto delegate = [&myhits](auto const &it, auto const & needle, unsigned const occErrors, bool const rev) {
           for (auto occ : getOccurrences(it.fwdIter)){
                hit me;
                me.occ = occ;
                me.errors = occErrors;
                myhits.push_back(me);

            }
        };

        auto delegateDirect = [&myhits](auto const & pos, DnaString const & needle, uint8_t const occErrors)
        {
            hit me;
            me.occ = pos;
            me.errors = occErrors;
            myhits.push_back(me);
        };


        DnaString/*auto*/ const & needle = infix(text, i, i + searchParams.length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);

        _optimalSearchScheme(delegate, delegateDirect, it, needle, scheme, EditDistance());

        std::sort(myhits.begin(), myhits.end(), occ_smaller);

/*
        if(i < 1){
            for(int j = 0; j < myhits.size(); ++j)
                std::cout << myhits[j].occ << "  "  << (int)myhits[j].errors << "\n";
        }*/

        myhits.erase(std::unique(myhits.begin(), myhits.end(), occ_similar<10>), myhits.end());
/*
        if(i < 1){
            for(int j = 0; j < myhits.size(); ++j)
                std::cout << myhits[j].occ << "  "  << (int)myhits[j].errors << "\n";
        }*/


        c[i] = myhits.size();

//         c[i] = hits;

//         std::cout << "Needle: " <<  needle << "\n";
//         std::cout << (int)hits << "\n";
    }

    resetLimits(indexText(index), c, searchParams.length);
}

