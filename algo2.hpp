using namespace seqan;

struct smallHit;

template <typename TIter, typename TOcc, typename TText>
inline void extendExact(TIter it, TOcc * hits, TText const & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
//     constexpr uint64_t max_val = std::numeric_limits<value_type>::max();

    if (b - a + 1 == length)
    {
        //Sven occ
        for (auto occ : getOccurrences(it)){
            smallHit hit(occ);
            hits[a-ab].push_back(hit);
        }
//         hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        return;
    }
    //if (b + 1 <= bb)
    //{
        auto it2 = it;
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            bool success = true;
            for (uint64_t i = b + 1; i <= b_new && success; ++i)
            {
                success = goDown(it2, text[i], Rev());
            }
            if (success)
                extendExact(it2, hits, text, length, a, b_new, ab, bb);
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) >> 1, (int64_t)0);
        for (int64_t i = a - 1; i >= a_new; --i)
        {
            if(!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact(it, hits, text, length, a_new, b, ab, bb);
    }
}

// forward
template <typename TIter, typename TOcc, typename TText, typename TDistanceTag>
inline void extend(TIter it, TOcc * hits, unsigned errors_left, TText const & text,
            unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb, // entire interval
            TDistanceTag const &
);

// TODO: remove text everywhere: auto & text = indexText(index(it));
template <typename TIter, typename TOcc, typename TText, typename TDistanceTag>
inline void approxSearch(TIter it, TOcc * hits, unsigned errors_left, TText const & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb, // entire interval
            uint64_t b_new,
            Rev const &,
            TDistanceTag const &
)
{
    if (b == b_new)
    {
        extend(it, hits, errors_left, text, length, a, b, ab, bb, TDistanceTag());
        return;
    }
    if (errors_left > 0)
    {
        // Insertion Sven
        if(std::is_same<TDistanceTag, EditDistance>::value)
            approxSearch(it, hits, errors_left - 1, text, length, a, b + 1, ab, bb, b_new, Rev(), TDistanceTag());

        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1]);
                approxSearch(it, hits, errors_left - delta, text, length, a, b + 1, ab, bb, b_new, Rev(), TDistanceTag());

                // Deletion Sven
                if(std::is_same<TDistanceTag, EditDistance>::value)
                    approxSearch(it, hits, errors_left - 1, text, length, a, b, ab, bb, b_new, Rev(), TDistanceTag());
            } while (goRight(it, Rev()));
        }
    }
    else
    {
        for (uint64_t i = b + 1; i <= b_new; ++i)
        {
            if (!goDown(it, text[i], Rev()))
                return;
        }
        extendExact(it, hits, text, length, a, b_new, ab, bb);
    }
}
template <typename TIter, typename TOcc, typename TText, typename TDistanceTag>
inline void approxSearch(TIter it, TOcc * hits, unsigned errors_left, TText const & text, unsigned const length,
                  uint64_t a, uint64_t b, // searched interval
                  uint64_t ab, uint64_t bb, // entire interval
                  int64_t a_new,
                  Fwd const & /*tag*/,
                  TDistanceTag const &
)
{
    if (a == a_new)
    {
        extend(it, hits, errors_left, text, length, a, b, ab, bb, TDistanceTag());
        return;
    }
    if (errors_left > 0)
    {
        //Insertion Sven
        if(std::is_same<TDistanceTag, EditDistance>::value)
            approxSearch(it, hits, errors_left - 1, text, length, a - 1, b, ab, bb, a_new, Fwd(), TDistanceTag());


        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1]);
                approxSearch(it, hits, errors_left - delta, text, length, a - 1, b, ab, bb, a_new, Fwd(), TDistanceTag());

                //Deletion Sven
                if(std::is_same<TDistanceTag, EditDistance>::value)
                    approxSearch(it, hits, errors_left - 1, text, length, a, b, ab, bb, a_new, Fwd(), TDistanceTag());
            } while (goRight(it, Fwd()));
        }
    }
    else
    {
        for (int64_t i = a - 1; i >= a_new; --i)
        {
            if (!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact(it, hits, text, length, a_new, b, ab, bb);
    }
}

template <typename TIter, typename TOcc, typename TText, typename TDistanceTag>
inline void extend(TIter it, TOcc * hits, unsigned errors_left, TText const & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb, // entire interval
    TDistanceTag const &
)
{
//     constexpr uint64_t max_val = std::numeric_limits<value_type>::max();

    if (errors_left == 0)
    {
        extendExact(it, hits, text, length, a, b, ab, bb);
        return;
    }
    if (b - a + 1 == length)
    {
        // Sven occ
        for (auto occ : getOccurrences(it)){
            smallHit hit(occ);
            hits[a-ab].push_back(hit);
        }
//         hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        return;
    }
    //if (b + 1 <= bb)
    //{
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            approxSearch(it, hits, errors_left, text, length,
                         a, b, // searched interval
                         ab, bb, // entire interval
                         b_new,
                         Rev(),
                         TDistanceTag()
            );
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) >> 1, (int64_t)0);
        approxSearch(it, hits, errors_left, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd(),
                     TDistanceTag()
        );
    }
}
/*
template <unsigned errors, typename TIndex, typename TText, typename TContainer>
inline void runAlgo2(TIndex & index, TText const & text, TContainer & c, SearchParams const & params)
{
    typedef typename TContainer::value_type value_type;

    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, params.overlap);

    uint64_t const textLength = length(text); // lengthSum() forwards to length() for a single string

    const uint64_t max_i = textLength - params.length + 1;
    const uint64_t step_size = params.length - params.overlap + 1;
   // #pragma omp parallel for schedule(dynamic, std::max(1ul, max_i/(step_size*params.threads*50))) num_threads(params.threads)
    for (uint64_t i = 0; i < max_i; i += step_size)
    {
        value_type hits[params.length - params.overlap + 1] = {};
        auto delegate = [&hits, i, textLength, &params, &text](auto it, auto const & , unsigned const errors_spent) {
            uint64_t const bb = std::min(textLength - 1, i + params.length - 1 + params.length - params.overlap);
            extend(it, hits, errors - errors_spent, text, params.length,
                i + params.length - params.overlap, i + params.length - 1, // searched interval
                i, bb // entire interval
            );
        };

        auto const & needle = infix(text, i + params.length - params.overlap, i + params.length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
        uint64_t max_pos = std::min(i + params.length - params.overlap, textLength - params.length);
        for (uint64_t j = i; j <= max_pos; ++j)
            c[j] = hits[j - i];
    }

    resetLimits(indexText(index), c, params.length);
}*/
