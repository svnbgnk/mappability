using namespace seqan;

template <typename TIter>
inline void extendExact(TIter it, auto hits, auto & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
    constexpr uint64_t max_val = (1 << 8) - 1;

    if (b - a + 1 == length)
    {
        for (auto i = it.fwdIter.vDesc.range.i1; i < it.fwdIter.vDesc.range.i2; ++i)
            hits[a-ab].insert(i);
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
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) >> 1, 0l);
        for (int64_t i = a - 1; i >= a_new; --i)
        {
            if(!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact(it, hits, text, length, a_new, b, ab, bb);
    }
}

// forward
template <typename TIter>
inline void extend(TIter it, auto hits, unsigned errors_left, auto & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb // entire interval
);

// TODO: remove text everywhere: auto & text = indexText(index(it));
template <typename TIter>
inline void approxSearch(TIter it, auto hits, unsigned errors_left, auto & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb, // entire interval
            uint64_t b_new,
            Rev const &
)
{
    if (b == b_new)
    {
        extend(it, hits, errors_left, text, length, a, b, ab, bb);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1]);
                approxSearch(it, hits, errors_left - delta, text, length, a, b + 1, ab, bb, b_new, Rev());
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
template <typename TIter>
inline void approxSearch(TIter it, auto hits, unsigned errors_left, auto & text, unsigned const length,
                  uint64_t a, uint64_t b, // searched interval
                  uint64_t ab, uint64_t bb, // entire interval
                  int64_t a_new,
                  Fwd const & /*tag*/
)
{
    if (a == a_new)
    {
        extend(it, hits, errors_left, text, length, a, b, ab, bb);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1]);
                approxSearch(it, hits, errors_left - delta, text, length, a - 1, b, ab, bb, a_new, Fwd());
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

template <typename TIter>
inline void extend(TIter it, auto hits, unsigned errors_left, auto & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
    constexpr uint64_t max_val = (1 << 8) - 1;

    if (errors_left == 0)
    {
        extendExact(it, hits, text, length, a, b, ab, bb);
        return;
    }
    if (b - a + 1 == length)
    {
        for (auto i = it.fwdIter.vDesc.range.i1; i < it.fwdIter.vDesc.range.i2; ++i)
            hits[a-ab].insert(i);
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
                         Rev()
            );
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) >> 1, 0l);
        approxSearch(it, hits, errors_left, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd()
        );
    }
}

template <unsigned errors, typename TIndex, typename TContainer>
inline void runAlgo2_approx(TIndex & index, auto const & text, unsigned const length, TContainer & c, unsigned const overlap, unsigned const threads, unsigned const threshold_min_hits)
{
    constexpr uint64_t max_val = (1 << 8) - 1;

    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, overlap);

    auto const & limits = stringSetLimits(indexText(index)); // don't call this in a parallel section because it might update limits if it hasn't done it before

    uint64_t const textLength = seqan::length(text); // lengthSum() forwards to length() for a single string

    const uint64_t max_i = textLength - length + 1;
    const uint64_t step_size = length - overlap + 1;
    #pragma omp parallel for schedule(guided) num_threads(threads)
    //#pragma omp parallel for schedule(dynamic, 1000000)
    // TODO: if we choose a multiple of step_size * |cyclic_rotations of int_vector| as chunks (not just chunksize), we dont need to worry about locking
    for (uint64_t i = 0; i < max_i; i += step_size)
    {
        std::vector<std::set<uint32_t>> hits(length - overlap + 1);
        auto delegate = [&hits, i, length, textLength, overlap, &text](auto it, auto const & /*read*/, unsigned const errors_spent) {
            uint64_t const bb = std::min(textLength - 1, i + length - 1 + length - overlap);
            extend(it, hits, errors - errors_spent, text, length,
                i + length - overlap, i + length - 1, // searched interval
                i, bb // entire interval
            );
        };

        auto const & needle = infix(text, i + length - overlap, i + length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());

        uint64_t max_pos = std::min(i + length - overlap, textLength - length);
        for (uint64_t j = i; j <= max_pos; ++j)
        {
            auto const & hit_sa_values = hits[j - i];

            auto no_hits = std::min((uint64_t) hit_sa_values.size(), max_val);
            if (hit_sa_values.size() >= threshold_min_hits /*&& */) // TODO: maybe antoher thread filled other spots? avoid locate
            {
                // assert(hit_sa_values.contains(it.fwdIter.index->sa[i]));
                for(auto h : hit_sa_values)
                {
                    // should even work with strings since stringSetLimits will return Nothing and posGlobalize is overloaded for Strings and Nothing
                    auto const text_pos = it.fwdIter.index->sa[h];
                    auto const text_global_pos = posGlobalize(text_pos, limits);
                    c[text_global_pos] = std::max(static_cast<uint32_t>(c[text_global_pos]), static_cast<uint32_t>(no_hits));
                }
            }
            else
            {
                c[j] = no_hits;
            }
        }




    }
}
