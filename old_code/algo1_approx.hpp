using namespace seqan;

template <unsigned errors, typename TIndex, typename TText, typename TContainer>
inline void runAlgo1_approx(TIndex & index, TText const & text, unsigned const length, TContainer & c, unsigned const threads, unsigned const threshold_min_hits)
{
    constexpr uint64_t max_val = (1 << 8) - 1;
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length);

    auto const & limits = stringSetLimits(indexText(index)); // don't call this in a parallel section because it might update limits if it hasn't done it before

    uint64_t const textLength = seqan::length(text);

    #pragma omp parallel for schedule(guided) num_threads(threads)
    for (uint64_t i = 0; i < textLength - length + 1; ++i)
    {
        if (c[i] == 0)
        {
            //unsigned hits = 0;
            std::set<uint32_t> hit_sa_values; // TODO: unordered set?
            auto delegate = [&hit_sa_values](auto const &it, auto const & /*read*/, unsigned const /*errors*/) {
                // if (hits + countOccurrences(it) <= max_val)
                //     hits += countOccurrences(it);
                // else
                //     hits = max_val;
                for (auto i = it.fwdIter.vDesc.range.i1; i < it.fwdIter.vDesc.range.i2; ++i)
                    hit_sa_values.insert(i);
            };

            auto const & needle = infix(text, i, i + length);
            // auto const & needle = extractNeedle(text, i, length);

            Iter<TIndex, VSTree<TopDown<> > > it(index);
            _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
            auto hits = std::min((uint64_t) hit_sa_values.size(), max_val);
            if (hit_sa_values.size() >= threshold_min_hits /*&& */) // TODO: maybe antoher thread filled other spots? avoid locate
            {
                // assert(hit_sa_values.contains(it.fwdIter.index->sa[i]));
                for(auto h : hit_sa_values)
                {
                    // should even work with strings since stringSetLimits will return Nothing and posGlobalize is overloaded for Strings and Nothing
                    auto const text_pos = it.fwdIter.index->sa[h];
                    auto const text_global_pos = posGlobalize(text_pos, limits);
                    c[text_global_pos] = std::max(static_cast<uint32_t>(c[text_global_pos]), static_cast<uint32_t>(hits));
                }
            }
            else
            {
                c[i] = hits;
                // std::cout << "Normal: c[" << i << "] = " << hit_sa_values.size() << '\n';
            }
        }
    }
}
