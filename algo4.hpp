using namespace seqan;

template <unsigned errors, typename TIndex, typename TContainer>
inline void runAlgo4(TIndex & index, auto const & text, unsigned const length, TContainer & c, unsigned const overlap, unsigned const threads)
{
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    auto const & limits = stringSetLimits(indexText(index));

    uint64_t const textLength = seqan::length(text); // lengthSum() forwards to length() for a single string

    unsigned count_skipped_windows = 0;
    unsigned count_forward_positions = 0;

    const uint64_t max_i = textLength - length + 1;
    const uint64_t step_size = length - overlap + 1;
    //#pragma omp parallel for schedule(guided) num_threads(threads)
    #pragma omp parallel for schedule(dynamic, max_i/(step_size*threads*50)) num_threads(threads)
    for (uint64_t i = 0; i < max_i; i += step_size)
    {
        uint64_t max_pos = std::min(i + length - overlap, textLength - length) + 1;

        // overlap is the length of the infix!

        uint64_t leading = 0, trailing = 0;
        for (uint64_t xx = i; c[xx] != 0 && xx < max_pos; ++xx)
            ++leading;
        for (uint64_t xx = max_pos - 1; c[xx] != 0 && xx >= i; --xx) // TODO: i could theoretically be 0 ... overflow because of unsigned value! but c[...] will be zero for i=0 (unless some really weired scheduling happens)
            ++trailing;

        if (trailing != max_pos - i)
        // doesn't work either: trailing != length - overlap + 1 because last interval might be smaller than length - overlap + 1
        // doesn't word: leading != max_pos - i. trailing is computed last and might have found the full range (while another thread writing) to be non-zero while leading didn't find a full range before!
        {
            uint64_t begin_pos = i + leading;
            uint64_t end_pos = max_pos - trailing; // excluding
            uint64_t new_overlap = length - (end_pos - begin_pos) + 1;

            auto scheme = OptimalSearchSchemes<0, errors>::VALUE; // TODO: move out as array
            _optimalSearchSchemeComputeFixedBlocklength(scheme, new_overlap); // only do when new_overlap != overlap

            TIter it_zero_errors[end_pos - begin_pos];
            unsigned hits[end_pos - begin_pos] = {};

            auto delegate = [&hits, &it_zero_errors, begin_pos, length, textLength, new_overlap, &text](auto it, auto const & /*read*/, unsigned const errors_spent) {
                uint64_t const bb = std::min(textLength - 1, begin_pos + length - 1 + length - new_overlap);
                if (errors_spent == 0)
                {
                    extend3<errors>(it, hits, it_zero_errors, errors - errors_spent, text, length,
                        begin_pos + length - new_overlap, begin_pos + length - 1, // searched interval
                        begin_pos, bb // entire interval
                    );
                }
                else
                {
                    extend(it, hits, errors - errors_spent, text, length,
                        begin_pos + length - new_overlap, begin_pos + length - 1, // searched interval
                        begin_pos, bb // entire interval
                    );
                }
            };

            auto const & needle = infix(text, begin_pos + length - new_overlap, begin_pos + length);
            TIter it(index);
            _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
            for (uint64_t j = begin_pos; j < end_pos; ++j)
            {
                if (countOccurrences(it_zero_errors[j - begin_pos]) > 1) // guaranteed to exist, since there has to be at least one match!
                {
                    if (c[j] == 0)
                        count_forward_positions += countOccurrences(it_zero_errors[j - begin_pos]) - 1;
                    for (auto const & occ : getOccurrences(it_zero_errors[j-begin_pos], Fwd()))
                    {
                        auto const occ_pos = posGlobalize(occ, limits);
                        c[occ_pos] = hits[j - begin_pos];
                    }
                }
                else
                {
                    c[j] = hits[j - begin_pos];
                }
            }
        }
        else
        {
            ++count_skipped_windows;
        }

    }
    std::cout << "Forwarded values: " << count_forward_positions << '\n';
    std::cout << "Skipped windows: " << count_skipped_windows << '\n';
}
