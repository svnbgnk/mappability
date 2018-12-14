using namespace seqan;

struct smallHit{
    Pair <uint16_t, uint32_t> occ;
//     uint8_t errors;
//     DnaString read;
    smallHit(Pair <uint16_t, uint32_t> & inOcc):
        occ(inOcc)
    {}
};

//why cant i sort occ in a custom way?

bool occ_s(const smallHit & x, const smallHit & y)
{
    if(getSeqNo(x.occ) == getSeqNo(y.occ))
        return getSeqOffset(x.occ) < getSeqOffset(y.occ);
    else
        return getSeqNo(x.occ) < getSeqNo(y.occ);

}

template<int disT>
bool occ_sim(const smallHit & x, const smallHit & y)
{
    return(getSeqNo(x.occ) == getSeqNo(y.occ) && getSeqOffset(x.occ) + disT >= getSeqOffset(y.occ) && getSeqOffset(x.occ) <= getSeqOffset(y.occ) + disT);
}

template <unsigned errors, typename TIndex, typename TText, typename TContainer>
inline void runAlgo4(TIndex & index, TText const & text, TContainer & c, SearchParams const & params)
{
    typedef typename TContainer::value_type value_type;

    constexpr uint64_t max_val = std::numeric_limits<value_type>::max();
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    auto const & limits = stringSetLimits(indexText(index));

    uint64_t const textLength = length(text); // lengthSum() forwards to length() for a single string

    const uint64_t max_i = textLength - params.length + 1;
    const uint64_t step_size = params.length - params.overlap + 1;
  //  #pragma omp parallel for schedule(dynamic, std::max(1ul, max_i/(step_size*params.threads*50))) num_threads(params.threads)
    for (uint64_t i = 0; i < max_i; i += step_size)
    {
        uint64_t max_pos = std::min(i + params.length - params.overlap, textLength - params.length) + 1;
//         std::cout << "I: " << i << "\t" << params.length - params.overlap << "\t" << "maxPos " << max_pos << "\n";

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
            uint64_t new_overlap = params.length - (end_pos - begin_pos) + 1;
//             std::cout << "bp " << begin_pos << "\t" << end_pos << "\n";

            auto scheme = OptimalSearchSchemes<0, errors>::VALUE; // TODO: move out as array
            _optimalSearchSchemeComputeFixedBlocklength(scheme, new_overlap); // only do when new_overlap != overlap

            TIter it_zero_errors[end_pos - begin_pos];
//             value_type hits[end_pos - begin_pos] = {};
            std::vector<smallHit> hits[end_pos - begin_pos] = {};

            auto delegate = [&hits, &it_zero_errors, begin_pos, &params, textLength, new_overlap, &text](
                TIter it, auto const & /*read*/, unsigned const errors_spent)
            {
                uint64_t const bb = std::min(textLength - 1, begin_pos + params.length - 1 + params.length - new_overlap);
                if (errors_spent == 0)
                {
                    extend3<errors>(it, hits, it_zero_errors, errors - errors_spent, text, params.length,
                        begin_pos + params.length - new_overlap, begin_pos + params.length - 1, // searched interval
                        begin_pos, bb // entire interval
                    );
                }
                else
                {
                    extend(it, hits, errors - errors_spent, text, params.length,
                        begin_pos + params.length - new_overlap, begin_pos + params.length - 1, // searched interval
                        begin_pos, bb // entire interval
                    );
                }
            };

//             std::cout << "old overlap: " << params.overlap << "\tnew_overlap: " << new_overlap << "\n";
            auto const & needle = infix(text, begin_pos + params.length - new_overlap, begin_pos + params.length);
//             std::cout << "Infix: " << begin_pos + params.length - new_overlap << "\t " << begin_pos + params.length << "\n";
            TIter it(index);
            if(params.indels)
                _optimalSearchScheme(delegate, it, needle, scheme, EditDistance());
            else
                _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
            for (uint64_t j = begin_pos; j < end_pos; ++j)
            {
                value_type cValue;
                if(params.indels){
                    std::vector<smallHit> & chits = hits[j - begin_pos];
                    std::sort(chits.begin(), chits.end(), occ_s);
                    chits.erase(std::unique(chits.begin(), chits.end(), occ_sim<15>), chits.end());
                    cValue = std::min((uint64_t) chits.size(), max_val);
                }
                else
                {
                    cValue = hits[j - begin_pos].size();
                }

                if (countOccurrences(it_zero_errors[j - begin_pos]) > 1) // guaranteed to exist, since there has to be at least one match!
                {
                    for (auto const & occ : getOccurrences(it_zero_errors[j-begin_pos], Fwd()))
                    {
                        auto const occ_pos = posGlobalize(occ, limits);
                        c[occ_pos] = cValue;//hits[j - begin_pos];
                    }
                }
                else
                {
                    c[j] = cValue;//hits[j - begin_pos];
                }
            }
        }
    }

    resetLimits(indexText(index), c, params.length);
}
