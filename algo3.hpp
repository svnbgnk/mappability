using namespace seqan;

template <unsigned max_errors, typename TIter>
inline void extendExact3(TIter it, unsigned * hits, TIter * it_zero_errors, auto & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
    constexpr uint64_t max_val = (1 << 8) - 1;

    if (b - a + 1 == length)
    {
        if (max_errors == 0){
            it_zero_errors[a-ab] = it;
            // std::cout << "set zero it for [" << a-ab << "]\n";
        } // TODO: könnte man für max_errors > 0 rausoptimieren
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
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
                extendExact3<max_errors>(it2, hits, it_zero_errors, text, length, a, b_new, ab, bb);
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
        extendExact3<max_errors>(it, hits, it_zero_errors, text, length, a_new, b, ab, bb);
    }
}

// forward
template <unsigned max_errors, typename TIter>
inline void extend3(TIter it, unsigned * hits, TIter * it_zero_errors, unsigned errors_left, auto & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb // entire interval
);

// TODO: remove text everywhere: auto & text = indexText(index(it));
template <unsigned max_errors, typename TIter>
inline void approxSearch3(TIter it, unsigned * hits, TIter * it_zero_errors, unsigned errors_left, auto & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb, // entire interval
            uint64_t b_new,
            Rev const &
)
{
    if (b == b_new)
    {
        extend3<max_errors>(it, hits, it_zero_errors, errors_left, text, length, a, b, ab, bb);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1]);
                approxSearch3<max_errors>(it, hits, it_zero_errors, errors_left - delta, text, length, a, b + 1, ab, bb, b_new, Rev());
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
        extendExact3<max_errors>(it, hits, it_zero_errors, text, length, a, b_new, ab, bb);
    }
}
template <unsigned max_errors, typename TIter>
inline void approxSearch3(TIter it, unsigned * hits, TIter * it_zero_errors, unsigned errors_left, auto & text, unsigned const length,
                  uint64_t a, uint64_t b, // searched interval
                  uint64_t ab, uint64_t bb, // entire interval
                  int64_t a_new,
                  Fwd const & /*tag*/
)
{
    if (a == a_new)
    {
        extend3<max_errors>(it, hits, it_zero_errors, errors_left, text, length, a, b, ab, bb);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1]);
                approxSearch3<max_errors>(it, hits, it_zero_errors, errors_left - delta, text, length, a - 1, b, ab, bb, a_new, Fwd());
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
        extendExact3<max_errors>(it, hits, it_zero_errors, text, length, a_new, b, ab, bb);
    }
}

template <unsigned max_errors, typename TIter>
inline void extend3(TIter it, unsigned * hits, TIter * it_zero_errors, unsigned errors_left, auto & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
    constexpr uint64_t max_val = (1 << 8) - 1;

    if (errors_left == 0)
    {
        extendExact3<max_errors>(it, hits, it_zero_errors, text, length, a, b, ab, bb);
        return;
    }
    if (b - a + 1 == length)
    {
        if (max_errors == errors_left)
            it_zero_errors[a-ab] = it;
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        return;
    }
    //if (b + 1 <= bb)
    //{
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            approxSearch3<max_errors>(it, hits, it_zero_errors, errors_left, text, length,
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
        approxSearch3<max_errors>(it, hits, it_zero_errors, errors_left, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd()
        );
    }
}

template <unsigned errors, typename TIndex, typename TContainer>
inline void runAlgo3(TIndex & index, auto const & text, unsigned const length, TContainer & c, unsigned const overlap, unsigned const threads)
{
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    auto const & limits = stringSetLimits(indexText(index));

    uint64_t const textLength = seqan::length(text); // lengthSum() forwards to length() for a single string

    unsigned count_skipped_windows = 0;
    unsigned count_forward_positions = 0;

    const uint64_t max_i = textLength - length + 1;
    const uint64_t step_size = length - overlap + 1;
    #pragma omp parallel for schedule(guided) num_threads(threads)
    //#pragma omp parallel for schedule(dynamic, 1000/*1000000*/) num_threads(threads)
    // TODO: if we choose a multiple of step_size * |cyclic_rotations of int_vector| as chunks (not just chunksize), we dont need to worry about locking
    for (uint64_t i = 0; i < max_i; i += step_size)
    {

	//if (omp_get_num_threads() < 4)
	//{
	//	std::cout << "x: " << omp_get_num_threads() << '\n';
	//}
        uint64_t max_pos = std::min(i + length - overlap, textLength - length) + 1;
        // std::cout << "max_pos + 1 = " << max_pos << '\n';

        // all are zero if (std::equal(c.begin() + i, c.begin() + max_pos, c.begin() + i) && c[i] == 0)
	
	// overlap is the length of the infix!
	//std::cout << overlap << '\n';

	uint64_t leading = 0, trailing = 0;
	
	for (uint64_t xx = i; c[xx] != 0 && xx < max_pos; ++xx)
		++leading;
	for (uint64_t xx = max_pos - 1; c[xx] != 0 && xx >= i; --xx) // TODO: i could theoretically be 0 ... overflow because of unsigned value! but c[...] will be zero for i=0 (unless some really weired scheduling happens)
		++trailing;

	/*if (trailing != 0 && leading != 0)
	{
		std::cout << "cut off " << leading << " and " << trailing << ": ";
		for (uint64_t xx = i; xx < max_pos; ++xx)
			std::cout << ((c[xx] != 0) ? '1' : '0') << ' ';
		std::cout << '\n';
	}*/

        //if (std::any_of(c.begin() + i, c.begin() + max_pos, [](auto value){ return value == 0; }))
	if (trailing != max_pos - i)
	// doesn't work either: trailing != length - overlap + 1 because last interval might be smaller than length - overlap + 1	
	// doesn't word: leading != max_pos - i. trailing is computed last and might have found the full range (while another thread writing) to be non-zero while leading didn't find a full range before!        
	{
		uint64_t begin_pos = i + leading;
		uint64_t end_pos = max_pos - trailing; // excluding
		uint64_t new_overlap = length - (end_pos - begin_pos) + 1;

	    	auto scheme = OptimalSearchSchemes<0, errors>::VALUE; // TODO: move out as array
	    	_optimalSearchSchemeComputeFixedBlocklength(scheme, new_overlap); // only do when new_overlap != overlap
		//if (std::any_of(c.begin() + i, c.begin() + max_pos, [](auto value){ return value != 0; }))
                //{
		//}

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

		//std::cout << "infix(" << (begin_pos + length - new_overlap) << ", " << (begin_pos + length) << ")\n";
            auto const & needle = infix(text, begin_pos + length - new_overlap, begin_pos + length);
            TIter it(index);
            _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
            for (uint64_t j = begin_pos; j < end_pos; ++j)
            {
                // std::cout << "access it_zero[" << j-i << "]\n";
                if (countOccurrences(it_zero_errors[j - begin_pos]) > 1) // guaranteed to exist, since there has to be at least one match!
                {
                    // ++count_forward_value;
                    if (c[j] == 0)
                        count_forward_positions += countOccurrences(it_zero_errors[j - begin_pos]) - 1;
                    for (auto const & occ : getOccurrences(it_zero_errors[j-begin_pos], Fwd()))
                    {
                        auto const occ_pos = posGlobalize(occ, limits);
                        c[occ_pos] = hits[j - begin_pos]; // tODO: occ pair ...
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
