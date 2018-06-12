#include "common.h"
#include <sdsl/bit_vectors.hpp>
#include <iterator> 
//#include <experimental/iterator>

template <unsigned errors, typename TDistance, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, StringSet<CharString> const & ids,
                CharString const & outputPath, unsigned const length, unsigned const chromosomeId, bool const digit)
{
    Iter<TIndex, VSTree<TopDown<> > > it(index);
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    // fill sheme with block length (according to permutation (+ cumulative)) and start position 
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length);

    uint64_t textLength = seqan::length(text);
    // TODO: is there an upper bound? are we interested whether a k-mer has 60.000 or 70.000 hits?
    std::vector<uint16_t> c; // TODO(cpockrandt): check whether this is space efficient. also memory mapping would be better?
    c.resize(textLength - length + 1);
    cout << mytime() << "Vector initialized (size: " << (textLength - length + 1) << ")." << endl;

    // TODO(cpockrandt): think about scheduling strategy
    #pragma omp parallel for schedule(dynamic, 1000000)
    for (uint64_t i = 0; i < textLength - length + 1; ++i)
    {
        unsigned hits = 0;
        // TODO(cpockrandt): for more than 2 errors we need to filter duplicates when counting or
        // choose disjunct search schemes. also we need to do this for EditDistance!
        auto delegate = [&hits](auto const &it, auto const & /*read*/, unsigned const /*errors*/) {
            hits += countOccurrences(it);
        };

        auto const & needle = infix(text, i, i + length);
        goRoot(it);
        _optimalSearchScheme(delegate, it, needle, scheme, TDistance());
        c[i] = hits;
    }
    cout << mytime() << "Done." << endl;

    std::ofstream outfile(toCString(outputPath) + to_string(chromosomeId), std::ios::out | std::ofstream::binary);
    if(!digit)
    {
        int counter = 0;
        vector<float> v(c.begin(), c.end());
        for (std::vector<float>::iterator it = v.begin() ; it != v.end(); ++it)
            if(*it == 0){
                *it = 1;
                counter++;
            }
            else{
                *it = static_cast<float>(1) / *it;
            }
        std::copy(v.begin(), v.end(), (std::ostream_iterator<float>(outfile), std::ostream_iterator<float>(outfile, " ")));
        cout << "Number of zeroes: " <<  counter << endl;
        cout << mytime() << "Final size: " << v.size() << endl;
    }else{
        std::copy(c.begin(), c.end(), std::ostream_iterator<int>(outfile));
    }
    
     
    
    outfile.close();
    cout << mytime() << "Saved to disk." << endl;
}

template <typename TDistance, typename TIndex, typename TText>
inline void run(TIndex /*const*/ & index, TText const & text, StringSet<CharString> const & ids,
                CharString const & outputPath, unsigned const errors, unsigned const length,
                unsigned const chromosomeId, bool const digit)
{
    switch (errors)
    {
        case 0: run<0, TDistance>(index, text, ids, outputPath, length, chromosomeId, digit);
                break;
        case 1: run<1, TDistance>(index, text, ids, outputPath, length, chromosomeId, digit);
                break;
        case 2: run<2, TDistance>(index, text, ids, outputPath, length, chromosomeId, digit);
                break;
        case 3: run<3, TDistance>(index, text, ids, outputPath, length, chromosomeId, digit);
               break;
        case 4: run<4, TDistance>(index, text, ids, outputPath, length, chromosomeId, digit);
               break;
        default: cout << "E = " << errors << " not yet supported." << endl;
                 exit(1);
    }
}

template <typename TChar, typename TAllocConfig, typename TDistance>
inline void run(StringSet<CharString>  & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const singleIndex, bool const digit)
{
    typedef String<TChar, TAllocConfig> TString;
    if (singleIndex)
    {
        Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
        open(index, toCString(indexPath), OPEN_RDONLY);

        // TODO(cpockrandt): replace with a ConcatView
        auto & text = indexText(index);
        typename Concatenator<StringSet<TString, Owner<ConcatDirect<> > >>::Type concatText = concat(text);

        cout << mytime() << "Index loaded." << endl;
        run<TDistance>(index, concatText, ids, outputPath, errors, length, 0, digit);
    }
    else
    {
        for (unsigned i = 0; i < seqan::length(ids); ++i)
        {
            std::string _indexPath = toCString(indexPath);
            _indexPath += "." + to_string(i);
            Index<TString, TIndexConfig> index;
            open(index, toCString(_indexPath), OPEN_RDONLY);
            auto & text = indexText(index);
            cout << mytime() << "Index of " << ids[i] << " loaded." << endl;
            run<TDistance>(index, text, ids, outputPath, errors, length, i, digit);
        }
    }
}

template <typename TChar, typename TAllocConfig>
inline void run(StringSet<CharString> & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const indels, bool const singleIndex, bool const digit)
{
    if (indels)
        run<TChar, TAllocConfig, EditDistance>(ids, indexPath, outputPath, errors, length, singleIndex, digit);
    else
        run<TChar, TAllocConfig, HammingDistance>(ids, indexPath, outputPath, errors, length, singleIndex, digit);
}

template <typename TChar>
inline void run(StringSet<CharString> & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const indels, bool const singleIndex,
                bool const mmap, bool const digit)
{
    if (mmap)
        run<TChar, MMap<> >(ids, indexPath, outputPath, errors, length, indels, singleIndex, digit);
    else
        run<TChar, Alloc<> >(ids, indexPath, outputPath, errors, length, indels, singleIndex, digit);
}

int main(int argc, char *argv[])
{
    // Argument parser
    ArgumentParser parser("Mappability");
    addDescription(parser,
        "App for calculating the mappability values. Only supports Dna4/Dna5 so far.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("i", "indels", "Turns on indels (EditDistance). "
        "If not selected, only mismatches will be considered."));

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("d", "digit",
        "Mappability vector now displays number of occurrences of the k-mers"));
    
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString indexPath, outputPath, _indexPath;
    unsigned errors, length;
    getOptionValue(errors, parser, "errors");
    getOptionValue(length, parser, "length");
    getOptionValue(indexPath, parser, "index");
    getOptionValue(outputPath, parser, "output");
    bool indels = isSet(parser, "indels");
    bool mmap = isSet(parser, "mmap");
    bool digit = isSet(parser, "digit");

    cout << mytime() << "Program started." << endl;

    bool singleIndex;
    CharString alphabet;
    StringSet<CharString> ids;

    _indexPath = indexPath;
    _indexPath += ".singleIndex";
    open(singleIndex, toCString(_indexPath));

    _indexPath = indexPath;
    _indexPath += ".alphabet";
    open(alphabet, toCString(_indexPath));

    _indexPath = indexPath;
    _indexPath += ".ids";
    open(ids, toCString(_indexPath), OPEN_RDONLY);

    if (alphabet == "dna4")
        run<Dna>(ids, indexPath, outputPath, errors, length, indels, singleIndex, mmap, digit);
    else
        run<Dna5>(ids, indexPath, outputPath, errors, length, indels, singleIndex, mmap, digit);
}
