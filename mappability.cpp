#include <vector>
#include <cstdint>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include <sdsl/int_vector.hpp>

// You can switch between different vector implementations. Consider that they have different thread safetyness!
// typedef sdsl::int_vector<8> TVector;
typedef std::vector<uint8_t> TVector; // thread-safe
// constexpr uint64_t max_val = (1 << 8) - 1; // not thread-safe

#include "common.h"

struct Options
{
    unsigned errors;
    bool mmap;
    bool indels;
    bool knut;
    bool knut2;
    bool singleIndex;
    seqan::CharString indexPath;
    seqan::CharString outputPath;
    seqan::CharString alphabet;
};

#include "algo2.hpp"
#include "algo3.hpp"
#include "algo4.hpp"

using namespace std;
using namespace seqan;

string get_output_path(Options const & opt, SearchParams const & searchParams, signed const chromosomeId)
{
    string output_path = toCString(opt.outputPath);
    output_path += "_" + to_string(opt.errors) + "_" + to_string(searchParams.length) + "_" + to_string(searchParams.overlap);
    if (opt.knut)
        output_path += "_knut";
    if (opt.knut2)
        output_path += "_knut2";
    if (chromosomeId >= 0)
        output_path += "-" + to_string(chromosomeId);
    return output_path;
}

template <typename T>
inline void save(vector<T> const & c, string const & output_path)
{
    ofstream outfile(output_path, ios::out | ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(TVector::value_type));
    outfile.close();

    // ofstream outfile(output_path, std::ios::out | std::ofstream::binary);
    // copy(c.begin(), c.end(), (std::ostream_iterator<uint8_t>(outfile), std::ostream_iterator<int>(outfile, " ")));
}

template <uint8_t width_t>
inline void save(sdsl::int_vector<width_t> const & c, string const & output_path)
{
    store_to_file(c, output_path);
}

template <typename TDistance, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, Options const & opt, SearchParams const & searchParams, signed const chromosomeId)
{
    TVector c(length(text) - searchParams.length + 1, 0);

    if (opt.knut)
    {
        switch (opt.errors)
        {
            case 0:  runAlgo3<0>(index, text, c, searchParams);
                     break;
            case 1:  runAlgo3<1>(index, text, c, searchParams);
                     break;
            case 2:  runAlgo3<2>(index, text, c, searchParams);
                     break;
            case 3:  runAlgo3<3>(index, text, c, searchParams);
                     break;
            case 4:  runAlgo3<4>(index, text, c, searchParams);
                     break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                     exit(1);
        }
    }
    else if (opt.knut2)
        {
            switch (opt.errors)
            {
                case 0:  runAlgo4<0>(index, text, c, searchParams);
                         break;
                case 1:  runAlgo4<1>(index, text, c, searchParams);
                         break;
                case 2:  runAlgo4<2>(index, text, c, searchParams);
                         break;
                case 3:  runAlgo4<3>(index, text, c, searchParams);
                         break;
                case 4:  runAlgo4<4>(index, text, c, searchParams);
                         break;
                default: cerr << "E = " << opt.errors << " not yet supported.\n";
                         exit(1);
            }
        }
    else
    {
        switch (opt.errors)
        {
            case 0:  runAlgo2<0>(index, text, c, searchParams);
                     break;
            case 1:  runAlgo2<1>(index, text, c, searchParams);
                     break;
            case 2:  runAlgo2<2>(index, text, c, searchParams);
                     break;
            case 3:  runAlgo2<3>(index, text, c, searchParams);
                     break;
            case 4:  runAlgo2<4>(index, text, c, searchParams);
                     break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                     exit(1);
        }
    }

    if (SearchParams::outputProgress)
        std::cout << '\r';
    std::cout << "Progress: 100.00%\n" << std::flush;
    // cout << mytime() << "Done.\n";
    cout.flush();

    string output_path = get_output_path(opt, searchParams, chromosomeId);
    save(c, output_path);
}

template <typename TChar, typename TAllocConfig, typename TDistance>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    typedef String<TChar, TAllocConfig> TString;
    if (opt.singleIndex)
    {
        typedef StringSet<TString, Owner<ConcatDirect<> > > TStringSet;
        TIndex<TStringSet> index;
        open(index, toCString(opt.indexPath), OPEN_RDONLY);
        auto const & text = indexText(index);
        run<TDistance>(index, text.concat, opt, searchParams, -1 /*no chromosomeId*/);
    }
    /*else
    {
        // load chromosome ids
        StringSet<CharString> ids;
        CharString _indexPath = opt.indexPath;
        _indexPath += ".ids";
        open(ids, toCString(_indexPath), OPEN_RDONLY);

        for (unsigned i = 0; i < seqan::length(ids); ++i)
        {
            std::string _indexPath = toCString(opt.indexPath);
            _indexPath += "." + to_string(i);
            TIndex<TString> index;
            open(index, toCString(_indexPath), OPEN_RDONLY);
            auto const & text = indexText(index);
            run<TDistance>(index, text, opt, searchParams, i);
        }
    }*/
}

template <typename TChar, typename TAllocConfig>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.indels) {
        cerr << "TODO: Indels are not yet supported.\n";
        exit(1);
        // run<TChar, TAllocConfig, EditDistance>(opt, searchParams);
    }
    else
        run<TChar, TAllocConfig, HammingDistance>(opt, searchParams);
}

template <typename TChar>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.mmap)
        run<TChar, MMap<> >(opt, searchParams);
    else
        run<TChar, Alloc<> >(opt, searchParams);
}

int main(int argc, char *argv[])
{
    // Argument parser
    ArgumentParser parser("Mappabilty");
    addDescription(parser,
        "App for calculating the mappability values. Only supports Dna4/Dna5 so far.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory (error number, length and overlap will be appended to the output file)", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("i", "indels", "Turns on indels (EditDistance). "
        "If not selected, only mismatches will be considered."));

    addOption(parser, ArgParseOption("z", "knut", "Turns on knuts trick."));
    addOption(parser, ArgParseOption("z2", "knut2", "Turns on knuts trick with considering leading/trailing non-zero values."));

    addOption(parser, ArgParseOption("o", "overlap", "Number of overlapping reads (o + 1 Strings will be searched at once beginning with their overlap region)", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "overlap");

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", omp_get_max_threads());

    // addOption(parser, ArgParseOption("x", "threshold", "Threshold for approximate calculation", ArgParseArgument::INTEGER, "INT"));
    // setDefaultValue(parser, "threshold", "7");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    Options opt;
    SearchParams searchParams;
    getOptionValue(opt.errors, parser, "errors");
    getOptionValue(opt.indexPath, parser, "index");
    getOptionValue(opt.outputPath, parser, "output");
    opt.mmap = isSet(parser, "mmap");
    opt.indels = isSet(parser, "indels");
    opt.knut = isSet(parser, "knut");
    opt.knut2 = isSet(parser, "knut2");

    getOptionValue(searchParams.length, parser, "length");
    getOptionValue(searchParams.threads, parser, "threads");
    getOptionValue(searchParams.overlap, parser, "overlap");

    searchParams.overlap = searchParams.length - searchParams.overlap;

    // if (isSet(parser, "threshold"))
    //     getOptionValue(opt.threshold, parser, "threshold");

    if (opt.knut && opt.knut2)
    {
        cerr << "ERROR: --knut and --knut2 are mutually exclusive.\n";
        exit(1);
    }

    if (searchParams.overlap > searchParams.length - opt.errors - 2)
    {
        cerr << "ERROR: overlap should be <= K - E - 2 (Common overlap has length K-O and will be split into E+2 parts).\n";
        exit(1);
    }

    CharString _indexPath;
    _indexPath = opt.indexPath;
    _indexPath += ".singleIndex";
    open(opt.singleIndex, toCString(_indexPath));

    _indexPath = opt.indexPath;
    _indexPath += ".alphabet";
    open(opt.alphabet, toCString(_indexPath));

    if (opt.alphabet == "dna4")
    {
        run<Dna>(opt, searchParams);
    }
    else
    {
        // run<Dna5>(opt, searchParams);
        cerr << "TODO: Dna5 alphabet has not been tested yet. Please do and remove this error message afterwards.\n";
        exit(1);
    }
}
