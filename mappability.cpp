#include <vector>
#include <cstdint>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include <sdsl/int_vector.hpp>

// You can switch between different vector implementations. Consider that they have different thread safetyness!
// typedef sdsl::int_vector<8> TVector;
typedef std::vector<uint8_t> TVector;
// constexpr uint64_t max_val = (1 << 8) - 1;

#include "common.h"

struct Options
{
    unsigned errors;
    unsigned length;
    unsigned overlap = 0;
    unsigned threshold = 0;
    unsigned threads;
    bool mmap;
    bool indels;
    bool singleIndex;
    seqan::CharString indexPath;
    seqan::CharString outputPath;
    seqan::CharString alphabet;
};

#include "algo1.hpp"
#include "algo1_approx.hpp"
#include "algo2.hpp"

using namespace std;
using namespace seqan;

string get_output_path(Options const & opt, signed const chromosomeId)
{
    string output_path = toCString(opt.outputPath);
    output_path += "_" + to_string(opt.errors) + "_" + to_string(opt.length) + "_" + to_string(opt.overlap);
    if (chromosomeId >= 0)
        output_path += "-" + to_string(chromosomeId);

    // for (unsigned i = 0; i < min(c.size(), 64ul); ++i)
    //     cout << (unsigned) c[i] << ' ';
    // cout << '\n';

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
inline void run(TIndex & index, TText const & text, Options const & opt, signed const chromosomeId)
{
    TVector c(seqan::length(text) - opt.length + 1, 0);

    // TODO: is there an upper bound? are we interested whether a k-mer has 60.000 or 70.000 hits?
    cout << mytime() << "Vector initialized (size: " << c.size() << ")." << endl;

    if (opt.overlap == 0 && opt.threshold > 0)
    {
        switch (opt.errors)
        {
            case 0: runAlgo1_approx<0>(index, text, opt.length, c, opt.threads, opt.threshold);
                    break;
            case 1: runAlgo1_approx<1>(index, text, opt.length, c, opt.threads, opt.threshold);
                    break;
            case 2: runAlgo1_approx<2>(index, text, opt.length, c, opt.threads, opt.threshold);
                    break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                     exit(1);
        }
    }
    else if (opt.overlap > 0 && opt.threshold == 0)
    {
        switch (opt.errors)
        {
            case 0: runAlgo2<0>(index, text, opt.length, c, opt.length - opt.overlap, opt.threads);
                    break;
            case 1: runAlgo2<1>(index, text, opt.length, c, opt.length - opt.overlap, opt.threads);
                    break;
            case 2: runAlgo2<2>(index, text, opt.length, c, opt.length - opt.overlap, opt.threads);
                    break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                     exit(1);
        }
    }
    else
    {
        cerr << "Overlap and Threshold are currently mutually exclusive.\n";
        exit(1);
    }
    
    //Sven count number of 0 events
    int counter = 0;
    for (TVector::iterator iter = c.begin() ; iter != c.end(); ++iter)
    {
        if(*iter == 0){
            ++counter;
            *iter = UINT_LEAST8_MAX;
        }

    }  
    cout << "Number of zeroes: " <<  counter << endl;
    cout << mytime() << "Done.\n";

    string output_path = get_output_path(opt, chromosomeId);
    save(c, output_path);

    cout << mytime() << "Saved to disk: " << output_path << '\n';
}

template <typename TChar, typename TAllocConfig, typename TDistance>
inline void run(Options const & opt)
{
    typedef String<TChar, TAllocConfig> TString;
    if (opt.singleIndex)
    {
        typedef StringSet<TString, Owner<ConcatDirect<> > > TStringSet;
        TIndex<TStringSet> index;
        open(index, toCString(opt.indexPath), OPEN_RDONLY);

        auto const & text = indexText(index);

        cout << mytime() << "Index loaded." << endl;
        run<TDistance>(index, text.concat, opt, -1 /*no chromosomeId*/);
    }
    else
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
            cout << mytime() << "Index of " << ids[i] << " loaded." << endl;
            auto const & text = indexText(index);
            run<TDistance>(index, text, opt, i);
        }
    }
}

template <typename TChar, typename TAllocConfig>
inline void run(Options const & opt)
{
    if (opt.indels) {
        cerr << "TODO: Indels are not yet supported.\n";
        exit(1);
        // run<TChar, TAllocConfig, EditDistance>(opt);
    }
    else
        run<TChar, TAllocConfig, HammingDistance>(opt);
}

template <typename TChar>
inline void run(Options const & opt)
{
    if (opt.mmap)
        run<TChar, MMap<> >(opt);
    else
        run<TChar, Alloc<> >(opt);
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

    addOption(parser, ArgParseOption("o", "overlap", "Length of overlap region (o + 1 Strings will be searched at once beginning with their overlap region)", ArgParseArgument::INTEGER, "INT"));
    // setRequired(parser, "overlap");

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", omp_get_max_threads());

    addOption(parser, ArgParseOption("x", "threshold", "Threshold for approximate calculation", ArgParseArgument::INTEGER, "INT"));
    // setDefaultValue(parser, "threshold", "7");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    Options opt;
    getOptionValue(opt.errors, parser, "errors");
    getOptionValue(opt.length, parser, "length");
    getOptionValue(opt.threads, parser, "threads");
    getOptionValue(opt.indexPath, parser, "index");
    getOptionValue(opt.outputPath, parser, "output");
    opt.mmap = isSet(parser, "mmap");
    opt.indels = isSet(parser, "indels");

    if (isSet(parser, "overlap"))
        getOptionValue(opt.overlap, parser, "overlap");
    if (isSet(parser, "threshold"))
        getOptionValue(opt.threshold, parser, "threshold");

    if (opt.overlap > opt.length - opt.errors - 2)
    {
        cerr << "ERROR: overlap should be <= K - E - 2\n";
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
        run<Dna>(opt);
    }
    else
    {
        // run<Dna5>(opt);
        cerr << "TODO: Dna5 alphabet has not been tested yet. Please do and remove this error message afterwards.\n";
        exit(1);
    }
}
