#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include "auxiliary.h"
#include "common_auxiliary.h"
// #include <seqan/index/find2_index_approx_itv.h>
// #include "find2_index_approx_extension.h"

#include "global.h"
#include <thread>         // std::this_thread::sleep_for

using namespace std;
using namespace seqan;

myGlobalParameters params;
int global;

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     TCondition & itvCondition,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     TDistanceTag const & /**/)
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(delegate,  delegateDirect, itvCondition, it, needle, scheme, TDistanceTag());
}

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag,
          typename TParallelTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     TCondition & itvCondition,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & /**/,
     TParallelTag const & /**/)
{
    typedef typename Iterator<StringSet<TNeedle, TStringSetSpec> const, Rooted>::Type TNeedleIt;
    typedef typename Reference<TNeedleIt>::Type                                       TNeedleRef;
    iterate(needles, [&](TNeedleIt const & needleIt)
    {
        TNeedleRef needle = value(needleIt);
        find<minErrors, maxErrors>(delegate, delegateDirect, itvCondition, index, needle, TDistanceTag());
    },
    Rooted(), TParallelTag());
}

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     TCondition & itvCondition,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & /**/)
{
    find<minErrors, maxErrors>(delegate, delegateDirect, itvCondition, index, needles, TDistanceTag(), Serial());
}

int main(int argc, char *argv[])
{
    ArgumentParser parser("Search");
    addDescription(parser, "App for searching reads. Only supports Dna4 so far.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");

    addOption(parser, ArgParseOption("B", "ibitvector", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");

    addOption(parser, ArgParseOption("R", "ireads", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory", ArgParseArgument::OUTPUT_FILE, "OUT"));
//     setRequired(parser, "output");


    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("E", "errors", "Max errors allowed during mapping", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "errors");

    addOption(parser, ArgParseOption("r", "r", "number of reads to test ", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("d", "default",
        "Test default with in Text Search"));

    addOption(parser, ArgParseOption("dd", "defaultT",
        "Test default with in Text Search"));

    addOption(parser, ArgParseOption("c", "ecompare",
        "Compare my Version and default version"));

    addOption(parser, ArgParseOption("T", "threshold", "Number of times a k-mer can occure (needed for compare)", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("p", "benchparams", "Which parameters set to select", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("n", "notmy",
        "Compare my Version and default version"));

    addOption(parser, ArgParseOption("su", "startuni",
        "Start Unidirectional"));


    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;


    CharString indexPath, bitvectorpath, readspath;
    string outputpath;
    //threshold is needed to check correctness
    int K, nerrors, benchparams, threshold = 0, r = 0;
    getOptionValue(indexPath, parser, "index");
    getOptionValue(bitvectorpath, parser, "ibitvector");
    getOptionValue(readspath, parser, "ireads");
    getOptionValue(outputpath, parser, "output");
    getOptionValue(K, parser, "length");
    getOptionValue(nerrors, parser, "errors");
    getOptionValue(r, parser, "r");
    getOptionValue(threshold, parser, "threshold");
    getOptionValue(benchparams, parser, "benchparams");
    bool mdefault = isSet(parser, "default");
    bool defaultT = isSet(parser, "defaultT");
    bool ecompare = isSet(parser, "ecompare");
    bool notmy = isSet(parser, "notmy");
    bool startuni = isSet(parser, "startuni");

    //load reads
    SeqFileIn seqFileIn(toCString(readspath));
    StringSet<CharString> ids;
    StringSet<DnaString> reads;
    readRecords(ids, reads, seqFileIn);
    int nreads = seqan::length(reads);
    cout << "Loaded reads: " << nreads << endl;
    if(r > nreads){
        cout << "not enought reads" << endl;
        exit(0);
    }
    if(r != 0){
        for(int i = 0; i < nreads - r; ++i)
            eraseBack(reads);
    }


    //load index
    cout << "Loading Index" << endl;
    MyIndex index;
    open(index, toCString(indexPath), OPEN_RDONLY);
    Iter<Index<TText, TIndexConfig>, VSTree<TopDown<> > > it(index);
    auto const & genome = indexText(index);
    cout << "Loaded Index. Size:" << seqan::length(index.fwd.sa) << endl;

/*
    // load bitvectors
    cout << "Loading bitvectors" << endl;
    vector<pair<TBitvector, TSupport>> bitvectors = loadBitvectors(bitvectorpath, K, nerrors);
    cout << "Bit vectors loaded. Number: " << bitvectors.size() << " Length: " << bitvectors[0].first.size() << endl;
*/

    // test OSS
    std::vector<hit> hitsDefault;
    {
    auto delegateDe = [&hitsDefault](auto & iter, DnaString const & needle, uint8_t const errors)
    {
        for (auto occ : getOccurrences(iter)){
            hit me;
            me.occ = occ;
            me.read = needle;
            me.errors = errors;
            me.rev = false;
            hitsDefault.push_back(me);
        }
    };
    cout << "Test default" << endl;
    auto start = std::chrono::high_resolution_clock::now();
    find<0, 2>(delegateDe, index, reads, HammingDistance());
    auto finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = finish - start;
    cout << "Default Version elapsed: " << elapsed.count() << "s" << endl;
    cout << "default Hits: " << hitsDefault.size() << endl;
    }

    // default with in text search
    std::vector<hit> hitsDe;
    std::vector<hit> dhitsDe;
    auto delegate2 = [&hitsDe](auto & iter, DnaString const & needle, uint8_t errors/*, bool const rev*/)
    {
        for (auto occ : getOccurrences(iter)){
            hit me;
            me.occ = occ;
            me.read = needle;
            me.errors = errors;
            me.rev = false;
            hitsDe.push_back(me);
        }
    };
    auto delegateDirect2 = [&dhitsDe](Pair<uint16_t, uint32_t> const & pos, DnaString const & needle, uint8_t const errors)
    {
        hit me;
        me.occ = pos;
        me.read = needle;
        me.errors = errors;
        me.rev = false;
        dhitsDe.push_back(me);
    };


    auto inTextSearchCondition = [](auto iter, uint32_t needleLeftPos, uint32_t needleRightPos, uint8_t errors, auto s, uint8_t const blockIndex)
    {
        return(iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1) * 5);
    };
    cout << "loaded lamda functions" << endl;
    auto start2 = std::chrono::high_resolution_clock::now();
    find<0, 2>(delegate2, delegateDirect2, inTextSearchCondition, index, reads, HammingDistance());
    auto finish2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish2 - start2;
    cout << "Default Version with DS: " << elapsed.count() << "s" << endl;
    cout << "default DS Hits: " << dhitsDe.size() << endl;
    cout << "default overall Hits: " << hitsDe.size() + dhitsDe.size() << endl;


    if(ecompare){
        for(uint32_t i = 0; i < dhitsDe.size(); ++i){
            hitsDe.push_back(dhitsDe[i]);
        }
        std::sort(hitsDe.begin(), hitsDe.end(), occ_smaller);
        //sort default hits
        std::sort(hitsDefault.begin(), hitsDefault.end(), occ_smaller);

//         hitsDe = print_readocc_sorted(hitsDefault, genome, true);
        cout << "Test if default and my version are the same: " << endl;
//     cout.setstate(std::ios_base::failbit); //TODO revert this
        vector<uint32_t> whitcount = compare(index, nerrors, threshold + 1, hitsDe, hitsDefault);
//     std::cout.clear();  //TODO revert this

        if(whitcount.size() == 0){
            cout << "MyVersion is still correct!" << endl;
        }else{
            cout << "Missed hits mappability" << endl;
        }
        cout << endl;
        cout << "M: " << endl;
        for(uint32_t i = 0; i < whitcount.size(); ++i)
            cout << whitcount[i] << endl;
        cout << endl;
    }


    return 0;

}

// ./testIts -I ../Data/chr13_index/index -B ../Data/chr13_tests/chr13_mappa_100/ -R ../Data/reads/chr13/1000k.fa -O ../Data/chr13_tests/chr13_mappa_100/ -K 100 -E 2 -r 10000 -c

// ./testIts -I ../Data/myhg3_index/index -B ../Data/myhg3_mappa/ -R ../Data/reads/myhg3/reads10k.fa -O ../Data/myhg3_mappa/ -K 100 -E 2 -r 10000 -c
