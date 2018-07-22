#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <seqan/arg_parse.h>
#include "auxiliary.h"
#include "common_auxiliary.h"
#include "find2_index_approx_extension.h"
#include "global.h"

using namespace std;
using namespace seqan;

myGlobalParameters params;
int global;

int main(int argc, char *argv[])
{
    ArgumentParser parser("Search");
    addDescription(parser, "App for searching reads. Only supports Dna4 so far.");
    
    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");
    
    addOption(parser, ArgParseOption("B", "ibitvector", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");
    
    addOption(parser, ArgParseOption("R", "ireads", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
//     setRequired(parser, "index");
    
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
    
    addOption(parser, ArgParseOption("p", "benchparams",
        "Compare my Version and default version"));
    
    addOption(parser, ArgParseOption("n", "notmy",
        "Compare my Version and default version"));
    
    addOption(parser, ArgParseOption("su", "startuni",
        "Start Unidirectional"));

    
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    
    
    CharString indexPath, bitvectorpath, readspath;
    string outputpath;
    int K, nerrors, threshold = 10, r = 0;
    getOptionValue(indexPath, parser, "index");
    getOptionValue(bitvectorpath, parser, "ibitvector");
    getOptionValue(readspath, parser, "ireads");
    getOptionValue(outputpath, parser, "output");
    getOptionValue(K, parser, "length");
    getOptionValue(nerrors, parser, "errors");
    getOptionValue(r, parser, "r");
    getOptionValue(threshold, parser, "threshold");
    bool mdefault = isSet(parser, "default");
    bool defaultT = isSet(parser, "defaultT");
    bool ecompare = isSet(parser, "ecompare");
    bool benchparams = isSet(parser, "benchparams");
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
//     auto iter = it.revIter;
//     print_genome(it, outputpath, 1); //TODO find solution for bidrectional iter test
//     print_genome(iter, outputpath, 1); //TODO use the solution

     // load bitvectors
    cout << "Loading bitvectors" << endl;
    vector<pair<TBitvector, TSupport>> bitvectors = loadBitvectors(bitvectorpath, K, nerrors);
    cout << "Bit vectors loaded. Number: " << bitvectors.size() << endl;
/*
    cout << "Test Index" << endl;
    // Test index
    vector<uint32_t> sl =  getSequencesLengths(it, bitvectors);
    cout << "Sequence Lengths: " << endl;
    for(int i = 0; i < sl.size(); ++i)
        cout << sl[i] << endl;
    cout << "count Sequences" << endl;
    cout << countSequences(index) << endl;
    cout << "count Sequences from iter" << endl;
    cout << countSequences(*it.fwdIter.index) << endl;
    
    
//     auto mylimits = stringSetLimits(genome);
    auto mylimits = getSeqLengths(index);
    cout << "print limits: " << endl;
    for(int i = 0; i < length(mylimits); ++i)
        cout << mylimits[i] << endl;
    cout << "end of mylimits" << endl;
    Pair<uint16_t, uint32_t> pos(2, 40000);
    auto pos2 = posGlobalize(pos, mylimits); //TODO use this?
    //     hitsv.reserve(countOccurrences(iter.fwdIter)); //TODO fix this
    cout << pos2 << endl;
    cout << endl;
    exit(0);
    */
    auto sll = getSeqLengths(index);
    for(int i = 0; i < countSequences(*it.fwdIter.index); ++i){
        if(sll[i + 1] != length(genome[i])){
            cout << "Wrong Lengths" << endl;
            cout << "sll:" << sll[i + 1] << "gen:" << length(genome[i]) << endl;
            exit(0);
        }
        cout << sll[i+1] << " ";
    }
    cout << endl;
    
    std::vector<hit> dhits;
    std::vector<hit> hits;
    auto delegate = [&hits](auto const & iter, DnaString const & needle, uint8_t errors, bool const rev)
    {
        for (auto occ : getOccurrences(iter)){
            hit me;
            me.occ = occ;
            me.read = needle;
            me.errors = errors;
            me.rev = rev;
            hits.push_back(me);
        }
    };
    auto delegateDirect = [&dhits](Pair<uint16_t, uint32_t> const & pos, DnaString const & needle, uint8_t const errors)
    {
        hit me;
        me.occ = pos;
        me.read = needle;
        me.errors = errors;
        me.rev = false;
        dhits.push_back(me);
    };
    
    if(startuni){
        params.startUnidirectional = true;
    }
    if(benchparams){
        params.normal.setbestnormal();
    }
    
    
    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;
    
    if(!notmy){
        cout << "Start My Search!" << endl;
        start = std::chrono::high_resolution_clock::now();
        find(0, nerrors, delegate, delegateDirect, index, reads, bitvectors);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        cout << "Finished My Search" << endl;
        
        auto scalc = std::chrono::high_resolution_clock::now();
        calcfwdPos(index, hits);
        auto ecalc = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedcalc = ecalc - scalc;
        cout << "Calc revPositions to forward positions: "<< elapsedcalc.count() << "s" << endl;
    }
    
    if(ecompare){
        for(int i = 0; i < dhits.size(); ++i){
            hits.push_back(dhits[i]);
        }
        std::sort(hits.begin(), hits.end(), occ_smaller);
        
        for(int i = 0; i < hits.size(); ++i){
            cout << "Errors: "<< (int)hits[i].errors;
            cout << "   "  << hits[i].occ << " " << hits[i].read << endl;
            cout << infix(genome[hits[i].occ.i1], hits[i].occ.i2, hits[i].occ.i2 + seqan::length(hits[i].read)) << endl;
        }
    }
        
    cout << "MyVersion elapsed: " << elapsed.count() << "s" << endl;
    cout << "normal Hits: (if compare this also includes dhits)" << hits.size() << endl;
    cout << "direct Hits: " << dhits.size() << endl;
    
    // Test default
    std::vector<hit> hitsDe;
    if(mdefault){
        auto delegateDe = [&hitsDe](auto & iter, DnaString const & needle, uint8_t const errors)
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
        cout << "Test default" << endl;
        start = std::chrono::high_resolution_clock::now();
        find(0, nerrors, delegateDe, index, reads, HammingDistance());
        finish = std::chrono::high_resolution_clock::now();
    
        elapsed = finish - start;
        cout << "Default Version elapsed: " << elapsed.count() << "s" << endl;
        cout << "default Hits: " << hitsDe.size() << endl;
    }
    
    // default with in text search
    if(defaultT){
        params.comp.directsearch_th = 5;
        std::vector<hit> hitsDe;
        std::vector<hit> dhitsDe;
        auto delegate2 = [&hitsDe](auto & iter, DnaString const & needle, uint8_t errors, bool const rev)
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
        auto start2 = std::chrono::high_resolution_clock::now();
        find(0, nerrors, delegate2, delegateDirect2, index, reads);
        auto finish2 = std::chrono::high_resolution_clock::now();
        elapsed = finish2 - start2;
        cout << "Default Version with DS: " << elapsed.count() << "s" << endl;
        cout << "default DS Hits: " << hitsDe.size() + dhitsDe.size() << endl;
    }
    
    
    if(ecompare){
        hitsDe = print_readocc_sorted(hitsDe, genome, true);
        cout << "Test if default and my version are the same: " << endl;
//     cout.setstate(std::ios_base::failbit); //TODO revert this
        vector<int> whitcount = compare(index, nerrors, threshold + 1, hits, hitsDe);
//     std::cout.clear();  //TODO revert this
    
        if(whitcount.size() == 0){
            cout << "MyVersion is still correct!" << endl;
        }else{
            cout << "Missed hits mappability" << endl;
        }
        cout << endl;
        cout << "M: " << endl;
        for(int i = 0; i < whitcount.size(); ++i)
            cout << whitcount[i] << endl;
        cout << endl;
    }
 
    return 0;
    
}
