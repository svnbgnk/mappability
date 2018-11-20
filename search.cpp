#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <seqan/arg_parse.h>
#include "auxiliary.h"
#include "common_auxiliary.h"
// #include "newStats.h"
#include "find2_index_approx_extension.h"
// #include "global.h"
#include <thread>         // std::this_thread::sleep_for

#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;

struct modusParameters{
public:
    bool nomappability;
    bool directsearch;
    bool compmappable;
    bool suspectunidirectional;

    bool testflipdensity;
    uint32_t step;
    uint32_t distancetoblockend;
    uint32_t directsearch_th;
    uint32_t directsearchblockoffset;
    float filter_th;
    float invflipdensity;
    uint32_t intervalsize;

    modusParameters(){
        setdefault();
    }

    void setdefault(){
        nomappability = true;
        directsearch = true;
        compmappable = true;
        suspectunidirectional = true;

        testflipdensity = true;
        //binaryNumber //has to be 2^x - 1 for fast modulo calculation
        step = 0b11;
        distancetoblockend = 2;

        directsearchblockoffset = 0;
        directsearch_th = 2;
        filter_th = 0.5;

        invflipdensity = 0.5;

        intervalsize = 3;
    }

    void print(){
        std::cout << "Cases Enabled: " << "\n";
        std::cout << nomappability << " " << directsearch << " " << compmappable << " " << suspectunidirectional << "\n";
        std::cout << "Params: " << "\n";

        std::cout << "step: " << step << "\n";
        std::cout << "distancetoblockend: " << distancetoblockend << "\n";
        std::cout << "directsearchblockoffset: " << directsearchblockoffset << "\n";
        std::cout << "directsearch_th: " << directsearch_th << "\n";
        std::cout << "filter_th: " << filter_th << "\n";
        std::cout << "invflipdensity: " << invflipdensity << "\n";
        std::cout << "intervalsize: " << intervalsize << "\n";
    }
};

// template <typename TTraits>
class OSSContext
{
public:

    //Parameters
    modusParameters normal;
    modusParameters comp;
    modusParameters uni;

    // Shared-memory read-write data.
    std::vector<hit> & hits;
    std::vector<hit> & dhits;
    std::vector<uint32_t> & readOccCount;
//     TMatches &          matches;

    // Shared-memory read-only data.
//     TContigSeqs const & contigSeqs;
    bool filterDelegate = true;

    OSSContext(std::vector<hit> & inhits,
               std::vector<hit> & indhits,
               std::vector<uint32_t> & inreadOccCount) :
        hits(inhits),
        dhits(indhits),
        readOccCount(inreadOccCount)
    {
        ;
    }

    void setdefault(){
        normal.setdefault();
        comp.setdefault();
        uni.setdefault();
    }

    void print(){
        std::cout << "Normal: ";
        normal.print();
        std::cout << "Comp: ";
        comp.print();
        std::cout << "Uni: ";
        uni.print();
    }

    template <size_t nbrBlocks>
    bool itvCondition(OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex,
                      uint32_t ivalOne)
    {
        return(ivalOne < (static_cast<int>(s.pi.size()) - blockIndex - 1 + normal.directsearchblockoffset) * normal.directsearch_th);
    }


    template <typename TText, typename TIndex, typename TIndexSpec,
              size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                      uint32_t const needleLeftPos,
                      uint32_t const needleRightPos,
                      uint8_t const errors,
                      OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex)
    {
        return(iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1) * 5);
    }

    template<size_t nbrBlocks>
    bool inBlockCheckMappabilityCondition(uint32_t needleLeftPos,
                                          uint32_t needleRightPos,
                                           OptimalSearch<nbrBlocks> const & s,
                                          uint8_t blockIndex)
    {
        uint32_t prevBlocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
        uint32_t nextBlocklength = s.blocklength[blockIndex];
        uint32_t step = (needleRightPos - needleLeftPos - 1);


        bool enoughDistanceToBlockEnds = step + normal.distancetoblockend < nextBlocklength && step - normal.distancetoblockend > prevBlocklength;
        return(((step & normal.step) == 0) && enoughDistanceToBlockEnds);
    }


};

StringSet<DnaString> createRCReads(StringSet<DnaString> & reads)
{
    StringSet<DnaString> rcReads;
    for(int i = 0; i < length(reads); ++i){
        //ModifiedString<ModifiedString<DnaString, ModComplementDna>, ModReverse> myModifier(reads[1]);
        DnaStringReverseComplement myModifier(reads[i]);
        appendValue(rcReads, myModifier);
    }
    return(rcReads);
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
    addOption(parser, ArgParseOption("rc", "rc",
        "Search on both strands"));
    addOption(parser, ArgParseOption("fr", "fr",
        "Create fastas with filtered reads in output directory"));

        addOption(parser, ArgParseOption("sp", "split",
        "split reads"));

    addOption(parser, ArgParseOption("c", "ecompare",
        "Compare my Version and default version"));

    addOption(parser, ArgParseOption("T", "threshold", "Number of times a k-mer can occure (needed for compare)", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("p", "benchparams", "Which parameters set to select", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("n", "notmy",
        "Compare my Version and default version"));

    addOption(parser, ArgParseOption("st", "stats",
        "Show stats for Default Search with ITV and Search with mappability"));

    addOption(parser, ArgParseOption("su", "startuni",
        "Start Unidirectional"));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;


    CharString indexPath, bitvectorpath, readspath;
    string outputpath;
    int K, nerrors, benchparams, threshold = 10, r = 0;
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
    bool rc = isSet(parser, "rc");
    bool fr = isSet(parser, "fr");
    bool split = isSet(parser, "sp");
    bool ecompare = isSet(parser, "ecompare");
    bool notmy = isSet(parser, "notmy");
    bool stats = isSet(parser, "stats");
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

    if(split){
        auto startT = std::chrono::high_resolution_clock::now();
        int readlength = length(reads[0]);
        for(int i = 0; i < length(reads); ++i){
//             appendValue(rcReads, infix(reads[i], readlength/2 + readlength));
            reads[i] = infix(reads[i], 0, 0 + readlength/2);
        }
        cout << "read length: ";
        cout << length(reads[0]) << endl;
        auto finishT = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedT = finishT - startT;
        cout << "splitted reads in " << elapsedT.count() << "s" << endl;
    }



    if(rc){
        //add reverse Complement of reads to search on both strands
        auto startT = std::chrono::high_resolution_clock::now();
        StringSet<DnaString> rcreads = createRCReads(reads);
        for(int i = 0; i < length(rcreads); i++){
            appendValue(reads, rcreads[i]);
        }
        auto finishT = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedT = finishT - startT;
        cout << "added reversed Reads in " << elapsedT.count() << "s" << endl;
    }

    //load index
    cout << "Loading Index" << endl;
    MyIndex index;
    open(index, toCString(indexPath), OPEN_RDONLY);
    //TODO load index with auto?
    Iter<Index<TText, TIndexConfig>, VSTree<TopDown<> > > it(index);
    auto const & genome = indexText(index);
    cout << "Loaded Index. Size:" << seqan::length(index.fwd.sa) << endl;

    // load bitvectors
    vector<pair<TBitvector, TSupport>> bitvectors;
    if(!notmy){
        cout << "Loading bitvectors" << endl;
        bitvectors = loadBitvectors(bitvectorpath, K, nerrors);
        cout << "Bit vectors loaded. Number: " << bitvectors.size() << " Length: " << bitvectors[0].first.size() << endl;
    }


    std::vector<hit> myhits;
    std::vector<hit> mydhits;
    std::vector<uint32_t> myreadOccCount;
    OSSContext myOSSContext(myhits, mydhits, myreadOccCount);

    auto delegate = [&myhits](auto const & iter, DnaString const & needle, uint32_t const needleId, uint8_t errors, bool const rev)
    {
        //NOTE have to get Occurrences from the forward iter otherwise filtering does not work properly
        for (auto occ : getOccurrences(iter)){
            hit me;
            me.occ = occ;
            me.read = needle;
            me.readId = needleId;
            me.errors = errors;
            me.rev = rev;
            myhits.push_back(me);
        }
    };
    auto delegateDirect = [&mydhits](Pair<uint16_t, uint32_t> const & pos, DnaString const & needle, uint32_t const needleId, uint8_t const errors)
    {
        hit me;
        me.occ = pos;
        me.read = needle;
        me.readId = needleId;
        me.errors = errors;
        me.rev = false;
        mydhits.push_back(me);
    };

//     std::this_thread::sleep_for (std::chrono::seconds(60));

    //TODO tranlate this for ossContext
    /*
    if(startuni){
        params.startUnidirectional = true;
    }
    switch(benchparams){
        case 1:
        {
            params.normal.setbestnormal();
            params.copyDirectsearchParamsfromNormal();
            break;
        }
        case 2:
        {
            params.normal.setbestnormalhg();
            params.copyDirectsearchParamsfromNormal();
            break;
        }
        case 3:
        {
            params.normal.setbestnormalhgE2();
            params.copyDirectsearchParamsfromNormal();
            break;
        }
        case 4:
        {
            params.normal.setbestnormalhgE3();
            params.copyDirectsearchParamsfromNormal();
            break;
        }
        case 5:
        {
            params.normal.setbestnormalhgE3();
            params.copyDirectsearchParamsfromNormal();
            params.startuni.setbestStartUnihgE3();
            break;
        }
        default:
        {
            break;
        }
    }
*/


    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;

    if(!notmy){
        cout << "Start My Search!" << endl;
        start = std::chrono::high_resolution_clock::now();
        find(0, nerrors, myOSSContext, delegate, delegateDirect, index, reads, bitvectors, HammingDistance());
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        cout << "Finished My Search" << endl;

        auto scalc = std::chrono::high_resolution_clock::now();
        calcfwdPos(index, myhits);
        auto ecalc = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedcalc = ecalc - scalc;
        cout << "Calc revPositions to forward positions: "<< elapsedcalc.count() << "s" << endl;
    }

    cout << "DirectHits: " << mydhits.size() << endl;

    if(ecompare){
        for(uint32_t i = 0; i < mydhits.size(); ++i){
            myhits.push_back(mydhits[i]);
        }
        std::sort(myhits.begin(), myhits.end(), occ_smaller);

        for(uint32_t i = 0; i < myhits.size(); ++i){
            cout << "Errors: "<< (uint32_t)myhits[i].errors;
            cout << "   "  << myhits[i].occ << " " << myhits[i].read << endl;
            cout << infix(genome[myhits[i].occ.i1], myhits[i].occ.i2, myhits[i].occ.i2 + seqan::length(myhits[i].read)) << endl;
        }
    }

    cout << "MyVersion elapsed: " << elapsed.count() << "s" << endl;
    cout << "normal Hits: (if compare this also includes mydhits)" << myhits.size() << endl;
    cout << "direct Hits: " << mydhits.size() << endl;

    // Test default
    //TODO change vector name of lambda function

    std::vector<hit> hitsDefault;
    if(mdefault){
        auto delegateDefault = [&hitsDefault](auto & iter, DnaString const & needle, uint8_t const errors)
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
        start = std::chrono::high_resolution_clock::now();
        find(0, nerrors, delegateDefault, index, reads, HammingDistance());
        finish = std::chrono::high_resolution_clock::now();

        elapsed = finish - start;
        cout << "Default Version elapsed: " << elapsed.count() << "s" << endl;
        cout << "default Hits: " << hitsDefault.size() << endl;
    }

/*
    // default with in text search
    if(defaultT){
        params.comp.directsearch_th = 5;
        params.comp.directsearchblockoffset = 0;
//         std::vector<hit> hitsDe;
//         std::vector<hit> dhitsDe;
        auto delegate2 = [](auto & iter, DnaString const & needle, uint8_t errors, bool const rev)
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
        auto delegateDirect2 = [](Pair<uint16_t, uint32_t> const & pos, DnaString const & needle, uint8_t const errors)
        {
            hit me;
            me.occ = pos;
            me.read = needle;
            me.errors = errors;
            me.rev = false;
            dhitsDe.push_back(me);
        };
        auto start2 = std::chrono::high_resolution_clock::now();
        find(0, nerrors, delegate2, delegateDirect2, index, reads, HammingDistance());
        auto finish2 = std::chrono::high_resolution_clock::now();
        elapsed = finish2 - start2;
        cout << "Default Version with DS: " << elapsed.count() << "s" << endl;
        cout << "default DS Hits: " << hitsDe.size() + dhitsDe.size() << endl;
    }
*/
/*
    if(!notmy){
        readOccurrences(reads, ids, outputpath, fr, rc, stats, notmy);
    }*/



    if(ecompare){
        hitsDefault = print_readocc_sorted(hitsDefault, genome, true);
        cout << "Test if default and my version are the same: " << endl;
        cout.setstate(std::ios_base::failbit);
        vector<uint32_t> whitcount = compare(index, nerrors, threshold + 1, myhits, hitsDefault);
        std::cout.clear();

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
//     params.print();

    return 0;

}
