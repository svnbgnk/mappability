#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <seqan/arg_parse.h>
#include "auxiliary.h"
#include "common_auxiliary.h"
#include "find2_index_approx_extension.h"
#include "global.h"
#include <thread>         // std::this_thread::sleep_for

#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;

myGlobalParameters params;
std::vector<hit> hits;
std::vector<hit> dhits;
std::vector<hit> hitsDe;
std::vector<hit> dhitsDe;
std::vector<uint32_t> readOccCount;
std::vector<uint32_t> readOccCountDeT;

vector<uint32_t> histogram(vector<uint32_t> & b , int const his_size, int const bucket_width)
{
    vector<uint32_t> hist(his_size, 0);
    auto it = b.begin();
    while(it != b.end()){
        if(*it < ((his_size) * bucket_width))
        {
            hist[floor(*it / bucket_width)] += 1;
        }
        else
        {
            ++hist[his_size - 1];
        }
        ++it;
    }
    return(hist);
}

/*
template <typename T>
inline void save(vector<T> const & c, string const & output_path)
{
    ofstream outfile(output_path, ios::out | ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(TVector::value_type));
    outfile.close();

    // ofstream outfile(output_path, std::ios::out | std::ofstream::binary);
    // copy(c.begin(), c.end(), (std::ostream_iterator<uint8_t>(outfile), std::ostream_iterator<int>(outfile, " ")));
}
*/


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

    addOption(parser, ArgParseOption("c", "ecompare",
        "Compare my Version and default version"));

    addOption(parser, ArgParseOption("T", "threshold", "Number of times a k-mer can occure (needed for compare)", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("p", "benchparams", "Which parameters set to select", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("n", "notmy",
        "Compare my Version and default version"));

    addOption(parser, ArgParseOption("su", "startuni",
        "Start Unidirectional"));;


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
    cout << "Loading bitvectors" << endl;
    vector<pair<TBitvector, TSupport>> bitvectors = loadBitvectors(bitvectorpath, K, nerrors);
    cout << "Bit vectors loaded. Number: " << bitvectors.size() << " Length: " << bitvectors[0].first.size() << endl;

//     std::vector<hit> dhits;
//     std::vector<hit> hits;
    auto delegate = [/*&hits*/](auto const & iter, DnaString const & needle, uint8_t errors, bool const rev)
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
    auto delegateDirect = [/*&dhits*/](Pair<uint16_t, uint32_t> const & pos, DnaString const & needle, uint8_t const errors)
    {
        hit me;
        me.occ = pos;
        me.read = needle;
        me.errors = errors;
        me.rev = false;
        dhits.push_back(me);
    };

//     std::this_thread::sleep_for (std::chrono::seconds(20));

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

    cout << "DirectHits: " << dhits.size() << endl;

    if(ecompare){
        for(uint32_t i = 0; i < dhits.size(); ++i){
            hits.push_back(dhits[i]);
        }
        std::sort(hits.begin(), hits.end(), occ_smaller);

        for(uint32_t i = 0; i < hits.size(); ++i){
            cout << "Errors: "<< (uint32_t)hits[i].errors;
            cout << "   "  << hits[i].occ << " " << hits[i].read << endl;
            cout << infix(genome[hits[i].occ.i1], hits[i].occ.i2, hits[i].occ.i2 + seqan::length(hits[i].read)) << endl;
        }
    }

    cout << "MyVersion elapsed: " << elapsed.count() << "s" << endl;
    cout << "normal Hits: (if compare this also includes dhits)" << hits.size() << endl;
    cout << "direct Hits: " << dhits.size() << endl;

    // Test default
    //TODO change vector name of lambda function
    /*
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
    */

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
        find(0, nerrors, delegate2, delegateDirect2, index, reads);
        auto finish2 = std::chrono::high_resolution_clock::now();
        elapsed = finish2 - start2;
        cout << "Default Version with DS: " << elapsed.count() << "s" << endl;
        cout << "default DS Hits: " << hitsDe.size() + dhitsDe.size() << endl;
    }


    int found = 0, foundD = 0;
    int notfound = 0, mymiss = 0, same = 0, nice = 0, verynice = 0;
    vector<uint8_t> mycase(length(reads), 255);




    if(rc){
        int nr = readOccCount.size()/2;
        cout << "Jump:" << nr << endl;
        for(int i = 0; i < readOccCountDeT.size()/2; ++i)
        {
            readOccCount[i] += readOccCount[i + nr];
            readOccCountDeT[i] += readOccCountDeT[i + nr];

        }
        readOccCount.erase(readOccCount.begin() + nr, readOccCount.end());
        readOccCountDeT.erase(readOccCountDeT.begin() + nr, readOccCountDeT.end());
    }


    cout << "Test " << readOccCount.size() << endl;
    for(int i = 0; i < readOccCount.size(); ++i)
    {
        found += readOccCount[i] > 0;
        foundD += readOccCountDeT[i] > 0;
//         cout << readOccCount[i] << " - " <<  readOccCountDeT[i] << endl;
        if(readOccCount[i] == readOccCountDeT[i]){
            if(readOccCount[i] == 0){
                ++notfound;
                mycase[i] = 0;
            }
            else
            {
                ++same;
                mycase[i] = 2;
            }
        }
        else
        {
            if(readOccCount[i] == 0){
                ++mymiss;
                mycase[i] = 1;
            }
            else
            {
                ++nice;
                if(static_cast<double>(readOccCountDeT[i]) / readOccCount[i] > 2)
                    ++verynice;

                mycase[i] = 3;
                if(readOccCount[i] > readOccCountDeT[i]){
                    cerr << "More occurrences with mappability" << endl;
                    exit(0);
                }
            }
        }

    }


    if(fr){
        //write filtered fastas
        SeqFileOut seqFileout0(toCString(outputpath + "/notfound.fa"));
        SeqFileOut seqFileout1(toCString(outputpath + "/mymiss.fa"));
        SeqFileOut seqFileout2(toCString(outputpath + "/same.fa"));
        SeqFileOut seqFileout3(toCString(outputpath + "/nice.fa"));


        for(int i = 0; i < length(reads); ++i)
        {
            switch(mycase[i])
            {
                case 0: writeRecord(seqFileout0, ids[i], reads[i]); break;
                case 1: writeRecord(seqFileout1, ids[i], reads[i]); break;
                case 2: writeRecord(seqFileout2, ids[i], reads[i]); break;
                case 3: writeRecord(seqFileout3, ids[i], reads[i]); break;
                default: break;
            }

        }

        close(seqFileout0);
        close(seqFileout1);
        close(seqFileout2);
        close(seqFileout3);
    }
    cout << "reads found with mappability: " << found << endl;
    cout << "reads found: " << foundD << endl;
    cout << "not found: " << notfound << endl;
    cout << "mymiss: " << mymiss << endl;
    cout << "same: " << same << endl;
    cout << "nice: " << nice << endl;
    cout << "thereof verynice: " << verynice << endl;



    // investigating the vectors
    int bucketSize = 10;
    int histSize = 10;
    vector<uint32_t> h = histogram(readOccCount, histSize, bucketSize);
    vector<uint32_t> hDeT = histogram(readOccCountDeT, histSize, bucketSize);


    cout << "Histogram buckets size " << bucketSize << ": " << endl;
    for(int i = 0; i < h.size(); ++i){
        cout << bucketSize*(i + 1) - 1 << "\t";
    }
    cout << endl;
    for(int i = 0; i < h.size(); ++i){
        cout << h[i] << "\t";
    }
    cout << endl;


    cout << "Histogram buckets size " << bucketSize << ": " << endl;
    for(int i = 0; i < hDeT.size(); ++i){
        cout << bucketSize*(i + 1) - 1 << "\t";
    }
    cout << endl;
    for(int i = 0; i < hDeT.size(); ++i){
        cout << hDeT[i] << "\t";
    }
    cout << endl;


    /*

    if(ecompare){
        hitsDe = print_readocc_sorted(hitsDe, genome, true);
        cout << "Test if default and my version are the same: " << endl;
//     cout.setstate(std::ios_base::failbit); //TODO revert this
        vector<uint32_t> whitcount = compare(index, nerrors, threshold + 1, hits, hitsDe);
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
//     params.print();

*/
    return 0;

}
