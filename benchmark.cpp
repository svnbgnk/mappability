#include "common.h"
#include "common_auxiliary.h"
#include "auxiliary.h"
#include "find2_index_approx_extension.h"
#include <chrono>
#include "global.h"

using namespace seqan;
using namespace std;

myGlobalParameters params;
int global;

template<typename Tdelegate, typename TdelegateD,
        typename TIndex,
        typename Treads,
        typename TBitvectors>
std::chrono::duration<double> callFunction(int nerrors,
                  vector<hit> & hits,
                  vector<hit> & dhits,
                  Tdelegate & delegate,
                  TdelegateD & delegateDirect,
                  TIndex & index,
                  Treads & reads,
                  TBitvectors & bitvectors)
{
    hits.clear();
    dhits.clear();
    auto start = std::chrono::high_resolution_clock::now();
//     cout.setstate(std::ios_base::failbit);
    find(0, nerrors, delegate, delegateDirect, index, reads, bitvectors);
//     std::cout.clear();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "MyVersion elapsed: " << elapsed.count() << "s" << endl;  

    
    for(int i = 0; i < dhits.size(); ++i){
        hits.push_back(dhits[i]);
    }
    calcfwdPos(index, bitvectors, hits);
    
    return(elapsed);
}

int main(int argc, char const ** argv)
{    
    ArgumentParser parser("benchmark");
    addDescription(parser, "App for search.cpp and finding optimal parameters.");
    
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
    
    
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    
    CharString indexPath, bitvectorpath, readspath;
    string outputpath;
    int K, nerrors, r = 0;
    getOptionValue(indexPath, parser, "index");
    getOptionValue(bitvectorpath, parser, "ibitvector");
    getOptionValue(readspath, parser, "ireads");
    getOptionValue(outputpath, parser, "output");
    getOptionValue(K, parser, "length");
    getOptionValue(nerrors, parser, "errors");
    getOptionValue(r, parser, "r");
    
    
    //load reads
    cout << "Loading reads" << endl;
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

    // load bitvectors
    cout << "Loading bitvectors" << endl;
    std::vector<pair<TBitvector, TSupport>> bitvectors = loadBitvectors(bitvectorpath, K, nerrors);
    cout << "Bit vectors loaded. Number: " << bitvectors.size() << endl;
    cout << "Start my Search" << endl;
    //start of loop
    
    //wrapp delegate into extension.h ?
    std::vector<hit> dhits;
    std::vector<hit> hits;
    auto delegate = [&hits](auto & iter, DnaString const & needle, uint8_t errors, bool const rev)
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
    auto delegateDirect = [&dhits](vector<Pair<uint16_t, uint32_t>> pos, DnaString const & needle, vector<uint8_t> errors)
    {
        for (int i = 0; i < pos.size(); ++i){
            hit me;
            me.occ = pos[i];
            me.read = needle;
            me.errors = errors[i];
            me.rev = false;
            dhits.push_back(me);
        }
    };
    
    
    std::chrono::duration<double> bestTime = params.terminateDuration;
    myGlobalParameters bestParams;
    
    params.clocking = true;
    
    //TODO modify all directSearch parameters otherwise it would be bad to go into another modus this they seem to have the largest influence
/*    
    //16 runs
    for(int i = 0; i < 16; ++i){
        params.wasStopped = false;
        params.normal.setCases(i);
        params.normal.printCases();
        auto time = callFunction(nerrors, hits, dhits, delegate, delegateDirect, index, reads, bitvectors);
        if(bestTime > time){
            bestTime = time;
            bestParams = params;
        }
    }
    
    */


/*   
    // 32 runs
//     bestTime = params.terminateDuration;
    params.normal.setdefault();
    
     
    int intervalsize = 1;
    while(intervalsize < 9){
        float invflipdensity = 0.9;//go lower too 0.1
        while(invflipdensity > 0){
            float filter_th = 0.9; // go lower too 0.1
            while(filter_th > 0){
                params.wasStopped = false;
                params.normal.intervalsize = intervalsize;
                params.normal.invflipdensity = invflipdensity;
                params.normal.filter_th = filter_th;
                cout << intervalsize << "\t" << invflipdensity << "\t" << filter_th << "\t";
                auto time = callFunction(nerrors, hits, dhits, delegate, delegateDirect, index, reads, bitvectors);
                if(bestTime > time){
                    bestTime = time;
                    bestParams = params;
                }
                filter_th -= 0.4;
            }
            
            invflipdensity -= 0.2;
        }
        intervalsize += 1;
    }
    
    
    // 108 runs
//     bestTime = params.terminateDuration;
    params.normal.setdefault();
    int directsearchblockoffset = 1;
    while(directsearchblockoffset < 5){
        int directsearch_th = 2; // go up 4
        while(directsearch_th < 6){
            int distancetoblockend = 1;//// go up to 4
            while(distancetoblockend < 4){
                int step = 2; // //test 8 // 16
                while(step < 9){
                    params.wasStopped = false;
                    params.normal.directsearchblockoffset = directsearchblockoffset;
                    params.normal.directsearch_th = directsearch_th;
                    params.normal.distancetoblockend = distancetoblockend;
                    params.normal.step = step;
                    cout << directsearchblockoffset << "\t" << directsearch_th << "\t" << distancetoblockend << "\t" << step << "\t";
                    auto time = callFunction(nerrors, hits, dhits, delegate, delegateDirect, index, reads, bitvectors);
                    if(bestTime > time){
                        bestTime = time;
                        bestParams = params;
                    }
                    if(step == 2)
                        step += 2;
                    else if(step == 4)
                        step += 4;
                    else
                        step += 8;
                }
                distancetoblockend += 1;
            }
            directsearch_th += 1;
        }
        directsearchblockoffset += 1;
    }
    
    
    */



    cout << "BestParams selected" << endl;
    bestParams.print();
    std::cout << "finished" << std::endl;
    
    
    return 0;
}