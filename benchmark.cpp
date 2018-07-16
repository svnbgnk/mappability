#include "common.h"
#include "common_auxiliary.h"
#include "auxiliary.h"
#include "find2_index_approx_extension.h"
#include <chrono>
#include "global.h"
#include <iostream>

using namespace seqan;
using namespace std;

myGlobalParameters params;
int global;

int main(int argc, char const ** argv)
{
//     params.print();
    
    global = 10;
    testglobal(); 
    
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
    std::vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> bit_vectors = loadBitvectors(bitvectorpath, K, nerrors);
    cout << "Bit vectors loaded. Number: " << bit_vectors.size() << endl;
    
    //start of loop
    
    //wrapp delegate into extension.h ?
    std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hits;
    std::vector<uint8_t> errors_v;
    auto delegate = [&hits, &errors_v](auto & iter, DnaString const & needle, uint8_t errors, bool const rev)
    {
        if(!rev)
        {
            for (auto occ : getOccurrences(iter)){
                hits.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
                errors_v.push_back(errors);
            }
        }else{
            auto const & rgenome = getUniIndexGenome(iter);
            for (auto occ : getOccurrences(iter)){
                cout << endl;
                 if(rev){
                    occ.i2 = seqan::length(rgenome[occ.i1]) - occ.i2 - length(needle);
                }
                hits.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
                errors_v.push_back(errors);
            }
        }
    };
    
    std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hitsD;
    std::vector<uint8_t> errors_vD;
    auto delegateDirect = [&hitsD, &errors_vD](std::vector<Pair<uint16_t, uint32_t>> poss, DnaString const & needle, std::vector<uint8_t> errors)
    {
        for (int i = 0; i < poss.size(); ++i){
            hitsD.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, poss[i]));
            errors_vD.push_back(errors[i]);
        }
    };
    
    
    auto start = std::chrono::high_resolution_clock::now();
//     cout.setstate(std::ios_base::failbit);
    find(0, nerrors, delegate, delegateDirect, index, reads, bit_vectors);
//     std::cout.clear();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "MyVersion elapsed: " << elapsed.count() << "s" << endl;
    
 
    std::cout << "finished" << std::endl;
    
    
    return 0;
}