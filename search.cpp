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
    cout << "Loaded Index. Size:" << seqan::length(index.fwd.sa) << endl;
    
//     auto iter = it.revIter;
//     print_genome(it, outputpath, 1); //TODO find solution for bidrectional iter test
//     print_genome(iter, outputpath, 1); //TODO use the solution

     // load bitvectors
    cout << "Loading bitvectors" << endl;
    vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> bit_vectors = loadBitvectors(bitvectorpath, K, nerrors);
    cout << "Bit vectors loaded. Number: " << bit_vectors.size() << endl;
  
    std::vector<readOcc> readOccs;
    std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hits;
    std::vector<uint8_t> errors_v;
    
    auto delegate = [&hits, &errors_v](auto & iter, DnaString const & needle, uint8_t errors, bool const rev)
    {
        cout << "delegate Call: " << endl;
        if(!rev ) // workaround check iter later
        {
//             cout << "Is bidirectional Iter" << endl;
            for (auto occ : getOccurrences(iter)){
                cout << "Occ: "  << occ.i2 << endl;
                //TODO replace with accual occurrence
                hits.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
                errors_v.push_back(errors);
            }
        }else{
            cout << "Is UniDirectional Iter" << endl;
            auto const & rgenome = getUniIndexGenome(iter);
            for (auto occ : getOccurrences(iter)){
                cout << endl;
                cout << "Occ: "  << occ.i2 << endl;
                 if(rev){
                    cout << "rev case" << endl;
                    occ.i2 = seqan::length(rgenome[occ.i1]) - occ.i2 - length(needle);
                    cout << "Occ after: "  << occ.i2 << endl;
                }
                hits.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
                errors_v.push_back(errors);
            }
        }
    };
    
      
    std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hitsD;
    std::vector<uint8_t> errors_vD;
    auto delegateDirect = [&hitsD, &errors_vD](vector<Pair<uint16_t, uint32_t>> poss, DnaString const & needle, vector<uint8_t> errors)
    {
        for (int i = 0; i < poss.size(); ++i){
            hitsD.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, poss[i]));
            errors_vD.push_back(errors[i]);
        }
    };
    
    params.startUnidirectional = false;
    
    cout << "Start My Search!" << endl;
    auto start = std::chrono::high_resolution_clock::now();
//     cout.setstate(std::ios_base::failbit);
    find(0, nerrors, delegate, delegateDirect, index, reads, bit_vectors);
//     std::cout.clear();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "Finished My Search" << endl;
    
    for(int i = 0; i < hits.size(); ++i){
        readOcc readOcc;
        readOcc.hit = hits[i];
        readOcc.errors = errors_v[i];
        readOccs.push_back(readOcc);
    }
    
    for(int i = 0; i < hitsD.size(); ++i){
        readOcc readOcc;
        readOcc.hit = hitsD[i];
        readOcc.errors = errors_vD[i];
        readOccs.push_back(readOcc);
    }
    
    std::sort(readOccs.begin(), readOccs.end(), occ_smaller);
    
    
    auto const & genome = indexText(index);
    
    cout << "normal Hits: " << hits.size() << endl;
    cout << "direct Hits: " << hitsD.size() << endl;
    for(int i = 0; i < readOccs.size(); ++i){
        cout << "Errors: "<< (int)readOccs[i].errors;
        cout << "   "  << readOccs[i].hit << endl;
        cout << infix(genome[readOccs[i].hit.i2.i1], readOccs[i].hit.i2.i2, readOccs[i].hit.i2.i2 + seqan::length(readOccs[i].hit.i1)) << endl;
        
    }

    
    cout << "Test default" << endl;
    std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hitsDe;
    std::vector<uint8_t> errors_vDe;
    auto delegateDe = [&hitsDe, &errors_vDe](auto & iter, DnaString const & needle, uint8_t errors)
    {
        for (auto occ : getOccurrences(iter)){
            hitsDe.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
            errors_vDe.push_back(errors);
        }
    };
    start = std::chrono::high_resolution_clock::now();
    find(0, nerrors, delegateDe, index, reads, HammingDistance());
    finish = std::chrono::high_resolution_clock::now();

    std::vector<readOcc> readOccsDe = print_readocc_sorted(hitsDe, errors_vDe, genome, true);
    
    cout << "MyVersion elapsed: " << elapsed.count() << "s" << endl;
    elapsed = finish - start;
    cout << "Default Version elapsed: " << elapsed.count() << "s" << endl;

//     int threshold = 11; 
    int threshold = 3; 
    cout << "Test if default and my version are the same: " << endl;
//     cout.setstate(std::ios_base::failbit); //TODO revert this
    vector<int> whitcount = compare(index, nerrors, threshold, readOccs, readOccsDe);
//     std::cout.clear();  //TODO revert this
    
    if(whitcount.size() == 0)
        cout << "MyVersion is still correct!" << endl;
    else{
        cout << "Missed hits mappability" << endl;
    }
    cout << endl;
    cout << "M: " << endl;
    for(int i = 0; i < whitcount.size(); ++i)
        cout << whitcount[i] << endl;
    cout << endl;
    
    return 0;
    
}


// inside    
//     auto scheme = OptimalSearchSchemes<0, 2>::VALUE;
//     print_search_scheme(scheme);    
//     _optimalSearchSchemeComputeFixedBlocklength(scheme, length(read));
//     print_search_scheme(scheme);
//     Iter<MyIndex, VSTree<TopDown<> > > it(index);
//     unsigned hits = 0;
//     auto delegate = [&hits](auto const &it, auto const & /*read*/, unsigned const /*errors*/) {
//         hits += countOccurrences(it);
//     };
//     auto const & needle = read;
//     goRoot(it);
//     _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
//     cout << hits << endl;