#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include "common.h"
#include "common_auxiliary.h"
#include "find2_index_approx_extension.h"
#include <sdsl/bit_vectors.hpp>
#include <chrono>

using namespace std;
using namespace seqan;

struct readOcc
{
    public:
    Pair<DnaString, Pair <unsigned, unsigned>> hit;
    uint8_t errors;
};

template <typename TText, typename TIndex, typename TIndexSpec>
void print_fullsa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
              bool const fwd)
{
    int size = seqan::length(iter.fwdIter.index->sa);
    uint32_t number_of_indeces = size - bitvectors[0].first.size();
    int noi = number_of_indeces;
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(!fwd){
        for(uint32_t i = 0; i < size; ++i){
            int seq = iter.revIter.index->sa[i].i1;
            int sa = iter.revIter.index->sa[i].i2;
            cout << i << "(" << seq << ")" << ": " << sa + sequenceLengths[seq] << endl;
        }
    }else{
        for(uint32_t i = 0; i < size; ++i){
            int seq = iter.fwdIter.index->sa[i].i1;
            int sa = iter.fwdIter.index->sa[i].i2;
            cout << i << "(" << seq << ")" << ": " << sa + sequenceLengths[seq] << endl;
        }
    }
}



template <typename TText, typename TIndex, typename TIndexSpec>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
              bool const fwd)
{
    int size = seqan::length(iter.fwdIter.index->sa);
    Pair<uint32_t, uint32_t> dirrange = (fwd) ? range(iter.fwdIter) : range(iter.revIter);
    uint32_t number_of_indeces = size - bitvectors[0].first.size();
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(!fwd){
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.revIter.index->sa[i].i1;
            int sa = iter.revIter.index->sa[i].i2;
            cout << i << ": " << iter.revIter.index->sa[i] << endl;
        }
    }else{
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.fwdIter.index->sa[i].i1;
            int sa = iter.fwdIter.index->sa[i].i2;
            cout << i << ": " << sa + sequenceLengths[seq] << endl;
        }
    }
}

bool compare(std::vector<readOcc> x, std::vector<readOcc> y){
    bool same = true;
    if(!(x.size() == y.size())){
        same = false;
        cout << "MyVersion has: " << x.size() << "hits while default version has: " << y.size() << endl;
    }
    
        int offset = 0;
        for(int i = 0; i < x.size(); ++i){
            if(x[i].hit.i2.i1 == y[i].hit.i2.i1 && x[i].hit.i2.i2 == y[i].hit.i2.i2)
                ;
            else{
                cout << "MyVersion has: " << x[i].hit.i2;
                cout << "while default version has: " << y[i].hit.i2 << endl;
                cout << y[i].hit.i1 << endl;
                same = false;
                ++offset;
            }
        }
    return(same);
}


bool occ_smaller(const readOcc& x, const readOcc& y)
{
    if(x.hit.i2.i1 == y.hit.i2.i1)
        return x.hit.i2.i2 < y.hit.i2.i2;
    else
        return x.hit.i2.i1 < y.hit.i2.i1;
}

void printPair(pair<uint32_t, uint32_t> p){
    cout << "<" << p.first << ", " << p.second << ">";
}

void printbit(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, Pair<uint8_t, Pair<uint32_t, uint32_t>> brange){
    sdsl::bit_vector const & rb = bitvectors[brange.i1].first;
    cout << "bitvector: " << (int)brange.i1 << " brange start: " << brange.i2.i1 << "  brange end: " << brange.i2.i2 << endl;
    for(int i = brange.i2.i1; i < brange.i2.i2; ++i)
        cout << i << " Bit: " << rb[i] << endl;
}

sdsl::bit_vector create_random_bit_v(int length){
//     srand (time(0));
    sdsl::bit_vector rv (length, 0);
    srand (time(0));
    for(sdsl::bit_vector::iterator it = rv.begin(); it != rv.end(); ++it){
        if((rand() % 5) == 1)
            *it = 1;
    }
    return(rv);
}

template <typename T> 
void printv(T a){
    for(int i = 0; i < a.size(); ++i){
        cout << static_cast<int> (a.at(i)) << ", ";
    }
    cout << endl;
}

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

template <size_t nbrBlocks, size_t N>
void print_search_scheme(std::array<OptimalSearch<nbrBlocks>, N> & searchsscheme){
    for(int i = 0; i < searchsscheme.size(); i++){
        cout << "Search sscheme: " << i << endl;
        cout << "Permutation: " << endl;
        printv(searchsscheme[i].pi);
        cout << "Lower bound: " << endl;
        printv(searchsscheme[i].l);
        cout << "Upper bound: " << endl;
        printv(searchsscheme[i].u);
        cout << "blockLengths: " << endl;
        printv(searchsscheme[i].blocklength);
        cout << "chronblockLengths: " << endl;
        printv(searchsscheme[i].chronBL);
        cout << "revchronblockLengths: " << endl;
        printv(searchsscheme[i].revChronBL);
        cout << "start Pos: " << endl;
        cout << searchsscheme[i].startPos << endl;
        cout << "minMax: " << endl;
        printv(searchsscheme[i].min);
        printv(searchsscheme[i].max);
        cout << "OneDirection" << endl << (int)searchsscheme[i].startUniDir << endl;
        cout << endl;
    }
}


template <typename TText, typename TIndex, typename TIndexSpec>
void print_genome(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                  string const & output_path, 
                  int chr)
{
    StringSet<DnaString> const & genome = indexText(*it.fwdIter.index);
    ofstream file(output_path, ios::out | ios::binary);
    for(int i = 0; i < chr; ++i){
        file << (">");
        file << to_string(i);
        file << ("\n");
        String<char> target;
        DnaString test = genome[i];
        move(target, test);
        file << target;
        file << ("\n");
    }
    file.close();
}


vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> loadBitvectors(CharString const bitvectorpath, const int K){
    vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> bit_vectors;
    if(file_exists(string("") + toCString(bitvectorpath) + "l_bit_vector_" + to_string(K) + "_shift_" + to_string(0)))
    {
    cout << "Load the following Bitvectors:" << endl;
    for(int i = 0; i < 10; ++i){
        string file_name = string("") + toCString(bitvectorpath) + "r_bit_vector_" + to_string(K) + "_shift_" + to_string(i);
        if(file_exists(file_name)){
            sdsl::bit_vector b;
            cout << "Filename: " << file_name << endl;
            load_from_file(b, file_name);
            sdsl::rank_support_v<> rb(& b);
            bit_vectors.push_back(make_pair(b, rb));
//              bit_vectors[i + 1].second.set_vector(&bit_vectors[i + 1].first);
         }
    }
    
    for(int i = 0; i < 10; ++i){
         string file_name = string("") + toCString(bitvectorpath) + "l_bit_vector_" + to_string(K) + "_shift_" + to_string(i);
         if(file_exists(file_name)){
             sdsl::bit_vector b;
             cout << "Filename: " << file_name << endl;
             load_from_file(b, file_name);
             sdsl::rank_support_v<> rb(& b);
             bit_vectors.push_back(make_pair(b, rb));
//              bit_vectors[i + 1].second.set_vector(&bit_vectors[i + 1].first);
         }
    }
    }else{
    string file_name = string("") + toCString(bitvectorpath) + "right_bit_vector_" + to_string(K);
    if(file_exists(file_name)){
        sdsl::bit_vector b;
        cout << "Filename: " << file_name << endl;
        load_from_file(b, file_name);
        sdsl::rank_support_v<> rb(& b);
        bit_vectors.push_back(make_pair(b, rb));
        }else{
            cerr << file_name << " not found" <<  "\n";
        }

    file_name = string("") + toCString(bitvectorpath) + "middle_bit_vector_" + to_string(K);
    if(file_exists(file_name)){
        sdsl::bit_vector b;
        cout << "Filename: " << file_name << endl;
        load_from_file(b, file_name);
        sdsl::rank_support_v<> rb(& b);
        bit_vectors.push_back(make_pair(b, rb));
    }else{
        cerr << file_name << " not found" <<  "\n";
    }
        
    file_name = string("") + toCString(bitvectorpath) + "left_bit_vector_" + to_string(K);
    if(file_exists(file_name)){
        sdsl::bit_vector b;
        cout << "Filename: " << file_name << endl;
        load_from_file(b, file_name);
        sdsl::rank_support_v<> rb(& b);
        bit_vectors.push_back(make_pair(b, rb));
    }else{
        cerr << file_name << " not found" <<  "\n";
    }
        
    if(bit_vectors.size() != 3){
        exit(0);
    }
        
    }
    return(bit_vectors);
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
    
    addOption(parser, ArgParseOption("r", "r", "number of reads to test ", ArgParseArgument::INTEGER, "INT"));

    
    
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    
    CharString indexPath, bitvectorpath, readspath;
    int K, r = 0;
    getOptionValue(indexPath, parser, "index");
    getOptionValue(bitvectorpath, parser, "ibitvector");
    getOptionValue(readspath, parser, "ireads");
    getOptionValue(K, parser, "length");
    getOptionValue(r, parser, "r");
    
    //load reads
    CharString seqFileName = readspath;
    SeqFileIn seqFileIn(toCString(seqFileName));
    StringSet<CharString> ids;
    StringSet<DnaString> reads;
    readRecords(ids, reads, seqFileIn);
    int nreads = seqan::length(reads);
    cout << "Loaded reads: " << nreads << endl;
    if(r >= nreads){
        cout << "not enought reads" << endl;
        exit(0);
    }
    if(r != 0){
        for(int i = 0; i < nreads - r; ++i)
            eraseBack(reads);
    }
    cout << "Loading Index" << endl;
    //load index
    typedef String<Dna, Alloc<>> TString;
    typedef Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> MyIndex;
    MyIndex index;      

    open(index, toCString(indexPath), OPEN_RDONLY);
    Iter<Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig>, VSTree<TopDown<> > > it(index);

//      open(index, toCString("/home/sven/devel/Data/hg38_test_index/index"), OPEN_RDONLY);
    cout << "Loaded Index. Size:" << seqan::length(index.fwd.sa) << endl;
    cout << "Loading bitvectors" << endl;
    
    // load bitvectors
    
    vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> bit_vectors = loadBitvectors(bitvectorpath, K);
    cout << "Bit vectors loaded. Number: " << bit_vectors.size() << endl;
    
    /*
    //Manuel reads
    
    StringSet<DnaString> reads;
//     String<Dna> read = "TGAGCGTAATTGTGTCGCGCGCACTGCCTGATGTCCGTGATGGGCTTAAGCCTGTCCATCGGCGCATTCTTCATGCGATGAATGAAATGGGACTTTTGTT";
    //                  12345678901234567890123456789012345678901234567890
//     appendValue(reads, "TGAGCGTAATTGTGTCGCGCGCACTGCCTGATGTCCGTGATGGGCTTAAGCGTGTCCATCGGCGCATTCTTCATGCGATGAATGAAATGGGACTTTTGTT");
//     appendValue(reads, "TGAGCGTAATTGTGTCGCGCGCACTGCCTGATGTCCGTGATGGGCTTAAGCCTGTCCATCGGCGCATTCTTCATGCGATGAATGAAATGGGACTTTTGTT"); //GA//GT
//     appendValue(reads, "TGAGCGTAATTGTGTCGCGCGCACTGCCTGATGTCCGTGTTGGGCTTAAGCCTGTCCATCGGCGCATTCTTCATGCGATGAATGAAATGGGACTTTTGTT");  //error pos 38 A-> T
//     appendValue(reads, "CGATCTTACTCGACTACCAGAACATGATGTGTCGACCGGTATTGAACCAGTCAGTATCATTGAAGAAATGCAGTGCTCTTATCTAGATTA"); // repeat
    appendValue(reads, "CGATCTTACTCGACTACCAGAACATGATGTGTCGACCGGTAT");// AT //short repeat start 5
//     appendValue(reads, "ACCAGAACATGATGTGTCGACCGGTATTGAACCAGTCAGT"); //short repeat start 20
//     appendValue(reads, "TGAGCGTAATTGTGTCGCGCGCACTG");
    */

  
    std::vector<readOcc> readOccs;
    std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hits;
    std::vector<uint8_t> errors_v;
//     std::vector<DnaString> reps;  
    auto delegate = [&hits, &errors_v](auto & iter, DnaString const & needle, uint8_t errors)
    {
        for (auto occ : getOccurrences(iter)){
            hits.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
            errors_v.push_back(errors);
        }
//         reps.push_back(representative(iter));
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
        
    cout << "Start My Search!" << endl;
    auto start = std::chrono::high_resolution_clock::now();
//     cout.setstate(std::ios_base::failbit);
    find<0, 2>(delegate, delegateDirect, index, reads, bit_vectors);
    //     std::cout.clear();
    cout << endl;
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "Finished elapsed: " << elapsed.count() << "s" << endl;
    
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
    
    cout << "normal Hits: " << hits.size() << endl;
    cout << "direct Hits: " << hitsD.size() << endl;
    for(int i = 0; i < readOccs.size(); ++i){
        cout << "Errors: "<< (int)readOccs[i].errors;
        cout << "   "  << readOccs[i].hit << endl;
        
    }

    
    std::vector<readOcc> readOccsDe;
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
    find<0, 2>(delegateDe, index, reads, HammingDistance());
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "Finished elapsed: " << elapsed.count() << "s" << endl;
    
    for(int i = 0; i < hitsDe.size(); ++i){
        readOcc readOcc;
        readOcc.hit = hitsDe[i];
        readOcc.errors = errors_vDe[i];
        readOccsDe.push_back(readOcc);
    }
    std::sort(readOccsDe.begin(), readOccsDe.end(), occ_smaller);
    
    cout << "Default Hits:" << hitsDe.size() << endl;
    
    for(int i = 0; i < readOccsDe.size(); ++i){
        cout << "Errors: "<< (int)readOccsDe[i].errors;
        cout << "   " << readOccsDe[i].hit << endl;
        
    }

    cout << "Test if default and my version are the same: " << endl;
    bool same = compare(readOccs, readOccsDe);
    cout << endl << same << endl;
    

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