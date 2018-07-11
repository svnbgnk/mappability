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
    if(fwd){
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.fwdIter.index->sa[i].i1;
            int sa = iter.fwdIter.index->sa[i].i2;
            cout << i << ": " << sa + sequenceLengths[seq] << endl;
        }
    }else{
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.revIter.index->sa[i].i1;
            int sa = iter.revIter.index->sa[i].i2;
            cout << i << ": " << sequenceLengths[seq + 1] - sa - 1 << endl;
        }

    }
}


bool occ_smaller(const readOcc& x, const readOcc& y)
{
    if(x.hit.i2.i1 == y.hit.i2.i1)
        return x.hit.i2.i2 < y.hit.i2.i2;
    else
        return x.hit.i2.i1 < y.hit.i2.i1;
}

std::vector<readOcc> print_readocc_sorted(std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hits, std::vector<uint8_t> errors_v, auto const & genome, bool const occEnabled)
{
    std::vector<readOcc> readOccs;
    for(int i = 0; i < hits.size(); ++i){
        readOcc readOcc;
        readOcc.hit = hits[i];
        readOcc.errors = errors_v[i];
        readOccs.push_back(readOcc);
    }
    std::sort(readOccs.begin(), readOccs.end(), occ_smaller);
    
    cout << "Default Hits:" << hits.size() << endl;
    
    for(int i = 0; i < readOccs.size(); ++i){
        cout << "Errors: "<< (int)readOccs[i].errors;
        cout << "   " << readOccs[i].hit << endl;
        if(occEnabled)
            cout << infix(genome[readOccs[i].hit.i2.i1], readOccs[i].hit.i2.i2, readOccs[i].hit.i2.i2 + 100) << endl;
        
    }
    return(readOccs);
}



template <size_t minErrors, size_t maxErrors,
          typename TText, typename TIndexSpec>
int testread(Index<TText, BidirectionalIndex<TIndexSpec> > & index,
              readOcc readOcc)
{
    auto const & genome = indexText(index);
    std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hits;
    std::vector<uint8_t> errors_v;
    auto delegate = [&hits, &errors_v](auto & iter, DnaString const & needle, uint8_t errors)
    {
        for (auto occ : getOccurrences(iter)){
            hits.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
            errors_v.push_back(errors);
        }
    };
    
    StringSet<DnaString> testocc;
    DnaString part = infix(genome[readOcc.hit.i2.i1], readOcc.hit.i2.i2, readOcc.hit.i2.i2 + 100);
    appendValue(testocc, part);
    cout << "Search occ: " << (int)readOcc.hit.i2.i2 << " which has seq: " << endl;
//     cout << part << endl; //TODO revert this

    find<minErrors, maxErrors>(delegate, index, testocc, HammingDistance());
//        print_readocc_sorted(hite, errors_v);
    cout << hits.size() << " hits!!!!!!!!!!" << endl;
    return(hits.size());
}

template <typename TText, typename TIndexSpec>
int testread(int minErrors, int maxErrors, Index<TText, BidirectionalIndex<TIndexSpec> > & index,
              readOcc readOcc){
    int nhits;
    switch (maxErrors)
    {
        case 1: nhits = testread<0, 1>(index, readOcc);
                break;
        case 2: nhits = testread<0, 2>(index, readOcc);
                break;
        case 3: nhits = testread<0, 3>(index, readOcc);
                break;
        default: cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
    return(nhits);
}





template <typename TText, typename TIndexSpec>
vector<int> compare(Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                    int errors,
                    int threshold,
                    std::vector<readOcc> x,
                    std::vector<readOcc> y)
{
    vector<int> wrongHitCount; 
    bool same2 = true;
    bool same = true;
    if(!(x.size() == y.size())){
        same2 = false;
        cout << "MyVersion has: " << x.size() << "hits while default version has: " << y.size() << " hits" << endl;
    }   
    int offset = 0;
    for(int i = 0; i + offset < y.size(); ++i){
        same = false;
        same = (i < x.size() && x[i].hit.i2.i1 == y[i + offset].hit.i2.i1 && x[i].hit.i2.i2 == y[i + offset].hit.i2.i2);
        while(!same && i + offset < y.size()){
            if(i < x.size())//TODO revert this
                cout << "MyVersion has: " << x[i].hit.i2 << " while " ; //TODO revert this
            cout << "default version has: " << y[i + offset].hit.i2 << endl;//TODO revert this
            int nhits = testread(0, errors, index, y[i + offset]);
            if(nhits < threshold){
                cout << "To few hits should have found this part!!!!" << endl;
                wrongHitCount.push_back(nhits);
//                 --offset;
//                 break;
//                 exit(0);
            }
            ++offset;
            same = (i < x.size() && x[i].hit.i2.i1 == y[i + offset].hit.i2.i1 && x[i].hit.i2.i2 == y[i + offset].hit.i2.i2);
        }
        if(i == x.size() && y.size() == i + offset){
            return(wrongHitCount);
        }
    }
    return(wrongHitCount);
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

//Code is now located in find2_index_approx_extension.h
/*
template <typename TIter>
struct isBidirectionalIter
{
     static constexpr bool VALUE = false;
};

template <typename TText, typename TIndex, typename TIndexSpec>
struct isBidirectionalIter<Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > >
{
     static constexpr bool VALUE = true;
};
*/


// template <typename TText, typename TConfig, typename TIndexSpec>
void print_genome(auto it,//Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > it,
                  string const & output_path, 
                  int chr)
{
    if(isBidirectionalIter<decltype(it)>::VALUE){
        cout << "Bidirectional iter" << endl;
    }else{
        cout << "UniDirectional Iter" << endl;
    }/*
    StringSet<DnaString> const & genome = indexText(*it.index);
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
    file.close();*/
}

/*
// template <typename TText, typename TIndex, typename TIndexSpec>
void print_genome(auto it, //Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                  string const & output_path, 
                  int chr)
{
    if(isBidirectionalIter<decltype(it)>::VALUE){
        cout << "Bidirectional Iter" << endl;
    }else{
        cout << "UniDirectional Iter" << endl;
    }
    /*
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
    file.close();*/
//}


vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> loadBitvectors(CharString const bitvectorpath, const int K, const int errors){
    vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> bit_vectors;
    if(file_exists(string("") + toCString(bitvectorpath) + "l_bit_vector_" + to_string(K) + "_" + to_string(errors) + "_shift_" + to_string(0)))
    {
    cout << "Load the following Bitvectors:" << endl;
    for(int i = 0; i < 10; ++i){
        string file_name = string("") + toCString(bitvectorpath) + "r_bit_vector_" + to_string(K) + "_" + to_string(errors) + "_shift_" + to_string(i);
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
         string file_name = string("") + toCString(bitvectorpath) + "l_bit_vector_" + to_string(K) + "_" + to_string(errors) + "_shift_" + to_string(i);
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
    string file_name = string("") + toCString(bitvectorpath) + "right_bit_vector_" + to_string(K) + "_" + to_string(errors);
    if(file_exists(file_name)){
        sdsl::bit_vector b;
        cout << "Filename: " << file_name << endl;
        load_from_file(b, file_name);
        sdsl::rank_support_v<> rb(& b);
        bit_vectors.push_back(make_pair(b, rb));
        }else{
            cerr << file_name << " not found" <<  "\n";
        }

    file_name = string("") + toCString(bitvectorpath) + "middle_bit_vector_" + to_string(K) + "_" + to_string(errors);
    if(file_exists(file_name)){
        sdsl::bit_vector b;
        cout << "Filename: " << file_name << endl;
        load_from_file(b, file_name);
        sdsl::rank_support_v<> rb(& b);
        bit_vectors.push_back(make_pair(b, rb));
    }else{
        cerr << file_name << " not found" <<  "\n";
    }
        
    file_name = string("") + toCString(bitvectorpath) + "left_bit_vector_" + to_string(K) + "_" + to_string(errors);
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
        cout << (string("") + toCString(bitvectorpath) + "l_bit_vector_" + to_string(K) + "_" + to_string(errors) + "_shift_" + to_string(0)) << endl;
        cout << "was not found in the first place maybe wrong k-mer length?" << endl;
        exit(0);
    }
        
    }
    return(bit_vectors);
}



template <typename TText, typename TConfig, typename TIndexSpec>
TText getUniIndexGenome(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > it)
{
    auto const & rgenome = indexText(*it.index);
    return(rgenome);
}

template <typename TText, typename TIndex, typename TIndexSpec>
TText getUniIndexGenome(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it){
    cout << "This thing should not be executed getUniIndexGenome" << endl;
    exit(0);
    auto const & rgenome = indexText(*it.revIter.index);
    return(rgenome);
}




template <typename TText, typename TIndex, typename TIndexSpec>
Pair<uint32_t, uint32_t> getUniRange(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it){
    cout << "This thing should not be executed getUniRange" << endl;
    exit(0);
    Pair<uint32_t, uint32_t> r = it.fwdIter.vDesc.range;
    return(r);
}

// template <typename TText, typename TSpec, typename TConfig>
template <typename TText, typename TConfig, typename TIndexSpec>
Pair<uint32_t, uint32_t> getUniRange(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > it){
    Pair<uint32_t, uint32_t> r = it.vDesc.range;
    return(r);
}


//TODO fix this
template <typename TText, typename TIndex, typename TIndexSpec>
String<unsigned> getUniSa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it){
    cout << "This thing should not be executed getUniSa" << endl;
    exit(0);
    auto const & sa = it.fwdIter.index->sa;
    return(sa);
}

// template <typename TText, typename TSpec, typename TConfig>
template <typename TText, typename TConfig, typename TIndexSpec>
String<unsigned> getUniSa(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > it){
    String<unsigned> const & sa = it.index->sa;
    return(sa);
}



// getUniRange


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
    CharString seqFileName = readspath;
    SeqFileIn seqFileIn(toCString(seqFileName)); //TODO clean this up
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
    cout << "Loading Index" << endl;
    //load index
    typedef String<Dna, Alloc<>> TString;
    typedef StringSet<TString, Owner<ConcatDirect<> > > TText;
    typedef Index<TText, TIndexConfig> MyIndex;
    MyIndex index;      

    open(index, toCString(indexPath), OPEN_RDONLY);
    Iter<Index<TText, TIndexConfig>, VSTree<TopDown<> > > it(index);
    auto iter = it.revIter;

//     print_genome(it, outputpath, 1); //TODO find solution for bidrectional iter test
//     print_genome(iter, outputpath, 1);
    
//      open(index, toCString("/home/sven/devel/Data/hg38_test_index/index"), OPEN_RDONLY);
    cout << "Loaded Index. Size:" << seqan::length(index.fwd.sa) << endl;
    cout << "Loading bitvectors" << endl;
    
    // load bitvectors
    
    vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> bit_vectors = loadBitvectors(bitvectorpath, K, nerrors);
    cout << "Bit vectors loaded. Number: " << bit_vectors.size() << endl;
    
    
    
    int size = seqan::length(it.fwdIter.index->sa);
    uint32_t number_of_indeces = size - bit_vectors[0].first.size();
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[it.fwdIter.index->sa[i].i1 + 1] = it.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    cout << "sequenceLengths: " << endl;
    for(int i = 0; i < sequenceLengths.size(); ++i)
        cout << sequenceLengths[i] << endl;
        
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
    
    //TODO correct hit also for only reverse index
    
/* 
    auto delegate = [&hits, &errors_v](auto & iter, DnaString const & needle, uint8_t errors,  bool const rev)
    {
        for (auto occ : getOccurrences(iter)){
            hits.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
            errors_v.push_back(errors);
        }
//         reps.push_back(representative(iter));
    };*/


    auto delegate = [&hits, &errors_v](auto & iter, DnaString const & needle, uint8_t errors, bool const rev)
    {
        cout << "delegate Call: " << endl;
        if(!rev /* workaround check iter later*/)
        {
//             cout << "Is bidirectional Iter" << endl;
            for (auto occ : getOccurrences(iter)){
                hits.push_back(Pair<DnaString, Pair <unsigned, unsigned>>(needle, occ));
                errors_v.push_back(errors);
            }
        }else{
//             cout << "Is UniDirectional Iter" << endl;
            auto const & rgenome = getUniIndexGenome(iter);
            for (auto occ : getOccurrences(iter)){
                cout << "Seq Size: " << endl;
                cout << seqan::length(rgenome[occ.i1]);
                cout << endl;
                cout << "Occ: "  << occ.i2 << endl;
                 if(rev){
                    cout << "rev case" << endl;
                    occ.i2 = seqan::length(rgenome[occ.i1]) - occ.i2 - length(needle);
                    cout << "Occ after: "  << occ.i2 << endl;
                }/*else{
                    cout << "fwd case" << endl;
                }*/
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
        cout << infix(genome[readOccs[i].hit.i2.i1], readOccs[i].hit.i2.i2, readOccs[i].hit.i2.i2 + 100) << endl;
        
    }

    
//     std::vector<readOcc> readOccsDe;
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

/*    
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
    }*/
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
    
 
    
    
//    appendValue(testocc, "TGAGCGTAATTGTGTCGCGCGCACTGCCTGACTTTTGTT");
/*
        StringSet<DnaString> testocc;
        DnaString part = infix(genome[readOcc.hit.i2.i1], readOcc.hit.i2.i2, readOcc.hit.i2.i2 + 100);
        appendValue(testocc, part);
        cout << "Test occ: " << i << readOcc.hit.i2.i2 << endl;
        hitsDe.clear();
        errors_vDe.clear();
        find<0, 2>(delegateDe, index, testocc, HammingDistance());
//         print_readocc_sorted(hitsDe, errors_vDe);
        cout << hitsDe.size() << endl;*/
    
    

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