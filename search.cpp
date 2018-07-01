#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include "common.h"
#include "common_auxiliary.h"
#include "find2_index_approx_extension.h"
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace seqan;


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

void printPair(pair<uint32_t, uint32_t> p){
    cout << "<" << p.first << ", " << p.second << ">";
}

void printbit(vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, Pair<uint8_t, Pair<uint32_t, uint32_t>> brange){
    sdsl::bit_vector const & rb = bitvectors[brange.i1].first;
    cout << "bitvector: " << (int)brange.i1 << " brange start: " << brange.i2.i1 << "  brange end: " << brange.i2.i2 << endl;
    for(int i = brange.i2.i1; i < brange.i2.i2; ++i)
        cout << i << " Bit: " << rb[i] << endl;
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

int main(int argc, char *argv[])
{
    //load index
    typedef String<Dna, Alloc<>> TString;
    typedef Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> MyIndex;
    MyIndex index;      
    
    
    if(argc == 4){
        open(index, toCString(argv[1]), OPEN_RDONLY);
        Iter<Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig>, VSTree<TopDown<> > > iter(index);
        int chr = stoi(argv[3]);
        print_genome(iter, static_cast<string>(argv[2]), chr);
        cout << "finished writing" << endl;
        exit(0);
        // /home/sven/devel/chr13.fa
    }else{
        open(index, toCString("/home/sven/devel/Data/ref_m_index/index"), OPEN_RDONLY);
        Iter<Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig>, VSTree<TopDown<> > > it(index);
    }

//      open(index, toCString("/home/sven/devel/Data/hg38_test_index/index"), OPEN_RDONLY);
     cout << "Loaded Index. Size:" << seqan::length(index.fwd.sa) << endl;
//     DnaString test = DnaString("ACCAGAACATGATGTGTCGACCGGTATTGAACCAGTCAGT");
    




    
    
//     cout << "; countSequences: " << seqan::countSequences(it.fwdIter.index) << endl;
//      for(int i = 0; i < 4; ++i)
//          cout << it.fwdIter.index->sa[i] << endl;
//     countSequences: 1
//     < 2 , 149 >
//     < 1 , 300 >
//     < 0 , 550 >
//     < 2 , 91 >

    
    // load bitvectors
    vector<pair<sdsl::bit_vector, sdsl::rank_support_v<>>> bit_vectors;
    
//     sdsl::bit_vector b1, b2;
//     load_from_file(b1, "/home/sven/devel/Data/mappability_ref_m.fa/r_bit_vector_100");
//     load_from_file(b2, "/home/sven/devel/Data/mappability_ref_m.fa/l_bit_vector_100");
//     sdsl::rank_support_v<> rb1 (& b1);
//     sdsl::rank_support_v<> rb2 (& b2);
//     bit_vectors.push_back(make_pair(b1, rb1));
//     bit_vectors[0].second.set_vector(&bit_vectors[0].first);
//         bit_vectors.push_back(make_pair(b2, rb2));
    for(int i = 0; i < 10; ++i){
         string file_name = toCString("/home/sven/devel/Data/ref_m_mappa/r_bit_vector_40_shift_") + to_string(i);
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
         string file_name = toCString("/home/sven/devel/Data/ref_m_mappa/l_bit_vector_40_shift_") + to_string(i);
         if(file_exists(file_name)){
             sdsl::bit_vector b;
             cout << "Filename: " << file_name << endl;
             load_from_file(b, file_name);
             sdsl::rank_support_v<> rb(& b);
             bit_vectors.push_back(make_pair(b, rb));
//              bit_vectors[i + 1].second.set_vector(&bit_vectors[i + 1].first);
         }
    }
    cout << "Bit vectors loaded. Size: " << bit_vectors.size() << endl;
    

//     sdsl::rank_support_v<> & rbt = bit_vectors[0].second;
//     rbt.set_vector(&bit_vectors[0].first);
//     cout << "Ranksupport Test: " << rbt(300) << endl;
     
//     String<Dna> read = "TATGGTGCTTAAATGCTCTTGGCTTTCTCCTGCCCACTTAAGGCCTGCCTGCAATTACAAGAGAAACCATTCATACTGGAAATGGTTGCTCTTTGCTGCT";
    
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
    
//     std::set<Pair<DnaString, unsigned> > hits;

    std::vector<Pair<DnaString, Pair <unsigned, unsigned>>> hits;
    std::vector<uint8_t> errors_v;
    std::vector<DnaString> reps;
    
//     auto scheme = OptimalSearchSchemes<0, 2>::VALUE;
//     print_search_scheme(scheme);
//     _optimalSearchSchemeSetMapParams(scheme);
//     cout << "Scheme used in the moment: " << endl;
//     print_search_scheme(scheme);
    
  
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
    
    find<0, 2>(delegate, delegateDirect, index, reads, bit_vectors);
    
//     find<0, 2>(delegate, index, reads, HammingDistance());
    cout << "Hits:" << endl;
    for (int i = 0; i < hits.size(); ++i){
        cout << hits[i] << endl;
        cout << "Errors: " << static_cast<int> (errors_v[i]) << endl;
    }
    cout << "Direct Hits:" << endl;
    for (int i = 0; i < hitsD.size(); ++i){
        cout << hitsD[i] << endl;
        cout << "Errors: " << static_cast<int> (errors_vD[i]) << endl;
    }
    
    // what i searched for
    cout << "Representatives: " << endl;
    for(int i = 0; i < reps.size(); ++i){
        cout << reps[i] << endl;
    }

    cout << "Hello!" << endl;
    
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