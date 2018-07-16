#ifndef AUXILLARY_H_
#define AUXILLARY_H_

#include "common.h"
#include <sdsl/bit_vectors.hpp>
using namespace seqan;

std::vector<std::pair<uint32_t, uint32_t>> getConsOnes(std::vector<std::pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, 
                                             Pair<uint8_t, Pair<uint32_t, uint32_t>> inside_bit_interval,
                                             int const intervalsize);

struct readOcc
{
    public:
    Pair<DnaString, Pair <unsigned, unsigned>> hit;
    uint8_t errors;
};

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


template <typename TText, typename TIndex, typename TIndexSpec>
void print_fullsa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              std::vector<std::pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
              bool const fwd)
{
    int size = seqan::length(iter.fwdIter.index->sa);
    uint32_t number_of_indeces = size - bitvectors[0].first.size();
    int noi = number_of_indeces;
    std::vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(!fwd){
        for(uint32_t i = 0; i < size; ++i){
            int seq = iter.revIter.index->sa[i].i1;
            int sa = iter.revIter.index->sa[i].i2;
            std::cout << i << "(" << seq << ")" << ": " << sa + sequenceLengths[seq] << "\n";
        }
    }else{
        for(uint32_t i = 0; i < size; ++i){
            int seq = iter.fwdIter.index->sa[i].i1;
            int sa = iter.fwdIter.index->sa[i].i2;
            std::cout << i << "(" << seq << ")" << ": " << sa + sequenceLengths[seq] << "\n";
        }
    }
}


template <typename TText, typename TIndex, typename TIndexSpec>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              int const number_of_indeces,
              bool const fwd)
{
    Pair<uint32_t, uint32_t> dirrange = (fwd) ? range(iter.fwdIter) : range(iter.revIter);
    std::vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(fwd){
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.fwdIter.index->sa[i].i1;
            int sa = iter.fwdIter.index->sa[i].i2;
            std::cout << i << ": " << sa + sequenceLengths[seq] << "\n";
        }
    }else{
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.revIter.index->sa[i].i1;
            int sa = iter.revIter.index->sa[i].i2;
            std::cout << i << ": " << sequenceLengths[seq + 1] - sa - 1 << "\n";
        }

    }
}

template <typename TText, typename TIndex, typename TIndexSpec>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              std::vector<std::pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors,
              bool const fwd)
{
    int size = seqan::length(iter.fwdIter.index->sa);
    Pair<uint32_t, uint32_t> dirrange = (fwd) ? range(iter.fwdIter) : range(iter.revIter);
    uint32_t number_of_indeces = size - bitvectors[0].first.size();
    std::vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(fwd){
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.fwdIter.index->sa[i].i1;
            int sa = iter.fwdIter.index->sa[i].i2;
            std::cout << i << ": " << sa + sequenceLengths[seq] << "\n";
        }
    }else{
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            int seq = iter.revIter.index->sa[i].i1;
            int sa = iter.revIter.index->sa[i].i2;
            std::cout << i << ": " << sequenceLengths[seq + 1] - sa - 1 << "\n";
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
    
    std::cout << "Default Hits:" << hits.size() << "\n";
    
    for(int i = 0; i < readOccs.size(); ++i){
        std::cout << "Errors: "<< (int)readOccs[i].errors;
        std::cout << "   " << readOccs[i].hit << "\n";
        if(occEnabled)
            std::cout << infix(genome[readOccs[i].hit.i2.i1], readOccs[i].hit.i2.i2, readOccs[i].hit.i2.i2 + seqan::length(readOccs[i].hit.i1)) << "\n";
        
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
    DnaString part = infix(genome[readOcc.hit.i2.i1], readOcc.hit.i2.i2, readOcc.hit.i2.i2 + seqan::length(readOcc.hit.i1));
    appendValue(testocc, part);
    std::cout << "Search occ: " << (int)readOcc.hit.i2.i2 << " which has seq: " << "\n";
    std::cout << part << "\n"; //TODO revert this

    find<minErrors, maxErrors>(delegate, index, testocc, HammingDistance());
//        print_readocc_sorted(hite, errors_v);
    std::cout << hits.size() << " hits!!!!!!!!!!" << "\n";
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
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                std::exit(1);
    }
    return(nhits);
}





template <typename TText, typename TIndexSpec>
std::vector<int> compare(Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                    int errors,
                    int threshold,
                    std::vector<readOcc> x,
                    std::vector<readOcc> y)
{
    std::vector<int> wrongHitCount; 
    bool same2 = true;
    bool same = true;
    if(!(x.size() == y.size())){
        same2 = false;
        std::cout << "MyVersion has: " << x.size() << "hits while default version has: " << y.size() << " hits" << "\n";
    }   
    int offset = 0;
    for(int i = 0; i + offset < y.size(); ++i){
        same = false;
        same = (i < x.size() && x[i].hit.i2.i1 == y[i + offset].hit.i2.i1 && x[i].hit.i2.i2 == y[i + offset].hit.i2.i2);
        while(!same && i + offset < y.size()){
            if(wrongHitCount.size() > 0)
                std::cout << "Something went wrong" << "\n";
            if(i < x.size())//TODO revert this
                std::cout << "MyVersion has: " << x[i].hit.i2 << " while " ; //TODO revert this
            std::cout << "default version has: " << y[i + offset].hit.i2 << "\n";//TODO revert this
            int nhits = testread(0, errors, index, y[i + offset]); //TODO  3 lines down
            if(nhits < threshold){      //TODO //3 lines down
                std::cout << "To few hits should have found this part!!!!" << "\n";
                wrongHitCount.push_back(nhits);
//                 --offset;
//                 break;
//                 std::exit(0);
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

void printPair(std::pair<uint32_t, uint32_t> p){
    std::cout << "<" << p.first << ", " << p.second << ">";
}
/*
void printbit(std::vector<std::pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, Pair<uint8_t, Pair<uint32_t, uint32_t>> brange){
    sdsl::bit_vector const & rb = bitvectors[brange.i1].first;
    std::cout << "bitvector: " << (int)brange.i1 << " brange start: " << brange.i2.i1 << "  brange end: " << brange.i2.i2 << "\n";
    for(int i = brange.i2.i1; i < brange.i2.i2; ++i)
        std::cout << i << " Bit: " << rb[i] << "\n";
}*/

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
        std::cout << static_cast<int> (a.at(i)) << ", ";
    }
    std::cout << "\n";
}

void printbit(std::vector<std::pair<sdsl::bit_vector, sdsl::rank_support_v<>>> & bitvectors, Pair<uint8_t, Pair<uint32_t, uint32_t>> brange){
    sdsl::bit_vector const & rb = bitvectors[brange.i1].first;
    std::cout << "bitvector: " << (int)brange.i1 << " brange start: " << brange.i2.i1 << "  brange end: " << brange.i2.i2 << "\n";
    for(int i = brange.i2.i1; i < brange.i2.i2; ++i)
        std::cout << i << " Bit: " << rb[i] << "\n";
}


inline bool file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

template <size_t nbrBlocks, size_t N>
void print_search_scheme(std::array<OptimalSearch<nbrBlocks>, N> & searchsscheme){
    for(int i = 0; i < searchsscheme.size(); i++){
        std::cout << "Search sscheme: " << i << "\n";
        std::cout << "Permutation: " << "\n";
        printv(searchsscheme[i].pi);
        std::cout << "Lower bound: " << "\n";
        printv(searchsscheme[i].l);
        std::cout << "Upper bound: " << "\n";
        printv(searchsscheme[i].u);
        std::cout << "blockLengths: " << "\n";
        printv(searchsscheme[i].blocklength);
        std::cout << "chronblockLengths: " << "\n";
        printv(searchsscheme[i].chronBL);
        std::cout << "revchronblockLengths: " << "\n";
        printv(searchsscheme[i].revChronBL);
        std::cout << "start Pos: " << "\n";
        std::cout << searchsscheme[i].startPos << "\n";
        std::cout << "minMax: " << "\n";
        printv(searchsscheme[i].min);
        printv(searchsscheme[i].max);
        std::cout << "OneDirection" << "\n" << (int)searchsscheme[i].startUniDir << "\n";
        std::cout << "\n";
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
                  std::string const & output_path, 
                  int chr)
{
    if(isBidirectionalIter<decltype(it)>::VALUE){
        std::cout << "Bidirectional iter" << "\n";
    }else{
        std::cout << "UniDirectional Iter" << "\n";
    }/*
    StringSet<DnaString> const & genome = indexText(*it.index);
    ofstream file(output_path, ios::out | ios::binary);
    for(int i = 0; i < chr; ++i){
        file << (">");
        file <<  std::to_string(i);
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
                  std::string const & output_path, 
                  int chr)
{
    if(isBidirectionalIter<decltype(it)>::VALUE){
        std::cout << "Bidirectional Iter" << "\n";
    }else{
        std::cout << "UniDirectional Iter" << "\n";
    }
    /*
    StringSet<DnaString> const & genome = indexText(*it.fwdIter.index);
    ofstream file(output_path, ios::out | ios::binary);
    for(int i = 0; i < chr; ++i){
        file << (">");
        file <<  std::to_string(i);
        file << ("\n");
        String<char> target;
        DnaString test = genome[i];
        move(target, test);
        file << target;
        file << ("\n");
    }
    file.close();*/
//}


std::vector<std::pair<sdsl::bit_vector, sdsl::rank_support_v<>>> loadBitvectors(CharString const bitvectorpath, const int K, const int errors){
    std::vector<std::pair<sdsl::bit_vector, sdsl::rank_support_v<>>> bit_vectors;
    if(file_exists(std::string("") + toCString(bitvectorpath) + "l_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors) + "_shift_" +  std::to_string(0)))
    {
    std::cout << "Load the following Bitvectors:" << "\n";
    for(int i = 0; i < 10; ++i){
        std::string file_name = std::string("") + toCString(bitvectorpath) + "r_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors) + "_shift_" +  std::to_string(i);
        if(file_exists(file_name)){
            sdsl::bit_vector b;
            std::cout << "Filename: " << file_name << "\n";
            load_from_file(b, file_name);
            sdsl::rank_support_v<> rb(& b);
            bit_vectors.push_back(std::make_pair(b, rb));
//              bit_vectors[i + 1].second.set_vector(&bit_vectors[i + 1].first);
         }
    }
    
    for(int i = 0; i < 10; ++i){
         std::string file_name = std::string("") + toCString(bitvectorpath) + "l_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors) + "_shift_" +  std::to_string(i);
         if(file_exists(file_name)){
             sdsl::bit_vector b;
             std::cout << "Filename: " << file_name << "\n";
             load_from_file(b, file_name);
             sdsl::rank_support_v<> rb(& b);
             bit_vectors.push_back(std::make_pair(b, rb));
//              bit_vectors[i + 1].second.set_vector(&bit_vectors[i + 1].first);
         }
    }
    }else{
    std::string file_name = std::string("") + toCString(bitvectorpath) + "right_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors);
    if(file_exists(file_name)){
        sdsl::bit_vector b;
        std::cout << "Filename: " << file_name << "\n";
        load_from_file(b, file_name);
        sdsl::rank_support_v<> rb(& b);
        bit_vectors.push_back(std::make_pair(b, rb));
        }else{
            std::cerr << file_name << " not found" <<  "\n";
        }

    file_name = std::string("") + toCString(bitvectorpath) + "middle_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors);
    if(file_exists(file_name)){
        sdsl::bit_vector b;
        std::cout << "Filename: " << file_name << "\n";
        load_from_file(b, file_name);
        sdsl::rank_support_v<> rb(& b);
        bit_vectors.push_back(std::make_pair(b, rb));
    }else{
        std::cerr << file_name << " not found" <<  "\n";
    }
        
    file_name = std::string("") + toCString(bitvectorpath) + "left_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors);
    if(file_exists(file_name)){
        sdsl::bit_vector b;
        std::cout << "Filename: " << file_name << "\n";
        load_from_file(b, file_name);
        sdsl::rank_support_v<> rb(& b);
        bit_vectors.push_back(std::make_pair(b, rb));
    }else{
        std::cerr << file_name << " not found" <<  "\n";
    }
        
    if(bit_vectors.size() != 3){
        std::cout << std::string("") + toCString(bitvectorpath) + "l_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors) + "_shift_" +  std::to_string(0) << "\n";
        std::cout << "was not found in the first place maybe wrong K or E parameter?" << "\n";
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
    std::cout << "This thing should not be executed getUniIndexGenome" << "\n";
    exit(0);
    auto const & rgenome = indexText(*it.revIter.index);
    return(rgenome);
}




template <typename TText, typename TIndex, typename TIndexSpec>
Pair<uint32_t, uint32_t> getUniRange(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it){
    std::cout << "This thing should not be executed getUniRange" << "\n";
    std::exit(0);
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
    std::cout << "This thing should not be executed getUniSa" << "\n";
    std::exit(0);
    auto const & sa = it.fwdIter.index->sa;
    return(sa);
}

// template <typename TText, typename TSpec, typename TConfig>
template <typename TText, typename TConfig, typename TIndexSpec>
String<unsigned> getUniSa(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > it){
    String<unsigned> const & sa = it.index->sa;
    return(sa);
}




#endif