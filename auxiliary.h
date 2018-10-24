#ifndef AUXILLARY_H_
#define AUXILLARY_H_

#include "common.h"
#include "common_auxiliary.h"

using namespace seqan;

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

template <typename TIndex>
auto getSeqLengths(TIndex & index){
//     auto mylimits = stringSetLimits(indexText(index)); //TODO test this
//     for(int i = 2; i < length(mylimits); ++i)
//         mylimits[i] -= mylimits[i - 1];
//     return(mylimits);
    auto const & genome = indexText(index);
    std::vector<uint32_t> sl;
    sl.push_back(0);
    for(uint32_t i = 0; i < countSequences(index)/*seqan::length(genome)*/; ++i)
        sl.push_back(seqan::length(genome[i]));
    return sl;
}


template <typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport>
std::vector<uint32_t> getSequencesLengths(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                std::vector<std::pair<TVector, TVSupport>> & bitvectors)
{
    uint32_t size = seqan::length(iter.fwdIter.index->sa);
    uint32_t number_of_indeces = size - bitvectors[0].first.size();

    std::cout << "Number of Indeces: " << number_of_indeces << "\n";
    std::vector<uint32_t> sequenceLengths(number_of_indeces + 1, 0);
    for(uint32_t i = 0; i < number_of_indeces; ++i){
        uint16_t seq = iter.fwdIter.index->sa[i].i1;
        uint32_t sa = iter.fwdIter.index->sa[i].i2;
        sequenceLengths[seq + 1] = sa;
        std::cout << "Saved length: " << sa << "\n";
        std::cout << "At position: " << seq + 1 << "\n"; //seq is 50%  448???????? debug
    }
    return sequenceLengths;
}

template<typename TIndex>
void calcfwdPos(TIndex & index,
                std::vector<hit> & hitsOutput,
                bool verbose = false)
{
    auto sl = getSeqLengths(index);

    for(uint32_t i = 0; i < hitsOutput.size(); ++i){
        if(hitsOutput[i].rev){
            if(verbose)
                std::cout << "before: " << hitsOutput[i].occ << "\n";
            hitsOutput[i].occ.i2 = sl[hitsOutput[i].occ.i1 + 1] - hitsOutput[i].occ.i2 - length(hitsOutput[i].read);
            if(verbose)
                std::cout << "after: " << hitsOutput[i].occ << "\n";
        }

    }
}

template <typename TText, typename TIndex, typename TIndexSpec>
void print_beginsa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              uint32_t number_of_indeces,
              CharString outputPath,
              bool const fwd)
{
    std::string name = (fwd) ? "start_fwd" : "start_rev";
    std::ofstream outfile(toCString(outputPath) + name, std::ios::out | std::ofstream::binary);
    uint32_t size = seqan::length(iter.fwdIter.index->sa);
    uint32_t noi = number_of_indeces;
    std::vector<uint32_t> sequenceLengths(number_of_indeces + 1, 0);
    for(uint32_t i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(uint32_t i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(!fwd){
        for(uint32_t i = 0; i < size && i < 1500; ++i){
            uint32_t seq = iter.revIter.index->sa[i].i1;
            uint32_t sa = iter.revIter.index->sa[i].i2;
            outfile << i << "\t(" << seq << ")\t" << sa << ":\t" << sa + sequenceLengths[seq] << "\n";
        }
    }else{
        for(uint32_t i = 0; i < size && i < 1500; ++i){
            uint32_t seq = iter.fwdIter.index->sa[i].i1;
            uint32_t sa = iter.fwdIter.index->sa[i].i2;
            outfile << i << "\t(" << seq << ")\t" << sa << ":\t" << sa + sequenceLengths[seq] << "\n";
        }
    }
    outfile.close();

}



template <typename TText, typename TIndex, typename TIndexSpec,
          typename TVector>
void print_beginsa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              std::vector<TVector> & bitvectors,
              CharString outputPath,
              bool const fwd)
{
    std::string name = (fwd) ? "start_fwd" : "start_rev";
    std::ofstream outfile(toCString(outputPath) + name, std::ios::out | std::ofstream::binary);
    uint32_t size = seqan::length(iter.fwdIter.index->sa);
    uint32_t number_of_indeces = size - bitvectors[0].size();
    uint32_t noi = number_of_indeces;
    std::vector<uint32_t> sequenceLengths(number_of_indeces + 1, 0);
    for(uint32_t i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(uint32_t i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(!fwd){
        for(uint32_t i = 0; i < size && i < 1500; ++i){
            uint32_t seq = iter.revIter.index->sa[i].i1;
            uint32_t sa = iter.revIter.index->sa[i].i2;
            outfile << i << "\t(" << seq << ")\t" << sa << ":\t" << sa + sequenceLengths[seq] << "\n";
        }
    }else{
        for(uint32_t i = 0; i < size && i < 1500; ++i){
            uint32_t seq = iter.fwdIter.index->sa[i].i1;
            uint32_t sa = iter.fwdIter.index->sa[i].i2;
            outfile << i << "\t(" << seq << ")\t" << sa << ":\t" << sa + sequenceLengths[seq] << "\n";
        }
    }
    outfile.close();

}


template <typename TText, typename TIndex, typename TIndexSpec>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              uint32_t const number_of_indeces,
              bool const fwd)
{
    Pair<uint32_t, uint32_t> dirrange = (fwd) ? range(iter.fwdIter) : range(iter.revIter);
    std::vector<int> sequenceLengths(number_of_indeces + 1, 0);
    for(uint32_t i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(uint32_t i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(fwd){
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            uint32_t seq = iter.fwdIter.index->sa[i].i1;
            uint32_t sa = iter.fwdIter.index->sa[i].i2;
            std::cout << i << ": " << sa + sequenceLengths[seq] << "\n";
        }
    }else{
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            uint32_t seq = iter.revIter.index->sa[i].i1;
            uint32_t sa = iter.revIter.index->sa[i].i2;
            std::cout << i << ": " << sequenceLengths[seq + 1] - sa - 1 << "\n";
        }

    }
}

template <typename TText, typename TIndex, typename TIndexSpec,
          typename TVector, typename TVSupport>
void print_sa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
              std::vector<std::pair<TVector, TVSupport>> & bitvectors,
              bool const fwd)
{
    uint32_t size = seqan::length(iter.fwdIter.index->sa);
    Pair<uint32_t, uint32_t> dirrange = (fwd) ? range(iter.fwdIter) : range(iter.revIter);
    uint32_t number_of_indeces = size - bitvectors[0].first.size();
    std::vector<uint32_t> sequenceLengths(number_of_indeces + 1, 0);
    for(uint32_t i = 0; i < number_of_indeces; ++i)
        sequenceLengths[iter.fwdIter.index->sa[i].i1 + 1] = iter.fwdIter.index->sa[i].i2;
        // cumulative sum seq
    for(uint32_t i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    if(fwd){
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            uint32_t seq = iter.fwdIter.index->sa[i].i1;
            uint32_t sa = iter.fwdIter.index->sa[i].i2;
            std::cout << i << ": " << sa + sequenceLengths[seq] << "\n";
        }
    }else{
        for(uint32_t i = dirrange.i1; i < dirrange.i2; ++i){
            uint32_t seq = iter.revIter.index->sa[i].i1;
            uint32_t sa = iter.revIter.index->sa[i].i2;
            std::cout << i << ": " << sequenceLengths[seq + 1] - sa - 1 << "\n";
        }

    }
}

bool occ_smaller(const hit & x, const hit & y)
{
    if(x.occ.i1 == y.occ.i1){
        if(x.occ.i2 == y.occ.i2)
            return x.errors < y.errors;
        else
            return x.occ.i2 < y.occ.i2;
    }
    else
    {
        return x.occ.i1 < y.occ.i1;
    }
}

bool occ_same(const hit & x, const hit & y)
{
    return(x.occ.i1 == y.occ.i1 && x.occ.i2 == y.occ.i2 && x.errors == y.errors);
}

bool occ_similar(const hit & x, const hit & y, uint32_t dis = 5)
{
    return(x.occ.i1 == y.occ.i1 && x.occ.i2 + dis >= y.occ.i2 && x.occ.i2 - dis <= y.occ.i2);
}

std::vector<hit> print_readocc_sorted(std::vector<hit> hits, auto const & genome, bool const editD, int nerrors, bool const occEnabled)
{
    if(!editD)
        nerrors = 0;
    std::sort(hits.begin(), hits.end(), occ_smaller);
    //TODO compare also needle to not remove identical reads
    if(editD){
        hits.erase(std::unique(hits.begin(), hits.end(), occ_same), hits.end());
        std::cout << "unique number: "  << hits.size() << "\n";
    }

    for(uint32_t i = 0; i < hits.size(); ++i){
        std::cout << "Errors: "<< (int)hits[i].errors;
        std::cout << "   " << hits[i].occ << " " << hits[i].read << "\n";
        if(occEnabled)
            std::cout << infix(genome[hits[i].occ.i1], hits[i].occ.i2 - nerrors, hits[i].occ.i2 + seqan::length(hits[i].read) + nerrors) << "\n";
    }
    return(hits);
}



template <size_t minErrors, size_t maxErrors,
          typename TText, typename TIndexSpec>
uint32_t testread(Index<TText, BidirectionalIndex<TIndexSpec> > & index,
              hit testhit, int const mErrors, bool const editD)
{
    auto const & genome = indexText(index);
    std::vector<hit> hits;
    auto delegate = [&hits](auto & iter, DnaString const & needle, uint8_t errors)
    {
        for (auto occ : getOccurrences(iter)){
            hit me;
            me.occ = occ;
            me.read = needle;
            me.errors = errors;
            me.rev = false;
            hits.push_back(me);
        }
    };

    StringSet<DnaString> testocc;

    DnaString part = infix(genome[testhit.occ.i1], testhit.occ.i2, testhit.occ.i2 + seqan::length(testhit.read));
    appendValue(testocc, part);

    bool common = testhit.occ.i2 > mErrors && testhit.occ.i2 + seqan::length(testhit.read) + mErrors < length(genome[testhit.occ.i1]);

    if(editD && !common){
        std::cout << "out of boundry skip!!" << "\n";
        return(666);
    }


    if(editD && common)
        part = infix(genome[testhit.occ.i1], testhit.occ.i2 - maxErrors, testhit.occ.i2 + seqan::length(testhit.read) + maxErrors);
    std::cout << "Search occ: " << (uint32_t)testhit.occ.i2 << " which has seq: " << "\n";
    std::cout << part << "\n"; //TODO revert this

    find<minErrors, maxErrors>(delegate, index, testocc, HammingDistance());
    std::cout << hits.size() << " hits!!!!!!!!!!" << "\n";
    return(hits.size());
}

template <typename TText, typename TIndexSpec>
uint32_t testread(int minErrors, int maxErrors, Index<TText, BidirectionalIndex<TIndexSpec> > & index,
              hit testhit, bool const editD){
    int nhits;
    switch (maxErrors)
    {
        case 1: nhits = testread<0, 1>(index, testhit, maxErrors, editD);
                break;
        case 2: nhits = testread<0, 2>(index, testhit, maxErrors, editD);
                break;
        case 3: nhits = testread<0, 3>(index, testhit, maxErrors, editD);
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                std::exit(1);
    }
    return(nhits);
}





template <typename TText, typename TIndexSpec>
std::vector<uint32_t> compare(Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                    int errors,
                    uint32_t threshold,
                    const bool editD,
                    std::vector<hit> x,
                    std::vector<hit> y)
{
    std::vector<uint32_t> wrongHitCount;
    bool same2 = true;
    bool same = true;
    if(!(x.size() == y.size())){
        same2 = false;
        std::cout << "MyVersion has: " << x.size() << "hits while default version has: " << y.size() << " hits" << "\n";
    }
    int offset = 0;
    for(uint32_t i = 0; i + offset < y.size(); ++i){
        same = i < x.size() && occ_same(x[i], y[i + offset])/*x[i].occ.i1 == y[i + offset].occ.i1 && x[i].occ.i2 == y[i + offset].occ.i2)*/;
        while(!same /*&& i + offset < y.size()*/){
            if(wrongHitCount.size() > 0)
                std::cout << "Something went wrong" << "\n";
            std::cout << "MyVersion has: " << x[i].occ.i2 << " while " ;
            std::cout << "default version has: " << y[i + offset].occ.i2 << "\n";
            uint32_t nhits = testread(0, errors, index, y[i + offset], editD);
            if(nhits < threshold){
                std::cout << "To few hits should have found this part!!!!" << "\n";
                wrongHitCount.push_back(nhits);
//                 --offset;
//                 break;
//                 std::exit(0);
            }
            ++offset;
            same = (i < x.size() && x[i].occ.i1 == y[i + offset].occ.i1 && x[i].occ.i2 == y[i + offset].occ.i2);
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

template<typename TVector, typename TVSupport>
void printbit(std::vector<std::pair<TVector, TVSupport>> & bitvectors, Pair<uint8_t, Pair<uint32_t, uint32_t>> brange){
    TVector const & rb = bitvectors[brange.i1].first;
    std::cout << "bitvector: " << (int)brange.i1 << " brange start: " << brange.i2.i1 << "  brange end: " << brange.i2.i2 << "\n";
    for(uint32_t i = brange.i2.i1; i < brange.i2.i2; ++i)
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
    auto const & genome = indexText(*it.index);
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
    auto const & genome = indexText(*it.fwdIter.index);
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


std::vector<std::pair<TBitvector, TSupport>> loadBitvectors(CharString const bitvectorpath, const int K, const int errors){
    std::vector<std::pair<TBitvector, TSupport>> bit_vectors;
    if(file_exists(std::string("") + toCString(bitvectorpath) + "l_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors) + "_shift_" +  std::to_string(0)))
    {
    std::cout << "Load the following Bitvectors:" << "\n";
    for(int i = 0; i < K + 2; ++i){
        std::string file_name = std::string("") + toCString(bitvectorpath) + "r_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors) + "_shift_" +  std::to_string(i);
        if(file_exists(file_name)){
            TBitvector b;
            std::cout << "Filename: " << file_name << "\n";
            load_from_file(b, file_name);
            TSupport rb(& b);
            bit_vectors.push_back(std::make_pair(b, rb));
         }
    }

    for(int i = 0; i < K + 2; ++i){
         std::string file_name = std::string("") + toCString(bitvectorpath) + "l_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors) + "_shift_" +  std::to_string(i);
         if(file_exists(file_name)){
             TBitvector b;
             std::cout << "Filename: " << file_name << "\n";
             load_from_file(b, file_name);
             TSupport rb(& b);
             bit_vectors.push_back(std::make_pair(b, rb));
         }
    }
    }else{
    std::cout << "Searching for: " << std::string("") + toCString(bitvectorpath) + "l_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors) + "_shift_" +  std::to_string(0) << "\n";
    std::string file_name = std::string("") + toCString(bitvectorpath) + "right_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors);
    if(file_exists(file_name)){
        TBitvector b;
        std::cout << "Filename: " << file_name << "\n";
        load_from_file(b, file_name);
        TSupport rb(& b);
        bit_vectors.push_back(std::make_pair(b, rb));
        }else{
            std::cerr << file_name << " not found" <<  "\n";
        }

    file_name = std::string("") + toCString(bitvectorpath) + "middle_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors);
    if(file_exists(file_name)){
        TBitvector b;
        std::cout << "Filename: " << file_name << "\n";
        load_from_file(b, file_name);
        TSupport rb(& b);
        bit_vectors.push_back(std::make_pair(b, rb));
    }else{
        std::cerr << file_name << " not found" <<  "\n";
    }

    file_name = std::string("") + toCString(bitvectorpath) + "left_bit_vector_" +  std::to_string(K) + "_" +  std::to_string(errors);
    if(file_exists(file_name)){
        TBitvector b;
        std::cout << "Filename: " << file_name << "\n";
        load_from_file(b, file_name);
        TSupport rb(& b);
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
auto getUniSa(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it){
    std::cout << "This thing should not be executed getUniSa" << "\n";
    std::exit(0);
    auto const & sa = it.fwdIter.index->sa;
    return(sa);
}

// template <typename TText, typename TSpec, typename TConfig>
template <typename TText, typename TConfig, typename TIndexSpec>
auto getUniSa(Iter<Index<TText, FMIndex<void, TConfig> >, VSTree<TopDown<TIndexSpec> > > it){
    auto const & sa = it.index->sa;
    return(sa);
}




#endif
