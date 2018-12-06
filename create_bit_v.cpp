#include <sdsl/bit_vectors.hpp>
#include <seqan/arg_parse.h>
#include "auxiliary.h"
#include "common.h"

using namespace std;
using namespace seqan;


struct bitvectors
{
    vector<bool> fwdd;
    vector<string> names;
    vector<sdsl::bit_vector> bv;
};

std::vector<int> getInt(std::string const& mappability_str)
{
  std::istringstream iss(mappability_str);
  return std::vector<int>{
    std::istream_iterator<int>(iss),
    std::istream_iterator<int>()
  };
}

vector<uint8_t> read(const string mappability_path){
    string mappability_str;

    vector<uint8_t> mappability_int;
    ifstream file(toCString(mappability_path), std::ios::binary);
    if (!file.eof() && !file.fail())
    {
        file.seekg(0, std::ios_base::end);
        std::streampos fileSize = file.tellg();
        mappability_int.resize(fileSize);
        file.seekg(0, std::ios_base::beg);
        file.read((char*)&mappability_int[0], fileSize);
        file.close();
        cout << "Load successful" << endl;
        }
    return(mappability_int);
}

auto loadBlockLengths(uint8_t se, uint32_t const len)
{
    vector<int> r;
    vector<int> l;
    switch (se)
    {
        case 0:
        {
            auto scheme = OptimalSearchSchemes<0, 0>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(int i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        case 1:
        {
            auto scheme = OptimalSearchSchemes<0, 1>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(int i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        case 2:
        {
            auto scheme = OptimalSearchSchemes<0, 2>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(int i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        case 3:
        {
            auto scheme = OptimalSearchSchemes<0, 3>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(int i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        case 4:
        {
            auto scheme = OptimalSearchSchemes<0, 4>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(int i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        default: std::cerr << "E = " << (int)se << " not yet supported.\n";
                exit(1);
    }
    return std::make_pair(r, l);
}

template<typename TVector, typename TElem>
bool checkForElem(TVector const & v, TElem const & e)
{
   for(int i = 0; i < v.size(); ++i){
       if(v[i].compare(e) == 0)
           return true;
   }
   return false;
}

bitvectors create_all_bit_vectors(const vector <uint8_t> & mappability,
                                  uint32_t const len, uint32_t const threshold, uint8_t const errors, uint8_t const strata){
    //TODO switch left and right in the moment they discribe in which direction the k-mere is
    bitvectors b;
    int e = errors;
    sdsl::bit_vector righti (mappability.size() + len - 1, 0);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 0);
    #pragma omp parallel for schedule(static)
    for(uint32_t i = 0; i < mappability.size(); ++i){
        lefti[i + len - 1] = (mappability[i] <= threshold);
        righti[i] = (mappability[i] <= threshold);
    }
    cout << "Finished Default Bit Vectors.  Length: " << righti.size() << endl;


    b.bv.push_back(righti);
    b.names.push_back("r_bit_vector_" + to_string(len) + "_shift_0");
    b.fwdd.push_back(true);
    b.bv.push_back(lefti);
    b.names.push_back("l_bit_vector_" + to_string(len) + "_shift_0");
    b.fwdd.push_back(false);

    std::cout << "\nAdditonal bitvectors besides shift_0\n";

    vector<int> shift_r;
    vector<int> shift_l;

    for(int s = 0; s <= e - strata; ++s){
        int se = e - s;
        std::cout << "Bitvectors for Scheme <" << se << ", " << se << ">: \n";
        auto blocklengths = loadBlockLengths(se, len);
        shift_r = blocklengths.first;
        shift_l = blocklengths.second;

        for(int i = 0; i < shift_r.size() - 1; ++i){
            uint32_t shift = shift_r[i];
            bool skip = checkForElem(b.names, ("left_anchored_bvector_" + to_string(len) + "_shift_" + to_string(shift)));
            if(skip){
                std::cout << "Bitvector already included" << "\n";
                continue;
            }

            b.names.push_back("left_anchored_bvector_" + to_string(len) + "_shift_" + to_string(shift));
            std::cout << b.names.back()<< /*"\t shift : " << shift <<*/ "\n";

            sdsl::bit_vector newright(mappability.size() + len - 1, 0);
            #pragma omp parallel for schedule(static)
            for(uint32_t j = 0; j < righti.size(); ++j){
                if(j >= shift)
                    newright[j] = righti[j - shift];
            }
            b.bv.push_back(newright);
            b.fwdd.push_back(true);
        }

        for(int i = 1; i < shift_l.size(); ++i){
            uint32_t shift = shift_l[i];
            bool skip = checkForElem(b.names, ("right_anchored_bvector_" + to_string(len) + "_shift_" + to_string(shift)));
            if(skip){
                std::cout << "Bitvector already included" << "\n";
                continue;
            }

            b.names.push_back("right_anchored_bvector_" + to_string(len) + "_shift_" + to_string(shift));
            std::cout << b.names.back() << /*"\t shift : " << shift <<*/ "\n";

            sdsl::bit_vector newleft(mappability.size() + len - 1, 0);
            #pragma omp parallel for schedule(static)
            for(uint32_t j = 0; j < righti.size(); ++j){
                if(j + shift < lefti.size() - 1)
                    newleft[j] = lefti[j + shift];
            }
            b.bv.push_back(newleft);
            b.fwdd.push_back(false);
        }
    }

    std::cout << "Number of Bitvectors: " << b.names.size() << "\n\n";

/*
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    auto s = scheme[0];

    uint8_t blocks = s.pi.size();
    if(errors != 0){
        for(uint32_t i = 0; i < blocks - 1; ++i){
            sdsl::bit_vector newright(mappability.size() + len - 1, 0);
            uint32_t shift = s.chronBL[i];
            cout << "r bitvector  name: " << to_string(shift) << endl;
            cout << "with shift: " << shift << endl;

            #pragma omp parallel for schedule(static)
            for(uint32_t j = 0; j < righti.size(); ++j){
                if(j >= shift)
                    newright[j] = righti[j - shift];
            }
            b.bv.push_back(newright);
            b.names.push_back("r_bit_vector_" + to_string(len) + "_shift_" + to_string(shift));
            b.fwdd.push_back(true);
        }

        for(uint32_t i = 1; i < blocks; ++i){
            sdsl::bit_vector newleft(mappability.size() + len - 1, 0);
            uint32_t shift = s.revChronBL[blocks - i];
            cout << "l bitvector  name: " << to_string(shift) << endl;
            cout << "with shift: " << shift << endl;
            #pragma omp parallel for schedule(static)
            for(uint32_t j = 0; j < righti.size(); ++j){
                if(j + shift < lefti.size() - 1)
                    newleft[j] = lefti[j + shift];
            }
            b.bv.push_back(newleft);
            b.names.push_back("l_bit_vector_" + to_string(len) + "_shift_" + to_string(shift));
            b.fwdd.push_back(false);
        }
    }*/


    return(b);
}


/*
bitvectors create_bit_vectors(const vector <uint8_t> & mappability, uint32_t const len, uint32_t const threshold, bool const bit3, uint8_t const errors){
    bitvectors result;
    if(bit3){
        switch (errors)
        {
            case 0: result = create_bit_vectors<0>(mappability, len, threshold);
                    break;
            case 1: result = create_bit_vectors<1>(mappability, len, threshold);
                    break;
            case 2: result = create_bit_vectors<2>(mappability, len, threshold);
                    break;
            case 3: result = create_bit_vectors<3>(mappability, len, threshold);
                    break;
            default: cerr << "E = " << errors << " not yet supported.\n";
                    exit(1);
        }
    }else{
        switch (errors)
        {
            case 0: result = create_all_bit_vectors<0>(mappability, len, threshold);
                    break;
            case 1: result = create_all_bit_vectors<1>(mappability, len, threshold);
                    break;
            case 2: result = create_all_bit_vectors<2>(mappability, len, threshold);
                    break;
            case 3: result = create_all_bit_vectors<3>(mappability, len, threshold);
                    break;
            default: cerr << "E = " << errors << " not yet supported.\n";
                    exit(1);
        }
    }
    return(result);
}*/

template <typename TChar, typename TAllocConfig>
void order_bit_vector(bitvectors & b, CharString const indexPath, uint32_t const threads)
{
    cout << mytime() << "Loading Index" << endl;
    typedef String<TChar, TAllocConfig> TString;
    Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
    open(index, toCString(indexPath), OPEN_RDONLY);
    uint32_t indexSize = seqan::length(index.fwd.sa);
    cout << mytime() << "Loaded Index. Size:" << indexSize << endl;
    vector<sdsl::bit_vector> bit_vectors_ordered (b.bv);

    uint32_t number_of_indeces = countSequences(index);

    std::vector<uint32_t> sequenceLengths = getSeqLengths(index);
    for(uint32_t i = 2; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += sequenceLengths[i - 1];

    cout << "Number of Sequences in index: " << countSequences(index) << endl;
    cout << mytime() << "Start sorting bitvectors" << endl;

    uint32_t mythreads;
    if(threads == 0)
        mythreads = omp_get_max_threads();
    else
        mythreads = threads;
    //dynamic since tasks can take different amount of time (SA Sampling) critical to guarantee thread safety

    uint8_t bsize = b.bv.size();

//     #pragma omp parallel for schedule(dynamic) num_threads(mythreads)
    #pragma omp parallel for schedule(static, (indexSize/(mythreads*100))) num_threads(mythreads)
    for (unsigned j = 0; j < indexSize - number_of_indeces; ++j)
    {
        // skip sentinels
        Pair<uint16_t, uint32_t> sa_f = index.fwd.sa[j + number_of_indeces];
        Pair<uint16_t, uint32_t> sa_r = index.rev.sa[j + number_of_indeces];
        uint32_t fpos = sa_f.i2 + sequenceLengths[sa_f.i1];
        uint32_t rpos = sequenceLengths[sa_r.i1 + 1] - sa_r.i2 - 1;
        vector<bool> values(bsize);

        for(uint32_t i = 0; i < bsize; ++i){
            if(b.fwdd[i]){
                values[i] = b.bv[i][fpos];
            }
            else
            {
                values[i] = b.bv[i][rpos];
            }
        }
        #pragma omp critical
        {
        for(uint32_t i = 0; i < bsize; ++i)
            bit_vectors_ordered[i][j] = values[i];
        }
    }
    b.bv = bit_vectors_ordered;
}



template <typename TChar>
void order_bit_vector(bitvectors & bit_vectors, CharString const indexPath, bool const mmap, uint32_t const threads)
{
    if(mmap)
        order_bit_vector<TChar, MMap<> >(bit_vectors, indexPath,threads);
    else
        order_bit_vector<TChar, Alloc<> >(bit_vectors, indexPath, threads);
}

void order_bit_vector(bitvectors & bit_vectors, CharString const indexPath, bool const mmap, CharString const alphabet, uint32_t const threads){
     if(alphabet == "dna4")
         order_bit_vector<Dna>(bit_vectors, indexPath, mmap, threads);
     else
         order_bit_vector<Dna5>(bit_vectors, indexPath, mmap, threads);
}


int main(int argc, char *argv[])
{

    ArgumentParser parser("Create bit vectors");
    addOption(parser, ArgParseOption("I", "map", "Path to the mappability file (including the number at the end)", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "map");

    addOption(parser, ArgParseOption("Idx", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("K", "length", "Length of reads that will be searched (for Hamming Distance read length == k-mere length)", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("T", "threshold", "Number of times a k-mer can occure and still be accepted as mappable", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "threshold");

    addOption(parser, ArgParseOption("E", "errors", "Max errors allowed during mapping", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "errors");

    addOption(parser, ArgParseOption("s", "strata", "Max errors allowed during mapping", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("d", "debug", "Also create chronical bit_vectors (for debugging)"));

    addOption(parser, ArgParseOption("S", "startSa", "Create first 1500 lines from SA array fwd and rev and then quit"));

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads used", ArgParseArgument::INTEGER, "INT"));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    //Retrieve input parameters
    CharString indexPath, _indexPath, outputPath;
    string mappability_path;
    unsigned len, threshold, errors, threads = 0;

    getOptionValue(mappability_path, parser, "map");
    getOptionValue(indexPath, parser, "index");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(len, parser, "length");
    getOptionValue(threshold, parser, "threshold");
    getOptionValue(errors, parser, "errors");
    unsigned strata = errors;
    getOptionValue(strata, parser, "strata");
    bool debug = isSet(parser, "debug");
    bool startSa = isSet(parser, "startSa");
    bool mmap = isSet(parser, "mmap");
    getOptionValue(threads, parser, "threads");

    StringSet<CharString> ids;
    CharString alphabet;

    _indexPath = indexPath;
    _indexPath += ".alphabet";
    open(alphabet, toCString(_indexPath));


    if(startSa){
        typedef String<Dna, Alloc<>> TString;
        Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
        open(index, toCString(indexPath), OPEN_RDONLY);
        Iter<Index<TText, TIndexConfig>, VSTree<TopDown<> > > it(index);
        print_beginsa(it, 500, outputPath, true);
        print_beginsa(it, 500, outputPath, false);
        exit(0);
    }

    if(!file_exists(mappability_path))
    {
        cout << "Cannot find mappability file" << endl;
        exit(0);
    }

    cout << mytime() << "Program start." << endl;
    cout << "Errors:" << errors << "\tStrata: " << strata << endl; //TODO comment
    vector<uint8_t> mappability = read(mappability_path);
    cout << mytime() << "Loaded Mappability vector. Size: " << mappability.size() << endl;

    bitvectors result = create_all_bit_vectors(mappability, len, threshold, errors, strata);

    cout << mytime() << "Finished bit vectors." << endl;

    if(debug)
    {
        for(uint32_t i = 0; i < result.bv.size(); ++i){
            sdsl::store_to_file(result.bv[i], toCString(outputPath) + result.names[i] + "_chron");
            std::ofstream outfile((toCString(outputPath) + result.names[i] + "_debug"), std::ios::out | std::ofstream::binary);
            std::copy(result.bv[i].begin(), result.bv[i].end(), std::ostream_iterator<bool>(outfile));
            outfile.close();
        }
    }

    //order bitvectors according to the forward or reverse suffix array
    order_bit_vector(result, indexPath, mmap, alphabet, threads);
    cout << mytime() << "Finished sorting" << endl;
    for(uint32_t i = 0; i < result.bv.size(); ++i){
        sdsl::store_to_file(result.bv[i], toCString(outputPath) + result.names[i]);
        if(debug){
            std::ofstream outfile((toCString(outputPath) + result.names[i] + "_osa_debug"), std::ios::out | std::ofstream::binary);
            std::copy(result.bv[i].begin(), result.bv[i].end(), std::ostream_iterator<bool>(outfile));
            outfile.close();
        }
    }

    cout << mytime() << "Finished saving bit vectors" << endl;


    return 0;
}
