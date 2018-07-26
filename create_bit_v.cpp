#include <sdsl/bit_vectors.hpp>
#include <seqan/arg_parse.h>
#include "auxiliary.h"
#include "common_auxiliary.h"

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

template <unsigned errors>
bitvectors create_all_bit_vectors(const vector <uint8_t> & mappability, uint32_t const len, uint32_t const threshold){
    bitvectors b;
    int e = errors;
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    auto s = scheme[0];
    sdsl::bit_vector righti (mappability.size() + len - 1, 0);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 0);
    #pragma omp parallel for schedule(static)
    for(uint32_t i = 0; i < mappability.size(); ++i){
        lefti[i + len - 1] = (mappability[i] <= threshold);
        righti[i] = (mappability[i] <= threshold);
    }
    cout << "Finished Default Bit Vectors.  Length: " << righti.size() << endl;


    b.bv.push_back(righti);
    b.names.push_back("r_bit_vector_" + to_string(len) + "_" + to_string(e) + "_shift_0");
    b.fwdd.push_back(true);

    uint8_t blocks = s.pi.size();
    if(errors != 0){
        for(uint32_t i = 0; i < blocks - 1; ++i){
            sdsl::bit_vector newright(mappability.size() + len - 1, 0); //TODO think 0 or 1 in edge cases
            uint32_t shift = s.chronBL[i];
            cout << "r bitvector  name: " << to_string(i + 1) << endl;
            cout << "with shift: " << shift << endl;

            #pragma omp parallel for schedule(static)
            for(uint32_t j = 0; j < righti.size(); ++j){
                if(j >= shift)
                    newright[j] = righti[j - shift];
            }
            b.bv.push_back(newright);
            b.names.push_back("r_bit_vector_" + to_string(len) + "_" + to_string(e) + "_shift_" + to_string(i + 1));
            b.fwdd.push_back(true);
        }

        for(uint32_t i = 1; i < blocks; ++i){
            sdsl::bit_vector newleft(mappability.size() + len - 1, 0);//TODO think 0 or 1 in edge cases
            uint32_t shift = s.revChronBL[blocks - i];
            cout << "l bitvector  name: " << to_string(i) << endl;
            cout << "with shift: " << shift << endl;
            #pragma omp parallel for schedule(static)
            for(uint32_t j = 0; j < righti.size(); ++j){
                if(j + shift < lefti.size() - 1)
                    newleft[j] = lefti[j + shift];
            }
            b.bv.push_back(newleft);
            b.names.push_back("l_bit_vector_" + to_string(len) + "_" + to_string(e) + "_shift_" + to_string(i));
            b.fwdd.push_back(false);
        }
    }

    b.bv.push_back(lefti);
    b.names.push_back("l_bit_vector_" + to_string(len) + "_" + to_string(e) + "_shift_0");
    b.fwdd.push_back(false);
    return(b);
}


template <unsigned errors>
bitvectors create_bit_vectors(const vector <uint8_t> & mappability, uint32_t const len, uint32_t const threshold){

    int e = errors;
    cout << "Create minimum amount of bitvectors" << endl;
    bitvectors b;
    sdsl::bit_vector righti (mappability.size() + len - 1, 0);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 0);
    #pragma omp parallel for schedule(static)
    for(uint32_t i = 0; i < mappability.size(); ++i){
        lefti[i + len - 1] = (mappability[i] <= threshold);
        righti[i] = (mappability[i] <= threshold);
    }
    cout << "Finished Default Bit Vectors.  Length: " << righti.size() << endl;
    vector<sdsl::bit_vector> bit_vectors;
    vector<string> names;


    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    for (auto & s : scheme){
        bool fwd = (s.pi[0] < s.pi[1]);
        uint32_t pos = s.pi[0];
        cout << pos << endl;
        cout << "Direction forward " << fwd << endl;
        uint8_t blocks = s.pi.size();
        if(pos == 1)
            {
                b.bv.push_back(righti);
                b.names.push_back("right_bit_vector_" + to_string(len) + "_" + to_string(e));
                b.fwdd.push_back(true);
                cout << "case1" << endl;
            }
        else if(pos == blocks)
            {
                b.bv.push_back(lefti);
                b.names.push_back("left_bit_vector_" + to_string(len) + "_" + to_string(e));
                b.fwdd.push_back(false);
                cout << "case2" << endl;
            }
        else
            {
                if(fwd){
                    sdsl::bit_vector newright(mappability.size() + len - 1, 0); //TODO think 0 or 1 in edge cases
                    uint32_t shift = s.chronBL[pos - 2];
                    cout << "shift r_bit for: " << shift << endl;
                    cout << "pos:  " << pos << endl;
                    for(uint32_t j = 0; j < righti.size(); ++j){
                        if(j >= shift)
                            newright[j] = righti[j - shift];
                    }
                    b.bv.push_back(newright);
                    b.names.push_back("middle_bit_vector_" + to_string(len) + "_" + to_string(e));
                    b.fwdd.push_back(true);
                }else{
                    //NOTE int32_t is not large enought to handle all the positions in the hg
                    sdsl::bit_vector newleft(mappability.size() + len - 1, 0);//TODO think 0 or 1 in edge cases
                    uint32_t shift = s.revChronBL[pos];
                    cout << "shift l_bit for: " << shift << endl;
                    cout << "pos:  " << pos << endl;
                    for(uint32_t j = 0; j < righti.size(); ++j){
                        if(j + shift < lefti.size() - 1)
                            newleft[j] = lefti[j + shift];
                    }
                    b.bv.push_back(newleft);
                    b.names.push_back("middle_bit_vector_" + to_string(len) + "_" + to_string(e));
                    b.fwdd.push_back(false);
                }
            }
    }
    return(b);
}


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
}

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
    //dynamic since tasks can take different amount (SA Sampling) ordered to guarantee thread safety

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

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers in the mappability vector", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("T", "threshold", "Number of times a k-mer can occure and still be accepted as mappable", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "threshold");

    addOption(parser, ArgParseOption("E", "errors", "Max errors allowed during mapping", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "errors");
    addOption(parser, ArgParseOption("d", "debug", "Also create chronical bit_vectors (for debugging)"));

    addOption(parser, ArgParseOption("s", "startSa", "Create first 1500 lines from SA array fwd and rev and then quit"));

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("min", "3bitversion",
        "Only create the required 3 bitvectors needed for acquiring mappability in non-unidirectional cases"));

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
    bool debug = isSet(parser, "debug");
    bool startSa = isSet(parser, "startSa");
    bool mmap = isSet(parser, "mmap");
    bool bit3 = isSet(parser, "3bitversion");
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
    vector<uint8_t> mappability = read(mappability_path);
    cout << mytime() << "Loaded Mappability vector. Size: " << mappability.size() << endl;


    bitvectors result = create_bit_vectors(mappability, len, threshold, bit3, errors);

    cout << mytime() << "Finished bit vectors." << endl;

    if(debug)
    {
        for(uint32_t i = 0; i < result.bv.size(); ++i){
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
