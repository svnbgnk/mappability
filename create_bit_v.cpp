#include "common.h"
#include <iostream>
#include <fstream>
#include <iterator> 
#include <sstream>
#include <tgmath.h>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace seqan;

std::vector<int> getInt(std::string const& mappability_str)
{
  std::istringstream iss(mappability_str);
  return std::vector<int>{ 
    std::istream_iterator<int>(iss),
    std::istream_iterator<int>()
  };
}

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

template <typename TChar, typename TAllocConfig>
void loadIndex(vector<sdsl::bit_vector> &bit_vectors, CharString const indexPath)
{    
    typedef String<TChar, TAllocConfig> TString;
    Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
    open(index, toCString(indexPath), OPEN_RDONLY);
    
    cout << "Index size:" << seqan::length(index.fwd.sa) << endl;
    vector<sdsl::bit_vector> bit_vectors_ordered (bit_vectors);    
    int number_of_indeces = seqan::length(index.fwd.sa) - bit_vectors[0].size();
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    cout << "Number of Indeces: " << number_of_indeces << endl;
    
    //sequenceLengths first value is 0
    for(int i = 0; i < number_of_indeces; ++i)
        sequenceLengths[getValueI1(index.fwd.sa[i]) + 1] = getValueI2(index.fwd.sa[i]);
    // cumulative sum seq
    for(int i = 1; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    // skip sentinels

    for (unsigned j = 0; j < seqan::length(index.fwd.sa) - number_of_indeces; ++j)
    {
        uint32_t sa_j = getValueI2(index.fwd.sa[j + number_of_indeces]);
        uint16_t seq = getValueI1(index.fwd.sa[j + number_of_indeces]);
        for(int i = 0; i < bit_vectors.size(); ++i){
            bit_vectors_ordered[i][sa_j + sequenceLengths[seq]] = bit_vectors[i][j];
        }
    }
    bit_vectors = bit_vectors_ordered;
}



template <typename TChar>
void loadIndex(vector<sdsl::bit_vector> &bit_vectors, CharString const indexPath, bool const mmap)
{
    if(mmap)
        loadIndex<TChar, MMap<> >(bit_vectors, indexPath);
    else
        loadIndex<TChar, Alloc<> >(bit_vectors, indexPath);
}

void order_bit_vector(vector<sdsl::bit_vector> &bit_vectors, CharString const indexPath, bool const mmap, CharString const alphabet){
     if(alphabet == "dna4")
         loadIndex<Dna>(bit_vectors, indexPath, mmap);
     else
         loadIndex<Dna5>(bit_vectors, indexPath, mmap);
}


/*
bool sdsl::operator[](int idx){
    
}
*/
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
    
    addOption(parser, ArgParseOption("T", "threshold", "Threshold for inverse frequency that gets accepted", ArgParseArgument::DOUBLE, "DOUBLE"));
    setRequired(parser, "threshold");
    
    addOption(parser, ArgParseOption("e", "errors", "Max errors allowed during mapping", ArgParseArgument::INTEGER, "INT"));
    
    addOption(parser, ArgParseOption("d", "debug", "Also create chronical bit_vectors (for debugging)"));
    
    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));
    
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    
    //Retrieve input parameters
    CharString indexPath, _indexPath, outputPath;
    string mappability_path;
    int len, errors = 0, readlength = 100, allowederrors = 2; 
    double threshold;
    
    getOptionValue(mappability_path, parser, "map");
    getOptionValue(indexPath, parser, "index");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(len, parser, "length");
    getOptionValue(threshold, parser, "threshold");
    getOptionValue(errors, parser, "errors");
    bool debug = isSet(parser, "debug");
    bool mmap = isSet(parser, "mmap");
    
    StringSet<CharString> ids;
    CharString alphabet;
    string mappability_str;
//     int shift = floor(static_cast<float>(readlength)/(allowederrors + 2));
    //ceiling for len
    
    _indexPath = indexPath;
    _indexPath += ".alphabet";
    open(alphabet, toCString(_indexPath));

    if(!file_exists(mappability_path))
    {
        cout << "Cannot find mappability file" << endl;
        exit(0);
    }
    ifstream myfile;
    myfile.open(mappability_path, std::ios::in | std::ofstream::binary);
    std::getline(myfile, mappability_str);
    vector<int> mappability_int = getInt(mappability_str);
    vector<float> mappability(mappability_int.begin(), mappability_int.end());
    myfile.close();
    for(vector<float>::iterator it = mappability.begin(); mappability.end() != it; ++it)
        *it = static_cast<float>(1) / *it;
    cout << "mappability size:" << mappability.size() << endl;    
    // TODO Merge both bit_vectors into for creation
    sdsl::bit_vector righti (mappability.size() + len - 1, 1);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 1);
        
    for(unsigned i = 0; i < mappability.size(); ++i){
        lefti[i + len - 1] = (mappability[i] >= threshold);
        righti[i] = (mappability[i] >= threshold);
    }
    cout << "Bit Vector Length: " << righti.size() << endl;
    vector<sdsl::bit_vector> bit_vectors;
    vector<string> names;
    bit_vectors.push_back(lefti);
    bit_vectors.push_back(righti);
    names.push_back("l_bit_vector_" + to_string(len));
    names.push_back("r_bit_vector_" + to_string(len));
    if(errors != 0){
        for(int i = 1; i < (errors + 2); ++i){
            sdsl::bit_vector newleft(mappability.size() + len - 1, 1);
            sdsl::bit_vector newright(mappability.size() + len - 1, 1);
            int shift = i * std::ceil(len / (errors + 2));
            for(int j = 0; j < righti.size(); ++j){
                if(j - shift >= 0)
                    newright[j] = righti[j - shift];
                if(j + shift < lefti.size() - 1)
                    newleft[j] = lefti[j + shift];
            }
            bit_vectors.push_back(newleft);
            bit_vectors.push_back(newright);
            names.push_back("l_bit_vector_" + to_string(len) + "_shift_" + to_string(shift));
            names.push_back("r_bit_vector_" + to_string(len) + "_shift_" + to_string(shift));
        }
    }
    if(debug)
    {
        sdsl::store_to_file(bit_vectors[0], toCString(outputPath) + names[0] + "_for_heatmap");
        for(int i = 0; i < bit_vectors.size(); ++i){
            std::ofstream outfile((toCString(outputPath) + names[i] + "_debug"), std::ios::out | std::ofstream::binary);
            std::copy(bit_vectors[i].begin(), bit_vectors[i].end(), std::ostream_iterator<bool>(outfile));
            outfile.close();
            
        }
    }
    //order in suffix array
    order_bit_vector(bit_vectors, indexPath, mmap, alphabet);
    for(int i = 0; i < bit_vectors.size(); ++i){
        sdsl::store_to_file(bit_vectors[i], toCString(outputPath) + names[i]);
//         std::ofstream outfile((toCString(outputPath) + names[i]), std::ios::out | std::ofstream::binary);
//         std::copy(bit_vectors[i].begin(), bit_vectors[i].end(), std::ostream_iterator<bool>(outfile));
//         outfile.close();
    }    
    return 0;
}
