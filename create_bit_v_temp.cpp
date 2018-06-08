#include "common.h"
#include <iostream>
#include <fstream>
#include <iterator> 
#include <sstream>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace seqan;

std::vector<double> getDouble(std::string const& mappability_str)
{
  std::istringstream iss(mappability_str);

  return std::vector<double>{ 
    std::istream_iterator<double>(iss),
    std::istream_iterator<double>()
  };
}


Pair<sdsl::bit_vector, sdsl::bit_vector> create_bit_vector_pair(sdsl::rank_support_v<> const &rank_left, sdsl::rank_support_v<> const &rank_right, int const len, int const blockL)
{
    sdsl::bit_vector right_len2 (rank_right.size(), 1);
    sdsl::bit_vector left_len2 (rank_left.size(), 1);
    
    // cant use unsigned because expect negative values (i - len)
    for(int i = 0; i < rank_left.size(); i++)
    {
        if(i - blockL > 0){
            left_len2[i] = rank_left(i + 1) - rank_left(i - blockL) == blockL + 1;
        }
        else
        {
//             cout << "rank " << rank_left(i + 1) << " " << i + 1 << " B:" << blockL << endl;
            left_len2[i] = rank_left(i + 1) == i + 1;
        }
        if(i + blockL + 1 <  rank_left.size()){
            right_len2[i] = rank_right(i + blockL + 1) - rank_right(i) == blockL + 1;
        }else{
            right_len2[i] = rank_right(rank_left.size()) - rank_right(i) == rank_left.size() - i;
        }
    }
    Pair<sdsl::bit_vector, sdsl::bit_vector> result(left_len2, right_len2);
    return (result);
}


vector<sdsl::bit_vector> create_bit_vectors(sdsl::rank_support_v<> const &rank_left, sdsl::rank_support_v<> const &rank_right, int const len, vector<int> const sizes){
    vector<sdsl::bit_vector> bit_vectors;
    for(int i = 0; i < sizes.size(); ++i)
    {
        Pair<sdsl::bit_vector, sdsl::bit_vector> result = create_bit_vector_pair(rank_left, rank_right, len, sizes[i]);
        bit_vectors.push_back(getValueI1(result));
        bit_vectors.push_back(getValueI2(result));
    }
    return bit_vectors;    
}


template <typename TChar, typename TAllocConfig>
void loadIndex(vector<sdsl::bit_vector> &bit_vectors, CharString const indexPath)
{    
    typedef String<TChar, TAllocConfig> TString;
    Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
    open(index, toCString(indexPath), OPEN_RDONLY);
    vector<sdsl::bit_vector> bit_vectors_ordered (bit_vectors);
    for(int i = 0; i < bit_vectors.size(); ++i){
        for (unsigned j = 0; j < seqan::length(index.fwd.sa) - 1; ++j)
        {
            uint32_t sa_j = getValueI2(index.fwd.sa[j + 1]);
            bit_vectors_ordered[i][sa_j] = bit_vectors[i][j];
        }
    }
    bit_vectors = bit_vectors_ordered;
}

template <typename TChar>
void loadIndex(vector<sdsl::bit_vector> &bit_vectors, CharString const indexPath, bool const mmap)
{
    if(mmap)
        return(loadIndex<TChar, MMap<> >(bit_vectors, indexPath));
    else
        return(loadIndex<TChar, Alloc<> >(bit_vectors, indexPath));
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
    addOption(parser, ArgParseOption("I", "map", "Path to the mappability file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "map");
    addOption(parser, ArgParseOption("Idx", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");
    addOption(parser, ArgParseOption("O", "output", "Path to output directory", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");
    addOption(parser, ArgParseOption("K", "length", "Length of k-mers in the mappability vector", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");
    addOption(parser, ArgParseOption("T", "threshold", "Threshold for inverse frequency that get accepted", ArgParseArgument::DOUBLE, "DOUBLE"));
    setRequired(parser, "threshold");
    addOption(parser, ArgParseOption("l", "level", "Length of k-mers in the mappability vector", ArgParseArgument::INTEGER, "INT"));
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
    int len, level = 2; 
    double threshold;
    
    getOptionValue(mappability_path, parser, "map");
    getOptionValue(indexPath, parser, "index");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(len, parser, "length");
    getOptionValue(threshold, parser, "threshold");
    bool mmap = isSet(parser, "mmap");
    
    StringSet<CharString> ids;
    CharString alphabet;
    string mappability_str;
    
    
    bool singleIndex;
    _indexPath = indexPath;
    _indexPath += ".singleIndex";
    open(singleIndex, toCString(_indexPath));
    
    _indexPath = indexPath;
    _indexPath += ".alphabet";
    open(alphabet, toCString(_indexPath));

    _indexPath = indexPath;
    _indexPath += ".ids";
    open(ids, toCString(_indexPath), OPEN_RDONLY);
    
//     template <typename TIndex> 
    
    /*
    if(alphabet == "dna4")
        auto index = loadIndex<Dna>(indexPath, mmap);
    else
        auto index = loadIndex<Dna5>(indexPath, mmap);
    */   
/*
    typedef String<Dna, Alloc<> > TString;
    Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
    open(index, toCString(indexPath), OPEN_RDONLY);

    */
    
    // The seqan way
//     open(mappability_str, toCString("../Data/30_mappability0"));
//     cout << mappability_str << endl;
    ifstream myfile;
    myfile.open(mappability_path, std::ios::in | std::ofstream::binary);
    std::getline(myfile, mappability_str);
    myfile.close();
    vector<double> mappability = getDouble(mappability_str);
    

    // for comma delimiter 
    /*
    std::istringstream iss(mappability_str);
    std::vector<std::string> mappability_v (std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>());
    class WordDelimitedByCommas : public std::string
    {};
    std::istream& operator>>(std::istream& is, WordDelimitedByComma& output)
    {
        std::getline(is, output, ',');
        return is;
    }
    */
//     for(unsigned i = 0; i < mappability.size(); i++)
//          cout << mappability[i];
    
    // merge both into one?
    /*
    sdsl::bit_vector rl (mappability.size() + 2*(length + 1), 1);
    for(unsigned i = 0; i < mappability.size(); i++){
        rl[i + length + 1] = (mappability[i] >= threshold);  
    }
    // position length + 1 for l 
    cout << endl;
    */

    sdsl::bit_vector righti (mappability.size() + len - 1, 1);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 1);

    //skip first value since we want pos discribe what happens to 30 position 
    for(unsigned i = 0; i < (mappability.size() - 1); ++i){
        lefti[i + len] = (mappability[i + 1] >= threshold);
        righti[i] = (mappability[i] >= threshold);
    }
    
    sdsl::rank_support_v<> rank_righti(&righti);
    sdsl::rank_support_v<> rank_lefti(&lefti);
    vector<int> sizes{0};
    vector<string> names{"l_bit_vector_1", "r_bit_vector_1"};
    for(int i = 0; i < level; ++i){
        int a = (i + 1) * len - 1;
//         cout << "i" << i << endl;
//         cout << a << endl;
        sizes.push_back(a);
        names.push_back("l_bit_vector_" + to_string(a + 1));
        names.push_back("r_bit_vector_" + to_string(a + 1));
    }
    
    vector<sdsl::bit_vector> bit_vectors = create_bit_vectors(rank_lefti, rank_righti, len, sizes);
  
    for(int i = 0; i < bit_vectors.size(); ++i){
        std::ofstream outfile((toCString(outputPath) + names[i]), std::ios::out | std::ofstream::binary);
        std::copy(bit_vectors[i].begin(), bit_vectors[i].end(), std::ostream_iterator<bool>(outfile));
        outfile.close();
    }
       
    //order in suffix array
    order_bit_vector(bit_vectors, indexPath, mmap, alphabet);
    for(int i = 0; i < bit_vectors.size(); ++i){
        std::ofstream outfile((toCString(outputPath) + names[i] + "_o"), std::ios::out | std::ofstream::binary);
        std::copy(bit_vectors[i].begin(), bit_vectors[i].end(), std::ostream_iterator<bool>(outfile));
        outfile.close();
    }
    
    return 0;
}

/*
Pair<sdsl::bit_vector, sdsl::bit_vector> create_bit_vector_pair(sdsl::rank_support_v<> const &rank_left, sdsl::rank_support_v<> const &rank_right, int const len, int const blockL)
{
    sdsl::bit_vector right_len2 (rank_right.size(), 1);
    sdsl::bit_vector left_len2 (rank_left.size(), 1);
    
    // cant use unsigned because expect negative values (i - len)
    for(int i = 0; i < rank_right.size(); i++)
    {
        //right bound
        if(rank_left.size() < i + len + 1)
        {
             left_len2[i] = rank_left(rank_left.size()) - rank_left(i - len * blockL) == rank_left.size() - i - len * blockL + 1;
        }
        else
        {
            // left bound
            if(i - len * blockL < 0)
            {
                  left_len2[i] = rank_left(i + len + 1) == i + len;
            }
            else
            {
                   left_len2[i] = rank_left(i + len + 1) - rank_left(i - len * blockL) == (blockL + 1) * len + 1;
            }
        }
        //left bound
        if(i < len + 1)
        {
              right_len2[i] = (i == 0) ? 0 : rank_right(i + len * blockL) == i + len * blockL;
        }
        else
        {
            // right bound
            if(rank_right.size() < i + len * blockL + 1){
                  right_len2[i] = rank_right(rank_right.size()) - rank_right(i - len) == rank_right.size() - i - len * blockL + 1;
            }
            else
            {
                right_len2[i] = rank_right(i + len * blockL + 1) - rank_right(i - len) == (blockL + 1) * len + 1;
            }
            
        }
    }
    Pair<sdsl::bit_vector, sdsl::bit_vector> result(left_len2, right_len2);
    return (result);
}
*/
