#include "common.h"
#include <iostream>
#include <fstream>
#include <iterator> 
#include <sstream>
#include <tgmath.h>
#include <sdsl/bit_vectors.hpp>

using namespace std;
// using namespace seqan;

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


/*


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

        // for (unsigned i = 0; i < 200; ++i)
        //     cout << mappability_int[i] << ' ';

        b.resize(mappability_int.size());
        for (unsigned i = 0; i < mappability_int.size(); ++i)
            b[i] = !(mappability_int[i] > threshold);
        cout << "Bit vector constructed" << endl;

*/
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
//     for(vector<float>::iterator it = mappability.begin(); mappability.end() != it; ++it)
//         *it = static_cast<float>(1) / *it;
    return(mappability_int); 
}


pair<vector<string>, vector<sdsl::bit_vector>> create_bit_vectors(const vector <uint8_t> & mappability, const int len, double threshold, const int errors){
    
    int th = round(1/threshold);    
    uint8_t blocks = errors + 2;
    uint32_t blocklength = len / blocks;
    uint8_t rest = len - blocks * blocklength;
    std::vector<uint32_t> blocklengths;
    std::vector<uint32_t> revBlocklengths(4, 0);
    for (uint8_t i = 0; i < blocks; ++i)
        blocklengths.push_back(blocklength + (i < rest));
    
    revBlocklengths[blocks - 1] = blocklengths[blocks - 1];    
    for(int8_t i = blocks - 2; i >= 0; --i)
        revBlocklengths[i] += blocklengths[i] + revBlocklengths[i + 1];
       
    for(uint8_t i = 1; i < blocks; ++i)
        blocklengths[i] += blocklengths[i - 1];
   
    //revBlocklength is not the reverse of blocklength since the blocklength values can differ!
/*    cout << "revBlockLengths" << endl;
    for(uint8_t i = 0; i < revBlocklengths.size(); ++i)
        cout << revBlocklengths[i] << endl;
//     
//     cout << "BlockLengths" << endl;
    for(uint8_t i = 0; i < blocklengths.size(); ++i)
        cout << blocklengths[i] << endl;
  */  
    sdsl::bit_vector righti (mappability.size() + len - 1, 0);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 0);
    #pragma omp parallel for schedule(static)   
    for(unsigned i = 0; i < mappability.size(); ++i){
        lefti[i + len - 1] = (mappability[i] <= th);
        righti[i] = (mappability[i] <= th);
//         cout << righti[i];
    }
//     cout << endl;
    cout << "Finished Default Bit Vectors.  Length: " << righti.size() << endl;
    vector<sdsl::bit_vector> bit_vectors;
    vector<string> names;
    
    bit_vectors.push_back(righti);
    names.push_back("r_bit_vector_" + to_string(len) + "_shift_0");

    
    if(errors != 0){
        for(int i = 0; i < blocks - 1; ++i){
            sdsl::bit_vector newright(mappability.size() + len - 1, 0); //TODO think 0 or 1 in edge cases
            int shift = blocklengths[i];
            cout << "shift for r_bit  " << shift << endl;
            cout << "name:  " << i + 1 << endl; 
            for(int j = 0; j < righti.size(); ++j){
                if(j - shift >= 0)
                    newright[j] = righti[j - shift];
            }
            bit_vectors.push_back(newright);
            names.push_back("r_bit_vector_" + to_string(len) + "_shift_" + to_string(i + 1));
        }
        
        for(int i = 1; i < blocks; ++i){
            sdsl::bit_vector newleft(mappability.size() + len - 1, 0);//TODO think 0 or 1 in edge cases
            int shift = revBlocklengths[i];
            cout << "shift for l_bit  " << shift << endl;
            cout << "name:  " << blocks - i << endl; 
            for(int j = 0; j < righti.size(); ++j){
                if(j + shift < lefti.size() - 1)
                    newleft[j] = lefti[j + shift];
            }
            bit_vectors.push_back(newleft);
            names.push_back("l_bit_vector_" + to_string(len) + "_shift_" + to_string(blocks - i));
        }
    }
    
    bit_vectors.push_back(lefti);
    names.push_back("l_bit_vector_" + to_string(len) + "_shift_0");

    return(std::make_pair(names, bit_vectors));
}


void print_SA(CharString const indexPath, vector<sdsl::bit_vector> &bit_vectors, CharString const & outputPath, bool const fwd){
    string name = "_SA_debug";
    string dir = (fwd) ? "fwd" : "rev";
    std::ofstream outfile((toCString(outputPath) + dir + name), std::ios::out | std::ofstream::binary);
    typedef String<Dna, Alloc<>> TString;
    Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
    open(index, toCString(indexPath), OPEN_RDONLY);
    
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
    uint32_t sa_j;
    uint16_t seq;
    for (unsigned j = 0; j < seqan::length(index.fwd.sa); ++j)
    {
        if(fwd){
            sa_j = index.fwd.sa[j].i2;
            seq = index.fwd.sa[j].i1;
        }else{
            sa_j = index.rev.sa[j].i2;
            seq = index.rev.sa[j].i1;
        }
        outfile << j << " " << "(" << seq << ", " << sa_j << "):\t" << sa_j + sequenceLengths[seq] << endl;
    }
    outfile.close();
//     for(int i = 0; i < 10; ++i){
//         cout << index.fwd.sa[i] << endl;
//     }
//     cout << "reverse:" << endl;
//        for(int i = 0; i < 10; ++i){
//         cout << index.rev.sa[i] << endl;
//     }
}


template <typename TChar, typename TAllocConfig>
void loadIndex(vector<sdsl::bit_vector> &bit_vectors, CharString const indexPath)
{
    cout << mytime() << "Loading Index" << endl;
    typedef String<TChar, TAllocConfig> TString;
    Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
    open(index, toCString(indexPath), OPEN_RDONLY);
    cout << mytime() << "Loaded Index. Size:" << seqan::length(index.fwd.sa) << endl;
    vector<sdsl::bit_vector> bit_vectors_ordered (bit_vectors);    
    int number_of_indeces = seqan::length(index.fwd.sa) - bit_vectors[0].size();
    vector<int> sequenceLengths(number_of_indeces + 1, 0);
    cout << "Number of Indeces: " << number_of_indeces << endl;
    
    int ssize = sequenceLengths.size();
    //sequenceLengths first value is 0
    for(int i = 0; i < ssize - 1; ++i)
        sequenceLengths[(index.fwd.sa[i]).i1 + 1] = index.fwd.sa[i].i2;  
    for(int i = 1; i < ssize; ++i)
        sequenceLengths[i] += (sequenceLengths[i - 1]);
    cout << "Sequence Lengths:" << endl;
    for(int i = 0; i < ssize; ++i){
        cout << sequenceLengths[i] << endl;        
    }

    // skip sentinels
    #pragma omp parallel for schedule(static)
    for (unsigned j = 0; j < seqan::length(index.fwd.sa) - number_of_indeces; ++j)
    {
        uint32_t sa_f = index.fwd.sa[j + number_of_indeces].i2;
        uint16_t seq_f = index.fwd.sa[j + number_of_indeces].i1;
        uint32_t sa_r = index.rev.sa[j + number_of_indeces].i2;
        uint16_t seq_r = index.rev.sa[j + number_of_indeces].i1;
        
//         cout << j << " < " << (sa_j + sequenceLengths[seq]) << endl;
        for(int i = 0; i < bit_vectors.size()/2; ++i){
            bit_vectors_ordered[i][j] = bit_vectors[i][sa_f + sequenceLengths[seq_f]];
//             cout << i << endl;
        }
        for(int i = bit_vectors.size()/2; i < bit_vectors.size(); ++i){
            
//             cout << i << endl;
//             int calc = sequenceLengths[seq_r + 1] - sa_r - 1;
//             cout << "Seq: " << seq_r << "  SA:" << sa_r << "   calc: " << calc << endl;
            bit_vectors_ordered[i][j] = bit_vectors[i][sequenceLengths[seq_r + 1] - sa_r - 1];
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
    
    addOption(parser, ArgParseOption("E", "errors", "Max errors allowed during mapping", ArgParseArgument::INTEGER, "INT"));
    
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
    int len, errors = 0; 
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
    
    _indexPath = indexPath;
    _indexPath += ".alphabet";
    open(alphabet, toCString(_indexPath));

    if(!file_exists(mappability_path))
    {
        cout << "Cannot find mappability file" << endl;
        exit(0);
    }
    
    cout << mytime() << "Program start." << endl;
    vector<uint8_t> mappability = read(mappability_path);
    cout << mytime() << "Loaded Mappability vector. Size: " << mappability.size() << endl;
//     for(int i = 0; i < mappability.size(); ++i)
//         cout << (int)mappability[i];
//     cout << endl;
    // TODO Merge both bit_vectors into for creation
    
    
    pair<vector<string>, vector<sdsl::bit_vector>> result = create_bit_vectors(mappability, len, threshold, errors);
    vector<string> names = result.first;
    vector<sdsl::bit_vector> bit_vectors = result.second;
    cout << mytime() << "Finished bit vectors." << endl;

    if(debug)
    {
        for(int i = 0; i < bit_vectors.size(); ++i){
            std::ofstream outfile((toCString(outputPath) + names[i] + "_debug"), std::ios::out | std::ofstream::binary);
            std::copy(bit_vectors[i].begin(), bit_vectors[i].end(), std::ostream_iterator<bool>(outfile));
            outfile.close();
            
        }
        print_SA(indexPath, bit_vectors, outputPath, true);
        print_SA(indexPath, bit_vectors, outputPath, false);
    
    }
    
    cout << "Start sorting" << endl;
    //order in suffix array
    order_bit_vector(bit_vectors, indexPath, mmap, alphabet);
    cout << mytime() << "Ordering (Suffix array) bit vectors" << endl;
    for(int i = 0; i < bit_vectors.size(); ++i){
        sdsl::store_to_file(bit_vectors[i], toCString(outputPath) + names[i]);
        if(debug){
            std::ofstream outfile((toCString(outputPath) + names[i] + "_osa_debug"), std::ios::out | std::ofstream::binary);
            std::copy(bit_vectors[i].begin(), bit_vectors[i].end(), std::ostream_iterator<bool>(outfile));
            outfile.close();
        }
    }   
    cout << mytime() << "Finished saving bit vectors" << endl;

    return 0;
}
