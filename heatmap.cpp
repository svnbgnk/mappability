// #include "common.h"
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <fstream>
#include <tgmath.h>

#include <seqan/arg_parse.h>

using namespace std;
using namespace sdsl;
using namespace seqan;

#define HEIGHT 1080
#define WIDTH 1920
#define HISTROGRAM_SIZE 16
#define BUCKET_WIDTH 15

vector<char> encode(sdsl::bit_vector const & bv, char const limit = SCHAR_MAX)
{
    auto it = bv.begin();
//     cout << static_cast<int>(limit) << endl;
    char counter;
    vector <char> bv_encoded;
    while(it != bv.end()){
//         cout << "it" << endl;
        counter = SCHAR_MIN;
        while(*it == 0 && it != bv.end() && counter < limit){
            ++counter;
//             cout << static_cast<int>(counter) << endl;
            ++it;
        }
        bv_encoded.push_back(counter);
        if(it == bv.end())
            return(bv_encoded);
        counter = SCHAR_MIN;
        while(*it == 1 && it != bv.end() && counter < limit){
            counter++;
            ++it;
        }
        bv_encoded.push_back(counter);
    }
    /*
    int *b_encode = new int [bv_encoded.size()];
    for(int i = 0; i < bv_encoded.size(); ++i)
        b_encode[i] = bv_encoded[i];
    std::pair<int * , int> result (b_encode, bv_encoded.size());*/
    return(bv_encoded);
}

void heatmap(sdsl::bit_vector & b, string output)
{
    int height = HEIGHT;
    int width = WIDTH;
    sdsl::rank_support_v<> rb(&b);
    if(b.size() < height * width - 1){
        cout << "changed image dimension to" << endl;
        height = round(pow(static_cast<float>(b.size()) / (8 * 16 / 9), 0.5));
        width = round(static_cast<float> (height) * (9 / 16));
        cout << "height: " << height << " width: " << width << endl;
    }
    cout << "Bases per pixel: "  << b.size()/(height * width) << endl;
    cout << "Number of non-mappability bits overall: " << rb.size() - rb(b.size()) << endl;
    int window = floor(static_cast<float>(b.size()) / (height * width));
    int pos = 0;

    ofstream img(output + ".ppm");
    img << "P3" << endl;
    img << width << " " << height << endl;
    img << "255" << endl;
    for(int i = 0; i < height * width; ++i)
    {
        float den = static_cast<float> (window - rb(pos + window) - rb(pos)) / window;
        int grey = static_cast<int>(round(den * 254));
        img << grey  << " " << grey << " " << grey << endl;
        pos += window;

    }
    img.close();
//     system("eog heatmap.ppm");
}

vector<int> histogram(sdsl::bit_vector const & b , int const his_size, int const bucket_width)
{
    vector<int> hist(his_size, 0);
    auto it = b.begin();
    int n = 0;
    while(it != b.end()){
        while(*it == 0 && it != b.end()){
            ++it;
            ++n;
        }
        if(n != 0){
            if(n < ((his_size) * bucket_width))
                hist[floor(n / bucket_width)] += 1;
        }
        n = 0;
        if(it != b.end())
            ++it;
    }
    return(hist);
}


vector<int> histogram(sdsl::bit_vector b ,const int his_size){
    vector<int> hist_lg(his_size, 0);
    sdsl::bit_vector::iterator it = b.begin();
    int n = 0;
    while(it != b.end()){
        while(*it == 0 && it != b.end()){
            ++it;
            ++n;
        }
        if(n != 0){
             if(n < static_cast<int>(pow(2, his_size))){
                hist_lg[static_cast<int>(ceil(log2(static_cast<float>(n))))] += 1;
             }
        }
        n = 0;
        if(it != b.end())
            ++it;
    }
    return(hist_lg);
}

template<typename T>
void print_vector(vector <T> & v){
    for(int i = 0; i < v.size(); ++i){
        cout << v[i] << " ";
    }
    cout << "\n";
}


// void dotest(vector<int> test){
//     for(int i = 0; i < test.size(); ++i){
//         sdsl::bit_vector rv = create_random_bit_v(test[i]);
//         /*
//         for(sdsl::bit_vector::iterator it = rv.begin(); it != rv.end(); ++it){
//             cout << *it;
//         }
//         cout << endl;
//         */
//         vector <char> rv_encoded = encode(rv);
//         vector <char> rv_encoded2 = encode(rv, -1);
//         vector <char> rv_encoded3 = encode(rv, -125);
//     //     std::pair <int *, int> result = encode(rv);
//     //     int * rv_encoded = result.first;
//     //     int rv_encoded_size = result.second;
//     //     sdsl::bit_vector::iterator it = rv.begin();
//
//         cout << "encoded" << endl;
//         cout << "length bit_v: " << rv.size() << endl;
//         cout << sizeof(rv[1]) << "??" << endl;
//         cout << "Overall size: " << rv.size()/8 << endl;
//         cout << "length encoded: " << rv_encoded.size() << endl;
//         cout << sizeof(rv_encoded[0]) << endl;
//         cout << "Overall size: " << rv_encoded.size() * sizeof(rv_encoded[0]) << endl;
//
//         float comp = ((float)rv.size()/8)/((float)rv_encoded.size()*sizeof(rv_encoded[0]));
//         float comp2 = ((float)rv.size()/8)/((float)rv_encoded2.size()/2);
//         float comp3 = ((float)rv.size()/8)/((float)rv_encoded3.size()/4);
//
//         cout << "Compression rate: " << comp << endl;
//         cout << "Compression rate2(4bit): " << comp2 << endl;
//         cout << "Compression rate2(2bit): " << comp3 << endl;
//         }
// }

std::vector<int> getInt(std::string const& mappability_str)
{
  std::istringstream iss(mappability_str);
  return std::vector<int>{
    std::istream_iterator<int>(iss),
    std::istream_iterator<int>()
  };
}

int main(int argc, char *argv[])
{
    // Argument Parser
    ArgumentParser parser("Heatmap Creation");
    addDescription(parser, "App for drawing a heatmap (highlighting repeat regions).");

    addOption(parser, ArgParseOption("M", "mappability-vector", "Text file with integers separated by a space character", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "mappability-vector");

    addOption(parser, ArgParseOption("O", "output", "Heatmap output file (*.ppm)", ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output");

    addOption(parser, ArgParseOption("T", "threshold", "Threshold for repeats (number of hits necessary to mark it as a repeat)", ArgParseArgument::INTEGER, "INT"));
	setRequired(parser, "threshold");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    unsigned threshold;
    CharString mappability_path, output;
    getOptionValue(mappability_path, parser, "mappability-vector");
    getOptionValue(output, parser, "output");
    getOptionValue(threshold, parser, "threshold");

    bit_vector b;
    {
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
            b[i] = !(mappability_int[i] <= threshold);
        cout << "Bit vector constructed" << endl;
    }

    heatmap(b, toCString(output));
    int his_size = HISTROGRAM_SIZE;
    int bucket_width = BUCKET_WIDTH;

    vector<int> hist = histogram(b, his_size, bucket_width);
    vector<int> hist_lg = histogram(b, his_size);
    print_vector(hist);
    cout << "log" << endl;
    print_vector(hist_lg);


    return 0;
}



// 10,000
// Compression rate: 0.249401
// Compression rate2(4bit): 0.498803
// Compression rate2(2bit): 0.774953

// 100,000

// Compression rate: 0.251869
// Compression rate2(4bit): 0.503738
// Compression rate2(2bit): 0.779873

// 1,000,000
// Compression rate: 0.250062
// Compression rate2(4bit): 0.500123
// Compression rate2(2bit): 0.778362


// 10,000,000
// Compression rate: 0.24991
// Compression rate2(4bit): 0.499819
// Compression rate2(2bit): 0.777673

// 100,000,000
//Compression rate: 0.250017
//Compression rate2(4bit): 0.500035
//Compression rate2(2bit): 0.777727
