#include "common.h"
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <climits>
using namespace std;
using namespace seqan;


sdsl::bit_vector create_random_bit_v(int length){
//     srand (time(0));
    sdsl::bit_vector rv (length, 1);
    srand (time(0));
    for(sdsl::bit_vector::iterator it = rv.begin(); it != rv.end(); ++it){
        int a = (rand() % 2);
        *it = a;
    }
    return(rv);
}




vector<char> encode(sdsl::bit_vector & bv){
    sdsl::bit_vector::iterator it = bv.begin();
    char limit = SCHAR_MAX; //pow(2, 8) - 1;
//     cout << static_cast<int>(limit) << endl;
//     char one = SCHAR_MIN + 1;
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


vector<char> encode(sdsl::bit_vector & bv, char limit){
    sdsl::bit_vector::iterator it = bv.begin();
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

void dotest(vector<int> test){
    for(int i = 0; i < test.size(); ++i){
        sdsl::bit_vector rv = create_random_bit_v(test[i]);
        /*
        for(sdsl::bit_vector::iterator it = rv.begin(); it != rv.end(); ++it){
            cout << *it;
        }
        cout << endl;
        */
        vector <char> rv_encoded = encode(rv);
        vector <char> rv_encoded2 = encode(rv, -1);
        vector <char> rv_encoded3 = encode(rv, -125);
    //     std::pair <int *, int> result = encode(rv);
    //     int * rv_encoded = result.first;
    //     int rv_encoded_size = result.second;
    //     sdsl::bit_vector::iterator it = rv.begin();
    
        cout << "encoded" << endl;    
        cout << "length bit_v: " << rv.size() << endl;
        cout << sizeof(rv[1]) << "??" << endl;
        cout << "Overall size: " << rv.size()/8 << endl;
        cout << "length encoded: " << rv_encoded.size() << endl;
        cout << sizeof(rv_encoded[0]) << endl;
        cout << "Overall size: " << rv_encoded.size() * sizeof(rv_encoded[0]) << endl;

        float comp = ((float)rv.size()/8)/((float)rv_encoded.size()*sizeof(rv_encoded[0]));
        float comp2 = ((float)rv.size()/8)/((float)rv_encoded2.size()/2);
        float comp3 = ((float)rv.size()/8)/((float)rv_encoded3.size()/4);
    
        cout << "Compression rate: " << comp << endl;
        cout << "Compression rate2(4bit): " << comp2 << endl;
        cout << "Compression rate2(2bit): " << comp3 << endl;  
        }
}


int main(int argc, char *argv[])
{
    vector<int> test = {10000, 100000, 1000000, 10000000, 100000000};
    dotest(test);
    
    

    //char test = -128; //127    
    
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

