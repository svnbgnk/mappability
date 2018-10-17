#ifndef NEW_STATS_H
#define NEW_STATS_H

#include "global.h"
using namespace std;
using namespace seqan;


std::vector<uint32_t> histogram(std::vector<uint32_t> & b , int const his_size, int const bucket_width)
{
    std::vector<uint32_t> hist(his_size, 0);
    auto it = b.begin();
    while(it != b.end()){
        if(*it < ((his_size) * bucket_width))
        {
            hist[floor(*it / bucket_width)] += 1;
        }
        else
        {
            ++hist[his_size - 1];
        }
        ++it;
    }
    return(hist);
}

void readOccurrences(StringSet<DnaString> const & reads, StringSet<CharString> & ids, string outputpath, bool const fr, bool const rc, bool const stats, bool const notmy)
{
        int found = 0, foundD = 0;
        int notfound = 0, mymiss = 0, same = 0, nice = 0, verynice = 0;
        std::vector<uint8_t> mycase(length(reads), 255);

        if(rc){
            int nr = readOccCountDeT.size()/2;
            for(int i = 0; i < readOccCountDeT.size()/2; ++i)
            {
                readOccCount[i] += readOccCount[i + nr];
                readOccCountDeT[i] += readOccCountDeT[i + nr];
            }
            readOccCount.erase(readOccCount.begin() + nr, readOccCount.end());
            readOccCountDeT.erase(readOccCountDeT.begin() + nr, readOccCountDeT.end());
            mycase.erase(mycase.begin() + nr, mycase.end());
        }
        if(fr || stats){
            for(int i = 0; i < readOccCountDeT.size(); ++i)
            {
                found += readOccCount[i] > 0;
                foundD += readOccCountDeT[i] > 0;
//                 cout << readOccCount[i] << " - " <<  readOccCountDeT[i] << std::endl;
                if(readOccCount[i] == readOccCountDeT[i]){
                    if(readOccCount[i] == 0){
                        ++notfound;
                        mycase[i] = 0;
                    }
                    else
                    {
                        ++same;
                        mycase[i] = 2;
                    }
                }
                else
                {
                    if(readOccCount[i] == 0){
                        ++mymiss;
                        mycase[i] = 1;
                    }
                    else
                    {
                        ++nice;
                        if(static_cast<double>(readOccCountDeT[i]) / readOccCount[i] > 2)
                            ++verynice;
                        mycase[i] = 3;
                        if(readOccCount[i] > readOccCountDeT[i]){
                            cerr << "More occurrences with mappability" << std::endl;
                            exit(0);
                        }
                    }
                }
                if(fr){
                    //write filtered fastas
                    SeqFileOut seqFileout0(toCString(outputpath + "/notfound.fa"));
                    SeqFileOut seqFileout1(toCString(outputpath + "/mymiss.fa"));
                    SeqFileOut seqFileout2(toCString(outputpath + "/same.fa"));
                    SeqFileOut seqFileout3(toCString(outputpath + "/nice.fa"));
                    for(int i = 0; i < mycase.size(); ++i)
                    {
                        switch(mycase[i])
                        {
                            case 0: writeRecord(seqFileout0, ids[i], reads[i]); break;
                            case 1: writeRecord(seqFileout1, ids[i], reads[i]); break;
                            case 2: writeRecord(seqFileout2, ids[i], reads[i]); break;
                            case 3: writeRecord(seqFileout3, ids[i], reads[i]); break;
                            default: break;
                        }
                    }
                    close(seqFileout0);
                    close(seqFileout1);
                    close(seqFileout2);
                    close(seqFileout3);
                }
            }
            std::cout << "reads found with mappability: " << found << std::endl;
            std::cout << "reads found without considering mappability: " << foundD << std::endl;
            std::cout << "not found: " << notfound << std::endl;
            std::cout << "mymiss: " << mymiss << std::endl;
            std::cout << "same: " << same << std::endl;
            std::cout << "nice: " << nice << std::endl;
            std::cout << "thereof verynice: " << verynice << std::endl;

            // investigating the vectors
            int bucketSize = 10;
            int histSize = 10;

            if(!notmy){
                std::vector<uint32_t> h = histogram(readOccCount, histSize, bucketSize);
                std::cout << "Histogram buckets size " << bucketSize << ": " << std::endl;
                for(int i = 0; i < h.size(); ++i){
                    std::cout << bucketSize*(i + 1) - 1 << "\t";
                }
                std::cout << std::endl;
                for(int i = 0; i < h.size(); ++i){
                    std::cout << h[i] << "\t";
                }
                std::cout << std::endl;
            }

            std::vector<uint32_t> hDeT = histogram(readOccCountDeT, histSize, bucketSize);
            cout << "Histogram buckets size " << bucketSize << ": " << std::endl;
            for(int i = 0; i < hDeT.size(); ++i){
                std::cout << bucketSize*(i + 1) - 1 << "\t";
            }
            std::cout << std::endl;
            for(int i = 0; i < hDeT.size(); ++i){
                std::cout << hDeT[i] << "\t";
            }
            std::cout << std::endl;
        }
}


#endif

