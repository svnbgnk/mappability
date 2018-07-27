#ifndef GLOBAL_H_
#define GLOBAL_H_


/*
normal:
all cases on
filter_th: 0.1 // 0.3
invflipdensity: 0.5
intervalsize = 2;
directsearchblockoffset: 2
directsearch_th: 3

step: 2 //4
distancetoblockend: 1  // 2

*/

struct majorCaseParameters{
private:
    uint32_t cases;
public:
    bool nomappability;
    bool directsearch;
    bool compmappable;
    bool suspectunidirectional;

    bool testflipdensity;
    uint32_t step;
    uint32_t distancetoblockend;
    uint32_t directsearch_th;
    uint32_t directsearchblockoffset;
    float filter_th;
    float invflipdensity;
    uint32_t intervalsize;

    majorCaseParameters(){
        setdefault();
    }

    void setbestnormal(){
        filter_th = 0.3;
        invflipdensity = 0.5;
        intervalsize = 2;
        directsearchblockoffset = 2;
        directsearch_th = 3;
        step = 2;
        distancetoblockend = 1;
    }


    void setbestnormalhg(){
        filter_th = 0.1; //around
        invflipdensity = 0.3; // lower or turned off
        intervalsize = 2; //1-3 or maybe 89 0.1 filter 109 invflipdensity turned off
        directsearchblockoffset = 4;
        directsearch_th = 4; // around
        step = 2;  // or 8
        distancetoblockend = 1; //????
    }

        void setbestnormalhgv2(){
        filter_th = 0.1;
        invflipdensity = 0.1;
        intervalsize = 3;
        directsearchblockoffset = 6;
        directsearch_th = 2;
        step = 2;
        distancetoblockend = 1;
    }

    void setdefault(){
        nomappability = true;
        directsearch = true;
        compmappable = true;
        suspectunidirectional = true;

        testflipdensity = true;
        //binaryNumber
        step = 0b11;
        distancetoblockend = 2;

        directsearchblockoffset = 0;
        directsearch_th = 2;
        filter_th = 0.5;

        invflipdensity = 0.5;

        intervalsize = 3;
    }

    void setCases(uint32_t a){
        if(a > 16){
            std::cerr << "Integer to large. Try again!!! (<16)" << "\n";
        }else{
            std::bitset<4> bit (a);
            nomappability = (bit[0]);
            directsearch = (bit[1]);
            compmappable = (bit[2]);
            suspectunidirectional = (bit[3]);
            cases = a;
        }
    }

    void printCases(){
        std::cout << std::bitset<4>(cases) << "\t";
    }

    void print(){
        std::cout << "Cases Enabled: " << "\n";
        std::cout << nomappability << " " << directsearch << " " << compmappable << " " << suspectunidirectional << "\n";
        std::cout << "Params: " << "\n";

        std::cout << "step: " << step << "\n";
        std::cout << "distancetoblockend: " << distancetoblockend << "\n";
        std::cout << "directsearchblockoffset: " << directsearchblockoffset << "\n";
        std::cout << "directsearch_th: " << directsearch_th << "\n";
        std::cout << "filter_th: " << filter_th << "\n";
        std::cout << "invflipdensity: " << invflipdensity << "\n";
        std::cout << "intervalsize: " << intervalsize << "\n";
    }
};


struct myGlobalParameters{
public:
    bool startUnidirectional;
    bool clocking;
    bool wasStopped;
//     std::chrono::duration<double> terminateDuration;
    std::chrono::seconds dtd = std::chrono::seconds{11};
    std::chrono::seconds terminateDuration;

    majorCaseParameters normal;
    majorCaseParameters comp;
    majorCaseParameters uni;
    majorCaseParameters startuni;

    myGlobalParameters(){
        startUnidirectional = false;
        clocking = false;
        wasStopped = false;
        setdefaultTime();
    }

    void setdefaultTime(){
//         std::chrono::microseconds(sec).count()
        terminateDuration = dtd;
    }

    void setdefault(){
        startUnidirectional = false;
        wasStopped = false;
        normal.setdefault();
        comp.setdefault();
        uni.setdefault();
        startuni.setdefault();
    }

    //INFO need to copy more parameters for startUni
    void copyDirectsearchParamsfromNormal(){
        comp.directsearchblockoffset = normal.directsearchblockoffset;
        comp.directsearch_th = normal.directsearch_th;

        uni.directsearchblockoffset = normal.directsearchblockoffset;
        uni.directsearch_th = normal.directsearch_th;

        startuni.directsearchblockoffset = normal.directsearchblockoffset;
        startuni.directsearch_th = normal.directsearch_th;
    }

    void print(){
        std::cout << "\n";
        std::cout << "Start Unidirectional: " << startUnidirectional << "\n";
        std::cout << "Was Stopped: " << wasStopped << "\n";
        std::cout << "Normal: ";
        normal.print();
        std::cout << "Comp: ";
        comp.print();
        std::cout << "Uni: ";
        uni.print();
        std::cout << "Startuni: ";
        startuni.print();
    }
};

struct checkTime{
public:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int steps = 1024;
    int k;
    bool stopped;

    checkTime(){
        k = 0;
        stopped = false;
        start = std::chrono::high_resolution_clock::now();
    }

    bool stopnow(std::chrono::seconds & d){
        if(stopped)
            return true;
        ++k;
        if(k % steps == 0)
        {
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed  = end - start;
            if(end - start > d){
                std::cout << "Stopped at:" << elapsed.count() << "s" << "\n";
                stopped = true;
                return(true);
            }
            else
            {
                return(false);
            }
        }
        else
        {
            return(false);
        }
    }
};

extern myGlobalParameters params;


/*
struct trackCases{
public:
    int returncodes[10];
};*/




#endif
