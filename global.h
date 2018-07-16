#ifndef GLOBAL_H_
#define GLOBAL_H_

struct majorCaseParameters{
    bool nomappability = true;
    bool directsearch = true;
    bool compmappable = true;
    bool suspectunidirectional = true;
    
    //binaryNumber
    int step = 0b11;
    int distancetoblockend = 2;
    
    int directsearch_th = 2;
    float filter_th = 0.5;
    
    float invflipdensity = 0.5;
    
    int intervalsize = 3;
    
    void print(){
        std::cout << "Cases Enabled: " << "\n";
        std::cout << nomappability << " " << directsearch << " " << compmappable << " " << suspectunidirectional << "\n";
        std::cout << "Params: " << "\n";
        std::cout << "step: " << step << "\n";
        std::cout << "distancetoblockend: " << distancetoblockend << "\n";
        std::cout << "directsearch_th: " << directsearch_th << "\n";
        std::cout << "filter_th: " << filter_th << "\n";
        std::cout << "invflipdensity: " << invflipdensity << "\n";
        std::cout << "intervalsize: " << intervalsize << "\n";
    }
};


struct myGlobalParameters{
public:
    bool startUnidirectional = false;
    majorCaseParameters normal;
    majorCaseParameters comp;
    majorCaseParameters uni;
    majorCaseParameters startuni;
    
    
    void print(){
        normal.print();
        comp.print();
        uni.print();
    }
};

extern myGlobalParameters params;
extern int global;


/*
struct trackCases{
public:
    int returncodes[10];
};*/




#endif
