#ifndef Analyzer_h
#define Analyzer_h

//#include "Analyzer_Linkdef.h"
#include "Library.h"

//TString fOptionWrite = "RECREATE";
//TString fOptionRead = "READ";

class Analyzer : public TObject {
    protected:
        TTree* tree;
        TFile* file;
        Long64_t nentries;
        void printProgress(int iEvt);
        int barWidth = 70;
        float progress;
        TString lineDivider = "================================================================================================";
    public:
    public:
        Analyzer();
        ~Analyzer();
        ClassDef(Analyzer,1);
        void getTree(TString treeName);
        void getFile(TString filePath);
        virtual int process();
        virtual void setup();
        virtual void end();
        virtual void initTree();
        virtual bool passSelection();
        void loop(TString filePath, TString treeName);
        void printInit(TString fileName);        
        void printEnd(); 

        TString outputDir;
};

ClassImp(Analyzer)

Analyzer::Analyzer(){

}

Analyzer::~Analyzer(){

}

void Analyzer::printInit(TString fileName){
    std::cout << lineDivider << std::endl;
    std::cout << "Processing "+fileName << std::endl;
    std::cout << "Number of entries " << nentries << std::endl;
}

void Analyzer::printEnd(){
    std::cout << std::endl;
    std::cout << lineDivider << std::endl;
}

bool Analyzer::passSelection(){
    return true;
};

void Analyzer::getTree(TString treeName){
    tree = (TTree*) file->Get(treeName);
}

void Analyzer::getFile(TString filePath) {
    file = TFile::Open(filePath,fOptionRead);
}

int Analyzer::process(){
    return 1; 
}

void Analyzer::setup(){
    
}

void Analyzer::end(){
    
}

void Analyzer::loop(TString filePath, TString treeName){

    getFile(filePath);

    getTree(treeName);

    initTree();
    
    setup();

    nentries = tree->GetEntries();

    nentries = tree->GetEntries();

    printInit(filePath);

    for(int iEvt=0; iEvt < nentries; iEvt++){
        tree->GetEntry(iEvt);

        printProgress(iEvt+1);

        if (!passSelection()) continue;
        
        process();
    };
    
    end();

    printEnd();
}

void Analyzer::initTree(){

}

void Analyzer::printProgress(int iEvt){
    //if(iEvt%100000==0) cout<<"Event "<<iEvt<<"/"<<tree->GetEntries()<<endl;
    if(iEvt%10000==0){
        
        progress = float(iEvt)/float(nentries);
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
    }
}

#endif
