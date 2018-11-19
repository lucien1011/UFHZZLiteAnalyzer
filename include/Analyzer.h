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
        void printProgress(int iEvt);
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

        TString outputDir;
};

ClassImp(Analyzer)

Analyzer::Analyzer(){

}

Analyzer::~Analyzer(){

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

    Long64_t nentries = tree->GetEntries();

    for(int iEvt=0; iEvt < nentries; iEvt++){
        tree->GetEntry(iEvt);

        printProgress(iEvt);

        if (!passSelection()) continue;
        
        process();
    };
    
    end();
}

void Analyzer::initTree(){

}

void Analyzer::printProgress(int iEvt){
    if(iEvt%100000==0) cout<<"Event "<<iEvt<<"/"<<tree->GetEntries()<<endl;
}

#endif
