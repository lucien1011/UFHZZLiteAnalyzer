
#ifndef FRWeightTreeProducer_h
#define FRWeightTreeProducer_h

#include "FRWeightTreeProducer_Linkdef.h"
#include "Analyzer.h"
#include "HZZTree.h"
#include "deltaR.h"

using namespace std;  

class FRWeightTreeProducer : public Analyzer 
{
    public:
        // static const double Zmass = 91.1876;
        const double Zmass = 91.1876;

        FRWeightTreeProducer();
        ~FRWeightTreeProducer();
        FRWeightTreeProducer(
                double m4lHighCut_in,
                double m4lLowCut_in,
                double mZ2HighCut_in,
                double mZ2LowCut_in,
                double mZ1HighCut_in,
                double mZ1LowCut_in,
                double isoCutEl_in,
                double isoCutMu_in,
                TString outputDir_in,
                TString outFileName_in,
                bool do_wrong_fc_in=false
                );
         FRWeightTreeProducer(
                TString outputDir_in,
                TString outFileName_in
                );

        int process();
        void setup();
        void end();
        void initTree();
        bool passSelection();

        double m4lHighCut = 9999999.;
        double m4lLowCut = 70.0;
        double mZ2HighCut=120.0;
        //double mZ2High=999999.;
        double mZ2LowCut=4.0;
        double mZ1HighCut=120.0;
        double mZ1LowCut=40.0;
        //double isoCutEl=999999.;
        double isoCutEl=0.35;
        double isoCutMu=0.35;
        double leadingPtCut=20.0; 
        double subleadingPtCut=10.0; 

        TString treeName = "Ana/passedEvents";
        TString outTreeName = "passedEvents";
        bool debug = false;
        bool do_wrong_fc = false;

        TString outFileName;
        TFile* outFile=0;
        TTree* outTree=0;
};

FRWeightTreeProducer::FRWeightTreeProducer(){

}

FRWeightTreeProducer::~FRWeightTreeProducer(){

}

bool FRWeightTreeProducer::passSelection(){
    return true;
}


void FRWeightTreeProducer::initTree(){
    setHZZTree(tree);
}

void FRWeightTreeProducer::setup(){
    outFile = TFile::Open(outputDir+outFileName,fOptionWrite);
    outTree = new TTree(outTreeName,outTreeName);

    initNewLiteTree(outTree);
}

void FRWeightTreeProducer::end(){
    outFile->cd();
    outTree->Write(outTreeName,TObject::kOverwrite);
    outFile->Close(); 
}

FRWeightTreeProducer::FRWeightTreeProducer(
                double m4lHighCut_in,
                double m4lLowCut_in,
                double mZ2HighCut_in,
                double mZ2LowCut_in,
                double mZ1HighCut_in,
                double mZ1LowCut_in,
                double isoCutEl_in,
                double isoCutMu_in,
                TString outputDir_in,
                TString outFileName_in,
                bool do_wrong_fc_in
                ){
    m4lHighCut  = m4lHighCut_in;
    m4lLowCut   = m4lLowCut_in;
    mZ2HighCut  = mZ2HighCut_in;
    mZ2LowCut   = mZ2LowCut_in;
    mZ1HighCut  = mZ1HighCut_in;
    mZ1LowCut   = mZ1LowCut_in;
    isoCutEl    = isoCutEl_in;
    isoCutMu    = isoCutMu_in;
    outputDir   = outputDir_in;
    outFileName = outFileName_in;
    do_wrong_fc = do_wrong_fc_in;
}

FRWeightTreeProducer::FRWeightTreeProducer(
                TString outputDir_in,
                TString outFileName_in
                ){
    outputDir   = outputDir_in;
    outFileName = outFileName_in;
}

int FRWeightTreeProducer::process(){
