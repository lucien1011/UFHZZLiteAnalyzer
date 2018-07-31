#ifndef ZZ4LConfig_h
#define ZZ4LConfig_h


using namespace std;  

bool isData;
bool bestCandMela=false;
bool redoEventSelection=true;
bool redoJets=true;
bool redoMela=false;
bool redoEbE=false;
float dxycut=0.5;
float dzcut=1.0;
double sip3dcut=4.0;
float lep_ptcut=5.0;
float vtxcut =0.0;
double m4lLowCut=105.0;
double m4lHighCut=140.0;
double mZ2Low=4.0;
//double mZ2Low=12.0;
double mZ2High=120.0;
double mZ1Low=40.0;
double mZ1High=120.0;
double isoCutEl=0.35;
double isoCutMu=0.35;
double leadingPtCut=20.0;
double subleadingPtCut=10.0;
double muPtCut = 5.0;
double elPtCut = 7.0;
double BTagCut=0.8484;
int job=-1;
int njobs=-1;

#endif
