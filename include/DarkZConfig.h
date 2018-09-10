#ifndef DarkZConfig_h
#define DarkZConfig_h


using namespace std;  

bool isData;
bool bestCandMela=false;
bool redoJets=false;
bool redoEbE=false;

// HZZ, dark photon, Jpsi or Upsilon
//double m4lLowCut=0.0;
//double m4lHighCut=9999999.;
double m4lLowCut=70.0;
double m4lHighCut=999999.;
//double m4lLowCut=105.0;
//double m4lHighCut=140.0;

// HZZ, dark photon, Jpsi or Upsilon
double mZ2High=120.0;
//double mZ2High=999999.;
double mZ2Low=4.0;
//double mZ2Low=12.0;
//double mZ2Low=0.;

// Z
double mZ1Low=40.0;
double mZ1High=120.0;

float dxycut=0.5;
float dzcut=1.0;
double sip3dcut=4.0;
float lep_ptcut=5.0;
float vtxcut =0.0;
//double isoCutEl=999999.;
double isoCutEl=0.35;
double isoCutMu=0.35;
double leadingPtCut=20.0;
double subleadingPtCut=10.0;
double BTagCut=0.8484;
int job=-1;
int njobs=-1;

#endif
