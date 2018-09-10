#ifndef UpsilonConfig_h
#define UpsilonConfig_h


using namespace std;  

bool isData;
bool bestCandMela=false;
bool redoJets=true;
bool redoEbE=false;

double m4lLowCut=0.0;
double m4lHighCut=9999999.;

double mZ2High=999999.;
double mZ2Low=0.;

double mZ1Low=8.96;
double mZ1High=9.96;
//double mZ1Low=0.;
//double mZ1High=20.;

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
