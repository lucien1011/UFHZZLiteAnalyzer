#ifndef ZXDrawConfig_h
#define ZXDrawConfig_h

#include "TString.h"

double lumi=41.4;
double Zmass = 91.1876;
double deltaZmass = 7.;
TString sumWeightPath = "Ana/sumWeights";
//TString elFilePath = "/home/lucien/UF-PyNTupleRunner/DarkZ/Data/FakeRate/fakeRates_el_v2.root";
//TString muFilePath = "/home/lucien/UF-PyNTupleRunner/DarkZ/Data/FakeRate/fakeRates_mu_v2.root";
TString elFilePath = "/home/lucien/Higgs/DarkZ/CMSSW_9_4_2/src/liteUFHZZ4LAnalyzer/fakeRate.root";
TString muFilePath = "/home/lucien/Higgs/DarkZ/CMSSW_9_4_2/src/liteUFHZZ4LAnalyzer/fakeRate.root";
int print_per_event=100000;
double isoCutEl=0.35;
double isoCutMu=0.35;

#endif
