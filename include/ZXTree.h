#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <map>
#include <utility>
#include <iterator>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"

#include "TSystem.h"
#include "TStyle.h"
#include "TPaveText.h"

#include "TPaveLabel.h"
#include "TLegend.h"

#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <fstream>
#include "TRandom3.h"

#include <algorithm>

#include "TColor.h"

using namespace std;  

////
// output tree
//float me_qqZZ_MCFM;  
bool passedFullSelection, passedZ4lSelection;
bool passedZ1LSelection;
bool passedZXCRSelection, passedZ4lZXCRSelection;
bool passSmartCut;
int nZXCRFailedLeptons;
int nFailedLeptonsZ2;
int finalState;
//
bool passed_4l_reco,passed_properLep_ID,passed_LepTightID_Z1,passed_Lep_Leading_subleading;
bool passed_Lep_overlaping, passed_QCD_cut, passed_smartcut, passed_LepIso_Z2, passed_LepTightID_Z2;
bool passed_mZ1_mZ2,passed_Lep_OSSF; 
bool passed_LepIso_Z1;
int Event_NLeptons;
//
bool passTrig;
float pTL1, etaL1;
float pTL2, etaL2;
float pTL3, etaL3;
float pTL4, etaL4;
float phiL1, deltaphiL13;
float phiL2, deltaphiL14;
float phiL3, deltaphiL23;
float phiL4, deltaphiL24;
float deltaphiZZ;
int idL1, idL2, idL3, idL4;

float mass4l, mass4lErr;
float mass3l;
float mass4lREFIT, mass4lErrREFIT;
float massZ1REFIT, massZ2REFIT;
float mass4mu, mass4e, mass2e2mu;
float pT4l;
float massZ1, massZ1_Z1L, massZ2;
int njets_pt30_eta4p7;
int njets_pt30_eta2p5;
float pTj1, etaj1;
float pTj2, etaj2;
float qgj1, qgj2;
float pTj1_2p5, pTj2_2p5;

float D_bkg_kin;
float D_bkg;
float Dgg10_VAMCFM;
float D_g4;
float D_g1g4;
float D_VBF;
float D_VBF1j;
float D_HadWH;
float D_HadZH;
float D_VBF_QG;
float D_VBF1j_QG;
float D_HadWH_QG;
float D_HadZH_QG;

int EventCat;
int nisoleptons, nbjets_pt30_eta4p7;
float met;

// input tree variables
std::string *triggersPassed;
ULong64_t Run, LumiSect, Event;
bool passedTrig;
bool passedFiducialSelection;
float dataMCWeight, genWeight, pileupWeight, crossSection, sumweight ;
float k_qqZZ_qcd_M,k_qqZZ_ewk,k_ggZZ;
float sumW;
int nVtx, nInt; //nPV
float me_qqZZ_MCFM;
std::vector<float>* lep_mass;
std::vector<float> *lep_pt; std::vector<float> *lep_eta; std::vector<float> *lep_phi;
std::vector<float>* lepFSR_mass;
std::vector<float> *lepFSR_pt; std::vector<float> *lepFSR_eta; std::vector<float> *lepFSR_phi;
std::vector<int> *lep_Hindex_stdvec;

std::vector<int> *lep_tightId;
std::vector<int> *lep_ecalDriven;
std::vector<int> *lep_id;
std::vector<int> *lep_Sip;
std::vector<float> *lep_dxy;
std::vector<float> *lep_dz;
int lep_Hindex[4];
std::vector<float> *lep_RelIso;
std::vector<float> *lep_RelIsoNoFSR;
std::vector<float> *lep_pterr;
std::vector<float> *lep_dataMC;

std::vector<float> *jet_mass;
std::vector<float> *jet_pt; std::vector<float> *jet_eta; std::vector<float> *jet_phi;
std::vector<int> *jet_iscleanH4l; std::vector<float> *jet_QGTagger; std::vector<float> *jet_csvv2;

std::vector<int> *fsrPhotons_lepindex;
std::vector<float> *fsrPhotons_pt; std::vector<float> *fsrPhotons_eta; std::vector<float> *fsrPhotons_phi;
std::vector<float> *fsrPhotons_pterr;
