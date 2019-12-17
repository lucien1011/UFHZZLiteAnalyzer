#include "TTree.h"

using namespace std;  

bool passedFullSelection, passedZ4lSelection;
bool passedZ1LSelection;
bool passedZXCRSelection, passedZ4lZXCRSelection;
bool passSmartCut;
int nZXCRFailedLeptons;
int finalState;
//
bool passed_4l_reco,passed_properLep_ID,passed_LepTightID_Z1,passed_Lep_Leading_subleading;
bool passed_Lep_overlaping, passed_QCD_cut, passed_smartcut, passed_LepIso_Z2, passed_LepTightID_Z2;
bool passed_mZ1_mZ2,passed_Lep_OSSF; 
bool passed_LepIso_Z1;
int Event_NLeptons;
//
bool passTrig;
float pTL1, etaL1, IsoL1, massL1;
float pTL2, etaL2, IsoL2, massL2;
float pTL3, etaL3, IsoL3, massL3;
float pTL4, etaL4;
float phiL1, deltaphiL13;
float phiL2, deltaphiL14;
float phiL3, deltaphiL23;
float phiL4, deltaphiL24;
float deltaphiZZ;
int idL1, idL2, idL3, idL4, MomIdL1, MomIdL2, MomIdL3, PDG_IdL1, PDG_IdL2, PDG_IdL3, MomMomIdL1, MomMomIdL2, MomMomIdL3;

float mass4l, mass4lErr;
float mass3l;
float mass4lREFIT, mass4lErrREFIT;
float massZ1REFIT, massZ2REFIT;
float mass4mu, mass4e, mass2e2mu;
float pT3l, pT4l;
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
float met, met_phi;

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

std::vector<float> *lep_matchedR03_PdgId;
std::vector<float> *lep_matchedR03_MomId;
std::vector<float> *lep_matchedR03_MomMomId;

std::vector<float> *jet_mass;
std::vector<float> *jet_pt; std::vector<float> *jet_eta; std::vector<float> *jet_phi;
std::vector<int> *jet_iscleanH4l; std::vector<float> *jet_QGTagger; std::vector<float> *jet_csvv2;

std::vector<int> *fsrPhotons_lepindex;
std::vector<float> *fsrPhotons_pt; std::vector<float> *fsrPhotons_eta; std::vector<float> *fsrPhotons_phi;
std::vector<float> *fsrPhotons_pterr;

void setHZZTree(TTree* tree){
    
    tree->SetBranchStatus("*",0);

    tree->SetBranchStatus("Run",1);
    tree->SetBranchStatus("LumiSect",1);
    tree->SetBranchStatus("Event",1);
    tree->SetBranchStatus("nVtx",1);
    tree->SetBranchStatus("nInt",1);
    tree->SetBranchStatus("genWeight",1);
    tree->SetBranchStatus("crossSection",1);
    tree->SetBranchStatus("pileupWeight",1);
    tree->SetBranchStatus("passedTrig",1);
    tree->SetBranchStatus("triggersPassed",1);
    tree->SetBranchStatus("passedFiducialSelection",1);
    tree->SetBranchStatus("lep_id",1);
    tree->SetBranchStatus("lep_Sip",1);
    tree->SetBranchStatus("lep_dxy",1);
    tree->SetBranchStatus("lep_dz",1);
    tree->SetBranchStatus("lep_tightId",1);
    tree->SetBranchStatus("lep_ecalDriven",1);
    tree->SetBranchStatus("lep_pt",1);
    tree->SetBranchStatus("lep_pterr",1);
    tree->SetBranchStatus("lep_eta",1);
    tree->SetBranchStatus("lep_phi",1);
    tree->SetBranchStatus("lep_mass",1);
    tree->SetBranchStatus("lep_RelIso",1);
    tree->SetBranchStatus("lep_RelIsoNoFSR",1);
    tree->SetBranchStatus("lepFSR_pt",1);
    tree->SetBranchStatus("lepFSR_eta",1);
    tree->SetBranchStatus("lepFSR_phi",1);
    tree->SetBranchStatus("lepFSR_mass",1);
    tree->SetBranchStatus("lep_matchedR03_PdgId",1);
    tree->SetBranchStatus("lep_matchedR03_MomId",1);
    tree->SetBranchStatus("lep_matchedR03_MomMomId",1);
    tree->SetBranchStatus("jet_pt",1);
    tree->SetBranchStatus("jet_eta",1);
    tree->SetBranchStatus("jet_phi",1);
    tree->SetBranchStatus("jet_mass",1);
    tree->SetBranchStatus("jet_iscleanH4l",1);
    tree->SetBranchStatus("jet_QGTagger",1);
    tree->SetBranchStatus("jet_csvv2",1);
    tree->SetBranchStatus("fsrPhotons_pt",1);
    tree->SetBranchStatus("fsrPhotons_pterr",1);
    tree->SetBranchStatus("fsrPhotons_eta",1);
    tree->SetBranchStatus("fsrPhotons_phi",1);
    tree->SetBranchStatus("fsrPhotons_lepindex",1);
    tree->SetBranchStatus("met",1);
    tree->SetBranchStatus("met_phi",1);
    tree->SetBranchStatus("me_qqZZ_MCFM",1);
    tree->SetBranchAddress("Run",&Run);
    tree->SetBranchAddress("LumiSect",&LumiSect);
    tree->SetBranchAddress("Event",&Event);
    tree->SetBranchAddress("nVtx",&nVtx);
    tree->SetBranchAddress("nInt",&nInt);
    tree->SetBranchAddress("genWeight",&genWeight);
    tree->SetBranchAddress("crossSection",&crossSection);
    tree->SetBranchAddress("pileupWeight",&pileupWeight);
    tree->SetBranchAddress("passedTrig",&passedTrig);
    tree->SetBranchAddress("triggersPassed",&triggersPassed);       
    tree->SetBranchAddress("lep_tightId", &lep_tightId);
    tree->SetBranchAddress("lep_ecalDriven", &lep_ecalDriven);
    tree->SetBranchAddress("passedFiducialSelection",&passedFiducialSelection);
    tree->SetBranchAddress("lep_id", &lep_id);
    tree->SetBranchAddress("lep_Sip",&lep_Sip);
    tree->SetBranchAddress("lep_dxy",&lep_dxy);  
    tree->SetBranchAddress("lep_dz",&lep_dz);
    tree->SetBranchAddress("lep_pt",&lep_pt);
    tree->SetBranchAddress("lep_pterr",&lep_pterr);
    tree->SetBranchAddress("lep_eta",&lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass);
    tree->SetBranchAddress("lep_RelIso",&lep_RelIso);
    tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
    tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt);
    tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta);
    tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi);
    tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass);
    tree->SetBranchAddress("lep_matchedR03_PdgId",&lep_matchedR03_PdgId);
    tree->SetBranchAddress("lep_matchedR03_MomId",&lep_matchedR03_MomId);
    tree->SetBranchAddress("lep_matchedR03_MomMomId",&lep_matchedR03_MomMomId);
    tree->SetBranchAddress("jet_pt",&jet_pt);
    tree->SetBranchAddress("jet_eta",&jet_eta);
    tree->SetBranchAddress("jet_phi",&jet_phi);
    tree->SetBranchAddress("jet_mass",&jet_mass);
    tree->SetBranchAddress("jet_iscleanH4l",&jet_iscleanH4l);
    tree->SetBranchAddress("jet_QGTagger",&jet_QGTagger);
    tree->SetBranchAddress("jet_csvv2",&jet_csvv2);
    tree->SetBranchAddress("fsrPhotons_pt",&fsrPhotons_pt);
    tree->SetBranchAddress("fsrPhotons_pterr",&fsrPhotons_pterr);
    tree->SetBranchAddress("fsrPhotons_eta",&fsrPhotons_eta);
    tree->SetBranchAddress("fsrPhotons_phi",&fsrPhotons_phi);
    tree->SetBranchAddress("fsrPhotons_lepindex",&fsrPhotons_lepindex);
    tree->SetBranchAddress("met",&met);
    tree->SetBranchAddress("met_phi",&met_phi);
    tree->SetBranchAddress("me_qqZZ_MCFM",&me_qqZZ_MCFM);
    // Event Selection

    tree->SetBranchStatus("lep_Hindex",1);
    tree->SetBranchStatus("passedZ4lSelection",1);
    tree->SetBranchStatus("passedZ1LSelection",1);
    tree->SetBranchStatus("passedFullSelection",1);
    tree->SetBranchStatus("passedZXCRSelection",1);
    tree->SetBranchStatus("nZXCRFailedLeptons",1);
    tree->SetBranchStatus("finalState",1);
    tree->SetBranchStatus("dataMCWeight",1);
    tree->SetBranchStatus("k_qqZZ_qcd_M",1);
    tree->SetBranchStatus("k_qqZZ_ewk",1);
    tree->SetBranchStatus("k_ggZZ",1);
    tree->SetBranchStatus("D_bkg_kin",1);
    tree->SetBranchStatus("njets_pt30_eta4p7",1);
    tree->SetBranchStatus("njets_pt30_eta2p5",1);
    tree->SetBranchStatus("nbjets_pt30_eta4p7",1);
    tree->SetBranchStatus("EventCat",1); 
   
    tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
    tree->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
    tree->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
    tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
    tree->SetBranchAddress("passedZXCRSelection",&passedZXCRSelection);
    tree->SetBranchAddress("nZXCRFailedLeptons",&nZXCRFailedLeptons);
    tree->SetBranchAddress("finalState",&finalState);
    tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
    tree->SetBranchAddress("k_qqZZ_qcd_M",&k_qqZZ_qcd_M);
    tree->SetBranchAddress("k_qqZZ_ewk",&k_qqZZ_ewk);
    tree->SetBranchAddress("k_ggZZ",&k_ggZZ);
    tree->SetBranchAddress("D_bkg_kin",&D_bkg_kin);
    tree->SetBranchAddress("njets_pt30_eta4p7",&njets_pt30_eta4p7);
    tree->SetBranchAddress("njets_pt30_eta2p5",&njets_pt30_eta2p5);
    tree->SetBranchAddress("nbjets_pt30_eta4p7",&nbjets_pt30_eta4p7);
    tree->SetBranchAddress("EventCat",&EventCat);
 
}

void initNewLiteTree(TTree* newtree){

    newtree->Branch("Run",&Run,"Run/l");
    newtree->Branch("Event",&Event,"Event/l");
    newtree->Branch("LumiSect",&LumiSect,"LumiSect/l");
    newtree->Branch("nVtx",&nVtx,"nVtx/I");
    newtree->Branch("passedTrig",&passedTrig,"passedTrig/O");
    newtree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
    newtree->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
    newtree->Branch("passedZXCRSelection",&passedZXCRSelection,"passedZXCRSelection/O");
    newtree->Branch("passSmartCut",&passSmartCut,"passSmartCut/O");
    newtree->Branch("nZXCRFailedLeptons",&nZXCRFailedLeptons,"nZXCRFailedLeptons/I");
    newtree->Branch("finalState",&finalState,"finalState/I");    
    newtree->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/F");
    newtree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
    newtree->Branch("genWeight",&genWeight,"genWeight/F");
    newtree->Branch("sumweight",&sumweight,"sumweight/F");
    newtree->Branch("crossSection",&crossSection,"crossSection/F");
    newtree->Branch("k_qqZZ_qcd_M",&k_qqZZ_qcd_M,"k_qqZZ_qcd_M/F");
    newtree->Branch("k_qqZZ_ewk",&k_qqZZ_ewk,"k_qqZZ_ewk/F");
    newtree->Branch("k_ggZZ",&k_ggZZ,"k_ggZZ/F");
    
    newtree->Branch("lep_id",&lep_id);
    newtree->Branch("lep_pt",&lep_pt);
    newtree->Branch("lep_eta",&lep_eta);
    newtree->Branch("lep_phi",&lep_phi);
    newtree->Branch("lep_mass",&lep_mass);
    newtree->Branch("lep_tightId",&lep_tightId);
    newtree->Branch("lep_RelIso",&lep_RelIso);
    newtree->Branch("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
    newtree->Branch("lep_Hindex",&lep_Hindex_stdvec);

    newtree->Branch("pTL1",&pTL1,"pTL1/F");
    newtree->Branch("pTL2",&pTL2,"pTL2/F");
    newtree->Branch("pTL3",&pTL3,"pTL3/F");
    newtree->Branch("pTL4",&pTL4,"pTL4/F");
    newtree->Branch("idL1",&idL1,"idL1/I");
    newtree->Branch("idL2",&idL2,"idL2/I");
    newtree->Branch("idL3",&idL3,"idL3/I");
    newtree->Branch("idL4",&idL4,"idL4/I");
    newtree->Branch("etaL1",&etaL1,"etaL1/F");
    newtree->Branch("etaL2",&etaL2,"etaL2/F");
    newtree->Branch("etaL3",&etaL3,"etaL3/F");
    newtree->Branch("etaL4",&etaL4,"etaL4/F");
    newtree->Branch("phiL1",&phiL1,"phiL1/F");
    newtree->Branch("phiL2",&phiL2,"phiL2/F");
    newtree->Branch("phiL3",&phiL3,"phiL3/F");
    newtree->Branch("phiL4",&phiL4,"phiL4/F");
    newtree->Branch("deltaphiL13",&deltaphiL13,"deltaphiL13/F");
    newtree->Branch("deltaphiL14",&deltaphiL14,"deltaphiL14/F");
    newtree->Branch("deltaphiL23",&deltaphiL23,"deltaphiL23/F");
    newtree->Branch("deltaphiL24",&deltaphiL24,"deltaphiL24/F");
    newtree->Branch("deltaphiZZ",&deltaphiZZ,"deltaphiZZ/F");
    newtree->Branch("mass4l",&mass4l,"mass4l/F");

    newtree->Branch("mass4mu",&mass4mu,"mass4mu/F");
    newtree->Branch("mass4e",&mass4e,"mass4e/F");
    newtree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/F");
    newtree->Branch("pT4l",&pT4l,"pT4l/F");
    newtree->Branch("massZ1",&massZ1,"massZ1/F");
    newtree->Branch("massZ2",&massZ2,"massZ2/F"); 
    newtree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/I");
    newtree->Branch("njets_pt30_eta2p5",&njets_pt30_eta2p5,"njets_pt30_eta2p5/I");
    newtree->Branch("met",&met,"met/F"); 
    newtree->Branch("pTj1",&pTj1,"pTj1/F");
    newtree->Branch("etaj1",&etaj1,"etaj1/F");
    newtree->Branch("qgj1",&qgj1,"qgj1/F");
    newtree->Branch("pTj2",&pTj2,"pTj2/F");
    newtree->Branch("etaj2",&etaj2,"etaj2/F");
    newtree->Branch("qgj2",&qgj2,"qgj1/F");

    newtree->Branch("pt_leadingjet_pt30_eta4p7",&pTj1,"pt_leadingjet_pt30_eta4p7/F");
    newtree->Branch("pt_leadingjet_pt30_eta2p5",&pTj1_2p5,"pt_leadingjet_pt30_eta2p5/F");
   
}

void initNewLiteTree_fakerate(TTree* newtree){

    newtree->Branch("Run",&Run,"Run/l");
    newtree->Branch("Event",&Event,"Event/l");
    newtree->Branch("LumiSect",&LumiSect,"LumiSect/l");
    newtree->Branch("nVtx",&nVtx,"nVtx/I");
    newtree->Branch("nZXCRFailedLeptons",&nZXCRFailedLeptons,"nZXCRFailedLeptons/I");
    newtree->Branch("finalState",&finalState,"finalState/I");    
    newtree->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/F");
    newtree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
    newtree->Branch("genWeight",&genWeight,"genWeight/F");
    newtree->Branch("sumweight",&sumweight,"sumweight/F");
    newtree->Branch("crossSection",&crossSection,"crossSection/F");
    newtree->Branch("k_qqZZ_qcd_M",&k_qqZZ_qcd_M,"k_qqZZ_qcd_M/F");
    newtree->Branch("k_qqZZ_ewk",&k_qqZZ_ewk,"k_qqZZ_ewk/F");
    newtree->Branch("k_ggZZ",&k_ggZZ,"k_ggZZ/F");
    
    newtree->Branch("lep_id",&lep_id);
    newtree->Branch("lep_pt",&lep_pt);
    newtree->Branch("lep_eta",&lep_eta);
    newtree->Branch("lep_phi",&lep_phi);
    newtree->Branch("lep_mass",&lep_mass);
    newtree->Branch("lep_tightId",&lep_tightId);
    newtree->Branch("lep_RelIso",&lep_RelIso);
    newtree->Branch("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
    newtree->Branch("lep_Hindex",&lep_Hindex_stdvec);
    newtree->Branch("lepFSR_pt",&lepFSR_pt);
    newtree->Branch("lepFSR_eta",&lepFSR_eta);
    newtree->Branch("lepFSR_phi",&lepFSR_phi);
    newtree->Branch("lepFSR_mass",&lepFSR_mass);
    newtree->Branch("lep_matchedR03_PdgId",&lep_matchedR03_PdgId);
    newtree->Branch("lep_matchedR03_MomId",&lep_matchedR03_MomId);
    newtree->Branch("lep_matchedR03_MomMomId",&lep_matchedR03_MomMomId);

    newtree->Branch("pTL1",&pTL1,"pTL1/F");
    newtree->Branch("pTL2",&pTL2,"pTL2/F");
    newtree->Branch("pTL3",&pTL3,"pTL3/F");
    newtree->Branch("pTL4",&pTL4,"pTL4/F");
    newtree->Branch("idL1",&idL1,"idL1/I");
    newtree->Branch("idL2",&idL2,"idL2/I");
    newtree->Branch("idL3",&idL3,"idL3/I");
    newtree->Branch("idL4",&idL4,"idL4/I");
    newtree->Branch("etaL1",&etaL1,"etaL1/F");
    newtree->Branch("etaL2",&etaL2,"etaL2/F");
    newtree->Branch("etaL3",&etaL3,"etaL3/F");
    newtree->Branch("etaL4",&etaL4,"etaL4/F");
    newtree->Branch("phiL1",&phiL1,"phiL1/F");
    newtree->Branch("phiL2",&phiL2,"phiL2/F");
    newtree->Branch("phiL3",&phiL3,"phiL3/F");
    newtree->Branch("phiL4",&phiL4,"phiL4/F");
    newtree->Branch("IsoL1",&IsoL1,"IsoL1/F");
    newtree->Branch("IsoL2",&IsoL2,"IsoL2/F");
    newtree->Branch("IsoL3",&IsoL3,"IsoL3/F");
    newtree->Branch("massL1",&massL1,"massL1/F");
    newtree->Branch("massL2",&massL2,"massL2/F");
    newtree->Branch("massL3",&massL3,"massL3/F");
    newtree->Branch("MomIdL1",&MomIdL1,"MomIdL1/I");
    newtree->Branch("MomIdL2",&MomIdL2,"MomIdL2/I");
    newtree->Branch("MomIdL3",&MomIdL3,"MomIdL3/I");
    newtree->Branch("PDG_IdL1",&PDG_IdL1,"PDG_IdL1/I");
    newtree->Branch("PDG_IdL2",&PDG_IdL2,"PDG_IdL2/I");
    newtree->Branch("PDG_IdL3",&PDG_IdL3,"PDG_IdL3/I");
    newtree->Branch("MomMomIdL1",&MomMomIdL1,"MomMomIdL1/I");
    newtree->Branch("MomMomIdL2",&MomMomIdL2,"MomMomIdL2/I");
    newtree->Branch("MomMomIdL3",&MomMomIdL3,"MomMomIdL3/I");

    newtree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/F");
    newtree->Branch("pT3l",&pT3l,"pT3l/F");
    newtree->Branch("massZ1",&massZ1,"massZ1/F");
    newtree->Branch("massZ2",&massZ2,"massZ2/F"); 
    newtree->Branch("met",&met,"met/F"); 
    newtree->Branch("met_phi",&met_phi,"met_phi/F");
    newtree->Branch("pTj1",&pTj1,"pTj1/F");
    newtree->Branch("etaj1",&etaj1,"etaj1/F");
    newtree->Branch("qgj1",&qgj1,"qgj1/F");
    newtree->Branch("pTj2",&pTj2,"pTj2/F");
    newtree->Branch("etaj2",&etaj2,"etaj2/F");
    newtree->Branch("qgj2",&qgj2,"qgj1/F");
    
}
