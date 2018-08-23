#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TChain.h>
#include <TLatex.h>
#include <TMath.h>
#include "TRandom.h"
#include <TFeldmanCousins.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include "Math/VectorUtil.h"
#include "TClonesArray.h"

// other includes
//#include "setTDRStyle.C"

using namespace std;

// Lumi
const double LUMI_INT = 36816; // /pb

// globals
const int DEBUG_LEVEL = 0; // 0 - ERRORs, 1 - ERRORs & WARNINGs, 2 - ERRORs & WARNINGs & INFO
const double CUT_ELPT = 7.;
const double CUT_MUPT = 5.;
const double CUT_MZLOW = 71.; 
const double CUT_MZHIGH = 111.;
const double CUT_MZ1LOW = 40.; 
const double CUT_MZ1HIGH = 120.;
const double CUT_MZ2LOW = 4.;
const double CUT_MZ2HIGH = 120.;
const double CUT_METHIGH = 100.;
const double CUT_M4LLOW = 105.;
const double CUT_M4LHIGH = 140.;
const double CUT_M4LLOW_FULL = 0.; 
const double CUT_M4LHIGH_FULL = 1000.;
const TString sPlotsStore = "plotsZX/";
const int SORT_EVENTS = false;
// print out and debugging settings
const int SILENT = true;
int printOutWidth = 12;
int printOutPrecision = 3;
long int lNEvents = 1;
//TString slimmedZXFileName="/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180820/SkimTree_DarkPhoton_ZX_Run2016Data_v1/Data_Run2016_noDuplicates.root";
TString slimmedZXFileName="/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180823/SkimTree_Data80X_HIG-16-041-ZXCRSelectionWithFlag_v3_liteHZZAna/Data_Run2016_noDuplicates_1.root";
//TString slimmedZXFileName="/raid/raid5/predragm/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Feb21/Data_ZX_Run2017-03Feb2017_slimmedZX.root";

double getFR(int lep_id, double lep_pt, double lep_eta, TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE);
void getEstimateZX(TString slimmedZXFileName, double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);
int getEstimatesFromCR(TTree* tree,
                  TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE,
                  TH1D* &h1D_m4l_SR_2P2F, TH1D* &h1D_m4l_SR_3P1F,
                  double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);

int estimateZX(){
    getEstimateZX(slimmedZXFileName);
    return -1;
}

//_______________________________________________________________________________________________________________________________________________
double getFR(int lep_id, float lep_pt, float lep_eta, TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE){
    if ((abs(lep_id) == 11) && (fabs(lep_eta) < 1.47)) return h1D_FRel_EB->GetBinContent(h1D_FRel_EB->FindBin(lep_pt));
    if ((abs(lep_id) == 11) && (fabs(lep_eta) > 1.47)) return h1D_FRel_EE->GetBinContent(h1D_FRel_EE->FindBin(lep_pt));
    if ((abs(lep_id) == 13) && (fabs(lep_eta) < 1.47)) return h1D_FRmu_EB->GetBinContent(h1D_FRmu_EB->FindBin(lep_pt));
    if ((abs(lep_id) == 13) && (fabs(lep_eta) > 1.47)) return h1D_FRmu_EE->GetBinContent(h1D_FRmu_EE->FindBin(lep_pt));
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
void getEstimateZX(TString slimmedZXFileName, double ptElCut, double ptMuCut, double mZ2Cut){

    // common properties
    Width_t lineWidth = 2;
    double leg_xl = 0.50, leg_xr = 0.90, leg_yb = 0.72, leg_yt = 0.90;

    // get the FR histograms and slimmed ZX tree
    TString elFilePath = "/home/lucien/UF-PyNTupleRunner/DarkZ/Data/FakeRate/fakeRates_el_v2.root";
    TString muFilePath = "/home/lucien/UF-PyNTupleRunner/DarkZ/Data/FakeRate/fakeRates_mu_v2.root";
    //TString elFilePath = "/home/lucien/UF-PyNTupleRunner/DarkZ/Data/FakeRate/fakeRates_el.root";
    //TString muFilePath = "/home/lucien/UF-PyNTupleRunner/DarkZ/Data/FakeRate/fakeRates_mu.root";
       
    TFile* elFile = new TFile(elFilePath,"READ");
    TH1D* h1D_FRel_EB = (TH1D*) elFile->Get("h1D_FRel_EB");
    TH1D* h1D_FRel_EE = (TH1D*) elFile->Get("h1D_FRel_EE");
    
    TFile* muFile = new TFile(muFilePath,"READ");
    TH1D* h1D_FRmu_EB = (TH1D*) muFile->Get("h1D_FRmu_EB");
    TH1D* h1D_FRmu_EE = (TH1D*) muFile->Get("h1D_FRmu_EE");

    TFile *f = TFile::Open(slimmedZXFileName);
    TTree* zxTree = (TTree*) f->Get("passedEvents");
    //TTree* zxTree = (TTree*) f->Get("selectedEvents");

    // define dummy histogram for CRs
    double var_plotHigh = 600.0; double var_plotLow = 50.0; double var_nBins = 110;
    TH1D* h1D_dummyCR = new TH1D("dummyCR", "dummyCR", var_nBins, var_plotLow, var_plotHigh);
    TString varAxLabel = "m_{4l} [GeV]";
    double binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;

    // define CR histograms
    TH1D* h1D_m4l_SR_2P2F = new TH1D("h1D_m4l_SR_2P2F","h1D_m4l_SR_2P2F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_2P2F->Sumw2();
    TH1D* h1D_m4l_SR_3P1F = new TH1D("h1D_m4l_SR_3P1F","h1D_m4l_SR_3P1F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_3P1F->Sumw2();

    // get histograms from tree AnaZX/passedEvents
    getEstimatesFromCR(zxTree, h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE, h1D_m4l_SR_2P2F, h1D_m4l_SR_3P1F, ptElCut, ptMuCut, mZ2Cut);

    // total contribution (3P1F - 2P2F)
    TH1D* h1D_m4l_SR_tot = (TH1D*) h1D_m4l_SR_3P1F->Clone("h1D_m4l_SR_tot");
    h1D_m4l_SR_tot->Add(h1D_m4l_SR_2P2F, -1);

    // store slimmed tree and histograms in .root file
    TString fileNameTag = "estimatesZX";
    TString fOption = "RECREATE";
    TFile* fTemplateTree = new TFile(sPlotsStore+"/"+fileNameTag+".root", fOption);
    fTemplateTree->cd();
    h1D_FRel_EB->SetName("h1D_FRel_EB"); h1D_FRel_EB->Write();
    h1D_FRel_EE->SetName("h1D_FRel_EE"); h1D_FRel_EE->Write();
    h1D_FRmu_EB->SetName("h1D_FRmu_EB"); h1D_FRmu_EB->Write();
    h1D_FRmu_EE->SetName("h1D_FRmu_EE"); h1D_FRmu_EE->Write();
    h1D_m4l_SR_2P2F->SetName("h1D_m4l_SR_2P2F"); h1D_m4l_SR_2P2F->Write();
    h1D_m4l_SR_3P1F->SetName("h1D_m4l_SR_3P1F"); h1D_m4l_SR_3P1F->Write();
    h1D_m4l_SR_tot->SetName("h1D_m4l_SR_tot"); h1D_m4l_SR_tot->Write();
    fTemplateTree->Close();
}


//_______________________________________________________________________________________________________________________________________________
int getEstimatesFromCR(TTree* tree,
              TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE,
              TH1D* &h1D_m4l_SR_2P2F, TH1D* &h1D_m4l_SR_3P1F,
              double ptElCut, double ptMuCut, double mZ2Cut){

    // define vars and branches
    float mass4l,massZ1,massZ2;
    float eventWeight, dataMCWeight, crossSection;
    ULong64_t Run, LumiSect, Event;
    int finalState, nVtx;
    bool passedFullSelection, passedZ1LSelection, passedZXCRSelection, passedZ4lSelection;
    vector<float> *lep_pt = 0; TBranch *b_lep_pt = 0;
    vector<float> *lep_eta = 0; TBranch *b_lep_eta = 0;
    vector<float> *lep_phi = 0; TBranch *b_lep_phi = 0;
    vector<float> *lep_mass = 0; TBranch *b_lep_mass = 0;
    //int lep_Hindex[4];
    vector<int>* lep_Hindex = 0;
    vector<int> *lep_id = 0; TBranch *b_lep_id = 0;
    vector<int> *lep_tightId = 0; TBranch *b_lep_tightId = 0;
    vector<float> *lep_RelIsoNoFSR = 0; TBranch *b_lep_RelIsoNoFSR = 0;

    // counters
    int nEvtPassedZXCRSelection = 0;
    int nEvt2P2FLeptons = 0;
    int nEvt3P1FLeptons = 0;
    int nFailedLeptonsZ2 = 0;

    // get branches
    tree->SetBranchAddress("Run",&Run);
    tree->SetBranchAddress("LumiSect",&LumiSect);
    tree->SetBranchAddress("Event",&Event);
    tree->SetBranchAddress("crossSection",&crossSection);
    tree->SetBranchAddress("eventWeight",&eventWeight);
    tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
    tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
    tree->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
    tree->SetBranchAddress("passedZXCRSelection",&passedZXCRSelection);
    tree->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
    tree->SetBranchAddress("nVtx",&nVtx);
    tree->SetBranchAddress("finalState",&finalState);
    tree->SetBranchAddress("mass4l",&mass4l);
    tree->SetBranchAddress("massZ1",&massZ1);
    tree->SetBranchAddress("massZ2",&massZ2);
    tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
    tree->SetBranchAddress("lep_pt",&lep_pt,&b_lep_pt);
    tree->SetBranchAddress("lep_eta",&lep_eta,&b_lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi,&b_lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass,&b_lep_mass);
    tree->SetBranchAddress("lep_id",&lep_id,&b_lep_id);
    tree->SetBranchAddress("lep_tightId",&lep_tightId,&b_lep_tightId);
    tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR,&b_lep_RelIsoNoFSR);

    Long64_t nentries = tree->GetEntries();
    cout << "nentries: " << nentries << endl;

    int *index = new int[nentries];
    if (SORT_EVENTS) {
        tree->Draw("Event","","goff");
        int nentries_int = static_cast<int>(nentries);
        TMath::Sort(nentries_int, tree->GetV1(), index, false);
    }
    for(int iEvt=0; iEvt < nentries; iEvt++){
        if (SORT_EVENTS) {
            tree->GetEntry(index[iEvt]);
        } else {
            tree->GetEntry(iEvt);
        }

        // weight
        //float weight = eventWeight*dataMCWeight*crossSection*LUMI_INT/lNEvents;
        float weight = 1.;
        
        //if (finalState!=4 && finalState!=3) continue;
        //if (finalState!=1) continue;
        //if (finalState!=2) continue;

        if (passedZXCRSelection) {
            nEvtPassedZXCRSelection++;
            
            int lep_tight[4];
            float lep_iso[4];
            int idL[4];
            float pTL[4];
            float etaL[4];
            float phiL[4];
            for(unsigned int k = 0; k <= 3; k++) {
                //lep_tight[k] = lep_tightId->at(lep_Hindex[k]);
                //lep_iso[k]= lep_RelIsoNoFSR->at(lep_Hindex[k]);
                //idL[k] = lep_id->at(lep_Hindex[k]);

                lep_tight[k] = lep_tightId->at(lep_Hindex->at(k));
                lep_iso[k]= lep_RelIsoNoFSR->at(lep_Hindex->at(k));
                idL[k] = lep_id->at(lep_Hindex->at(k));
                TLorentzVector lep;
                lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex->at(k)),lep_eta->at(lep_Hindex->at(k)),lep_phi->at(lep_Hindex->at(k)),lep_mass->at(lep_Hindex->at(k)));
                //lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex[k]),lep_eta->at(lep_Hindex[k]),lep_phi->at(lep_Hindex[k]),lep_mass->at(lep_Hindex[k]));
                pTL[k]  = lep.Pt();
                etaL[k] = lep.Eta();
                phiL[k] = lep.Phi();
            }

            // count the failed leptons
            nFailedLeptonsZ2 = !(lep_tight[2] && ((abs(idL[2])==11 && lep_iso[2]<0.35) || (abs(idL[2])==13 && lep_iso[2]<0.35))) + !(lep_tight[3] && ((abs(idL[3])==11 && lep_iso[3]<0.35) || (abs(idL[3])==13 && lep_iso[3]<0.35))); 
            
            if (nFailedLeptonsZ2 == 1) {
                nEvt3P1FLeptons++;
                float fr3 = getFR(idL[2], pTL[2], etaL[2], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
                float fr4 = getFR(idL[3], pTL[3], etaL[3], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
                float fr = (!(lep_tight[2] && ((abs(idL[2])==11 && lep_iso[2]<0.35) || (abs(idL[2])==13 && lep_iso[2]<0.35))))*(fr3/(1-fr3)) +
                            (!(lep_tight[3] && ((abs(idL[3])==11 && lep_iso[3]<0.35) || (abs(idL[3])==13 && lep_iso[3]<0.35))))*(fr4/(1-fr4));
                h1D_m4l_SR_3P1F->Fill(mass4l, weight * fr);
            }
            if (nFailedLeptonsZ2 == 2){
                nEvt2P2FLeptons++;
                float fr3 = getFR(idL[2], pTL[2], etaL[2], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
                float fr4 = getFR(idL[3], pTL[3], etaL[3], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
                float fr = (fr3/(1-fr3)) * (fr4/(1-fr4));

                h1D_m4l_SR_2P2F->Fill(mass4l, weight * fr);
            }
        }
    }

    cout  << setw(printOutWidth) << "selected events in ZX CRs:      " << nEvtPassedZXCRSelection << endl;
    cout  << setw(printOutWidth) << "selected events in 2P2F region: " << nEvt2P2FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 3P1F region: " << nEvt3P1FLeptons << endl;

    cout  << setw(printOutWidth) << "int. contributions from 2P2F region: " << h1D_m4l_SR_2P2F->Integral(0,h1D_m4l_SR_2P2F->GetNbinsX()+1) << endl;
    cout  << setw(printOutWidth) << "int. contributions from 3P1F region: " << h1D_m4l_SR_3P1F->Integral(0,h1D_m4l_SR_3P1F->GetNbinsX()+1) << endl;
    cout  << setw(printOutWidth) << "Prediction from ZX CR: " << h1D_m4l_SR_3P1F->Integral(0,h1D_m4l_SR_3P1F->GetNbinsX()+1) - h1D_m4l_SR_2P2F->Integral(0,h1D_m4l_SR_2P2F->GetNbinsX()+1) << endl;

    return 0;
}
