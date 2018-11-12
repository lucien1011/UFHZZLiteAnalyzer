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

TString t2_prefix="root://cmsio5.rc.ufl.edu/";
//TString file_path="/store/user/t2/users/klo/Higgs/DarkZ/NTuples/Z1LSkim/Run2016/Data_Run2016_noDuplicates.root";
TString file_path="/store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/skimZ1L_Data_Run2017-17Nov2017_noDuplicates.root";
TString slimmedZXFileName=t2_prefix+file_path;
//TString slimmedZXFileName="/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180820/SkimTree_DarkPhoton_ZX_Run2016Data_v1/Data_Run2016_noDuplicates.root";
//TString slimmedZXFileName="/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180823/SkimTree_Data80X_HIG-16-041-ZXCRSelect//ionWithFlag_v3_liteHZZAna/Data_Run2016_noDuplicates_1.root";
//TString slimmedZXFileName="/raid/raid5/predragm/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Feb21/Data_ZX_Run2017-03Feb2017_slimmedZX.root";
//TString slimmedZXFileName="root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/SingleDoubleMuon_Run2016-03Feb2017.root";
//TString slimmedZXFileName="root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/SingleDoubleMuon_Run2016-03Feb2017.root";

TString outputPath = "Data/fakeRate2017.root";

double isoCutEl=999999.;
//double isoCutEl=0.35;
double isoCutMu=0.35;

void getEstimateFR( TTree* tree, 
                    TH1D* &h1D_FRel_EB,   TH1D* &h1D_FRel_EE,   TH1D* &h1D_FRmu_EB,   TH1D* &h1D_FRmu_EE,
                    TH1D* &h1D_FRel_EB_d, TH1D* &h1D_FRel_EE_d, TH1D* &h1D_FRmu_EB_d, TH1D* &h1D_FRmu_EE_d,
                    double ptElCut, double ptMuCut, double mZ2Cut);
const int SORT_EVENTS = false;

int estimateFR(){
    
    TFile *f = TFile::Open(slimmedZXFileName);
    TTree* zxTree = (TTree*) f->Get("passedEvents");
    //TTree* zxTree = (TTree*) f->Get("selectedEvents");
    //TTree* zxTree = (TTree*) f->Get("Ana/passedEvents");
    
    // define dummy histogram for FRs
    double var_plotHigh = 80.0, var_plotLow = 0.0; int var_nBins = 8;
    TH1D* h1D_dummy = new TH1D("dummy", "dummy", var_nBins, var_plotLow, var_plotHigh);
    TString varAxLabel = "p_{T} [GeV]";
    double binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
    TString sUnit = (varAxLabel.Contains(" [GeV]"))?" GeV":" ";
    TString sBinWidth = TString::Format("%.1f",binWidth) + sUnit;

    // define FR numerator histograms
    TH1D* h1D_FRel_EB   = new TH1D("h1D_FRel_EB","h1D_FRel_EB",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_EB->Sumw2();
    TH1D* h1D_FRel_EE   = new TH1D("h1D_FRel_EE","h1D_FRel_EE",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_EE->Sumw2();
    TH1D* h1D_FRmu_EB   = new TH1D("h1D_FRmu_EB","h1D_FRmu_EB",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_EB->Sumw2();
    TH1D* h1D_FRmu_EE   = new TH1D("h1D_FRmu_EE","h1D_FRmu_EE",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_EE->Sumw2();

    // define FR denominator histograms
    TH1D* h1D_FRel_EB_d   = new TH1D("h1D_FRel_EB_d","h1D_FRel_EB_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_EB_d->Sumw2();
    TH1D* h1D_FRel_EE_d   = new TH1D("h1D_FRel_EE_d","h1D_FRel_EE_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_EE_d->Sumw2();
    TH1D* h1D_FRmu_EB_d   = new TH1D("h1D_FRmu_EB_d","h1D_FRmu_EB_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_EB_d->Sumw2();
    TH1D* h1D_FRmu_EE_d   = new TH1D("h1D_FRmu_EE_d","h1D_FRmu_EE_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_EE_d->Sumw2();

    float ptElCut = 7.;
    float ptMuCut = 5.;
    float mZ2Cut = 12.;

    getEstimateFR(zxTree,
              h1D_FRel_EB,   h1D_FRel_EE,   h1D_FRmu_EB,   h1D_FRmu_EE,
              h1D_FRel_EB_d, h1D_FRel_EE_d, h1D_FRmu_EB_d, h1D_FRmu_EE_d,
              ptElCut, ptMuCut, mZ2Cut);
   
    // divide hists to get the fake rates
    h1D_FRel_EB->Divide(h1D_FRel_EB_d);
    h1D_FRel_EE->Divide(h1D_FRel_EE_d);
    h1D_FRmu_EB->Divide(h1D_FRmu_EB_d);
    h1D_FRmu_EE->Divide(h1D_FRmu_EE_d);

    TString fOption = "RECREATE";
    TFile* fTemplateTree = new TFile(outputPath, fOption);
    fTemplateTree->cd();
    
    h1D_FRel_EB->SetName("h1D_FRel_EB"); h1D_FRel_EB->Write();
    h1D_FRel_EE->SetName("h1D_FRel_EE"); h1D_FRel_EE->Write();
    h1D_FRmu_EB->SetName("h1D_FRmu_EB"); h1D_FRmu_EB->Write();
    h1D_FRmu_EE->SetName("h1D_FRmu_EE"); h1D_FRmu_EE->Write();
    
    fTemplateTree->Close();

    return -1;
}

//_______________________________________________________________________________________________________________________________________________
void getEstimateFR( TTree* tree, 
                    TH1D* &h1D_FRel_EB,   TH1D* &h1D_FRel_EE,   TH1D* &h1D_FRmu_EB,   TH1D* &h1D_FRmu_EE,
                    TH1D* &h1D_FRel_EB_d, TH1D* &h1D_FRel_EE_d, TH1D* &h1D_FRmu_EB_d, TH1D* &h1D_FRmu_EE_d,
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
    vector<float> *lepFSR_pt = 0; TBranch *b_lepFSR_pt = 0;
    vector<float> *lepFSR_eta = 0; TBranch *b_lepFSR_eta = 0;
    vector<float> *lepFSR_phi = 0; TBranch *b_lepFSR_phi = 0;
    vector<float> *lepFSR_mass = 0; TBranch *b_lepFSR_mass = 0;
    int lep_Hindex[4];
    //vector<int>* lep_Hindex = 0;
    vector<int> *lep_id = 0; TBranch *b_lep_id = 0;
    vector<int> *lep_tightId = 0; TBranch *b_lep_tightId = 0;
    vector<float> *lep_RelIsoNoFSR = 0; TBranch *b_lep_RelIsoNoFSR = 0;
    float met = 0.;

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
    tree->SetBranchAddress("met",&met);
    tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
    tree->SetBranchAddress("lep_pt",&lep_pt,&b_lep_pt);
    tree->SetBranchAddress("lep_eta",&lep_eta,&b_lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi,&b_lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass,&b_lep_mass);
    tree->SetBranchAddress("lep_id",&lep_id,&b_lep_id);
    tree->SetBranchAddress("lep_tightId",&lep_tightId,&b_lep_tightId);
    tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR,&b_lep_RelIsoNoFSR);
    tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt,&b_lepFSR_pt);
    tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta,&b_lepFSR_eta);
    tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi,&b_lepFSR_phi);
    tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass,&b_lepFSR_mass);

    Long64_t nentries = tree->GetEntries();
    cout << "nentries: " << nentries << endl;
    
    int nEvtPassedZ1LSelection = 0;

    const double pdg_massZ1 = 91.1876;

    for(int iEvt=0; iEvt < nentries; iEvt++){
        tree->GetEntry(iEvt);

        if(iEvt%100000==0) cout<<"Event "<<iEvt<<"/"<<tree->GetEntries()<<endl;
        
        float weight = 1.;

        if (!passedZ1LSelection) continue;

        TLorentzVector lepTmp0, lepTmp1, lepTmp2;
        lepTmp0.SetPtEtaPhiM(lepFSR_pt->at(0),lepFSR_eta->at(0),lepFSR_phi->at(0),lepFSR_mass->at(0));
        lepTmp1.SetPtEtaPhiM(lepFSR_pt->at(1),lepFSR_eta->at(1),lepFSR_phi->at(1),lepFSR_mass->at(1));
        lepTmp2.SetPtEtaPhiM(lepFSR_pt->at(2),lepFSR_eta->at(2),lepFSR_phi->at(2),lepFSR_mass->at(2));
        int idTmp0=lep_id->at(0); int idTmp1=lep_id->at(1); int idTmp2=lep_id->at(2);
        double massZ1_01=0; double massZ1_12=0; double massZ1_02=0;
        if (lep_tightId->at(0) && ((lep_RelIsoNoFSR->at(0)<isoCutEl && abs(lep_id->at(0))==11) || (lep_RelIsoNoFSR->at(0)<isoCutMu && abs(lep_id->at(0))==13)) && lep_tightId->at(1) && ((lep_RelIsoNoFSR->at(1)<isoCutEl && abs(lep_id->at(1))==11) || (lep_RelIsoNoFSR->at(1)<isoCutMu && abs(lep_id->at(1))==13)) && (idTmp0+idTmp1)==0)
            massZ1_01=(lepTmp0+lepTmp1).M();
        if (lep_tightId->at(0) && ((lep_RelIsoNoFSR->at(0)<isoCutEl && abs(lep_id->at(0))==11) || (lep_RelIsoNoFSR->at(0)<isoCutMu && abs(lep_id->at(0))==13)) && lep_tightId->at(2) && ((lep_RelIsoNoFSR->at(2)<isoCutEl && abs(lep_id->at(2))==11) || (lep_RelIsoNoFSR->at(2)<isoCutMu && abs(lep_id->at(2))==13)) && (idTmp0+idTmp2)==0)
            massZ1_02=(lepTmp0+lepTmp2).M();
        if (lep_tightId->at(1) && ((lep_RelIsoNoFSR->at(1)<isoCutEl && abs(lep_id->at(1))==11) || (lep_RelIsoNoFSR->at(1)<isoCutMu && abs(lep_id->at(1))==13)) && lep_tightId->at(2) && ((lep_RelIsoNoFSR->at(2)<isoCutEl && abs(lep_id->at(2))==11) || (lep_RelIsoNoFSR->at(2)<isoCutMu && abs(lep_id->at(2))==13)) && (idTmp1+idTmp2)==0)
            massZ1_12=(lepTmp1+lepTmp2).M();

        if ((abs(pdg_massZ1-massZ1_01) < TMath::Min(abs(pdg_massZ1-massZ1_02), abs(pdg_massZ1-massZ1_12)))){
            lep_Hindex[2] = 2;
            massZ1 = massZ1_01;
        }else if ((abs(pdg_massZ1-massZ1_02) < abs(pdg_massZ1-massZ1_12))){
            lep_Hindex[2] = 1;
            massZ1 = massZ1_02;
        }else{
            lep_Hindex[2] = 0;
            massZ1 = massZ1_12;
        }

        if (!(massZ1 > 4.)) continue;
        if (!(met < 25.)) continue;
        if (!(abs(massZ1-pdg_massZ1) < 7.)) continue;

        nEvtPassedZ1LSelection++;

        int lep_tight = lep_tightId->at(lep_Hindex[2]);
        float lep_iso = lep_RelIsoNoFSR->at(lep_Hindex[2]);
        int idL3 = lep_id->at(lep_Hindex[2]);
//        TLorentzVector *lep = (TLorentzVector*) lep_p4->At(lep_Hindex[2]);
        TLorentzVector lep;
        lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex[2]),lep_eta->at(lep_Hindex[2]),lep_phi->at(lep_Hindex[2]),lep_mass->at(lep_Hindex[2]));
        float pTL3  = lep.Pt();
        float etaL3 = lep.Eta();

        // fill 1D hists
        if ((abs(idL3) == 11) && (fabs(etaL3) < 1.497)) {
            h1D_FRel_EB_d->Fill(pTL3, weight);
            if (lep_tight && (lep_iso<isoCutEl)) {
                h1D_FRel_EB->Fill(pTL3, weight);
            }
        }
        if ((abs(idL3) == 11) && (fabs(etaL3) > 1.497)) {
            h1D_FRel_EE_d->Fill(pTL3, weight);
            if (lep_tight && (lep_iso<isoCutEl)) {
                h1D_FRel_EE->Fill(pTL3, weight);
            }
        }
        if ((abs(idL3) == 13) && (fabs(etaL3) < 1.2)) {
            h1D_FRmu_EB_d->Fill(pTL3, weight);
            if (lep_tight && (lep_iso<isoCutMu)) {
                h1D_FRmu_EB->Fill(pTL3, weight);
            }
        }
        if ((abs(idL3) == 13) && (fabs(etaL3) > 1.2)) {
            h1D_FRmu_EE_d->Fill(pTL3, weight);
            if (lep_tight && (lep_iso<isoCutMu)) {
                h1D_FRmu_EE->Fill(pTL3, weight);
            }
        }
    }
}
