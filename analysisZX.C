// C++ includes
#include <iostream>
#include <iomanip>
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
#include "setTDRStyle.C"

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
const double CUT_MZ2LOW = 12;
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

// declarations
TString analysisZ1LZXCR(TString processFileName, double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);
int getHistFR(TChain* tree, TTree* TT,
              TH1D* &h1D_FRel_EB,   TH1D* &h1D_FRel_EE,   TH1D* &h1D_FRmu_EB,   TH1D* &h1D_FRmu_EE,
              TH1D* &h1D_FRel_EB_d, TH1D* &h1D_FRel_EE_d, TH1D* &h1D_FRmu_EB_d, TH1D* &h1D_FRmu_EE_d,
              TH1D* &h1D_m4l_2P2F, TH1D* &h1D_m4l_3P1F, TH1D* &h1D_m4l_4P,
              TH1D* &h1D_FRel_mZ,   TH1D* &h1D_FRel_mZ_d, TH1D* &h1D_FRmu_mZ,   TH1D* &h1D_FRmu_mZ_d,
              double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);
double getFR(int lep_id, double lep_pt, double lep_eta, TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE);
void estimateZX(TString processFileName, double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);
int getEstimatesFromCR(TTree* tree,
                  TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE,
                  TH1D* &h1D_m4l_SR_2P2F, TH1D* &h1D_m4l_SR_3P1F,
                  double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);
void analysisInit();
void cmsPreliminary(TCanvas* &c);
void setCavasAndStyles(TString canvasName, TCanvas* &c, TString stat = "", double leftMaring = 0.15, double rightMaring = 0.05, double bottomMaring = 0.15, double topMaring = 0.05);
int setHistProperties(TH1D* &hist, Width_t lineWidth, Style_t lineStyle, Color_t lineColor, Style_t fillStyle=0, Color_t fillColor=0, TString xAxisTitle = "skip", TString yAxisTitle = "skip");
int setHistProperties2D(TH2D* &hist, TString xAxisTitle = "skip", TString yAxisTitle = "skip");
int setLegendProperties(TLegend* &leg, TString sHeader = "skip", Style_t fillStyle=0, Color_t fillColor=0);
int normaliseHist2D(TH2D* &h2D, double norm = 1.);
int fillEmptyBinsHist2D(TH2D* &h2D, double floor);
int normaliseHist1D(TH1D* &h1D, double norm = 1.);
int fillEmptyBinsHist1D(TH1D* &h1D, double floor);

//_______________________________________________________________________________________________________________________________________________
// temporary solution
long int lNEvents = 1;
void setNEvents(TString processFileName){
//    if(processFileName=="Data_Run2015D_ZX.root") lNEvents = 1;
//    if(processFileName=="DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_slimmedZX.root") lNEvents = 21843370;
//    if(processFileName=="DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_slimmedZX.root") lNEvents = 19259106.0;
//    if(processFileName=="TTTo2L2Nu_13TeV-powheg_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_slimmedZX.root") lNEvents = 4997000;
//    if(processFileName=="WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_slimmedZX.root") lNEvents = 8136410;
//    if(processFileName=="ZZTo4L_13TeV_powheg_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_slimmedZX.root") lNEvents = 6652512;
    //    if(processFileName=="DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_slimmedZX.root") lNEvents = 21843370;
    //    if(processFileName=="DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_slimmedZX.root") lNEvents = 19259106.0;
    //    if(processFileName=="TTTo2L2Nu_13TeV-powheg_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_slimmedZX.root") lNEvents = 4997000;
    //    if(processFileName=="WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_slimmedZX.root") lNEvents = 8136410;
    //    if(processFileName=="ZZTo4L_13TeV_powheg_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_slimmedZX.root") lNEvents = 6652512;
}

//_______________________________________________________________________________________________________________________________________________
void analysisZX(TString processFileName){
    analysisInit();

    // temp solution
    setNEvents(processFileName);

    // analyze Z+1L region and compute fake rates
    cout << "[sample: " << processFileName <<"]" << endl;
    cout << "[analyzing Z+1L & Z+OS control regions]" << "["<<CUT_M4LLOW<<" < m4l < "<<CUT_M4LHIGH<<"]" << endl;
    TString slimmedZXFileName = analysisZ1LZXCR(processFileName);
    cout << "[estimating Z+X in signal region]" << "["<<CUT_M4LLOW<<" < m4l < "<<CUT_M4LHIGH<<"]" << endl;
    estimateZX(slimmedZXFileName);
}

TString analysisZ1LZXCR(TString processFileName, double ptElCut, double ptMuCut, double mZ2Cut){

    // define the chain
    TChain* zxTree = 0;
    // prepare the chain and number of MC events in simulation
    if (!(processFileName.Contains("_Run201"))) {
        zxTree = new TChain("Ana/passedEvents");
        TFile* file = new TFile(processFileName,"READ");
        TH1D* hSEvents = (TH1D*) file->Get("Ana/sumWeights");
        lNEvents = (long int)hSEvents->GetBinContent(1);
        file->Close();
        cout  << setw(printOutWidth) << "   MC - weighted number of events: " << lNEvents << endl;
    } else {
        zxTree = new TChain("Ana/passedEvents");
        cout  << setw(printOutWidth) << "   Data - number of events: " << lNEvents << endl;
    }
    zxTree->Add(processFileName);

    // common properties
    Width_t lineWidth = 2;
    double leg_xl = 0.50, leg_xr = 0.90, leg_yb = 0.72, leg_yt = 0.90;

    // define dummy histogram for FRs
    double var_plotHigh = 50.0, var_plotLow = 0.0; int var_nBins = 10;
    TH1D* h1D_dummy = new TH1D("dummy", "dummy", var_nBins, var_plotLow, var_plotHigh);
    TString varAxLabel = "p_{T} [GeV]";
    double binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
    TString sUnit = (varAxLabel.Contains(" [GeV]"))?" GeV":" ";
    TString sBinWidth = TString::Format("%.1f",binWidth) + sUnit;
    setHistProperties(h1D_dummy,1,1,kBlack,0,0,varAxLabel,"Misidentification rate");

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

    // define dummy histogram for CRs
    var_plotHigh = 600.0; var_plotLow = 50.0; var_nBins = 110;
    TH1D* h1D_dummyCR = new TH1D("dummyCR", "dummyCR", var_nBins, var_plotLow, var_plotHigh);
    varAxLabel = "m_{4l} [GeV]";
    binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
    sUnit = (varAxLabel.Contains(" [GeV]"))?" GeV":" ";
    sBinWidth = TString::Format("%.1f",binWidth) + sUnit;
    setHistProperties(h1D_dummyCR,1,1,kBlack,0,0,varAxLabel,"Events/"+sBinWidth);

    // define CR histograms
    TH1D* h1D_m4l_2P2F = new TH1D("h1D_m4l_2P2F","h1D_m4l_2P2F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_2P2F->Sumw2();
    TH1D* h1D_m4l_3P1F = new TH1D("h1D_m4l_3P1F","h1D_m4l_3P1F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_3P1F->Sumw2();
    TH1D* h1D_m4l_4P   = new TH1D("h1D_m4l_4P","h1D_m4l_4P",    var_nBins, var_plotLow, var_plotHigh); h1D_m4l_4P->Sumw2();

    // define dummy histogram for FR vs. mZ
    var_plotHigh = 120.0, var_plotLow = 60.0; var_nBins = 120;
    TH1D* h1D_dummy_mZ = new TH1D("dummy_mZ", "dummy_mZ", var_nBins, var_plotLow, var_plotHigh);
    varAxLabel = "m_{Z} [GeV]";
    binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
    sUnit = (varAxLabel.Contains(" [GeV]"))?" GeV":" ";
    sBinWidth = TString::Format("%.1f",binWidth) + sUnit;
    setHistProperties(h1D_dummy_mZ,1,1,kBlack,0,0,varAxLabel,"Misidentification rate");

    // define FR numerator/denominator histograms
    TH1D* h1D_FRel_mZ   = new TH1D("h1D_FRel_mZ",  "h1D_FRel_mZ",  var_nBins, var_plotLow, var_plotHigh); h1D_FRel_mZ->Sumw2();
    TH1D* h1D_FRel_mZ_d = new TH1D("h1D_FRel_mZ_d","h1D_FRel_mZ_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_mZ_d->Sumw2();
    TH1D* h1D_FRmu_mZ   = new TH1D("h1D_FRmu_mZ",  "h1D_FRmu_mZ",  var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_mZ->Sumw2();
    TH1D* h1D_FRmu_mZ_d = new TH1D("h1D_FRmu_mZ_d","h1D_FRmu_mZ_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_mZ_d->Sumw2();

    // tree for a set of variables, for selected events
    TTree* slimmedTreeZX = new TTree("selectedEvents","selectedEvents");

    // get histograms from tree AnaZX/passedEvents
    getHistFR(zxTree, slimmedTreeZX,
              h1D_FRel_EB,   h1D_FRel_EE,   h1D_FRmu_EB,   h1D_FRmu_EE,
              h1D_FRel_EB_d, h1D_FRel_EE_d, h1D_FRmu_EB_d, h1D_FRmu_EE_d,
              h1D_m4l_2P2F,  h1D_m4l_3P1F,  h1D_m4l_4P,
              h1D_FRel_mZ,   h1D_FRel_mZ_d, h1D_FRmu_mZ,   h1D_FRmu_mZ_d,
              ptElCut, ptMuCut, mZ2Cut);

    // divide hists to get the fake rates
    h1D_FRel_EB->Divide(h1D_FRel_EB_d);
    h1D_FRel_EE->Divide(h1D_FRel_EE_d);
    h1D_FRmu_EB->Divide(h1D_FRmu_EB_d);
    h1D_FRmu_EE->Divide(h1D_FRmu_EE_d);

    // plot the histograms
    // setup environment & canvas
    TCanvas *c1;
    setCavasAndStyles("c1",c1,"");
    TLegend* leg = 0;

    // Fake rates
    // electrons FR
    h1D_dummy_mZ->SetMaximum(4.0*h1D_FRel_mZ->GetMaximum());
    h1D_dummy_mZ->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRel_mZ,lineWidth,1,kRed-7); h1D_FRel_mZ->Draw("same E1 goff"); leg->AddEntry(h1D_FRel_mZ, "electrons, numerator","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_electrons_mZ.pdf");
    h1D_dummy_mZ->SetMaximum(4.0*h1D_FRel_mZ_d->GetMaximum());
    h1D_dummy_mZ->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRel_mZ_d,lineWidth,1,kRed-7); h1D_FRel_mZ_d->Draw("same E1 goff"); leg->AddEntry(h1D_FRel_mZ_d, "electrons, denominator","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_electrons_mZ_d.pdf");
    // muons FR
    h1D_dummy_mZ->SetMaximum(4.0*h1D_FRmu_mZ->GetMaximum());
    h1D_dummy_mZ->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRmu_mZ,lineWidth,1,kRed-7); h1D_FRmu_mZ->Draw("same E1 goff"); leg->AddEntry(h1D_FRmu_mZ, "muons, numerator","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_muons_mZ.pdf");
    h1D_dummy_mZ->SetMaximum(4.0*h1D_FRmu_mZ_d->GetMaximum());
    h1D_dummy_mZ->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRmu_mZ_d,lineWidth,1,kRed-7); h1D_FRmu_mZ_d->Draw("same E1 goff"); leg->AddEntry(h1D_FRmu_mZ_d, "muons, denominator","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_muons_mZ_d.pdf");

    // Fake rates vs. mZ
    // electrons FR
    h1D_dummy->SetMaximum(4.0*h1D_FRel_EE->GetMaximum());
    h1D_dummy->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRel_EB,lineWidth,1,kRed-7); h1D_FRel_EB->Draw("same E1 goff"); leg->AddEntry(h1D_FRel_EB, "electrons, EB","L");
    setHistProperties(h1D_FRel_EE,lineWidth,1,kBlue-7); h1D_FRel_EE->Draw("same E1 goff"); leg->AddEntry(h1D_FRel_EE, "electrons, EE","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_electrons.pdf");
    // muons FR
    h1D_dummy->SetMaximum(4.0*h1D_FRmu_EE->GetMaximum());
    h1D_dummy->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRmu_EB,lineWidth,1,kRed-7); h1D_FRmu_EB->Draw("same E1 goff"); leg->AddEntry(h1D_FRmu_EB, "muons, EB","L");
    setHistProperties(h1D_FRmu_EE,lineWidth,1,kBlue-7); h1D_FRmu_EE->Draw("same E1 goff"); leg->AddEntry(h1D_FRmu_EE, "muons, EE","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_muons.pdf");

    // Control Regions
    // CR 2P2F
    h1D_dummyCR->SetMaximum(2.0*h1D_m4l_2P2F->GetMaximum());
    h1D_dummyCR->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Region 2P+2F");
    setHistProperties(h1D_m4l_2P2F,lineWidth,1,kRed-7); h1D_m4l_2P2F->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_2P2F, "4l","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/CR_2P2F.pdf");
    // CR 3P1F
    h1D_dummyCR->SetMaximum(2.0*h1D_m4l_3P1F->GetMaximum());
    h1D_dummyCR->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Region 3P+1F");
    setHistProperties(h1D_m4l_3P1F,lineWidth,1,kRed-7); h1D_m4l_3P1F->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_3P1F, "4l","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/CR_3P1F.pdf");
    // CR 4P
    h1D_dummyCR->SetMaximum(2.0*h1D_m4l_4P->GetMaximum());
    h1D_dummyCR->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Signal Region");
    setHistProperties(h1D_m4l_4P,lineWidth,1,kRed-7); h1D_m4l_4P->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_4P, "4l","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/CR_4P.pdf");

    // store slimmed tree and histograms in .root file
    // get the file name by stripping off the path
    TObjArray* tad = (TObjArray*) processFileName.Tokenize("/");
    TString fNameExt = ((TObjString*) tad->At(tad->GetEntries() - 1))->GetString();
    if (fNameExt.Contains(".root")){
        TObjArray* taf = (TObjArray*) fNameExt.Tokenize("."); // assumes "." is not used in the base file name
        fNameExt = (((TObjString*) taf->At(taf->GetEntries() - 2))->GetString()).Data();
    }
    TString fileNameTag = /*sPlotsStore+"/"+*/fNameExt+"_slimmedZX.root";
    TString fOption = "RECREATE";
    TFile* fTemplateTree = new TFile(fileNameTag, fOption);
    fTemplateTree->cd();
    slimmedTreeZX->Write();
    h1D_FRel_EB->SetName("h1D_FRel_EB"); h1D_FRel_EB->Write();
    h1D_FRel_EE->SetName("h1D_FRel_EE"); h1D_FRel_EE->Write();
    h1D_FRmu_EB->SetName("h1D_FRmu_EB"); h1D_FRmu_EB->Write();
    h1D_FRmu_EE->SetName("h1D_FRmu_EE"); h1D_FRmu_EE->Write();
    h1D_m4l_2P2F->SetName("h1D_m4l_2P2F"); h1D_m4l_2P2F->Write();
    h1D_m4l_3P1F->SetName("h1D_m4l_3P1F"); h1D_m4l_3P1F->Write();

    h1D_FRel_mZ->SetName("h1D_FRel_mZ"); h1D_FRel_mZ->Write();
    h1D_FRel_mZ_d->SetName("h1D_FRel_mZ_d"); h1D_FRel_mZ_d->Write();
    h1D_FRmu_mZ->SetName("h1D_FRmu_mZ"); h1D_FRmu_mZ->Write();
    h1D_FRmu_mZ_d->SetName("h1D_FRmu_mZ_d"); h1D_FRmu_mZ_d->Write();

    fTemplateTree->Close();
    slimmedTreeZX->Delete();
    zxTree->Delete();
    
    return fileNameTag;
}

//_______________________________________________________________________________________________________________________________________________
int getHistFR(TChain* tree, TTree* TT,
              TH1D* &h1D_FRel_EB,   TH1D* &h1D_FRel_EE,   TH1D* &h1D_FRmu_EB,   TH1D* &h1D_FRmu_EE,
              TH1D* &h1D_FRel_EB_d, TH1D* &h1D_FRel_EE_d, TH1D* &h1D_FRmu_EB_d, TH1D* &h1D_FRmu_EE_d,
              TH1D* &h1D_m4l_2P2F,  TH1D* &h1D_m4l_3P1F,  TH1D* &h1D_m4l_4P,
              TH1D* &h1D_FRel_mZ,   TH1D* &h1D_FRel_mZ_d, TH1D* &h1D_FRmu_mZ,   TH1D* &h1D_FRmu_mZ_d,
              double ptElCut, double ptMuCut, double mZ2Cut){

    // define vars and branches
    float mass3l, mass4l, massZ1, massZ2, met, met_phi, mass4lREFIT, mass4lErr, mass4lErrREFIT;
    double D_bkg_kin, Djet_VAJHU, D_bkg;
    float eventWeight, dataMCWeight, crossSection, genWeight, pileupWeight;
    long int Run, LumiSect, Event;
    int finalState, nVtx;
    bool passedFullSelection, passedZ1LSelection, passedZXCRSelection, passedZ4lSelection, passedZ4lZXCRSelection;
//    TClonesArray *lep_p4 = new TClonesArray("TLorentzVector", 10);
//    TClonesArray *met_p4 = new TClonesArray("TLorentzVector", 1);
    vector<float> *lep_pt = 0; TBranch *b_lep_pt = 0;
    vector<float> *lep_eta = 0; TBranch *b_lep_eta = 0;
    vector<float> *lep_phi = 0; TBranch *b_lep_phi = 0;
    vector<float> *lep_mass = 0; TBranch *b_lep_mass = 0;
    int lep_Hindex[4];
    vector<float> *lepFSR_pt = 0; TBranch *b_lepFSR_pt = 0;
    vector<float> *lepFSR_eta = 0; TBranch *b_lepFSR_eta = 0;
    vector<float> *lepFSR_phi = 0; TBranch *b_lepFSR_phi = 0;
    vector<float> *lepFSR_mass = 0; TBranch *b_lepFSR_mass = 0;
    vector<int> *lep_id = 0; TBranch *b_lep_id = 0;
    vector<int> *lep_tightId = 0; TBranch *b_lep_tightId = 0;
    //vector<float> *lep_RelIso = 0; TBranch *b_lep_RelIso = 0;
    vector<float> *lep_RelIsoNoFSR = 0; TBranch *b_lep_RelIsoNoFSR = 0;
//    vector<int> *lep_matchedR03_PdgId = 0; TBranch *b_lep_matchedR03_PdgId = 0;
//    vector<int> *lep_matchedR03_MomId = 0; TBranch *b_lep_matchedR03_MomId = 0;
//    vector<int> *lep_matchedR03_MomMomId = 0; TBranch *b_lep_matchedR03_MomMomId = 0;

    // counters
    int nEvtPassedZ1LSelection = 0;
    int nEvtPassedZXCRSelection = 0;
    int nEvtZ1WithFailedLeptons = 0;
    int nEvt1P3FLeptons = 0;
    int nEvt2P2FLeptons = 0;
    int nEvt3P1FLeptons = 0;
    int nEvt4PLeptons = 0;

    // get branches
    tree->SetBranchAddress("Run",&Run);
    tree->SetBranchAddress("LumiSect",&LumiSect);
    tree->SetBranchAddress("Event",&Event);
    tree->SetBranchAddress("crossSection",&crossSection);
    tree->SetBranchAddress("eventWeight",&eventWeight);
    tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
    tree->SetBranchAddress("pileupWeight",&pileupWeight);
    tree->SetBranchAddress("genWeight",&genWeight);
    tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
    tree->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
    tree->SetBranchAddress("passedZXCRSelection",&passedZXCRSelection);
    tree->SetBranchAddress("passedZ4lZXCRSelection",&passedZ4lZXCRSelection);
    tree->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
    tree->SetBranchAddress("nVtx",&nVtx);
    tree->SetBranchAddress("finalState",&finalState);
    tree->SetBranchAddress("mass4l",&mass4l);
    tree->SetBranchAddress("massZ1",&massZ1);
    tree->SetBranchAddress("massZ2",&massZ2);
    tree->SetBranchAddress("mass3l",&mass3l);
    tree->SetBranchAddress("mass4lREFIT",&mass4lREFIT);
    tree->SetBranchAddress("mass4lErr",&mass4lErr);
    tree->SetBranchAddress("mass4lErrREFIT",&mass4lErrREFIT);
    tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
    //    tree->GetBranch("lep_p4")->SetAutoDelete(kFALSE);
//    tree->SetBranchAddress("lep_p4",&lep_p4);
    tree->SetBranchAddress("lep_pt",&lep_pt,&b_lep_pt);
    tree->SetBranchAddress("lep_eta",&lep_eta,&b_lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi,&b_lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass,&b_lep_mass);
    tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt,&b_lepFSR_pt);
    tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta,&b_lepFSR_eta);
    tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi,&b_lepFSR_phi);
    tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass,&b_lepFSR_mass);
//    tree->SetBranchAddress("met_p4",&met_p4);
    tree->SetBranchAddress("met",&met);
    tree->SetBranchAddress("met_phi",&met_phi);
    tree->SetBranchAddress("lep_id",&lep_id,&b_lep_id);
    tree->SetBranchAddress("lep_tightId",&lep_tightId,&b_lep_tightId);
    tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR,&b_lep_RelIsoNoFSR);
//    tree->SetBranchAddress("lep_matchedR03_PdgId",&lep_matchedR03_PdgId,&b_lep_matchedR03_PdgId);
//    tree->SetBranchAddress("lep_matchedR03_MomId",&lep_matchedR03_MomId,&b_lep_matchedR03_MomId);
//    tree->SetBranchAddress("lep_matchedR03_MomMomId",&lep_matchedR03_MomMomId,&b_lep_matchedR03_MomMomId);
    tree->SetBranchAddress("D_bkg_kin",&D_bkg_kin);
//    tree->SetBranchAddress("Djet_VAJHU",&Djet_VAJHU);
    tree->SetBranchAddress("D_bkg",&D_bkg);

    // prepare the new tree branches
    int nFailedLeptonsZ1, nFailedLeptonsZ2, nFailedLeptons;
    TT->Branch("Run",&Run,"Run/l");
    TT->Branch("LumiSect",&LumiSect,"LumiSect/l");
    TT->Branch("Event",&Event,"Event/l");
    TT->Branch("crossSection",&crossSection,"crossSection/F");
    TT->Branch("eventWeight",&eventWeight,"eventWeight/F");
    TT->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/F");
    TT->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
    TT->Branch("genWeight",&genWeight,"genWeight/F");
    TT->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
    TT->Branch("passedZ1LSelection",&passedZ1LSelection,"passedZ1LSelection/O");
    TT->Branch("passedZXCRSelection",&passedZXCRSelection,"passedZXCRSelection/O");
    TT->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
    TT->Branch("passedZ4lZXCRSelection",&passedZ4lZXCRSelection,"passedZ4lZXCRSelection/O");
    TT->Branch("nVtx",&nVtx,"nVtx/I");
    TT->Branch("finalState",&finalState,"finalState/I");
    TT->Branch("mass4l",&mass4l,"mass4l/F");
    TT->Branch("massZ1",&massZ1,"massZ1/F");
    TT->Branch("massZ2",&massZ2,"massZ2/F");
    TT->Branch("mass3l",&mass3l,"mass3l/F");
    TT->Branch("mass4lREFIT",&mass4lREFIT,"mass4lREFIT/F");
    TT->Branch("mass4lErr",&mass4lErr,"mass4lErr/F");
    TT->Branch("mass4lErrREFIT",&mass4lErrREFIT,"mass4lErrREFIT/F");

//    TT->Branch("lep_p4","TClonesArray", &lep_p4, 128000, 0);
    TT->Branch("lep_pt",&lep_pt);
    TT->Branch("lep_eta",&lep_eta);
    TT->Branch("lep_phi",&lep_phi);
    TT->Branch("lep_mass",&lep_mass);
    TT->Branch("lep_Hindex",&lep_Hindex,"lep_Hindex[4]/I");
    //    TT->Branch("lep_genindex",&lep_genindex);
    TT->Branch("lep_id",&lep_id);
    TT->Branch("lep_tightId",&lep_tightId);
    TT->Branch("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
//    TT->Branch("lep_matchedR03_PdgId",&lep_matchedR03_PdgId);
//    TT->Branch("lep_matchedR03_MomId",&lep_matchedR03_MomId);
//    TT->Branch("lep_matchedR03_MomMomId",&lep_matchedR03_MomMomId);
    TT->Branch("nFailedLeptonsZ1",&nFailedLeptonsZ1,"nFailedLeptonsZ1/I");
    TT->Branch("nFailedLeptonsZ2",&nFailedLeptonsZ2,"nFailedLeptonsZ2/I");
    TT->Branch("nFailedLeptons",&nFailedLeptons,"nFailedLeptons/I");
//    TT->Branch("met_p4","TClonesArray", &met_p4, 128000, 0);
    TT->Branch("met",&met,"met/F");
    TT->Branch("met_phi",&met_phi,"met_phi/F");
    TT->Branch("D_bkg_kin",&D_bkg_kin,"D_bkg_kin/D");
//    TT->Branch("Djet_VAJHU",&Djet_VAJHU,"Djet_VAJHU/D");
    TT->Branch("D_bkg",&D_bkg,"D_bkg/D");

    // fill histograms
    Long64_t nentries = tree->GetEntries();
    Long64_t nentriesOffset = 0;
//    cout << "nentries: " << nentries << endl;
//    nentries = tree->GetEntries("passedZ1LSelection");
//    cout << "nentries: " << nentries << endl;
//    nentries = tree->GetEntries("passedZXCRSelection");
//    cout << "nentries: " << nentries << endl;
//    nentries = 20000;
//    nentriesOffset = nentries/2;
    cout << "nentries: " << nentries << ", nentriesOffset: " << nentriesOffset << endl;

    if(!SILENT) {
        cout << setw(printOutWidth) << "nEntries: " << nentries << setw(printOutWidth) << endl;
        cout << "Event" << ":" << "Run" << ":" << "LumiSect" << ":" << "mass3l" << ":" << "massZ1" << ":" << "weight" << endl;
    }
    // sort tree in "Event" number, ascending
    int *index = new int[nentries];
    if (SORT_EVENTS) { // not sure if it works when nentriesOffset is non zero
        tree->Draw("Event","","goff");
        int nentries_int = static_cast<int>(nentries);
        TMath::Sort(nentries_int, tree->GetV1(), index, false);
    }
    for(int iEvt=0+nentriesOffset; iEvt < nentries; iEvt++){
        if (SORT_EVENTS) {
            tree->GetEntry(index[iEvt]);	//index[iEvt]);}// take sorted entries
        } else {
            tree->GetEntry(iEvt);
        }// take unsroted entries

        // print out
        if(iEvt%5000000==0) cout << "   event: " << iEvt << "/" << nentries << endl;

        // weight
//        float weight = eventWeight*dataMCWeight*crossSection*LUMI_INT/lNEvents;
        float weight = genWeight*dataMCWeight*crossSection/lNEvents; // weigth per 1 \pb of int. lumi

        double pdg_massZ1 = 91.1876;
        if (passedZ1LSelection) {
            //        cout << "   nEvtPassedZ1LSelection: " << nEvtPassedZ1LSelection << "/" << iEvt << endl;
            //        cout << Event << ":" << Run << ":" << LumiSect << ":" << massZ1 << ":" << passedZ1LSelection << ":" << lep_Hindex << ":" << b_lep_id << ":" << lep_p4 << endl;
            //        cout << "iEvt: " << iEvt << ", lep_Hindex[2]: " << lep_Hindex[2] << ", lep_id->size():" << lep_id->size() << ", lep_tightId->size():" << lep_tightId->size() << endl;

            //----------------------
            // temporary fix
            TLorentzVector lepTmp0, lepTmp1, lepTmp2;
            lepTmp0.SetPtEtaPhiM(lepFSR_pt->at(0),lepFSR_eta->at(0),lepFSR_phi->at(0),lepFSR_mass->at(0));
            lepTmp1.SetPtEtaPhiM(lepFSR_pt->at(1),lepFSR_eta->at(1),lepFSR_phi->at(1),lepFSR_mass->at(1));
            lepTmp2.SetPtEtaPhiM(lepFSR_pt->at(2),lepFSR_eta->at(2),lepFSR_phi->at(2),lepFSR_mass->at(2));
            int idTmp0=lep_id->at(0); int idTmp1=lep_id->at(1); int idTmp2=lep_id->at(2);
            double massZ1_01=0; double massZ1_12=0; double massZ1_02=0;
            if (lep_tightId->at(0) && lep_RelIsoNoFSR->at(0)<0.35 && lep_tightId->at(1) && lep_RelIsoNoFSR->at(1)<0.35 && (idTmp0+idTmp1)==0)
                massZ1_01=(lepTmp0+lepTmp1).M();
            if (lep_tightId->at(0) && lep_RelIsoNoFSR->at(0)<0.35 && lep_tightId->at(2) && lep_RelIsoNoFSR->at(2)<0.35 && (idTmp0+idTmp2)==0)
                massZ1_02=(lepTmp0+lepTmp2).M();
            if (lep_tightId->at(1) && lep_RelIsoNoFSR->at(1)<0.35 && lep_tightId->at(2) && lep_RelIsoNoFSR->at(2)<0.35 && (idTmp1+idTmp2)==0)
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
            //----------------------

//            if (lep_tightId->at(0) && lep_RelIsoNoFSR->at(0)<0.35 && lep_tightId->at(1) && lep_RelIsoNoFSR->at(1)<0.35 && lep_tightId->at(2) && lep_RelIsoNoFSR->at(2)<0.35) {
//                cout << "\n  abs(pdg_massZ1-massZ1_01): " << abs(pdg_massZ1-massZ1_01) << "\n  abs(pdg_massZ1-massZ1_02): " << abs(pdg_massZ1-massZ1_02) << "\n  abs(pdg_massZ1-massZ1_12): " << abs(pdg_massZ1-massZ1_12) <<endl;
//                cout << "\n  id/iso 0: " << (lep_tightId->at(0) && lep_RelIsoNoFSR->at(0)<0.35) << "\n  id/iso 1: " << (lep_tightId->at(1) && lep_RelIsoNoFSR->at(1)<0.35) << "\n  id/iso 2: " << (lep_tightId->at(2) && lep_RelIsoNoFSR->at(2)<0.35) <<endl;
//                cout << "\n  massZ1: " << massZ1 << ", massZ1_01: " << massZ1_01 << ", massZ1_12: " << massZ1_12 << ", massZ1_02: " << massZ1_02 <<endl;
//                cout <<   "  lep_Hindex[2]: " << lep_Hindex[2] <<endl;
//            }

            // check conditions on met and mZ1 in order to suprees true leptons
//            if ((met > 25.0) || (abs(pdg_massZ1-massZ1)>7.0))  continue;
            
            nEvtPassedZ1LSelection++;
            // get properties of the probe lepton (3rd)
            int lep_tight = lep_tightId->at(lep_Hindex[2]);
            float lep_iso = lep_RelIsoNoFSR->at(lep_Hindex[2]);
            int idL3 = lep_id->at(lep_Hindex[2]);
//            TLorentzVector *lep = (TLorentzVector*) lep_p4->At(lep_Hindex[2]);
            TLorentzVector lep;
            lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex[2]),lep_eta->at(lep_Hindex[2]),lep_phi->at(lep_Hindex[2]),lep_mass->at(lep_Hindex[2]));
            float pTL3  = lep.Pt();
            float etaL3 = lep.Eta();

            // fill 1D hists
            if ((abs(idL3) == 11) && (fabs(etaL3) < 1.497)) {
                h1D_FRel_EB_d->Fill(pTL3, weight);
                h1D_FRel_mZ_d->Fill(massZ1, weight);
                if (lep_tight && (lep_iso<0.35)) {h1D_FRel_EB->Fill(pTL3, weight);h1D_FRel_mZ->Fill(massZ1, weight);}
            }
            if ((abs(idL3) == 11) && (fabs(etaL3) > 1.497)) {
                h1D_FRel_EE_d->Fill(pTL3, weight);
                h1D_FRel_mZ_d->Fill(massZ1, weight);
                if (lep_tight && (lep_iso<0.35)) {h1D_FRel_EE->Fill(pTL3, weight);h1D_FRel_mZ->Fill(massZ1, weight);}
            }
            if ((abs(idL3) == 13) && (fabs(etaL3) < 1.2)) {
                h1D_FRmu_EB_d->Fill(pTL3, weight);
                h1D_FRmu_mZ_d->Fill(massZ1, weight);
                if (lep_tight && (lep_iso<0.35)) {h1D_FRmu_EB->Fill(pTL3, weight);h1D_FRmu_mZ->Fill(massZ1, weight);}
            }
            if ((abs(idL3) == 13) && (fabs(etaL3) > 1.2)) {
                h1D_FRmu_EE_d->Fill(pTL3, weight);
                h1D_FRmu_mZ_d->Fill(massZ1, weight);
                if (lep_tight && (lep_iso<0.35)) {h1D_FRmu_EE->Fill(pTL3, weight);h1D_FRmu_mZ->Fill(massZ1, weight);}
            }
            TT->Fill();
        } else if (passedZ4lZXCRSelection) { // passedZ4lZXCRSelection is without trigger requirement

            nEvtPassedZXCRSelection++;
            //        cout << "   nEvtPassedZ1LSelection: " << nEvtPassedZ1LSelection << "/" << iEvt << endl;
            //        cout << Event << ":" << Run << ":" << LumiSect << ":" << massZ1 << ":" << passedZ1LSelection << ":" << lep_Hindex << ":" << b_lep_id << ":" << lep_p4 << endl;
            //        cout << "iEvt: " << iEvt << ", lep_Hindex[2]: " << lep_Hindex[2] << ", lep_id->size():" << lep_id->size() << ", lep_tightId->size():" << lep_tightId->size() << endl;

            // get properties of the non-Z1 leptons (3rd, 4th)
            int lep_tight[4];
            float lep_iso[4];
            int idL[4];
            for(unsigned int k = 0; k <= 3; k++) {
                lep_tight[k] = lep_tightId->at(lep_Hindex[k]);
                lep_iso[k]= lep_RelIsoNoFSR->at(lep_Hindex[k]);
                idL[k] = lep_id->at(lep_Hindex[k]);
            }

            // count the failed leptons
            nFailedLeptonsZ1 = !(lep_tight[0] && ((abs(idL[0])==11 && lep_iso[0]<0.35) || (abs(idL[0])==13 && lep_iso[0]<0.35))) +
                                !(lep_tight[1] && ((abs(idL[1])==11 && lep_iso[1]<0.35) || (abs(idL[1])==13 && lep_iso[1]<0.35)));
            nFailedLeptonsZ2 = !(lep_tight[2] && ((abs(idL[2])==11 && lep_iso[2]<0.35) || (abs(idL[2])==13 && lep_iso[2]<0.35))) +
                                !(lep_tight[3] && ((abs(idL[3])==11 && lep_iso[3]<0.35) || (abs(idL[3])==13 && lep_iso[3]<0.35)));
            nFailedLeptons   = nFailedLeptonsZ1 + nFailedLeptonsZ2;

            // fill 1D hists
            if (nFailedLeptons == 0) {
                nEvt4PLeptons++;
                h1D_m4l_4P->Fill(mass4l, weight);
            }
            if (nFailedLeptons == 1) {
                nEvt3P1FLeptons++;
                h1D_m4l_3P1F->Fill(mass4l, weight);
            }
            if (nFailedLeptons == 2){
                nEvt2P2FLeptons++;
                h1D_m4l_2P2F->Fill(mass4l, weight);
            }
            if (nFailedLeptons == 3){
                nEvt1P3FLeptons++;
//                h1D_m4l_2P2F->Fill(mass4l, weight);
            }
            if (nFailedLeptonsZ1) nEvtZ1WithFailedLeptons++;
            TT->Fill();
        }

        if(!SILENT) {
            cout << fixed;
            cout.precision(printOutPrecision);
            cout << Event << ":" << Run << ":" << LumiSect << ":" << mass3l << ":" << massZ1 << ":" << weight << endl;
        }
    }

    //    if(!SILENT) {
    cout  << setw(printOutWidth) << "selected events in Z+1L region: " << nEvtPassedZ1LSelection << endl;
    cout  << setw(printOutWidth) << "selected events in ZX CRs:      " << nEvtPassedZXCRSelection << endl;
    cout  << setw(printOutWidth) << "selected events in 1P3F region: " << nEvt1P3FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 2P2F region: " << nEvt2P2FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 3P1F region: " << nEvt3P1FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 4P region:   " << nEvt4PLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in ZX CRs with failing Z1 leptons: " << nEvtZ1WithFailedLeptons << endl;
    //    }
    
    return 0;
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
int getEstimatesFromCR(TTree* tree,
              TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE,
              TH1D* &h1D_m4l_SR_2P2F, TH1D* &h1D_m4l_SR_3P1F,
              double ptElCut, double ptMuCut, double mZ2Cut){

    // define vars and branches
    float mass4l;
    float eventWeight, dataMCWeight, crossSection;
//    long int Run, LumiSect, Event;
    ULong64_t Run, LumiSect, Event;
    int finalState, nVtx;
    bool passedFullSelection, passedZ1LSelection, passedZXCRSelection, passedZ4lSelection;
//    TClonesArray *lep_p4 = new TClonesArray("TLorentzVector", 10);
    vector<float> *lep_pt = 0; TBranch *b_lep_pt = 0;
    vector<float> *lep_eta = 0; TBranch *b_lep_eta = 0;
    vector<float> *lep_phi = 0; TBranch *b_lep_phi = 0;
    vector<float> *lep_mass = 0; TBranch *b_lep_mass = 0;
    int lep_Hindex[4];
    vector<int> *lep_id = 0; TBranch *b_lep_id = 0;
    vector<int> *lep_tightId = 0; TBranch *b_lep_tightId = 0;
    //vector<float> *lep_RelIso = 0; TBranch *b_lep_RelIso = 0;
    vector<float> *lep_RelIsoNoFSR = 0; TBranch *b_lep_RelIsoNoFSR = 0;
//    vector<int> *lep_matchedR03_PdgId = 0; TBranch *b_lep_matchedR03_PdgId = 0;
//    vector<int> *lep_matchedR03_MomId = 0; TBranch *b_lep_matchedR03_MomId = 0;
//    vector<int> *lep_matchedR03_MomMomId = 0; TBranch *b_lep_matchedR03_MomMomId = 0;

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
    //    tree->SetBranchAddress("mass3l",&mass3l);
    tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
    //    tree->GetBranch("lep_p4")->SetAutoDelete(kFALSE);
//    tree->SetBranchAddress("lep_p4",&lep_p4);
    tree->SetBranchAddress("lep_pt",&lep_pt,&b_lep_pt);
    tree->SetBranchAddress("lep_eta",&lep_eta,&b_lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi,&b_lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass,&b_lep_mass);
    tree->SetBranchAddress("lep_id",&lep_id,&b_lep_id);
    tree->SetBranchAddress("lep_tightId",&lep_tightId,&b_lep_tightId);
    //tree->SetBranchAddress("lep_RelIso",&lep_RelIso,&b_lep_RelIso);
    tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR,&b_lep_RelIsoNoFSR);

    // fill histograms
    Long64_t nentries = tree->GetEntries();
    cout << "nentries: " << nentries << endl;

    // sort tree in "Event" number, ascending
    int *index = new int[nentries];
    if (SORT_EVENTS) {
        tree->Draw("Event","","goff");
        int nentries_int = static_cast<int>(nentries);
        TMath::Sort(nentries_int, tree->GetV1(), index, false);
    }
    for(int iEvt=0; iEvt < nentries; iEvt++){
        if (SORT_EVENTS) {
            tree->GetEntry(index[iEvt]);	//index[iEvt]);}// take sorted entries
        } else {
            tree->GetEntry(iEvt);
        }// take unsroted entries

        // weight
        float weight = eventWeight*dataMCWeight*crossSection*LUMI_INT/lNEvents;

        if (passedZXCRSelection) {
            nEvtPassedZXCRSelection++;
            //        cout << "   nEvtPassedZ1LSelection: " << nEvtPassedZ1LSelection << "/" << iEvt << endl;
            //        cout << Event << ":" << Run << ":" << LumiSect << ":" << massZ1 << ":" << passedZ1LSelection << ":" << lep_Hindex << ":" << b_lep_id << ":" << lep_p4 << endl;
            //        cout << "iEvt: " << iEvt << ", lep_Hindex[2]: " << lep_Hindex[2] << ", lep_id->size():" << lep_id->size() << ", lep_tightId->size():" << lep_tightId->size() << endl;

            // get properties of the non-Z1 leptons (3rd, 4th)
            int lep_tight[4];
            float lep_iso[4];
            int idL[4];
            float pTL[4];
            float etaL[4];
            float phiL[4];
            for(unsigned int k = 0; k <= 3; k++) {
                lep_tight[k] = lep_tightId->at(lep_Hindex[k]);
                lep_iso[k]= lep_RelIsoNoFSR->at(lep_Hindex[k]);
                idL[k] = lep_id->at(lep_Hindex[k]);
//                TLorentzVector *lep = (TLorentzVector*) lep_p4->At(lep_Hindex[k]);
                TLorentzVector lep;
                lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex[k]),lep_eta->at(lep_Hindex[k]),lep_phi->at(lep_Hindex[k]),lep_mass->at(lep_Hindex[k]));
                pTL[k]  = lep.Pt();
                etaL[k] = lep.Eta();
                phiL[k] = lep.Phi();
            }

            // count the failed leptons
//            nFailedLeptonsZ1 = !(lep_tight[0] && ((abs(idL[0])==11 && lep_iso[0]<0.35) || (abs(idL[0])==13 && lep_iso[0]<0.35))) +
//                                   !(lep_tight[1] && ((abs(idL[1])==11 && lep_iso[1]<0.35) || (abs(idL[1])==13 && lep_iso[1]<0.35)));
            nFailedLeptonsZ2 = !(lep_tight[2] && ((abs(idL[2])==11 && lep_iso[2]<0.35) || (abs(idL[2])==13 && lep_iso[2]<0.35))) +
                                   !(lep_tight[3] && ((abs(idL[3])==11 && lep_iso[3]<0.35) || (abs(idL[3])==13 && lep_iso[3]<0.35)));
//            nFailedLeptons   = nFailedLeptonsZ1 + nFailedLeptonsZ2;
//            nFailedIDOnlyLeptonsZ2 = !(lep_tight[2]) + !(lep_tight[3]);

//            float Riso = 0.3;
//            float dR34 = pow( (etaL[2] - etaL[3])*(etaL[2] - etaL[3]) + (phiL[2] - phiL[3])*(phiL[2] - phiL[3]), 0.5);

            // prepare the estimates
            if (nFailedLeptonsZ2 == 1) {
                nEvt3P1FLeptons++;
                float fr3 = getFR(idL[2], pTL[2], etaL[2], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
                float fr4 = getFR(idL[3], pTL[3], etaL[3], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
                float fr = (!(lep_tight[2] && ((abs(idL[2])==11 && lep_iso[2]<0.35) || (abs(idL[2])==13 && lep_iso[2]<0.35))))*(fr3/(1-fr3)) +
                            (!(lep_tight[3] && ((abs(idL[3])==11 && lep_iso[3]<0.35) || (abs(idL[3])==13 && lep_iso[3]<0.35))))*(fr4/(1-fr4));
//                if (dR34<2*Riso) fr = 0;
                h1D_m4l_SR_3P1F->Fill(mass4l, weight * fr);
            }
            if (nFailedLeptonsZ2 == 2){
                nEvt2P2FLeptons++;
                float fr3 = getFR(idL[2], pTL[2], etaL[2], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
                float fr4 = getFR(idL[3], pTL[3], etaL[3], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
                float fr = (fr3/(1-fr3)) * (fr4/(1-fr4));
//                if (dR34<2*Riso) {
//                    nEvt2P2FLeptonsIsoID_lt06 += 2;
//                    nEvt2P2FLeptonsIDonly_lt06 += nFailedIDOnlyLeptonsZ2;
//                    float frN = fr3    *fr4    *(dR34/(2*Riso))*(dR34/(2*Riso)) + (pow(fr3    *fr4    ,0.5))*(1 - (dR34/(2*Riso)));
//                    float frD = (1-fr3)*(1-fr4)*(dR34/(2*Riso))*(dR34/(2*Riso)) + (pow((1-fr3)*(1-fr4),0.5))*(1 - (dR34/(2*Riso)));
//                    fr = frN/frD;
//                }
                h1D_m4l_SR_2P2F->Fill(mass4l, weight * fr);
            }
        }
    }

    cout  << setw(printOutWidth) << "selected events in ZX CRs:      " << nEvtPassedZXCRSelection << endl;
    cout  << setw(printOutWidth) << "selected events in 2P2F region: " << nEvt2P2FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 3P1F region: " << nEvt3P1FLeptons << endl;

    cout  << setw(printOutWidth) << "int. contributions from 2P2F region: " << h1D_m4l_SR_2P2F->Integral() << endl;
    cout  << setw(printOutWidth) << "int. contributions from 3P1F region: " << h1D_m4l_SR_3P1F->Integral() << endl;

    return 0;
}

//_______________________________________________________________________________________________________________________________________________
void estimateZX(TString slimmedZXFileName, double ptElCut, double ptMuCut, double mZ2Cut){

    // common properties
    Width_t lineWidth = 2;
    double leg_xl = 0.50, leg_xr = 0.90, leg_yb = 0.72, leg_yt = 0.90;

    // get the FR histograms and slimmed ZX tree
    TFile *f = TFile::Open(slimmedZXFileName);
    TTree* zxTree = (TTree*) f->Get("selectedEvents");
    TH1D* h1D_FRel_EB = (TH1D*) f->Get("h1D_FRel_EB");
    TH1D* h1D_FRel_EE = (TH1D*) f->Get("h1D_FRel_EE");
    TH1D* h1D_FRmu_EB = (TH1D*) f->Get("h1D_FRmu_EB");
    TH1D* h1D_FRmu_EE = (TH1D*) f->Get("h1D_FRmu_EE");

    // define dummy histogram for CRs
    double var_plotHigh = 600.0; double var_plotLow = 50.0; double var_nBins = 110;
    TH1D* h1D_dummyCR = new TH1D("dummyCR", "dummyCR", var_nBins, var_plotLow, var_plotHigh);
    TString varAxLabel = "m_{4l} [GeV]";
    double binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
    TString sUnit = (varAxLabel.Contains(" [GeV]"))?" GeV":" ";
    TString sBinWidth = TString::Format("%.1f",binWidth) + sUnit;
    setHistProperties(h1D_dummyCR,1,1,kBlack,0,0,varAxLabel,"Events/"+sBinWidth);

    // define CR histograms
    TH1D* h1D_m4l_SR_2P2F = new TH1D("h1D_m4l_SR_2P2F","h1D_m4l_SR_2P2F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_2P2F->Sumw2();
    TH1D* h1D_m4l_SR_3P1F = new TH1D("h1D_m4l_SR_3P1F","h1D_m4l_SR_3P1F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_3P1F->Sumw2();

    // get histograms from tree AnaZX/passedEvents
    getEstimatesFromCR(zxTree, h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE, h1D_m4l_SR_2P2F, h1D_m4l_SR_3P1F, ptElCut, ptMuCut, mZ2Cut);

    // total contribution (3P1F - 2P2F)
    TH1D* h1D_m4l_SR_tot = (TH1D*) h1D_m4l_SR_3P1F->Clone("h1D_m4l_SR_tot");
    h1D_m4l_SR_tot->Add(h1D_m4l_SR_2P2F, -1);

    // plot the histograms
    // setup environment & canvas
    TCanvas *c1;
    setCavasAndStyles("c1",c1,"");
    TLegend* leg = 0;

    // Estimates from Control Regions
    // SR from 2P2F
    h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_2P2F->GetMaximum());
    h1D_dummyCR->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 2P+2F");
    setHistProperties(h1D_m4l_SR_2P2F,lineWidth,1,kRed-7); h1D_m4l_SR_2P2F->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_2P2F, "4l","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_2P2F.pdf");
    // SR form 3P1F
    h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_3P1F->GetMaximum());
    h1D_dummyCR->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 3P+1F");
    setHistProperties(h1D_m4l_SR_3P1F,lineWidth,1,kRed-7); h1D_m4l_SR_3P1F->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_3P1F, "4l","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_3P1F.pdf");
    // SR total
    h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_tot->GetMaximum());
    h1D_dummyCR->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Total Z+X");
    setHistProperties(h1D_m4l_SR_tot,lineWidth,1,kRed-7); h1D_m4l_SR_tot->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_tot, "4l","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_total.pdf");

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
void getFRDependance(TString slimmedZXFileName){

    enum    lrun                 {r254231,  r256868,  r257613,  r257804,  r258157,  r258159, MAX_LRUN};
    TString runName [MAX_LRUN] = {"254231", "256868", "257613", "257804", "258157", "258159"};
    // equal lumi (~100-120 \pb) run-stop-points: 254231... 256868... 257613... 257804... 258157... 258159

    // get the FR histograms and slimmed ZX tree
    TFile *f = TFile::Open(slimmedZXFileName);
    TTree* zxTree = (TTree*) f->Get("selectedEvents");

    // cuts
    TString tsCutCommonZ1L = "(passedZ1LSelection && (lep_pt[lep_Hindex[2]]<100))";
    TString tsCutEl =        "(abs(lep_id[lep_Hindex[2]])==11)";
    TString tsCutElPassed =  "(lep_tightId[lep_Hindex[2]] && lep_RelIsoNoFSR[lep_Hindex[2]]<0.35)";
    TString tsCutMu =        "(abs(lep_id[lep_Hindex[2]])==13)";
    TString tsCutMuPassed =  "(lep_tightId[lep_Hindex[2]] && lep_RelIsoNoFSR[lep_Hindex[2]]<0.35)";
    TString tsCutEl_d = "(" + tsCutCommonZ1L + " && " + tsCutEl + ")";
    TString tsCutEl_n = "(" + tsCutCommonZ1L + " && " + tsCutEl + " && " + tsCutElPassed + ")";
    TString tsCutMu_d = "(" + tsCutCommonZ1L + " && " + tsCutMu + ")";
    TString tsCutMu_n = "(" + tsCutCommonZ1L + " && " + tsCutMu + " && " + tsCutMuPassed + ")";

    // loop
    for(unsigned int iRun = 1; iRun<MAX_LRUN; iRun++){
        TString cutEl_d = tsCutEl_d + " && Run>" + runName[iRun-1] + " && Run<=" + runName[iRun];
        TString cutEl_n = tsCutEl_n + " && Run>" + runName[iRun-1] + " && Run<=" + runName[iRun];
        TString cutMu_d = tsCutMu_d + " && Run>" + runName[iRun-1] + " && Run<=" + runName[iRun];
        TString cutMu_n = tsCutMu_n + " && Run>" + runName[iRun-1] + " && Run<=" + runName[iRun];

        double nFakesEl_d = zxTree->GetEntries("(" + cutEl_d + ")*eventWeight*dataMCWeight");
        double nFakesEl_n = zxTree->GetEntries("(" + cutEl_n + ")*eventWeight*dataMCWeight");
        double nFakesMu_d = zxTree->GetEntries("(" + cutMu_d + ")*eventWeight*dataMCWeight");
        double nFakesMu_n = zxTree->GetEntries("(" + cutMu_n + ")*eventWeight*dataMCWeight");

//        cout << "(" + cutEl_d + ")*eventWeight*dataMCWeight" << endl;
//        cout << "(" + cutEl_n + ")*eventWeight*dataMCWeight" << endl;
//        cout << "(" + cutMu_d + ")*eventWeight*dataMCWeight" << endl;
//        cout << "(" + cutMu_n + ")*eventWeight*dataMCWeight" << endl;

        cout << "run range [" << runName[iRun-1] << ":" << runName[iRun]  << "], nFakesEl_d: " << nFakesEl_d
            << ", nFakesEl_n: " << nFakesEl_n << ", nFakesMu_d: " << nFakesMu_d << ", nFakesMu_n: " << nFakesMu_n << endl;

        cout << "run range [" << runName[iRun-1] << ":" << runName[iRun]  << "], FR el: " << (double) nFakesEl_n/nFakesEl_d
            << ", Fr mu: " << (double) nFakesMu_n/nFakesMu_d << endl;
    }
}

//_______________________________________________________________________________________________________________________________________________
void cmsPreliminary(TCanvas* &c){
    c->cd();

    TLatex *CMSPrelim = new TLatex();
    CMSPrelim->SetNDC(kTRUE);

    CMSPrelim->SetTextSize(0.5*c->GetTopMargin());
    CMSPrelim->SetTextFont(42);
    CMSPrelim->SetTextAlign(31); // align right
    CMSPrelim->DrawLatex(0.93, 0.96,TString::Format("%.1f",LUMI_INT/1000) + " fb^{-1} at #sqrt{s} = 13 TeV");

    //    CMSPrelim->SetTextSize(0.9*c->GetTopMargin());
    //    CMSPrelim->SetTextFont(62);
    //    CMSPrelim->SetTextAlign(11); // align right
    //    CMSPrelim->DrawLatex(0.27, 0.85, "CMS");
    //
    //    CMSPrelim->SetTextSize(0.7*c->GetTopMargin());
    //    CMSPrelim->SetTextFont(52);
    //    CMSPrelim->SetTextAlign(11);
    //    CMSPrelim->DrawLatex(0.25, 0.8, "Preliminary");
}

//_______________________________________________________________________________________________________________________________________________
void analysisInit() {
    gErrorIgnoreLevel = kWarning;
    gErrorIgnoreLevel = kError;
}

//_______________________________________________________________________________________________________________________________________________
void setCavasAndStyles(TString canvasName, TCanvas* &c, TString stat, double leftMaring, double rightMaring, double bottomMaring, double topMaring){
    // setup environment
//    gROOT->ProcessLine(".L setTDRStyle.C");
    setTDRStyle();
    // setup canvas
    c = new TCanvas(canvasName,"myPlots",0,0,800,600);
    c->cd(1); c->SetLogy(0);
    gStyle->SetOptStat(stat);
    gStyle->SetPalette(1);

    c->GetPad(0)->SetRightMargin(rightMaring);
    c->GetPad(0)->SetLeftMargin(leftMaring);
    c->GetPad(0)->SetTopMargin(topMaring);
    c->GetPad(0)->SetBottomMargin(bottomMaring);
}

//_______________________________________________________________________________________________________________________________________________
int setHistProperties(TH1D* &hist, Width_t lineWidth, Style_t lineStyle, Color_t lineColor, Style_t fillStyle, Color_t fillColor, TString xAxisTitle, TString yAxisTitle){
    if (!hist) return -1;
    // line
    hist->SetLineWidth(lineWidth);
    hist->SetLineStyle(lineStyle);
    hist->SetLineColor(lineColor);
    // fill
    hist->SetFillStyle(fillStyle);
    hist->SetFillColor(fillColor);
    // divisions, offsets, sizes
    hist->GetXaxis()->SetNdivisions(510);
    hist->GetYaxis()->SetNdivisions(510);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    // titles
    if (xAxisTitle!="skip") hist->GetXaxis()->SetTitle(xAxisTitle);
    if (yAxisTitle!="skip") hist->GetYaxis()->SetTitle(yAxisTitle);
    // return
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int setHistProperties2D(TH2D* &hist, TString xAxisTitle, TString yAxisTitle){
    if (!hist) return -1;
    // divisions, offsets, sizes
    hist->GetXaxis()->SetNdivisions(510);
    hist->GetYaxis()->SetNdivisions(510);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    // titles
    if (xAxisTitle!="skip") hist->GetXaxis()->SetTitle(xAxisTitle);
    if (yAxisTitle!="skip") hist->GetYaxis()->SetTitle(yAxisTitle);
    // return
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int setLegendProperties(TLegend* &leg, TString sHeader, Style_t fillStyle, Color_t fillColor){
    // sanity-check
    if (!leg) return -1;
    // titles
    if (sHeader!="skip") leg->SetHeader(sHeader);;
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    // return
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int normaliseHist2D(TH2D* &h2D, double norm){
    if (h2D->Integral()==0) return -1;
    h2D->Scale(norm/h2D->Integral());

    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int fillEmptyBinsHist2D(TH2D* &h2D, double floor) {
    int nXbins=h2D->GetNbinsX();
    int nYbins=h2D->GetNbinsY();
    for(int i=1; i<=nXbins; i++){
        for(int j=1; j<=nYbins; j++){
            h2D->SetBinContent(i,j,h2D->GetBinContent(i,j)+floor);
        }
    }

    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int normaliseHist1D(TH1D* &h1D, double norm){
    if (h1D->Integral()==0) return -1;
    h1D->Scale(norm/h1D->Integral());

    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int fillEmptyBinsHist1D(TH1D* &h1D, double floor) {
    int nXbins=h1D->GetNbinsX();
    for(int i=1; i<=nXbins; i++){
        h1D->SetBinContent(i,h1D->GetBinContent(i)+floor);
    }
    
    return 0;
}

double yields13Comb[3], uncert13Comblow[3], uncert13Combhigh[3];
//_______________________________________________________________________________________________________________________________________________
void plotCombinedResults(TString strOption = "13TeV") {
    //    gROOT->ProcessLine(".L setTDRStyle.C");
    //	setTDRStyle();
    TCanvas *c1 = new TCanvas("c1","myPlots",0,0,700,500);
    c1->cd(1);
    c1->SetLogy(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c1->GetPad(0)->SetRightMargin(0.04);
    c1->GetPad(0)->SetLeftMargin(0.07);
    c1->SetTopMargin(0.06);

    // yields
    double yields13A[3] = {3.1, 2.0, 3.1};
    double yields13AA[3] = {1.7, 2.3, 3.4};
    // errors
    double uncert13A[3] = {0.35*yields13A[0], 0.35*yields13A[1], 0.35*yields13A[2]};
    double uncert13AA[3] = {0.45*yields13AA[0], 0.30*yields13AA[1], 0.35*yields13AA[2]};
    // comb
    yields13Comb[0] = ( yields13A[0]/(uncert13A[0]*uncert13A[0]) +  yields13AA[0]/(uncert13AA[0]*uncert13AA[0]) ) /  ( 1/(uncert13A[0]*uncert13A[0]) +  1/(uncert13AA[0]*uncert13AA[0]) ) ;
    yields13Comb[1] = ( yields13A[1]/(uncert13A[1]*uncert13A[1]) +  yields13AA[1]/(uncert13AA[1]*uncert13AA[1]) ) /  ( 1/(uncert13A[1]*uncert13A[1]) +  1/(uncert13AA[1]*uncert13AA[1]) ) ;
    yields13Comb[2] = ( yields13A[2]/(uncert13A[2]*uncert13A[2]) +  yields13AA[2]/(uncert13AA[2]*uncert13AA[2]) ) /  ( 1/(uncert13A[2]*uncert13A[2]) +  1/(uncert13AA[2]*uncert13AA[2]) ) ;
    uncert13Comblow[0] = ((yields13A[0] - uncert13A[0])<(yields13AA[0] - uncert13AA[0]))?(yields13Comb[0] - yields13A[0] + uncert13A[0]):(yields13Comb[0] - yields13AA[0] + uncert13AA[0]);
    uncert13Comblow[1] = ((yields13A[1] - uncert13A[1])<(yields13AA[1] - uncert13AA[1]))?(yields13Comb[1] - yields13A[1] + uncert13A[1]):(yields13Comb[1] - yields13AA[1] + uncert13AA[1]);
    uncert13Comblow[2] = ((yields13A[2] - uncert13A[2])<(yields13AA[2] - uncert13AA[2]))?(yields13Comb[2] - yields13A[2] + uncert13A[2]):(yields13Comb[2] - yields13AA[2] + uncert13AA[2]);
    uncert13Combhigh[0] = ((yields13A[0] + uncert13A[0])>(yields13AA[0] + uncert13AA[0]))?(yields13A[0] + uncert13A[0] - yields13Comb[0]):(yields13AA[0] + uncert13AA[0] - yields13Comb[0]);
    uncert13Combhigh[1] = ((yields13A[1] + uncert13A[1])>(yields13AA[1] + uncert13AA[1]))?(yields13A[1] + uncert13A[1] - yields13Comb[1]):(yields13AA[1] + uncert13AA[1] - yields13Comb[1]);
    uncert13Combhigh[2] = ((yields13A[2] + uncert13A[2])>(yields13AA[2] + uncert13AA[2]))?(yields13A[2] + uncert13A[2] - yields13Comb[2]):(yields13AA[2] + uncert13AA[2] - yields13Comb[2]);

    // graphs
    double xDummy[3] = {0.80, 2.5, 3.55};
    double errxDummy[3] = {0, 0, 0};
    double x13A[3] = {0.9, 1.9, 2.9};
    double x13AA[3] = {1.1, 2.1, 3.1};
    double x13Comb[3] = {1, 2, 3};
    double errx[3] = {0, 0, 0};
    TGraphErrors* grDummy = new TGraphErrors(6, xDummy, errxDummy);
    TGraphAsymmErrors* gr13A    = new TGraphAsymmErrors(3, x13A, yields13A, errx, errx, uncert13A, uncert13A);
    TGraphAsymmErrors* gr13AA   = new TGraphAsymmErrors(3, x13AA, yields13AA, errx, errx, uncert13AA, uncert13AA);
    TGraphAsymmErrors* gr13Comb = new TGraphAsymmErrors(3, x13Comb, yields13Comb, errx, errx, uncert13Comblow, uncert13Combhigh);

    grDummy->SetTitle("");
    grDummy->GetYaxis()->SetTitle("Expected yield");
    grDummy->GetXaxis()->SetTitle("4e                                  4#mu                                 2e2#mu                         ");
    grDummy->GetXaxis()->SetNdivisions(000);
    grDummy->GetXaxis()->SetTitleOffset(0.8);
    grDummy->GetYaxis()->SetNdivisions(510);
    grDummy->GetYaxis()->SetTitleOffset(0.9);
    grDummy->GetYaxis()->SetRangeUser(0,10);

    // graph properties - 13 TeV
    gr13A->SetMarkerColor(kBlue-2);
    gr13A->SetMarkerStyle(20);
    gr13A->SetMarkerSize(1.3);
    gr13A->SetLineColor(kBlue-2);
    gr13A->SetLineStyle(1);

    gr13AA->SetMarkerColor(kCyan-2);
    gr13AA->SetMarkerStyle(22);
    gr13AA->SetMarkerSize(1.3);
    gr13AA->SetLineColor(kCyan-2);
    gr13AA->SetLineStyle(1);

    gr13Comb->SetMarkerColor(kRed-6);
    gr13Comb->SetMarkerStyle(21);
    gr13Comb->SetMarkerSize(1.5);
    gr13Comb->SetLineColor(kRed+-6);

    grDummy->Draw("AP");
    gr13A->Draw("P same");
    gr13AA->Draw("P same");
    gr13Comb->Draw("P same");

    TLegend* leg = new TLegend(0.50,0.65,0.95,0.90);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetHeader("Reducbile background estimate");
    leg->AddEntry(gr13A, "OS method","pl");
    leg->AddEntry(gr13AA, "SS method","pl");
    leg->AddEntry(gr13Comb, "Combined estimate","pl");
    leg->Draw();

    cmsPreliminary(c1);
    c1->SaveAs(sPlotsStore+"/CombinedPrediction_"+strOption+".pdf");

    cout.precision(3);
    cout << "4e comb:    " << yields13Comb[0] << " + " << uncert13Combhigh[0] << " - " << uncert13Comblow[0] <<  endl;
    cout << "4mu comb:   " << yields13Comb[1] << " + " << uncert13Combhigh[1] << " - " << uncert13Comblow[1] <<  endl;
    cout << "2e2mu comb: " << yields13Comb[2] << " + " << uncert13Combhigh[2] << " - " << uncert13Comblow[2] <<  endl;
}

//_______________________________________________________________________________________________________________________________________________
void plotShapes(TString strOption = "13TeV") {
    //    gROOT->ProcessLine(".L setTDRStyle.C");
    //	setTDRStyle();
    TCanvas *c1 = new TCanvas("c1","myPlots",0,0,700,500);
    c1->cd(1);
    c1->SetLogy(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c1->GetPad(0)->SetRightMargin(0.04);
    c1->GetPad(0)->SetLeftMargin(0.09);
    c1->SetTopMargin(0.06);

    // graphs
    double xDummy[3] = {0, 125, 800};
    double errxDummy[3] = {0, 0, 0};
    TGraphErrors* grDummy = new TGraphErrors(3, xDummy, errxDummy);

    // shapes - 4e
    double rangeX1 = 70;
    double rangeX2 = 800;
    TF1 *f4eAA = new TF1("f4eAA", "[0]*(TMath::Landau(x, [1], [2]))", 70, 800);
    f4eAA->SetParameters(1.4, 157.3, 26.0);f4eAA->SetParameter(0,f4eAA->GetParameter(0)/f4eAA->Integral(rangeX1,rangeX2));
    TF1 *f4eA = new TF1("f4eA","landau( 0 )*(1 + exp( pol1(3)))",50,1000);
    f4eA->SetParameters(0.004265,151.2,36.6,7.06,-0.00497);f4eA->SetParameter(0,f4eA->GetParameter(0)/f4eA->Integral(rangeX1,rangeX2));
    // shapes - 4mu
    TF1 *f4muAA = new TF1("f4muAA", "[0]*(TMath::Landau(x, [1], [2]))", 70, 800);
    f4muAA->SetParameters(2.0, 134.4, 26.6);f4muAA->SetParameter(0,f4muAA->GetParameter(0)/f4muAA->Integral(rangeX1,rangeX2));
    TF1 *f4muA = new TF1("f4muA","landau( 0 )",50,1100);
    f4muA->SetParameters(1.852,134.76,22.23);f4muA->SetParameter(0,f4muA->GetParameter(0)/f4muA->Integral(rangeX1,rangeX2));
    // shapes - 2e2mu
    TF1 *f2e2muAA = new TF1("f2e2muAA", "[0]*(TMath::Landau(x, [1], [2])) + [3]*(TMath::Landau(x, [4], [5]))", 70, 800);
    f2e2muAA->SetParameters(1.7, 134.1, 26.1, 1.6, 155.9, 27.2); double int2e2muAA = f2e2muAA->Integral(rangeX1,rangeX2);
    f2e2muAA->SetParameter(0,f2e2muAA->GetParameter(0)/int2e2muAA); f2e2muAA->SetParameter(3,f2e2muAA->GetParameter(3)/int2e2muAA);
    TF1 *f2e2muA = new TF1("f2e2muA","landau(0)",50,900);
    f2e2muA->SetParameters(2.65,143.7,24.1);f2e2muA->SetParameter(0,f2e2muA->GetParameter(0)/f2e2muA->Integral(rangeX1,rangeX2));
    // shapes - combined
    TF1 *f4eC = new TF1("f4eC", "landau( 0 )*(1 + exp( pol1(3))) + [5]*(TMath::Landau(x, [6], [7]))", 70, 800);
    f4eC->SetParameters(0.004265,151.2,36.6,7.06,-0.00497, 1.4, 157.3, 26.0);double int4eC = f4eC->Integral(rangeX1,rangeX2);
    f4eC->SetParameter(0,f4eC->GetParameter(0)/int4eC); f4eC->SetParameter(5,f4eC->GetParameter(5)/int4eC);
    TF1 *f4muC = new TF1("f4muC","landau( 0 )",70,800);
    f4muC->SetParameters(1.9,134.6,24.4);f4muC->SetParameter(0,f4muC->GetParameter(0)/f4muC->Integral(rangeX1,rangeX2));
    TF1 *f2e2muC = new TF1("f2e2muC","landau(0)",70,800);
    f2e2muC->SetParameters(2.65,144.5,25.3);f2e2muC->SetParameter(0,f2e2muC->GetParameter(0)/f2e2muC->Integral(rangeX1,rangeX2));

    cout.precision(4);
    cout << "f4eC->SetParameters(" << f4eC->GetParameter(0) << ", " << f4eC->GetParameter(1) << ", " << f4eC->GetParameter(2) << ", " << f4eC->GetParameter(3) << ", " << f4eC->GetParameter(4) << ", " << f4eC->GetParameter(5) << ", " << f4eC->GetParameter(6) << ", " << f4eC->GetParameter(7) << ");" << endl;
    cout << "f4muC->SetParameters(" << f4muC->GetParameter(0) << ", " << f4muC->GetParameter(1) << ", " << f4muC->GetParameter(2) << ");" << endl;
    cout << "f2e2muC->SetParameters(" << f2e2muC->GetParameter(0) << ", " << f2e2muC->GetParameter(1) << ", " << f2e2muC->GetParameter(2) << ");" << endl;

//fsum = new TF1(”fsum”,”landau( 0 )*(1 + exp( pol1(3)))+landau(5)+landau(8) ”,50,1100); SetParameters(0.00426549,151.251,36.6082,7.06054,-0.00497032,1.85194,134.765,22.2293,2.65044,143.701,24.071);

    // dummy
    grDummy->SetTitle("");
    grDummy->GetYaxis()->SetTitle("Normalised yield");
    grDummy->GetXaxis()->SetTitle("m_{4l} [GeV]");
    grDummy->GetXaxis()->SetNdivisions(510);
    grDummy->GetXaxis()->SetTitleOffset(1.2);
    grDummy->GetXaxis()->SetRangeUser(rangeX1,rangeX2);
    grDummy->GetYaxis()->SetNdivisions(510);
    grDummy->GetYaxis()->SetTitleOffset(1.3);
    grDummy->GetYaxis()->SetRangeUser(0,0.012);

    // shape properties - 13 TeV
    f4muA->SetLineColor(kBlue-2);
    f4muA->SetLineStyle(2);
    f4muAA->SetLineColor(kBlue-2);
    f4muAA->SetLineStyle(3);

    f4eA->SetLineColor(kGreen-7);
    f4eA->SetLineStyle(2);
    f4eAA->SetLineColor(kGreen-7);
    f4eAA->SetLineStyle(3);

    f2e2muA->SetLineColor(kRed-7);
    f2e2muA->SetLineStyle(2);
    f2e2muAA->SetLineColor(kRed-7);
    f2e2muAA->SetLineStyle(3);

    f4muC->SetLineColor(kBlue-2);
    f4muC->SetLineStyle(1);
    f4eC->SetLineColor(kGreen-7);
    f4eC->SetLineStyle(1);
    f2e2muC->SetLineColor(kRed-7);
    f2e2muC->SetLineStyle(1);

//    // draw
//    grDummy->Draw("A P");
//    f4muA->Draw("P same");
//    f4muAA->Draw("P same");
//    f4eA->Draw("P same");
//    f4eAA->Draw("P same");
//    f2e2muA->Draw("P same");
//    f2e2muAA->Draw("P same");
//
//    f4muC->Draw("P same");
//    f4eC->Draw("P same");
//    f2e2muC->Draw("P same");
//
//    TLegend* leg = new TLegend(0.60,0.50,0.99,0.90);
//    leg->SetFillColor(0);
//    leg->SetFillStyle(0);
//    leg->SetBorderSize(0);
//    leg->SetHeader("Reducibile background shape");
//    leg->AddEntry(f4eA,     "OS method, 4e",    "pl");
//    leg->AddEntry(f4eAA,    "SS method, 4e",    "pl");
//    leg->AddEntry(f4eC,     "Combined,  4e",    "pl");
//    leg->AddEntry(f4muA,    "OS method, 4#mu",  "pl");
//    leg->AddEntry(f4muAA,   "SS method, 4#mu",  "pl");
//    leg->AddEntry(f4muC,    "Combined,  4#mu",  "pl");
//    leg->AddEntry(f2e2muA,  "OS method, 2e2#mu","pl");
//    leg->AddEntry(f2e2muAA, "SS method, 2e2#mu","pl");
//    leg->AddEntry(f2e2muC,  "Combined,  2e2#mu","pl");
//    leg->Draw();
//
//    cmsPreliminary(c1);
//    c1->SaveAs(sPlotsStore+"/Shapes_"+strOption+".pdf");
//
//    grDummy->GetXaxis()->SetRangeUser(100,250);

    // draw
    grDummy->Draw("A P");
    f4muA->Draw("P same");
    f4muAA->Draw("P same");
    f4eA->Draw("P same");
    f4eAA->Draw("P same");
    f2e2muA->Draw("P same");
    f2e2muAA->Draw("P same");

    f4muC->Draw("P same");
    f4eC->Draw("P same");
    f2e2muC->Draw("P same");

    TLegend* leg = new TLegend(0.60,0.50,0.99,0.90);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetHeader("Reducibile background shape");
    leg->AddEntry(f4eA,     "OS method, 4e",    "pl");
    leg->AddEntry(f4eAA,    "SS method, 4e",    "pl");
    leg->AddEntry(f4eC,     "Combined,  4e",    "pl");
    leg->AddEntry(f4muA,    "OS method, 4#mu",  "pl");
    leg->AddEntry(f4muAA,   "SS method, 4#mu",  "pl");
    leg->AddEntry(f4muC,    "Combined,  4#mu",  "pl");
    leg->AddEntry(f2e2muA,  "OS method, 2e2#mu","pl");
    leg->AddEntry(f2e2muAA, "SS method, 2e2#mu","pl");
    leg->AddEntry(f2e2muC,  "Combined,  2e2#mu","pl");
    //    leg->AddEntry(gr13Comb, "Combined estimate","pl");
    leg->Draw();

    cmsPreliminary(c1);
    c1->SaveAs(sPlotsStore+"/Shapes_"+strOption+"_zoomed.pdf");

    cout.precision(3);
    cout << "[m4l > 70 Gev]" << endl;
    cout << "4e comb:    " << f4eC->Integral(70,800)*yields13Comb[0] <<  endl;
    cout << "4mu comb:   " << f4muC->Integral(70,800)*yields13Comb[1] <<  endl;
    cout << "2e2mu comb: " << f2e2muC->Integral(70,800)*yields13Comb[2] <<  endl;

    cout << "[104 < m4l < 140 GeV]" << endl;
    cout << "4e comb:    " << f4eC->Integral(104,140)*yields13Comb[0] <<  endl;
    cout << "4mu comb:   " << f4muC->Integral(104,140)*yields13Comb[1] <<  endl;
    cout << "2e2mu comb: " << f2e2muC->Integral(104,140)*yields13Comb[2] <<  endl;

    cout << "[118 < m4l < 132 GeV]" << endl;
    cout << "4e comb:    " << f4eC->Integral(118,132)*yields13Comb[0] <<  endl;
    cout << "4mu comb:   " << f4muC->Integral(118,132)*yields13Comb[1] <<  endl;
    cout << "2e2mu comb: " << f2e2muC->Integral(118,132)*yields13Comb[2] <<  endl;

}



