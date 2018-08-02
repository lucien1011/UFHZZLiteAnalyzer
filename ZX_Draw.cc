#include "deltaPhi.h"
#include "ZZ4LAnalysisTree.h"
#include "LeptonEfficiency.h"
#include "PileupWeight.h"
#include "Helper.h"
#include "ZXDrawConfig.h"

using namespace std;

TTree* GetTree(TFile* infile);
void DrawZ1LPlot(TTree* tree);

TString filename;
bool debug;
TString mode;

int main(int argc, char *argv[])
{    
     
  debug = false;     

  if(argc == 0) return -1;

  mode = argv[1];

  if (mode=="DrawZ1LPlot") {
    filename = argv[2];
    TString sumWeightFilename = argv[3];
    TString outfilename = argv[4];

    TFile* sumWeightfile = TFile::Open(sumWeightFilename+".root");
    TH1D* sumWeightHist = (TH1D*) sumWeightfile->Get(sumWeightPath);
    sumweight = sumWeightHist->Integral();

    TFile* infile = TFile::Open(filename+".root");
    TTree* tree = GetTree(infile);
    if (!tree) return -1;

    TFile* tmpFile =  new TFile(outfilename+".root","RECREATE");
    tmpFile->cd();
    DrawZ1LPlot(tree);
    tmpFile->Write();
    tmpFile->Close();
  };

}

TTree* GetTree(TFile* infile)
{ 
  TTree* tree;
  tree = (TTree*) infile->Get("Ana/passedEvents");
  if(!tree) tree = (TTree*) infile->Get("passedEvents");
  if(!tree) tree = (TTree*) infile->Get("selectedEvents");
  return tree;
}


void DrawZ1LPlot(TTree* tree){
    
    ZZ4LAnalysisTree::setAddresses(tree, filename);

    int firstevt=0; int lastevt=tree->GetEntries();

    TH1D* h1D_loose_mu_pt = new TH1D("h1D_loose_mu_pt","",10,0.,200.);
    TH1D* h1D_tight_mu_pt = new TH1D("h1D_tight_mu_pt","",10,0.,200.);
    TH1D* h1D_loose_ele_pt = new TH1D("h1D_loose_ele_pt","",10,0.,200.);
    TH1D* h1D_tight_ele_pt = new TH1D("h1D_tight_ele_pt","",10,0.,200.);

    h1D_loose_mu_pt->Sumw2();
    h1D_tight_mu_pt->Sumw2();
    h1D_loose_ele_pt->Sumw2();
    h1D_tight_ele_pt->Sumw2();
   
    for(int evt=0; evt < tree->GetEntries(); evt++) { //event loop
    
        if (evt<firstevt) continue;
        if (evt>lastevt) continue;
       
        if(evt%1000==0) cout<<"Event "<<evt<<"/"<<tree->GetEntries()<<endl;
        tree->GetEntry(evt);

        if ( abs(Zmass-massZ1)<deltaZmass) continue;

        float weight = dataMCWeight*pileupWeight*lumi/sumweight;

        if ( abs( (*lep_id)[lep_Hindex[2]] )==13 && ( !((*lep_tightId)[lep_Hindex[2]]) || (*lep_RelIsoNoFSR)[lep_Hindex[2]]>0.35 ) ) {
            h1D_loose_mu_pt->Fill((*lep_pt)[lep_Hindex[2]],weight);
        } else if ( abs( (*lep_id)[lep_Hindex[2]] )==11 && !((*lep_tightId)[lep_Hindex[2]]) ) {
            h1D_loose_ele_pt->Fill((*lep_pt)[lep_Hindex[2]],weight);
        };
        
        if ( abs( (*lep_id)[lep_Hindex[2]] )==13 && (*lep_tightId)[lep_Hindex[2]] && (*lep_RelIsoNoFSR)[lep_Hindex[2]]<0.35 ) {
            h1D_tight_mu_pt->Fill((*lep_pt)[lep_Hindex[2]],weight);
        } else if ( abs( (*lep_id)[lep_Hindex[2]] )==11  && (*lep_tightId)[lep_Hindex[2]] ) {
            h1D_tight_ele_pt->Fill((*lep_pt)[lep_Hindex[2]],weight);
        };

    }; // Event loop    
    
    TH1D* fr_mu_pt = (TH1D*) h1D_tight_mu_pt->Clone("FakeRate_mu_pt");
    TH1D* fr_ele_pt = (TH1D*) h1D_tight_ele_pt->Clone("FakeRate_ele_pt");
    fr_mu_pt->Divide(h1D_loose_mu_pt);
    fr_ele_pt->Divide(h1D_loose_ele_pt);
}
