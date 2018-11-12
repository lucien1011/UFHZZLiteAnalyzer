#include "deltaR.h"
#include "ZXTree.h"
#include "Helper.h"
#include "ZXWeightConfig.h"
#include <math.h>

using namespace std;  

TTree* GetTree(TFile* infile);

void MakeFRWeight(
        TTree* tree, TTree* newtree, bool isData, TString filename, 
        TString elFilePath, TString muFilePath, 
        TString inputType
        );

//float getFR(float lepPt, TH1D* hist);
double getFR(int lep_id, double lep_pt, double lep_eta, TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE);
void SetNewTree(TTree* tree);
void SetPedjaTree(TTree* tree);

TString filename;
bool debug;
TString mode;
bool isData;
float FRWeightProd;
float FRWeightSum;
float FRWeightL3;
float FRWeightL4;
float FRWeightProd_UniIso;
float FRWeightProd_AsymIso;

int main(int argc, char *argv[])
{    
     
  debug = false;     

  if(argc == 0) return -1;

  mode = argv[1];

  if (mode=="MakeFRWeight") { 
    filename = argv[2];
    TString sumWeightFilename = argv[3];
    TString outfilename = argv[4];
    isData = atof(argv[5]) > 0;
    TString inputType = argv[6];
    TString elFilePath = argv[7];
    TString muFilePath = argv[8];

    if (!isData) {
        TFile* sumWeightfile = TFile::Open(sumWeightFilename+".root");
        TH1D* sumWeightHist = (TH1D*) sumWeightfile->Get(sumWeightPath);
        sumweight = sumWeightHist->Integral();
    };

    TFile* infile = TFile::Open(filename+".root");
    TFile* tmpFile =  new TFile(outfilename+".root","RECREATE");
    
    TTree* tree = GetTree(infile);
    if (!tree) return -1;

    TTree* newtree = tree->CloneTree(0);
    newtree->CopyAddresses(tree);

    MakeFRWeight(tree,newtree,isData,filename,elFilePath,muFilePath,inputType);

    tmpFile->cd();
    if (inputType=="liteHZZ") {
        newtree->Write("passedEvents",TObject::kOverwrite);
    } else if (inputType=="PedjaSkim") {
        newtree->Write("selectedEvents",TObject::kOverwrite);
    };

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

//float getFR(float lepPt, TH1D* hist)
//{
//    int ibin = hist->GetXaxis()->FindFixBin(lepPt);
//    return hist->GetBinContent(ibin);
//}

double getFR(int lep_id, float lep_pt, float lep_eta, TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE){
    if ((abs(lep_id) == 11) && (fabs(lep_eta) < 1.47)) return h1D_FRel_EB->GetBinContent(h1D_FRel_EB->FindBin(lep_pt));
    if ((abs(lep_id) == 11) && (fabs(lep_eta) > 1.47)) return h1D_FRel_EE->GetBinContent(h1D_FRel_EE->FindBin(lep_pt));
    if ((abs(lep_id) == 13) && (fabs(lep_eta) < 1.47)) return h1D_FRmu_EB->GetBinContent(h1D_FRmu_EB->FindBin(lep_pt));
    if ((abs(lep_id) == 13) && (fabs(lep_eta) > 1.47)) return h1D_FRmu_EE->GetBinContent(h1D_FRmu_EE->FindBin(lep_pt));
    return 0;
}

void MakeFRWeight(
        TTree* tree, 
        TTree* newtree, 
        bool isData, 
        TString filename, 
        TString elFilePath,
        TString muFilePath,
        TString inputType="liteHZZ"
        ){
    
    // ZZ4LAnalysisTree::setAddresses(tree, filename);
    if (inputType=="liteHZZ") {
        SetNewTree(tree);
    } else if (inputType=="PedjaSkim") {
        SetPedjaTree(tree);
    };

    int firstevt=0; int lastevt=tree->GetEntries();

    newtree->Branch("FRWeightProd",&FRWeightProd,"FRWeightProd/F");
    newtree->Branch("FRWeightSum",&FRWeightSum,"FRWeightSum/F");
    newtree->Branch("FRWeightL3",&FRWeightL3,"FRWeightSumCorrL3/F");
    newtree->Branch("FRWeightL4",&FRWeightL4,"FRWeightSumCorrL4/F");
    newtree->Branch("FRWeightProd_UniIso",&FRWeightProd_UniIso,"FRWeightProd_UniIso/F");
    newtree->Branch("FRWeightProd_AsymIso",&FRWeightProd_AsymIso,"FRWeightProd_AsymIso/F");
    newtree->Branch("nFailedLeptonsZ2",&nFailedLeptonsZ2,"nFailedLeptonsZ2/I");
    
    TFile* elFile = new TFile(elFilePath,"READ");
    TH1D* h1D_FRel_EB = (TH1D*) elFile->Get("h1D_FRel_EB");
    TH1D* h1D_FRel_EE = (TH1D*) elFile->Get("h1D_FRel_EE");
    
    TFile* muFile = new TFile(muFilePath,"READ");
    TH1D* h1D_FRmu_EB = (TH1D*) muFile->Get("h1D_FRmu_EB");
    TH1D* h1D_FRmu_EE = (TH1D*) muFile->Get("h1D_FRmu_EE");

    for(int evt=0; evt < tree->GetEntries(); evt++) { //event loop
    
        if (evt<firstevt) continue;
        if (evt>lastevt) continue;
       
        if(evt%print_per_event==0) cout<<"Event "<<evt<<"/"<<tree->GetEntries()<<endl;
        tree->GetEntry(evt);

        if (!passedZXCRSelection) continue;

        //nFailedLeptonsZ2 = 0;
        FRWeightProd=-1.;
        FRWeightSum=-1.;
        // get properties of the non-Z1 leptons (3rd, 4th)
        int lep_tight[4];
        float lep_iso[4];
        int idL[4];
        float pTL[4];
        float etaL[4];
        float phiL[4];
        std::vector<int> index_vec;
        for(unsigned int k = 0; k <= 3; k++) {
            TLorentzVector lep;
            if (inputType=="liteHZZ") {
                lep_tight[k] = lep_tightId->at(lep_Hindex_stdvec->at(k));
                lep_iso[k]= lep_RelIsoNoFSR->at(lep_Hindex_stdvec->at(k));
                idL[k] = lep_id->at(lep_Hindex_stdvec->at(k));
                //TLorentzVector *lep = (TLorentzVector*) lep_p4->At((*lep_Hindex_stdvec)[k]);
                //lep.SetPtEtaPhiM(lep_pt->at((*lep_Hindex_stdvec)[k]),lep_eta->at((*lep_Hindex_stdvec)[k]),lep_phi->at((*lep_Hindex_stdvec)[k]),lep_mass->at((*lep_Hindex_stdvec)[k]));
                lep.SetPtEtaPhiM(
                        lep_pt->at(lep_Hindex_stdvec->at(k)),
                        lep_eta->at(lep_Hindex_stdvec->at(k)),
                        lep_phi->at(lep_Hindex_stdvec->at(k)),
                        lep_mass->at(lep_Hindex_stdvec->at(k))
                        );
            } else if (inputType=="PedjaSkim") {
                lep_tight[k] = lep_tightId->at(lep_Hindex[k]);
                lep_iso[k]= lep_RelIsoNoFSR->at(lep_Hindex[k]);
                idL[k] = lep_id->at(lep_Hindex[k]);
                //TLorentzVector *lep = (TLorentzVector*) lep_p4->At((*lep_Hindex_stdvec)[k]);
                //lep.SetPtEtaPhiM(lep_pt->at((*lep_Hindex_stdvec)[k]),lep_eta->at((*lep_Hindex_stdvec)[k]),lep_phi->at((*lep_Hindex_stdvec)[k]),lep_mass->at((*lep_Hindex_stdvec)[k]));
                lep.SetPtEtaPhiM(
                        lep_pt->at(lep_Hindex[k]),
                        lep_eta->at(lep_Hindex[k]),
                        lep_phi->at(lep_Hindex[k]),
                        lep_mass->at(lep_Hindex[k])
                        );
            };

            pTL[k]  = lep.Pt();
            etaL[k] = lep.Eta();
            phiL[k] = lep.Phi();
            if ( !(abs(idL[k])==11 && lep_tight[k] && lep_iso[k]<isoCutEl) && !(abs(idL[k])==13 && lep_tight[k] && lep_iso[k]<isoCutMu)  ) {
                index_vec.push_back(k);
                //nFailedLeptonsZ2++;
            }
        };

        if (inputType=="liteHZZ") {
            nFailedLeptonsZ2 = !(lep_tight[2] && ((abs(idL[2])==11 && lep_iso[2]<isoCutEl) || (abs(idL[2])==13 && lep_iso[2]<isoCutMu))) + !(lep_tight[3] && ((abs(idL[3])==11 && lep_iso[3]<isoCutEl) || (abs(idL[3])==13 && lep_iso[3]<isoCutMu)));
        };

        if (nFailedLeptonsZ2 == 1) {
            float fr3 = getFR(idL[2], pTL[2], etaL[2], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            float fr4 = getFR(idL[3], pTL[3], etaL[3], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            float fr = (!(lep_tight[2] && ((abs(idL[2])==11 && lep_iso[2]<isoCutEl) || (abs(idL[2])==13 && lep_iso[2]<isoCutMu))))*(fr3/(1-fr3)) +
                        (!(lep_tight[3] && ((abs(idL[3])==11 && lep_iso[3]<isoCutEl) || (abs(idL[3])==13 && lep_iso[3]<isoCutMu))))*(fr4/(1-fr4));
            
            //float fr3 = getFR(idL[index_vec.at(0)], pTL[index_vec.at(0)], etaL[index_vec.at(0)], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            FRWeightProd=fr;
            FRWeightSum=fr;
            FRWeightL3=fr3;
            FRWeightL4=fr4;
            if (lep_iso[2] > lep_iso[3]) {
                FRWeightProd_AsymIso=fr3;
            } else {

                FRWeightProd_AsymIso=fr4;
            }
            FRWeightProd_UniIso=fr;

            newtree->Fill();
        }
        if (nFailedLeptonsZ2 == 2){
            float fr3 = getFR(idL[2], pTL[2], etaL[2], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            float fr4 = getFR(idL[3], pTL[3], etaL[3], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            //float fr3 = getFR(idL[index_vec.at(0)], pTL[index_vec.at(0)], etaL[index_vec.at(0)], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            //float fr4 = getFR(idL[index_vec.at(1)], pTL[index_vec.at(1)], etaL[index_vec.at(1)], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            float fr = (fr3/(1-fr3)) * (fr4/(1-fr4)); 
            FRWeightProd=fr;
            FRWeightSum=fr3/(1-fr3)+fr4/(1-fr4);
            FRWeightL3=fr3/(1-fr3);
            FRWeightL4=fr4/(1-fr4);
            if (lep_iso[2] > lep_iso[3]) {
                FRWeightProd_AsymIso=fr3;
            } else {

                FRWeightProd_AsymIso=fr4;
            }
            float dr2 = deltaR2(etaL[2],phiL[2],etaL[3],phiL[3]);
            float drIso = 0.3;
            FRWeightProd_UniIso = ( fr3*fr4*dr2/4./drIso/drIso + sqrt(fr3*fr4)*(1.-sqrt(dr2)/2./drIso) ) / ( (1.-fr3)*(1.-fr4)*dr2/4./drIso/drIso + sqrt((1.-fr3)*(1.-fr4))*(1.-sqrt(dr2)/2./drIso) );  

            newtree->Fill();
        }

    }; // Event loop

    elFile->Close();
    muFile->Close();
};

void SetPedjaTree(TTree* tree){

   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus("Run",1);
   tree->SetBranchStatus("LumiSect",1);
   tree->SetBranchStatus("Event",1);
   tree->SetBranchStatus("nVtx",1);
   tree->SetBranchStatus("passedFullSelection",1);
   tree->SetBranchStatus("passedZ4lSelection",1);
   tree->SetBranchStatus("passedZXCRSelection",1);
   tree->SetBranchStatus("finalState",1);
   tree->SetBranchStatus("dataMCWeight",1);
   tree->SetBranchStatus("pileupWeight",1);
   tree->SetBranchStatus("genWeight",1);
   tree->SetBranchStatus("sumWeight",1);
   tree->SetBranchStatus("crossSection",1);
   tree->SetBranchStatus("lep_id",1);
   tree->SetBranchStatus("lep_pt",1);
   tree->SetBranchStatus("lep_eta",1);
   tree->SetBranchStatus("lep_phi",1);
   tree->SetBranchStatus("lep_mass",1);
   tree->SetBranchStatus("lep_tightId",1);
   tree->SetBranchStatus("lep_RelIsoNoFSR",1);
   tree->SetBranchStatus("lep_Hindex",1);
   tree->SetBranchStatus("deltaphiZZ",1);
   tree->SetBranchStatus("mass4l",1);
   tree->SetBranchStatus("massZ1",1);
   tree->SetBranchStatus("massZ2",1);
   tree->SetBranchStatus("met",1);
   tree->SetBranchStatus("nFailedLeptonsZ2",1);

   tree->SetBranchAddress("Run", &Run);
   tree->SetBranchAddress("Event", &Event);
   tree->SetBranchAddress("LumiSect", &LumiSect);
   tree->SetBranchAddress("nVtx", &nVtx);
   tree->SetBranchAddress("passedTrig", &passedTrig);
   tree->SetBranchAddress("passedFullSelection", &passedFullSelection);
   tree->SetBranchAddress("passedZ4lSelection", &passedZ4lSelection);
   tree->SetBranchAddress("passedZXCRSelection", &passedZXCRSelection);
   tree->SetBranchAddress("finalState", &finalState);
   tree->SetBranchAddress("dataMCWeight", &dataMCWeight);
   tree->SetBranchAddress("pileupWeight", &pileupWeight);
   tree->SetBranchAddress("genWeight", &genWeight);
   tree->SetBranchAddress("sumweight", &sumweight);
   tree->SetBranchAddress("crossSection", &crossSection);
   tree->SetBranchAddress("lep_id", &lep_id);
   tree->SetBranchAddress("lep_pt", &lep_pt);
   tree->SetBranchAddress("lep_eta", &lep_eta);
   tree->SetBranchAddress("lep_phi", &lep_phi);
   tree->SetBranchAddress("lep_mass", &lep_mass);
   tree->SetBranchAddress("lep_tightId", &lep_tightId);
   tree->SetBranchAddress("lep_RelIsoNoFSR", &lep_RelIsoNoFSR);
   tree->SetBranchAddress("lep_Hindex", &lep_Hindex);
   tree->SetBranchAddress("deltaphiZZ", &deltaphiZZ);
   tree->SetBranchAddress("mass4l", &mass4l);
   tree->SetBranchAddress("massZ1", &massZ1);
   tree->SetBranchAddress("massZ2", &massZ2);
   tree->SetBranchAddress("met", &met);
   tree->SetBranchAddress("nFailedLeptonsZ2", &nFailedLeptonsZ2);

}

void SetNewTree(TTree* tree){

   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus("Run",1);
   tree->SetBranchStatus("LumiSect",1);
   tree->SetBranchStatus("Event",1);
   tree->SetBranchStatus("nVtx",1);
   tree->SetBranchStatus("passedTrig",1);
   tree->SetBranchStatus("passedFullSelection",1);
   tree->SetBranchStatus("passedZ4lSelection",1);
   tree->SetBranchStatus("passedZXCRSelection",1);
   tree->SetBranchStatus("passSmartCut",1);
   tree->SetBranchStatus("nZXCRFailedLeptons",1);
   tree->SetBranchStatus("finalState",1);
   tree->SetBranchStatus("dataMCWeight",1);
   tree->SetBranchStatus("pileupWeight",1);
   tree->SetBranchStatus("genWeight",1);
   tree->SetBranchStatus("sumWeight",1);
   tree->SetBranchStatus("crossSection",1);
   tree->SetBranchStatus("k_qqZZ_qcd_M",1);
   tree->SetBranchStatus("k_qqZZ_ewk",1);
   tree->SetBranchStatus("k_ggZZ",1);
   tree->SetBranchStatus("lep_id",1);
   tree->SetBranchStatus("lep_pt",1);
   tree->SetBranchStatus("lep_eta",1);
   tree->SetBranchStatus("lep_phi",1);
   tree->SetBranchStatus("lep_mass",1);
   tree->SetBranchStatus("lep_tightId",1);
   tree->SetBranchStatus("lep_RelIso",1);
   tree->SetBranchStatus("lep_RelIsoNoFSR",1);
   tree->SetBranchStatus("lep_Hindex",1);
   tree->SetBranchStatus("pTL1",1);
   tree->SetBranchStatus("pTL2",1);
   tree->SetBranchStatus("pTL3",1);
   tree->SetBranchStatus("pTL4",1);
   tree->SetBranchStatus("idL1",1);
   tree->SetBranchStatus("idL2",1);
   tree->SetBranchStatus("idL3",1);
   tree->SetBranchStatus("idL4",1);
   tree->SetBranchStatus("etaL1",1);
   tree->SetBranchStatus("etaL2",1);
   tree->SetBranchStatus("etaL3",1);
   tree->SetBranchStatus("etaL4",1);
   tree->SetBranchStatus("phiL1",1);
   tree->SetBranchStatus("phiL2",1);
   tree->SetBranchStatus("phiL3",1);
   tree->SetBranchStatus("phiL4",1);
   tree->SetBranchStatus("deltaphiL13",1);
   tree->SetBranchStatus("deltaphiL14",1);
   tree->SetBranchStatus("deltaphiL23",1);
   tree->SetBranchStatus("deltaphiL24",1);
   tree->SetBranchStatus("deltaphiZZ",1);
   tree->SetBranchStatus("mass4l",1);
   tree->SetBranchStatus("mass4mu",1);
   tree->SetBranchStatus("mass4e",1);
   tree->SetBranchStatus("mass2e2mu",1);
   tree->SetBranchStatus("pT4l",1);
   tree->SetBranchStatus("massZ1",1);
   tree->SetBranchStatus("massZ2",1);
   tree->SetBranchStatus("njets_pt30_eta4p7",1);
   tree->SetBranchStatus("njets_pt30_eta2p5",1);
   tree->SetBranchStatus("met",1);
   tree->SetBranchStatus("pTj1",1);
   tree->SetBranchStatus("etaj1",1);
   tree->SetBranchStatus("qgj1",1);
   tree->SetBranchStatus("pTj2",1);
   tree->SetBranchStatus("etaj2",1);
   tree->SetBranchStatus("qgj2",1);
   tree->SetBranchStatus("pt_leadingjet_pt30_eta4p7",1);
   tree->SetBranchStatus("pt_leadingjet_pt30_eta2p5",1);

   tree->SetBranchAddress("Run", &Run);
   tree->SetBranchAddress("Event", &Event);
   tree->SetBranchAddress("LumiSect", &LumiSect);
   tree->SetBranchAddress("nVtx", &nVtx);
   tree->SetBranchAddress("passedTrig", &passedTrig);
   tree->SetBranchAddress("passedFullSelection", &passedFullSelection);
   tree->SetBranchAddress("passedZ4lSelection", &passedZ4lSelection);
   tree->SetBranchAddress("passedZXCRSelection", &passedZXCRSelection);
   tree->SetBranchAddress("passSmartCut", &passSmartCut);
   tree->SetBranchAddress("nZXCRFailedLeptons", &nZXCRFailedLeptons);
   tree->SetBranchAddress("finalState", &finalState);
   tree->SetBranchAddress("dataMCWeight", &dataMCWeight);
   tree->SetBranchAddress("pileupWeight", &pileupWeight);
   tree->SetBranchAddress("genWeight", &genWeight);
   tree->SetBranchAddress("sumweight", &sumweight);
   tree->SetBranchAddress("crossSection", &crossSection);
   tree->SetBranchAddress("k_qqZZ_qcd_M", &k_qqZZ_qcd_M);
   tree->SetBranchAddress("k_qqZZ_ewk", &k_qqZZ_ewk);
   tree->SetBranchAddress("k_ggZZ", &k_ggZZ);
   tree->SetBranchAddress("lep_id", &lep_id);
   tree->SetBranchAddress("lep_pt", &lep_pt);
   tree->SetBranchAddress("lep_eta", &lep_eta);
   tree->SetBranchAddress("lep_phi", &lep_phi);
   tree->SetBranchAddress("lep_mass", &lep_mass);
   tree->SetBranchAddress("lep_tightId", &lep_tightId);
   tree->SetBranchAddress("lep_RelIso", &lep_RelIso);
   tree->SetBranchAddress("lep_RelIsoNoFSR", &lep_RelIsoNoFSR);
   tree->SetBranchAddress("lep_Hindex", &lep_Hindex_stdvec);
   tree->SetBranchAddress("pTL1", &pTL1);
   tree->SetBranchAddress("pTL2", &pTL2);
   tree->SetBranchAddress("pTL3", &pTL3);
   tree->SetBranchAddress("pTL4", &pTL4);
   tree->SetBranchAddress("idL1", &idL1);
   tree->SetBranchAddress("idL2", &idL2);
   tree->SetBranchAddress("idL3", &idL3);
   tree->SetBranchAddress("idL4", &idL4);
   tree->SetBranchAddress("etaL1", &etaL1);
   tree->SetBranchAddress("etaL2", &etaL2);
   tree->SetBranchAddress("etaL3", &etaL3);
   tree->SetBranchAddress("etaL4", &etaL4);
   tree->SetBranchAddress("phiL1", &phiL1);
   tree->SetBranchAddress("phiL2", &phiL2);
   tree->SetBranchAddress("phiL3", &phiL3);
   tree->SetBranchAddress("phiL4", &phiL4);
   tree->SetBranchAddress("deltaphiL13", &deltaphiL13);
   tree->SetBranchAddress("deltaphiL14", &deltaphiL14);
   tree->SetBranchAddress("deltaphiL23", &deltaphiL23);
   tree->SetBranchAddress("deltaphiL24", &deltaphiL24);
   tree->SetBranchAddress("deltaphiZZ", &deltaphiZZ);
   tree->SetBranchAddress("mass4l", &mass4l);
   tree->SetBranchAddress("mass4mu", &mass4mu);
   tree->SetBranchAddress("mass4e", &mass4e);
   tree->SetBranchAddress("mass2e2mu", &mass2e2mu);
   tree->SetBranchAddress("pT4l", &pT4l);
   tree->SetBranchAddress("massZ1", &massZ1);
   tree->SetBranchAddress("massZ2", &massZ2);
   tree->SetBranchAddress("njets_pt30_eta4p7", &njets_pt30_eta4p7);
   tree->SetBranchAddress("njets_pt30_eta2p5", &njets_pt30_eta2p5);
   tree->SetBranchAddress("met", &met);
   tree->SetBranchAddress("pTj1", &pTj1);
   tree->SetBranchAddress("etaj1", &etaj1);
   tree->SetBranchAddress("qgj1", &qgj1);
   tree->SetBranchAddress("pTj2", &pTj2);
   tree->SetBranchAddress("etaj2", &etaj2);
   tree->SetBranchAddress("qgj2", &qgj2);
   //tree->SetBranchAddress("pt_leadingjet_pt30_eta4p7", &pt_leadingjet_pt30_eta4p7);
   //tree->SetBranchAddress("pt_leadingjet_pt30_eta2p5", &pt_leadingjet_pt30_eta2p5);

}
