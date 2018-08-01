#include "ZXConfig.h" 
#include "deltaPhi.h"
#include "ZZ4LAnalysisTree.h"
#include "LeptonEfficiency.h"
#include "PileupWeight.h"
#include "EbECorrection.h"
#include "Helper.h"
#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "ZZMatrixElement/MELA/interface/TUtil.hh"
//#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"
#include "KinZfitter/KinZfitter.h"
using namespace std;
/////////////////////


void SetNewTree(TTree* newtree);
void ReadTree(TTree* tree, TTree* & newtree, TString filename);
void SkimTree(TTree* tree, TTree* & newtree, TString filename);

TString filename;
bool debug;

Mela* mela = new Mela(13.0,125.0,TVar::SILENT);

//KalmanMuonCalibrator *kalmanMuonCalibrator;

//KinZfitter *kinZfitter;

int main(int argc, char *argv[])
{    
     
  debug = false;     

  if(argc > 6)  {
      cout<<argv[0]<<" filename "<<argv[1]<<" outfile "<<argv[2]<<" isData "<<argv[3]<<endl;
      return -1;
    }

  //gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc530/pluginUFHZZAnalysisRun2UFHZZ4LAna.so");  

  //mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  /////////////////////

  filename = argv[1];
  TString outfilename = argv[2];
 
  if(atof(argv[3])>0) { isData = true; }
  else {isData = false; }
  if(atof(argv[4])>0) { 
      job=strtol(argv[4], NULL, 10); njobs=strtol(argv[5], NULL, 10);
      cout<<"job "<<job<<" of "<<njobs<<endl;
      outfilename+="_";
      outfilename+=argv[4];
  }

  TFile* infile = TFile::Open(filename+".root");
  TTree* tree;
  tree = (TTree*) infile->Get("Ana/passedEvents");
  if(!tree) tree = (TTree*) infile->Get("passedEvents");
  if(!tree) tree = (TTree*) infile->Get("selectedEvents");
  if(!tree) { cout<<"ERROR could not find the tree for "<<filename<<endl; return -1;}

  TH1F* sumWeights = (TH1F*) infile->Get("Ana/sumWeights");
  sumW = 1.0;
  if(sumWeights) sumW = sumWeights->GetBinContent(1); 

  cout<<"sumW is "<<sumW<<endl;

  // read tree     

  TString name = outfilename;
  TFile* tmpFile =  new TFile(name+".root","RECREATE");
  TTree* newtree = new TTree("passedEvents","passedEvents");

  if(debug)cout<<"start setting new tree "<<endl;

  SetNewTree(newtree);

  if(debug)cout<<"start reading tree "<<endl;

  //ReadTree(tree, newtree, filename);
  SkimTree(tree, newtree, filename);

  if(debug)cout<<"end reading tree"<<endl;

  tmpFile->cd();

  newtree->Write("passedEvents",TObject::kOverwrite);
  tmpFile->Close(); 

}

void SkimTree(TTree* tree, TTree* & newtree, TString filename){

    ZZ4LAnalysisTree::setAddresses(tree, filename);
    
    int firstevt=0; int lastevt=tree->GetEntries();
    if (job>0) {
        firstevt = tree->GetEntries()*(job-1)/njobs;
        lastevt = tree->GetEntries()*(job)/njobs-1;
    }
    
    for(int evt=0; evt < tree->GetEntries(); evt++) { //event loop
        
        if (evt<firstevt) continue;
        if (evt>lastevt) continue;
       
        if(evt%1000==0) cout<<"Event "<<evt<<"/"<<tree->GetEntries()<<endl;
        tree->GetEntry(evt);

        if (!passedZ1LSelection) continue; 

        passTrig=false;
        if (isData) {
            bool passSingleElectronTrig = false;
            bool passSingleMuonTrig = false;
            bool passDoubleMuonTrig = false;
            bool passMuonEGTrig = false;
            bool passDoubleEGTrig = false;
            if (strstr((*triggersPassed).c_str(),"HLT_Ele35_WPTight_Gsf_v")) {
                passSingleElectronTrig=true;
            } else if (strstr((*triggersPassed).c_str(),"HLT_Ele38_WPTight_Gsf_v")) {  
                passSingleElectronTrig=true;
            } else if (strstr((*triggersPassed).c_str(),"HLT_Ele40_WPTight_Gsf_v")) {
                passSingleElectronTrig=true;
            };

            if (strstr((*triggersPassed).c_str(),"HLT_IsoMu27")) {
                passSingleMuonTrig=true;
            };
            
            // double mu
            if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v")) passDoubleMuonTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v")) passDoubleMuonTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_TripleMu_12_10_5_v")) passDoubleMuonTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_TripleMu_10_5_5_DZ_v")) passDoubleMuonTrig=true;

            // double eg
            if (strstr((*triggersPassed).c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v")) passDoubleEGTrig=true;
           
            // muon eg 
            if (strstr((*triggersPassed).c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passMuonEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_DoubleEle33_CaloIdL_MW_v")) passMuonEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passMuonEGTrig=true;
        
            if (strstr(filename,"SingleElectron")) {
                if (passSingleElectronTrig) passTrig=true;
            } else if (strstr(filename,"SingleMuon")) {
                if (!passSingleElectronTrig && passSingleMuonTrig) passTrig=true;
            } else if (strstr(filename,"DoubleMuon")) {
                if (!passSingleElectronTrig && !passSingleMuonTrig && passDoubleMuonTrig) passTrig=true;
            } else if (strstr(filename,"DoubleEG")) {
                if (!passSingleElectronTrig && !passSingleMuonTrig && !passDoubleMuonTrig && passDoubleEGTrig) passTrig=true;
            } else if (strstr(filename,"MuonEG")) {
                if (!passSingleElectronTrig && !passSingleMuonTrig && !passDoubleMuonTrig && !passDoubleEGTrig && passMuonEGTrig) passTrig=true;
            };
        } else {
            passTrig = true;
            //pileupWeight = float(puweight(nInt));
            //sumweight += pileupWeight*genWeight;
        }

        if (!passTrig) continue;

        newtree->Fill();
    }
}

void ReadTree(TTree* tree, TTree* & newtree, TString filename){
    
    ZZ4LAnalysisTree::setAddresses(tree, filename);

    float npass = 0.0;
    float sumweight = 0.0;

    if(debug) cout<<"start looping"<<endl;

    int firstevt=0; int lastevt=tree->GetEntries();
    if (job>0) {
        firstevt = tree->GetEntries()*(job-1)/njobs;
        lastevt = tree->GetEntries()*(job)/njobs-1;
    }
    
    for(int evt=0; evt < tree->GetEntries(); evt++) { //event loop
    
        if (evt<firstevt) continue;
        if (evt>lastevt) continue;
       
        if(evt%1000==0) cout<<"Event "<<evt<<"/"<<tree->GetEntries()<<endl;
        tree->GetEntry(evt);
        passTrig=false;
        if (isData) {
            bool passSingleElectronTrig = false;
            bool passSingleMuonTrig = false;
            bool passDoubleMuonTrig = false;
            bool passMuonEGTrig = false;
            bool passDoubleEGTrig = false;
            if (strstr((*triggersPassed).c_str(),"HLT_Ele35_WPTight_Gsf_v")) {
                passSingleElectronTrig=true;
            } else if (strstr((*triggersPassed).c_str(),"HLT_Ele38_WPTight_Gsf_v")) {  
                passSingleElectronTrig=true;
            } else if (strstr((*triggersPassed).c_str(),"HLT_Ele40_WPTight_Gsf_v")) {
                passSingleElectronTrig=true;
            };

            if (strstr((*triggersPassed).c_str(),"HLT_IsoMu27")) {
                passSingleMuonTrig=true;
            };
            
            // double mu
            if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v")) passDoubleMuonTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v")) passDoubleMuonTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_TripleMu_12_10_5_v")) passDoubleMuonTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_TripleMu_10_5_5_DZ_v")) passDoubleMuonTrig=true;

            // double eg
            if (strstr((*triggersPassed).c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")) passDoubleEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v")) passDoubleEGTrig=true;
           
            // muon eg 
            if (strstr((*triggersPassed).c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passMuonEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_DoubleEle33_CaloIdL_MW_v")) passMuonEGTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passMuonEGTrig=true;
        
            if (strstr(filename,"SingleElectron")) {
                if (passSingleElectronTrig) passTrig=true;
            } else if (strstr(filename,"SingleMuon")) {
                if (!passSingleElectronTrig && passSingleMuonTrig) passTrig=true;
            } else if (strstr(filename,"DoubleMuon")) {
                if (!passSingleElectronTrig && !passSingleMuonTrig && passDoubleMuonTrig) passTrig=true;
            } else if (strstr(filename,"DoubleEG")) {
                if (!passSingleElectronTrig && !passSingleMuonTrig && !passDoubleMuonTrig && passDoubleEGTrig) passTrig=true;
            } else if (strstr(filename,"MuonEG")) {
                if (!passSingleElectronTrig && !passSingleMuonTrig && !passDoubleMuonTrig && !passDoubleEGTrig && passMuonEGTrig) passTrig=true;
            };
        } else {
            passTrig = true;
            //pileupWeight = float(puweight(nInt));
            sumweight += pileupWeight*genWeight;
        }

        if (!passTrig) continue;

        const double Zmass = 91.1876;
  
        unsigned int Nlep = (*lepFSR_pt).size();
        if (debug) cout<<Nlep<<" leptons in total"<<endl;
        if( Nlep != 3 ) return;
  
        // First, make all Z candidates including any FSR photons
        int n_Zs=0;
        vector<int> Z_Z1L_lepindex1;
        vector<int> Z_Z1L_lepindex2;  

        bool foundZ1LCandidate=false;

        for(unsigned int i=0; i<Nlep; i++){
            for(unsigned int j=i+1; j<Nlep; j++){
        
                // same flavor opposite charge
                if(((*lep_id)[i]+(*lep_id)[j])!=0) continue;
        
                TLorentzVector li, lj;
                li.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
                lj.SetPtEtaPhiM((*lep_pt)[j],(*lep_eta)[j],(*lep_phi)[j],(*lep_mass)[j]);
        
                TLorentzVector lifsr, ljfsr;
                lifsr.SetPtEtaPhiM((*lepFSR_pt)[i],(*lepFSR_eta)[i],(*lepFSR_phi)[i],(*lepFSR_mass)[i]);
                ljfsr.SetPtEtaPhiM((*lepFSR_pt)[j],(*lepFSR_eta)[j],(*lepFSR_phi)[j],(*lepFSR_mass)[j]);
        
                TLorentzVector liljfsr = lifsr+ljfsr;
         
                TLorentzVector Z, Z_noFSR;
                Z = lifsr+ljfsr;
                Z_noFSR = li+lj;
        
                if (debug) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;
        
                if (Z.M()>0.0) {
                    n_Zs++;
                    Z_Z1L_lepindex1.push_back(i);
                    Z_Z1L_lepindex2.push_back(j);
                    if (debug) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
                }
        
            } // lep i
        } // lep j

        bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
        for(unsigned int i =0; i<Nlep; i++) {
            if((*lep_id)[i]==-13) Nmm = Nmm+1;
            if((*lep_id)[i]==13) Nmp = Nmp+1;
        }
        for(unsigned int i =0; i<Nlep; i++) {
            if((*lep_id)[i]==-11) Nem = Nem+1;
            if((*lep_id)[i]==11) Nep = Nep+1;
        }

        if(Nmm>=1 && Nmp>=1) properLep_ID = true; //2mu + x
        if(Nem>=1 && Nep>=1) properLep_ID = true; //2e + x

        // proper charge flavor combination for Z + 1L
        if(!properLep_ID) return;

        if (debug) cout<<"found three leptons"<<endl;

        // Consider all Z candidates
        double minZ1DeltaM=9999.9;
        for (int i=0; i<n_Zs; i++) {

            int i1 = Z_Z1L_lepindex1[i]; int i2 = Z_Z1L_lepindex2[i];
            int j1 = 3 - i1 - i2; // index of the third lepton (check if this works)

            TLorentzVector lep_i1, lep_i2, lep_j1;
            lep_i1.SetPtEtaPhiM((*lepFSR_pt)[i1],(*lepFSR_eta)[i1],(*lepFSR_phi)[i1],(*lepFSR_mass)[i1]);
            lep_i2.SetPtEtaPhiM((*lepFSR_pt)[i2],(*lepFSR_eta)[i2],(*lepFSR_phi)[i2],(*lepFSR_mass)[i2]);
            lep_j1.SetPtEtaPhiM((*lepFSR_pt)[j1],(*lepFSR_eta)[j1],(*lepFSR_phi)[j1],(*lepFSR_mass)[j1]);

            TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr;
            lep_i1_nofsr.SetPtEtaPhiM((*lep_pt)[i1],(*lep_eta)[i1],(*lep_phi)[i1],(*lep_mass)[i1]);
            lep_i2_nofsr.SetPtEtaPhiM((*lep_pt)[i2],(*lep_eta)[i2],(*lep_phi)[i2],(*lep_mass)[i2]);
            lep_j1_nofsr.SetPtEtaPhiM((*lep_pt)[j1],(*lep_eta)[j1],(*lep_phi)[j1],(*lep_mass)[j1]);

            TLorentzVector Zi;
            Zi = lep_i1+lep_i2;
            //Zi.SetPtEtaPhiM(Z_Z1L_pt[i],Z_Z1L_eta[i],Z_Z1L_phi[i],Z_Z1L_mass[i]);

            if (debug) {cout<<"Z candidate Zi->M() "<<Zi.M()<<endl;}

            TLorentzVector Z1 = Zi;
            double Z1DeltaM = abs(Zi.M()-Zmass);
            int Z1_lepindex[2] = {0,0};
            if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
            else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }

            // Check Leading and Subleading pt Cut
            vector<double> allPt;
            allPt.push_back(lep_i1.Pt()); allPt.push_back(lep_i2.Pt());
            std::sort(allPt.begin(), allPt.end());
            if (debug) cout<<" leading pt: "<<allPt[1]<<" cut: "<<leadingPtCut<<" subleadingPt: "<<allPt[0]<<" cut: "<<subleadingPtCut<<endl;
            if (allPt[1]<leadingPtCut || allPt[0]<subleadingPtCut ) continue;

            // Check dR(li,lj)>0.02 for any i,j
            vector<double> alldR;
            alldR.push_back(lep_i1.DeltaR(lep_i2));
            alldR.push_back(lep_i1.DeltaR(lep_j1));
            alldR.push_back(lep_i2.DeltaR(lep_j1));
            if (debug) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
            if (*min_element(alldR.begin(),alldR.end())<0.02) continue;

            // Check M(l+,l-)>4.0 GeV for any OS pair
            // Do not include FSR photons
            vector<double> allM;
            TLorentzVector i1i2;
            i1i2 = (lep_i1_nofsr)+(lep_i2_nofsr); allM.push_back(i1i2.M());
            if ((*lep_id)[i1]*(*lep_id)[j1]<0) {
                TLorentzVector i1j1;
                i1j1 = (lep_i1_nofsr)+(lep_j1_nofsr); allM.push_back(i1j1.M());
            } else {
                TLorentzVector i2j1;
                i2j1 = (lep_i2_nofsr)+(lep_j1_nofsr); allM.push_back(i2j1.M());
            }
            if (debug) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
            if (*min_element(allM.begin(),allM.end())<4.0) continue;

            // Check isolation cut (without FSR ) for Z1 leptons
            //if (*(lep_RelIsoNoFSR)[Z1_lepindex[0]]>((abs(*(lep_id)[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
            //if (*(lep_RelIsoNoFSR)[Z1_lepindex[1]]>((abs(*(lep_id)[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
            
            if ((*lep_RelIsoNoFSR)[Z1_lepindex[0]]>((abs((*lep_id)[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue;
            if ((*lep_RelIsoNoFSR)[Z1_lepindex[1]]>((abs((*lep_id)[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue;

            // Check tight ID cut for Z1 leptons
            if (!((*lep_tightId)[Z1_lepindex[0]])) continue; // checking tight lepton ID
            if (!((*lep_tightId)[Z1_lepindex[1]])) continue; // checking tight lepton ID
            
            if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) ) continue;

            if (debug) cout<<"good Z1L candidate, Z1DeltaM: "<<Z1DeltaM<<" minZ1DeltaM: "<<minZ1DeltaM<<endl;

            // Check if this candidate has the best Z1 and highest scalar sum of Z2 lepton pt

            if ( Z1DeltaM<=minZ1DeltaM ) {

                minZ1DeltaM = Z1DeltaM;

                TLorentzVector Z1L;
                Z1L = Z1+lep_j1;

                massZ1_Z1L = Z1.M();
                mass3l = Z1L.M();

                lep_Hindex[0] = Z1_lepindex[0];
                lep_Hindex[1] = Z1_lepindex[1];
                lep_Hindex[2] = j1;

                if (debug) cout<<" new best Z1L candidate: massZ1: "<<massZ1<<" (mass3l: "<<mass3l<<")"<<endl;
                foundZ1LCandidate=true;

            }
        }
    } // Event loop    
}

void SetNewTree(TTree* newtree){

    newtree->Branch("Run",&Run,"Run/l");
    newtree->Branch("Event",&Event,"Event/l");
    newtree->Branch("LumiSect",&LumiSect,"LumiSect/l");
    newtree->Branch("nVtx",&nVtx,"nVtx/I");
    newtree->Branch("passedTrig",&passedTrig,"passedTrig/O");
    newtree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
    newtree->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
    newtree->Branch("passedZ1LSelection",&passedZ1LSelection,"passedZ1LSelection/O");
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
    //newtree->Branch("mass4lErr",&mass4lErr,"mass4lErr/F");
    //newtree->Branch("mass4lREFIT",&mass4lREFIT,"mass4lREFIT/F");
    //newtree->Branch("mass4lErrREFIT",&mass4lErrREFIT,"mass4lErrREFIT/F");

    newtree->Branch("mass4mu",&mass4mu,"mass4mu/F");
    newtree->Branch("mass4e",&mass4e,"mass4e/F");
    newtree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/F");
    newtree->Branch("pT4l",&pT4l,"pT4l/F");
    newtree->Branch("massZ1",&massZ1,"massZ1/F");
    newtree->Branch("massZ2",&massZ2,"massZ2/F"); 
    newtree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/I");
    newtree->Branch("njets_pt30_eta2p5",&njets_pt30_eta2p5,"njets_pt30_eta2p5/I");
    newtree->Branch("met",&met,"met/F"); 
    newtree->Branch("me_qqZZ_MCFM",&me_qqZZ_MCFM,"me_qqZZ_MCFM/F");
    newtree->Branch("pTj1",&pTj1,"pTj1/F");
    newtree->Branch("etaj1",&etaj1,"etaj1/F");
    newtree->Branch("qgj1",&qgj1,"qgj1/F");
    newtree->Branch("pTj2",&pTj2,"pTj2/F");
    newtree->Branch("etaj2",&etaj2,"etaj2/F");
    newtree->Branch("qgj2",&qgj2,"qgj1/F");

    newtree->Branch("pt_leadingjet_pt30_eta4p7",&pTj1,"pt_leadingjet_pt30_eta4p7/F");
    newtree->Branch("pt_leadingjet_pt30_eta2p5",&pTj1_2p5,"pt_leadingjet_pt30_eta2p5/F");

    newtree->Branch("D_bkg_kin", &D_bkg_kin, "D_bkg_kin/F");
    newtree->Branch("D_bkg", &D_bkg, "D_bkg/F");
    newtree->Branch("Dgg10_VAMCFM", &Dgg10_VAMCFM, "Dgg10_VAMCFM/F");
    newtree->Branch("D_g4", &D_g4, "D_g4/F");
    newtree->Branch("D_VBF", &D_VBF, "D_VBF/F");
    newtree->Branch("D_VBF1j",&D_VBF1j,"D_VBF1j/F");
    newtree->Branch("D_HadWH",&D_HadWH,"D_HadWH/F");
    newtree->Branch("D_HadZH",&D_HadZH,"D_HadZH/F");
    newtree->Branch("D_VBF_QG",&D_VBF_QG,"D_VBF_QG/F");
    newtree->Branch("D_VBF1j_QG",&D_VBF1j_QG,"D_VBF1j_QG/F");
    newtree->Branch("D_HadWH_QG",&D_HadWH_QG,"D_HadWH_QG/F");
    newtree->Branch("D_HadZH_QG",&D_HadZH_QG,"D_HadZH_QG/F");
    newtree->Branch("EventCat",&EventCat,"EventCat/I");
    newtree->Branch("massZ1_Z1L",&massZ1_Z1L,"massZ1_Z1L/F");
    newtree->Branch("mass3l",&mass3l,"mass3l/F");
 
}


