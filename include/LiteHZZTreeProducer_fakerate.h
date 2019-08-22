#ifndef LiteHZZTreeProducer_fakerate_h
#define LiteHZZTreeProducer_fakerate_h

#include "LiteHZZTreeProducer_fakerate_Linkdef.h"
#include "Analyzer.h"
#include "HZZTree_fakerate.h"
#include "deltaR.h"
#include <cmath>

using namespace std;  

class LiteHZZTreeProducer_fakerate : public Analyzer 
{
    public:
        // static const double Zmass = 91.1876;
        const double Zmass = 91.1876;

        LiteHZZTreeProducer_fakerate();
        ~LiteHZZTreeProducer_fakerate();
        LiteHZZTreeProducer_fakerate(
                double m4lHighCut_in,
                double m4lLowCut_in,
                double mZ2HighCut_in,
                double mZ2LowCut_in,
                double mZ1HighCut_in,
                double mZ1LowCut_in,
                double isoCutEl_in,
                double isoCutMu_in,
                TString outputDir_in,
                TString outFileName_in,
                bool do_wrong_fc_in=false
                );
         LiteHZZTreeProducer_fakerate(
                TString outputDir_in,
                TString outFileName_in
                );

        int process();
        void setup();
        void end();
        void initTree();
        bool passSelection();
        void setDebugMode(bool debug_in);
        void sortedArray(double x, double y, double z, double sortarray[3]);

        double m4lHighCut = 9999999.;
        double m4lLowCut = 70.0;
        double mZ2HighCut=120.0;
        //double mZ2High=999999.;
        double mZ2LowCut=4.0;
        double mZ1HighCut=120.0;
        double mZ1LowCut=40.0;
        //double isoCutEl=999999.;
        double isoCutEl=0.35;
        double isoCutMu=0.35;
        double leadingPtCut=20.0; 
        double subleadingPtCut=10.0; 

        TString treeName = "Ana/passedEvents";
        TString outTreeName = "passedEvents";
        bool debug = false;
        bool do_wrong_fc = false;

        TString outFileName;
        TFile* outFile=0;
        TTree* outTree=0;
};

LiteHZZTreeProducer_fakerate::LiteHZZTreeProducer_fakerate(){

}

LiteHZZTreeProducer_fakerate::~LiteHZZTreeProducer_fakerate(){

}

bool LiteHZZTreeProducer_fakerate::passSelection(){
    return true;
}


void LiteHZZTreeProducer_fakerate::initTree(){
    setHZZTree(tree);
}

void LiteHZZTreeProducer_fakerate::setDebugMode(bool debug_in){
    debug = debug_in;
}

void LiteHZZTreeProducer_fakerate::setup(){
    outFile = TFile::Open(outputDir+outFileName,fOptionWrite);
    outTree = new TTree(outTreeName,outTreeName);

    initNewLiteTree_fakerate(outTree);
}

void LiteHZZTreeProducer_fakerate::end(){
    outFile->cd();
    outTree->Write(outTreeName,TObject::kOverwrite);
    outFile->Close(); 
}

LiteHZZTreeProducer_fakerate::LiteHZZTreeProducer_fakerate(
                double m4lHighCut_in,
                double m4lLowCut_in,
                double mZ2HighCut_in,
                double mZ2LowCut_in,
                double mZ1HighCut_in,
                double mZ1LowCut_in,
                double isoCutEl_in,
                double isoCutMu_in,
                TString outputDir_in,
                TString outFileName_in,
                bool do_wrong_fc_in
                ){
    m4lHighCut  = m4lHighCut_in;
    m4lLowCut   = m4lLowCut_in;
    mZ2HighCut  = mZ2HighCut_in;
    mZ2LowCut   = mZ2LowCut_in;
    mZ1HighCut  = mZ1HighCut_in;
    mZ1LowCut   = mZ1LowCut_in;
    isoCutEl    = isoCutEl_in;
    isoCutMu    = isoCutMu_in;
    outputDir   = outputDir_in;
    outFileName = outFileName_in;
    do_wrong_fc = do_wrong_fc_in;
}

LiteHZZTreeProducer_fakerate::LiteHZZTreeProducer_fakerate(
                TString outputDir_in,
                TString outFileName_in
                ){
    outputDir   = outputDir_in;
    outFileName = outFileName_in;
}

void LiteHZZTreeProducer_fakerate::sortedArray(double x, double y, double z, double sortarray[3])
{
    double max_v = max(x,max(y,z));
    double min_v = min(x,min(y,z));
    double mid_v = -1;
    if (x != max_v && x != min_v) mid_v = x;
    if (y != max_v && y != min_v) mid_v = y;
    if (z != max_v && z != min_v) mid_v = z;

    sortarray[0] = max_v;
    sortarray[1] = mid_v;
    sortarray[2] = min_v;
    //return sortarray;
}

int LiteHZZTreeProducer_fakerate::process(){

    int Nlep = (*lep_id).size();
    int nTightLep = 0;
    int nLooseLep = 0;
    vector<int> tightIsoLepIndex;
    vector<int> looseIsoLepIndex;

    if((*lep_id).size()<3 || (*lep_pt).size()<3) return -1;

    if(passedTrig == 0) return -1; 

    for (int iLep = 0; iLep < Nlep; iLep++) {
        if ((*lep_tightId)[iLep] == 1 && (*lep_RelIso)[iLep] < 0.35 && (abs((*lep_id)[iLep]) == 13 || abs((*lep_id)[iLep]) == 11) && nTightLep < 2){
            nTightLep++;
            tightIsoLepIndex.push_back(iLep);
        }
        else if ((*lep_tightId)[iLep] == 1 && (abs((*lep_id)[iLep]) == 13 || abs((*lep_id)[iLep]) == 11)){
            nLooseLep++;
            looseIsoLepIndex.push_back(iLep);
        }
    }

    if (!(nTightLep == 2 && nLooseLep == 1)) return -1;

    int index1 = tightIsoLepIndex[0];
    int index2 = tightIsoLepIndex[1];
    int index3 = looseIsoLepIndex[0];

    if ( (*lep_id)[index1] + (*lep_id)[index2] != 0 ) return -1;

    if ((*lep_Sip)[index1] > 3 || (*lep_Sip)[index2] > 3 || (*lep_Sip)[index3] > 3) return -1;

    double pTs[3]; sortedArray((*lep_pt)[index1], (*lep_pt)[index2], (*lep_pt)[index3], pTs);
    if (pTs[0] < 20 || pTs[1] < 10 || pTs[2] < 5) return -1;

    double isos[3]; sortedArray((*lep_RelIso)[index1], (*lep_RelIso)[index2], (*lep_RelIso)[index3], isos);
    double sips[3]; sortedArray((*lep_Sip)[index1], (*lep_Sip)[index2], (*lep_Sip)[index3], sips);

    TLorentzVector Lep1,Lep2,Lep3;
    Lep1.SetPtEtaPhiM((*lep_pt)[index1],(*lep_eta)[index1],(*lep_phi)[index1],(*lep_mass)[index1]); 
    Lep2.SetPtEtaPhiM((*lep_pt)[index2],(*lep_eta)[index2],(*lep_phi)[index2],(*lep_mass)[index2]);
    Lep3.SetPtEtaPhiM((*lep_pt)[index3],(*lep_eta)[index3],(*lep_phi)[index3],(*lep_mass)[index3]);

    TLorentzVector Leps = Lep1+Lep2+Lep3;
    double LepsPhi = Leps.Phi(); double LepsPt = Leps.Pt();
    //double mT = sqrt(2*(met)*LepsPt*(1-cos(deltaPhi(LepsPhi, met_phi) ) ) );

    vector<double> dRs; double tmpDr = -1;
    vector<TLorentzVector> mlls;
    if ((*lep_id)[index1] != (*lep_id)[index2]) {
        mlls.push_back(Lep1 + Lep2);
        dRs.push_back(deltaR(Lep1.Eta(),Lep1.Phi(),Lep2.Eta(),Lep2.Phi()));

    } else {tmpDr = deltaR(Lep1.Eta(),Lep1.Phi(),Lep2.Eta(),Lep2.Phi());}

    if ((*lep_id)[index1] != (*lep_id)[index3]) {
        mlls.push_back(Lep1 + Lep3);
        dRs.push_back(deltaR(Lep1.Eta(),Lep1.Phi(),Lep3.Eta(),Lep3.Phi()));
  
    } else {tmpDr = deltaR(Lep1.Eta(),Lep1.Phi(),Lep3.Eta(),Lep3.Phi());}

    if ((*lep_id)[index2] != (*lep_id)[index3]) {
        mlls.push_back(Lep2 + Lep3);
        dRs.push_back(deltaR(Lep2.Eta(),Lep2.Phi(),Lep3.Eta(),Lep3.Phi()));

    } else {tmpDr = deltaR(Lep2.Eta(),Lep2.Phi(),Lep3.Eta(),Lep3.Phi());}

    dRs.push_back(tmpDr);

    TLorentzVector Z = Lep1+Lep2;
    if (!(Z.M() > 86 && Z.M() < 96)) return kTRUE;

    // 2 OSSF Pairs
    /*bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
    for(unsigned int i =0; i<Nlep; i++) {
        if((*lep_id)[i]==-13) Nmm = Nmm+1;
        if((*lep_id)[i]==13) Nmp = Nmp+1;
    };
    for(unsigned int i =0; i<Nlep; i++) {
        if((*lep_id)[i]==-11) Nem = Nem+1;
        if((*lep_id)[i]==11) Nep = Nep+1;
    };

    if(Nmm>=2 && Nmp>=2) properLep_ID = true; //4mu
    if(Nem>=2 && Nep>=2) properLep_ID = true; //4e
    if(Nmm>0 && Nmp>0 && Nem>0 && Nep>0) properLep_ID = true; //2e2mu 


    int n_Zs=0;
    vector<int> Z_lepindex1;
    vector<int> Z_lepindex2;
    vector<float> Z_pt, Z_eta, Z_phi, Z_mass;
    lep_Hindex_stdvec->clear();
   
    int properZCand=0; 
    for(unsigned int i=0; i<Nlep; i++){
        for(unsigned int j=i+1; j<Nlep; j++){
            // same flavor opposite charge
            if (!do_wrong_fc && ((*lep_id)[i]+(*lep_id)[j])!=0) continue;
            if (do_wrong_fc && ((*lep_id)[i]+(*lep_id)[j])==0) properZCand++;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
            lj.SetPtEtaPhiM((*lep_pt)[j],(*lep_eta)[j],(*lep_phi)[j],(*lep_mass)[j]);
            
            TLorentzVector lifsr, ljfsr;
            lifsr.SetPtEtaPhiM((*lepFSR_pt)[i],(*lepFSR_eta)[i],(*lepFSR_phi)[i],(*lepFSR_mass)[i]);
            ljfsr.SetPtEtaPhiM((*lepFSR_pt)[j],(*lepFSR_eta)[j],(*lepFSR_phi)[j],(*lepFSR_mass)[j]);
            
            TLorentzVector liljfsr = lifsr+ljfsr;
            
            if (debug) {
                cout<<"OSSF pair: i="<<i<<" id1="<<(*lep_id)[i]<<" j="<<j<<" id2="<<(*lep_id)[j]<<" pt1: "
                    <<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;    
            }
            
            TLorentzVector Z, Z_noFSR;
            Z = lifsr+ljfsr;
            Z_noFSR = li+lj;
            
            if (debug) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2LowCut<<endl;
            
            if (Z.M()>0.0) {
                n_Zs++;
                Z_pt.push_back(Z.Pt());
                Z_eta.push_back(Z.Eta());
                Z_phi.push_back(Z.Phi());
                Z_mass.push_back(Z.M());
                Z_lepindex1.push_back(i);
                Z_lepindex2.push_back(j);
                if (debug) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
            }
            
        } // lep i
    } // lep j*/

    // Consider all ZZ candidates
    /*TLorentzVector Z1Vec, Z2Vec, HVec;
    double minZ1DeltaM_SR=9999.9; double minZ1DeltaM_CR=99999.9;
    double maxZ2SumPt_SR=0.0; double maxZ2SumPt_CR=0.0;
    double max_D_bkg_kin_SR=0.0; double max_D_bkg_kin_CR=0.0;
    bool foundSRCandidate=false;
    
    passedZ4lSelection=false;
    passedFullSelection=false;*/
    
    /*vector<int> lep_Hindex, Z_Hindex;
    vector<int> Z_Hindex;
    for (int i=0; i<4; i++) {
        if (i<2) Z_Hindex.push_back(-1);
        //lep_Hindex.push_back(-1);
        lep_Hindex[i]=-1;
    }*/
    
    /*for (int i=0; i<n_Zs; i++) {
        for (int j=i+1; j<n_Zs; j++) {
            
            int i1 = Z_lepindex1[i]; int i2 = Z_lepindex2[i];                            
            int j1 = Z_lepindex1[j]; int j2 = Z_lepindex2[j];                            
            
            if (i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;
            
            TLorentzVector lep_i1, lep_i2, lep_j1, lep_j2;
            lep_i1.SetPtEtaPhiM((*lepFSR_pt)[i1],(*lepFSR_eta)[i1],(*lepFSR_phi)[i1],(*lepFSR_mass)[i1]);
            lep_i2.SetPtEtaPhiM((*lepFSR_pt)[i2],(*lepFSR_eta)[i2],(*lepFSR_phi)[i2],(*lepFSR_mass)[i2]);
            lep_j1.SetPtEtaPhiM((*lepFSR_pt)[j1],(*lepFSR_eta)[j1],(*lepFSR_phi)[j1],(*lepFSR_mass)[j1]);
            lep_j2.SetPtEtaPhiM((*lepFSR_pt)[j2],(*lepFSR_eta)[j2],(*lepFSR_phi)[j2],(*lepFSR_mass)[j2]);
            
            TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr, lep_j2_nofsr;
            lep_i1_nofsr.SetPtEtaPhiM((*lep_pt)[i1],(*lep_eta)[i1],(*lep_phi)[i1],(*lep_mass)[i1]);
            lep_i2_nofsr.SetPtEtaPhiM((*lep_pt)[i2],(*lep_eta)[i2],(*lep_phi)[i2],(*lep_mass)[i2]);
            lep_j1_nofsr.SetPtEtaPhiM((*lep_pt)[j1],(*lep_eta)[j1],(*lep_phi)[j1],(*lep_mass)[j1]);
            lep_j2_nofsr.SetPtEtaPhiM((*lep_pt)[j2],(*lep_eta)[j2],(*lep_phi)[j2],(*lep_mass)[j2]);
            
            TLorentzVector Zi, Zj;
            Zi.SetPtEtaPhiM(Z_pt[i],Z_eta[i],Z_phi[i],Z_mass[i]);
            Zj.SetPtEtaPhiM(Z_pt[j],Z_eta[j],Z_phi[j],Z_mass[j]);
            
            if (debug) {cout<<"ZZ candidate Zi->M() "<<Zi.M()<<" Zj->M() "<<Zj.M()<<endl;}
            
            TLorentzVector Z1, Z2;
            int Z1index, Z2index;
            int Z1_lepindex[2] = {0,0};
            int Z2_lepindex[2] = {0,0};
            double Z1DeltaM, Z2SumPt;
            
            if (abs(Zi.M()-Zmass)<abs(Zj.M()-Zmass)) { 
                Z1index = i; Z2index = j;
                Z1 = Zi; Z2 = Zj;                 
                if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
                else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }                
                if (lep_j1.Pt()>lep_j2.Pt()) { Z2_lepindex[0] = j1;  Z2_lepindex[1] = j2; } 
                else { Z2_lepindex[0] = j2;  Z2_lepindex[1] = j1; }                
                Z1DeltaM = abs(Zi.M()-Zmass); 
                Z2SumPt = lep_j1.Pt()+lep_j2.Pt();
            }
            else { 
                Z1index = j; Z2index = i;
                Z1 = Zj; Z2 = Zi; 
                if (lep_j1.Pt()>lep_j2.Pt()) { Z1_lepindex[0] = j1;  Z1_lepindex[1] = j2; }
                else { Z1_lepindex[0] = j2;  Z1_lepindex[1] = j1; }
                if (lep_i1.Pt()>lep_i2.Pt()) { Z2_lepindex[0] = i1;  Z2_lepindex[1] = i2; }
                else { Z2_lepindex[0] = i2;  Z2_lepindex[1] = i1; }
                Z1DeltaM = abs(Zj.M()-Zmass); 
                Z2SumPt = lep_i1.Pt()+lep_i2.Pt();
            }
            
            // Check isolation cut (without FSR ) for Z1 leptons
            if (debug) {cout << "RelIso NoFSR Z1: "<< (*lep_RelIsoNoFSR)[Z1_lepindex[0]] << " " << (*lep_RelIsoNoFSR)[Z1_lepindex[1]] << endl;}
            if ((*lep_RelIsoNoFSR)[Z1_lepindex[0]]>((abs((*lep_id)[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue;
            if ((*lep_RelIsoNoFSR)[Z1_lepindex[1]]>((abs((*lep_id)[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue;
            // Check tight ID cut for Z1 leptons
            if (debug) {cout << "Tight ID Z1: "<< (*lep_tightId)[Z1_lepindex[0]] << " " << (*lep_tightId)[Z1_lepindex[1]] << endl;}
            if (!((*lep_tightId)[Z1_lepindex[0]])) continue;
            if (!((*lep_tightId)[Z1_lepindex[1]])) continue;
            
            // Check Leading and Subleading pt Cut
            vector<double> allPt;
            allPt.push_back(lep_i1_nofsr.Pt()); allPt.push_back(lep_i2_nofsr.Pt());
            allPt.push_back(lep_j1_nofsr.Pt()); allPt.push_back(lep_j2_nofsr.Pt());
            std::sort(allPt.begin(), allPt.end());
            if (debug) cout<<" leading pt: "<<allPt[3]<<" cut: "<<leadingPtCut
                           <<" subleadingPt: "<<allPt[2]<<" cut: "<<subleadingPtCut<<endl;
            if (allPt[3]<leadingPtCut || allPt[2]<subleadingPtCut ) continue;
            
            // Check dR(li,lj)>0.02 for any i,j
            vector<double> alldR;
            alldR.push_back(lep_i1_nofsr.DeltaR(lep_i2_nofsr));
            alldR.push_back(lep_i1_nofsr.DeltaR(lep_j1_nofsr));
            alldR.push_back(lep_i1_nofsr.DeltaR(lep_j2_nofsr));
            alldR.push_back(lep_i2_nofsr.DeltaR(lep_j1_nofsr));
            alldR.push_back(lep_i2_nofsr.DeltaR(lep_j2_nofsr));
            alldR.push_back(lep_j1_nofsr.DeltaR(lep_j2_nofsr));
            if (debug) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
            if (*min_element(alldR.begin(),alldR.end())<0.02) continue;
            
            // Check M(l+,l-)>4.0 GeV for any OS pair
            vector<double> allM;
            TLorentzVector i1i2;
            i1i2 = (lep_i1_nofsr)+(lep_i2_nofsr); allM.push_back(i1i2.M());
            TLorentzVector j1j2;
            j1j2 = (lep_j1_nofsr)+(lep_j2_nofsr); allM.push_back(j1j2.M());

            if ((*lep_id)[i1]*(*lep_id)[j1]<0) {
                TLorentzVector i1j1;
                i1j1 = (lep_i1_nofsr)+(lep_j1_nofsr); allM.push_back(i1j1.M());
                TLorentzVector i2j2;
                i2j2 = (lep_i2_nofsr)+(lep_j2_nofsr); allM.push_back(i2j2.M());
            } else {
                TLorentzVector i1j2;
                i1j2 = (lep_i1_nofsr)+(lep_j2_nofsr); allM.push_back(i1j2.M());
                TLorentzVector i2j1;
                i2j1 = (lep_i2_nofsr)+(lep_j1_nofsr); allM.push_back(i2j1.M());
            }
            if (debug) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
            if (*min_element(allM.begin(),allM.end())<4.0) { continue;}
            // Do not include FSR photons
            //  if (*min_element(allM.begin(),allM.end())<0.1) { continue;}
            // Check the "smart cut": !( |mZa-mZ| < |mZ1-mZ| && mZb<12)
            // only for 4mu or 4e ZZ candidates
           // 
            bool passSmartCut=true;
            if ( abs((*lep_id)[i1])==abs((*lep_id)[j1])) {
                TLorentzVector Za, Zb;
                if ((*lep_id)[i1]==(*lep_id)[j1]) {                  
                    Za = (lep_i1)+(lep_j2);
                    Zb = (lep_i2)+(lep_j1);                    
                } else {
                    Za = (lep_i1)+(lep_j1);
                    Zb = (lep_i2)+(lep_j2);
                }                
                if ( abs(Za.M()-Zmass)<abs(Zb.M()-Zmass) ) {
                    if (debug) cout<<"abs(Za.M()-Zmass)-abs(Z1.M()-Zmass): "
                                   <<abs(Za.M()-Zmass)-abs(Z1.M()-Zmass)<<" Zb.M(): "<<Zb.M()<<endl;
                    if ( abs(Za.M()-Zmass)<abs(Z1.M()-Zmass) && Zb.M()<mZ2LowCut ) passSmartCut=false;
                }
                else {
                    if (debug) cout<<"abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass): "
                                   <<abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass)<<" Za.M(): "<<Za.M()<<endl;
                    if ( abs(Zb.M()-Zmass)<abs(Z1.M()-Zmass) && Za.M()<mZ2LowCut ) passSmartCut=false;
                }
                
            }
            
            if (!passSmartCut) continue; //Fix me 
            if (debug) cout<<" massZ1: "<<Z1.M()<<" massZ2: "<<Z2.M()<<endl;
            if (Z1.M() < mZ1LowCut) continue;
            if (Z1.M() > mZ1HighCut) continue;
            if (Z2.M() < mZ2LowCut) continue;
            if (Z2.M() > mZ2HighCut) continue; 
            if (debug) cout<<" pass Z mass cuts"<<endl;
            
            
            // Signal region if Z2 leptons are both tight ID Iso
            bool signalRegion=true;
            if ((*lep_RelIsoNoFSR)[Z2_lepindex[0]]>((abs((*lep_id)[Z2_lepindex[0]])==11) ? isoCutEl : isoCutMu)) signalRegion=false;
            if ((*lep_RelIsoNoFSR)[Z2_lepindex[1]]>((abs((*lep_id)[Z2_lepindex[1]])==11) ? isoCutEl : isoCutMu)) signalRegion=false;
            if (!((*lep_tightId)[Z2_lepindex[0]])) signalRegion=false; // checking tight lepton ID
            if (!((*lep_tightId)[Z2_lepindex[1]])) signalRegion=false; // checking tight lepton ID
            
            // Check if this candidate has the highest D_bkg_kin
            vector<TLorentzVector> P4s;
            P4s.clear();
            vector<int> tmpIDs;
            tmpIDs.clear();
            
            if (Z1_lepindex[0] == i1) {
                P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                if (Z2_lepindex[0] == j1) {
                    P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                } else {
                    P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                }
            } else if (Z1_lepindex[0] == i2) {
                P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                if (Z2_lepindex[0] == j1) {
                    P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                } else {
                    P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                }
            } else if (Z1_lepindex[0] == j1) {
                P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                if (Z2_lepindex[0] == i1) {
                    P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                } else {
                    P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                }
            } else if (Z1_lepindex[0] == j2) {
                P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                if (Z2_lepindex[0] == i1) {
                    P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                } else {
                    P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                }
            }
            
            tmpIDs.push_back((*lep_id)[Z1_lepindex[0]]); tmpIDs.push_back((*lep_id)[Z1_lepindex[1]]);
            tmpIDs.push_back((*lep_id)[Z2_lepindex[0]]); tmpIDs.push_back((*lep_id)[Z2_lepindex[1]]);

            bool same4l=false;
            bool foundZ11=false; bool foundZ12=false; bool foundZ21=false; bool foundZ22=false;
            for(int l = 0; l < 4; l++){
                if (lep_Hindex[l]==Z1_lepindex[0]) foundZ11 = true;
                if (lep_Hindex[l]==Z1_lepindex[1]) foundZ12 = true;
                if (lep_Hindex[l]==Z2_lepindex[0]) foundZ21 = true;
                if (lep_Hindex[l]==Z2_lepindex[1]) foundZ22 = true;
            }
            same4l = (foundZ11 && foundZ12 && foundZ21 && foundZ22);
            
            if (signalRegion) { // Signal Region has priority
                
                if (!foundSRCandidate) same4l=false;
                

                if (  Z1DeltaM<=minZ1DeltaM_SR ) { 
                    
                    minZ1DeltaM_SR = Z1DeltaM;
                    
                    if (Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_SR) continue;

                    Z_Hindex[0] = Z1index;
                    lep_Hindex[0] = Z1_lepindex[0];
                    lep_Hindex[1] = Z1_lepindex[1];
                
                    maxZ2SumPt_SR = Z2SumPt;    
                    Z_Hindex[1] = Z2index;
                    lep_Hindex[2] = Z2_lepindex[0];
                    lep_Hindex[3] = Z2_lepindex[1];
                    
                    Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                    massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();
                    
                    if (debug) cout<<" new best candidate SR: mass4l: "<<HVec.M()<<endl;
                    if ((HVec.M()>m4lLowCut)&&(HVec.M()<m4lHighCut))  {
                        foundHiggsCandidate=true;                    
                        foundSRCandidate=true;
                    }
                }
            } else if (!foundSRCandidate) { // Control regions get second priority

                if (  Z1DeltaM<=minZ1DeltaM_CR ) {                 

                    minZ1DeltaM_CR = Z1DeltaM;
                    
                    if (Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_CR) continue;

                    Z_Hindex[0] = Z1index;
                    lep_Hindex[0] = Z1_lepindex[0];
                    lep_Hindex[1] = Z1_lepindex[1];
                    
                    maxZ2SumPt_CR = Z2SumPt;
                    Z_Hindex[1] = Z2index;
                    lep_Hindex[2] = Z2_lepindex[0];
                    lep_Hindex[3] = Z2_lepindex[1];

                    Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                    massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();

                    if (debug) cout<<" new best candidate CR: mass4l: "<<HVec.M()<<endl;
                    if (HVec.M()>m4lLowCut&&HVec.M()<m4lHighCut) foundHiggsCandidate=true;                    
                }
            } 

            if (debug) cout<<"Z_Hindex[0]: "<<Z_Hindex[0]<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                           <<"Z_Hindex[1]: "<<Z_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;
            
        } // Zj
    } // Zi*/   

    /*TLorentzVector Lep1, Lep2, Lep3, Lep4,  Jet1, Jet2, Jet1_2p5, Jet2_2p5;            
    TLorentzVector nullFourVector(0, 0, 0, 0);                 
    Lep1.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[0]],(*lepFSR_eta)[lep_Hindex[0]],(*lepFSR_phi)[lep_Hindex[0]],(*lepFSR_mass)[lep_Hindex[0]]);
    Lep2.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[1]],(*lepFSR_eta)[lep_Hindex[1]],(*lepFSR_phi)[lep_Hindex[1]],(*lepFSR_mass)[lep_Hindex[1]]);
    Lep3.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[2]],(*lepFSR_eta)[lep_Hindex[2]],(*lepFSR_phi)[lep_Hindex[2]],(*lepFSR_mass)[lep_Hindex[2]]);
    Lep4.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[3]],(*lepFSR_eta)[lep_Hindex[3]],(*lepFSR_phi)[lep_Hindex[3]],(*lepFSR_mass)[lep_Hindex[3]]);*/

    /*TLorentzVector Lep1,Lep2,Lep3;
    Lep1.SetPtEtaPhiM(lep_pt[index1],lep_eta[index1],lep_phi[index1],lep_mass[index1]); 
    Lep2.SetPtEtaPhiM(lep_pt[index2],lep_eta[index2],lep_phi[index2],lep_mass[index2]);
    Lep3.SetPtEtaPhiM(lep_pt[index3],lep_eta[index3],lep_phi[index3],lep_mass[index3]);*/

    /*int jet1index=-1, jet2index=-1;
    float jet1pt=0.0, jet2pt=0.0;
    int jet1index2p5=-1, jet2index2p5=-1;
    float jet1pt2p5=0.0, jet2pt2p5=0.0;*/

    /*nZXCRFailedLeptons = 0;
    for(unsigned int i = 0; i <= 3; i++) {
        if (!(abs((*lep_id)[lep_Hindex[i]])==11 && ((*lep_tightId)[lep_Hindex[i]] && (*lep_RelIsoNoFSR)[lep_Hindex[i]]<isoCutEl)) &&
            !(abs((*lep_id)[lep_Hindex[i]])==13 && ((*lep_tightId)[lep_Hindex[i]] && (*lep_RelIsoNoFSR)[lep_Hindex[i]]<isoCutMu))){
            nZXCRFailedLeptons++;
        }
    }*/

    //jet1index=(*jet_iscleanH4l)[0]; jet2index=(*jet_iscleanH4l)[1];

    idL1 = (*lep_id)[index1]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta();
    idL2 = (*lep_id)[index2]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta();
    idL3 = (*lep_id)[index3]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta();       
    //idL4 = (*lep_id)[lep_Hindex[3]]; pTL4 = Lep4.Pt(); etaL4 = Lep4.Eta();
    phiL1 = Lep1.Phi();//deltaphiL13 = deltaPhi (Lep1.Phi(), Lep3.Phi()); 
    phiL2 = Lep2.Phi();//deltaphiL14 = deltaPhi (Lep1.Phi(), Lep4.Phi());
    phiL3 = Lep3.Phi();//deltaphiL23 = deltaPhi (Lep2.Phi(), Lep3.Phi());
    //phiL4 = Lep4.Phi();//deltaphiL24 = deltaPhi (Lep2.Phi(), Lep4.Phi());
    //deltaphiZZ = deltaPhi((Lep1+Lep2).Phi(), (Lep3+Lep4).Phi());
    vector<TLorentzVector> P4s; vector<int> tmpIDs;             
    P4s.push_back(Lep1); P4s.push_back(Lep2);
    P4s.push_back(Lep3); //P4s.push_back(Lep4);
    tmpIDs.push_back(idL1); tmpIDs.push_back(idL2);
    tmpIDs.push_back(idL3); //tmpIDs.push_back(idL4);
    IsoL1 = (*lep_RelIso)[index1]; IsoL2 = (*lep_RelIso)[index2];
    IsoL3 = (*lep_RelIso)[index3];
    //lep_Hindex_stdvec->push_back(lep_Hindex[0]);
    //lep_Hindex_stdvec->push_back(lep_Hindex[1]);
    //lep_Hindex_stdvec->push_back(lep_Hindex[2]);
    //lep_Hindex_stdvec->push_back(lep_Hindex[3]);

    //TLorentzVector threeleps = Lep1+Lep2+Lep3;

    massZ1 = Z.M();
    pT3l = Leps.Pt();
            
    //cout<<"fill tree"<<endl;
    //cout<<endl;
    outTree->Fill();   
return 1;
}

#endif
