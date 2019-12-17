#ifndef Wto3l_mem_controlregion_h
#define Wto3l_mem_controlregion_h

#include "Wto3l_mem_controlregion_Linkdef.h"
#include "Analyzer.h"
#include "HZZTree_fakerate.h"
#include "deltaR.h"
#include <cmath>

using namespace std;  

class Wto3l_mem_controlregion : public Analyzer 
{
    public:
        // static const double Zmass = 91.1876;
        const double Zmass = 91.1876;

        Wto3l_mem_controlregion();
        ~Wto3l_mem_controlregion();
        Wto3l_mem_controlregion(
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
         Wto3l_mem_controlregion(
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

Wto3l_mem_controlregion::Wto3l_mem_controlregion(){

}

Wto3l_mem_controlregion::~Wto3l_mem_controlregion(){

}

bool Wto3l_mem_controlregion::passSelection(){
    return true;
}


void Wto3l_mem_controlregion::initTree(){
    setHZZTree(tree);
}

void Wto3l_mem_controlregion::setDebugMode(bool debug_in){
    debug = debug_in;
}

void Wto3l_mem_controlregion::setup(){
    outFile = TFile::Open(outputDir+outFileName,fOptionWrite);
    outTree = new TTree(outTreeName,outTreeName);

    initNewLiteTree_fakerate(outTree);
}

void Wto3l_mem_controlregion::end(){
    outFile->cd();
    outTree->Write(outTreeName,TObject::kOverwrite);
    outFile->Close(); 
}

Wto3l_mem_controlregion::Wto3l_mem_controlregion(
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

Wto3l_mem_controlregion::Wto3l_mem_controlregion(
                TString outputDir_in,
                TString outFileName_in
                ){
    outputDir   = outputDir_in;
    outFileName = outFileName_in;
}

void Wto3l_mem_controlregion::sortedArray(double x, double y, double z, double sortarray[3])
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

int Wto3l_mem_controlregion::process(){

    int Nlep = (*lep_id).size();
    int nTightLep = 0;
    int nLooseLep = 0;
    vector<int> tightIsoLepIndex;
    vector<int> looseIsoLepIndex;

    if((*lep_id).size()<3 || (*lep_pt).size()<3) return -1;

    if(passedTrig == 0) return -1; 

    for (int iLep = 0; iLep < Nlep; iLep++) {
        if ((*lep_tightId)[iLep] == 1 && (*lep_RelIso)[iLep] < 0.35 && nTightLep < 2){
            nTightLep++;
            tightIsoLepIndex.push_back(iLep);
        }
        else if ((*lep_tightId)[iLep] == 1 && abs((*lep_id)[iLep]) == 13){
            nLooseLep++;
            looseIsoLepIndex.push_back(iLep);
        }
    }

    if (!(nTightLep == 2 && nLooseLep == 1)) return -1;

    int index1 = tightIsoLepIndex[0];
    int index2 = tightIsoLepIndex[1];
    int index3 = looseIsoLepIndex[0];

    //if ( (*lep_id)[index1] + (*lep_id)[index2] != 0 ) return -1;
    if ( ((*lep_id)[index1] > 0 && (*lep_id)[index2] > 0) || ((*lep_id)[index1] < 0 && (*lep_id)[index2] < 0) ) return -1;

    if ((*lep_Sip)[index1] > 3 || (*lep_Sip)[index2] > 3 || (*lep_Sip)[index3] > 3) return -1;

    double pTs[3]; sortedArray((*lep_pt)[index1], (*lep_pt)[index2], (*lep_pt)[index3], pTs);
    if (pTs[0] < 20 || pTs[1] < 10 || pTs[2] < 5) return -1;

    double isos[3]; sortedArray((*lep_RelIso)[index1], (*lep_RelIso)[index2], (*lep_RelIso)[index3], isos);
    double sips[3]; sortedArray((*lep_Sip)[index1], (*lep_Sip)[index2], (*lep_Sip)[index3], sips);

    int tmp = 0;
    TLorentzVector Lep1,Lep2,Lep3;
    Lep3.SetPtEtaPhiM((*lep_pt)[index3],(*lep_eta)[index3],(*lep_phi)[index3],(*lep_mass)[index3]);
    if ((*lep_pt)[index1] >= (*lep_pt)[index2]){
        Lep1.SetPtEtaPhiM((*lep_pt)[index1],(*lep_eta)[index1],(*lep_phi)[index1],(*lep_mass)[index1]); 
        Lep2.SetPtEtaPhiM((*lep_pt)[index2],(*lep_eta)[index2],(*lep_phi)[index2],(*lep_mass)[index2]);
    }
    if ((*lep_pt)[index1] < (*lep_pt)[index2]){
        tmp = index1;
        index1 = index2;
        index2 = tmp;
        Lep1.SetPtEtaPhiM((*lep_pt)[index1],(*lep_eta)[index1],(*lep_phi)[index1],(*lep_mass)[index1]); 
        Lep2.SetPtEtaPhiM((*lep_pt)[index2],(*lep_eta)[index2],(*lep_phi)[index2],(*lep_mass)[index2]);
    }

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
    //if (!(Z.M() > 86 && Z.M() < 96)) return kTRUE;

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
    massL1 = (*lep_mass)[index1]; massL1 = (*lep_mass)[index2];
    massL3 = (*lep_mass)[index3];
    MomIdL1 = (*lep_matchedR03_MomId)[index1];  MomIdL2 = (*lep_matchedR03_MomId)[index2];  MomIdL3 = (*lep_matchedR03_MomId)[index3];
    PDG_IdL1 = (*lep_matchedR03_PdgId)[index1];  PDG_IdL2 = (*lep_matchedR03_PdgId)[index2];  PDG_IdL3 = (*lep_matchedR03_PdgId)[index3];
    MomMomIdL1 = (*lep_matchedR03_MomMomId)[index1];  MomMomIdL2 = (*lep_matchedR03_MomMomId)[index2];  MomMomIdL3 = (*lep_matchedR03_MomMomId)[index3];
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
