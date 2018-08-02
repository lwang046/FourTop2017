#include "TROOT.h"
#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TNtuple.h>
#include <TLorentzVector.h>


#include <vector>
#include <iostream>

/*
bool CheckValue(ROOT::TTreeReaderArray& value) {
   if (value->GetSetupStatus() < 0) {
      std::cerr << "Error " << value->GetSetupStatus()
                << "setting up reader for " << value->GetBranchName() << '\n';
      return false;
   }
   return true;
}
*/
using namespace std;

int main() {
   TFile *file = TFile::Open("DoubleMuon_Run2017B_NANOAOD.root");
   TTreeReader reader("Events", file);

   TTreeReaderValue<unsigned int> nMuon(reader,"nMuon");
   TTreeReaderArray<float> Muon_pt(reader,"Muon_pt");
   TTreeReaderArray<float> Muon_eta(reader,"Muon_eta");
   TTreeReaderArray<float> Muon_phi(reader,"Muon_phi");
   TTreeReaderArray<float> Muon_mass(reader,"Muon_mass");
   TTreeReaderArray<int> Muon_charge(reader,"Muon_charge");
   //TTreeReaderArray<bool> Muon_softID(reader,"Muon_softId");
   TTreeReaderArray<float> Muon_RelIso(reader,"Muon_pfRelIso04_all");
   TTreeReaderArray<int> Muon_jetIdx(reader,"Muon_jetIdx");
   TTreeReaderArray<int> Muon_pdfId(reader,"Muon_pdgId");

   TTreeReaderValue<unsigned int> nElectron(reader,"nElectron");
   TTreeReaderArray<float> Electron_pt(reader,"Electron_pt");
   TTreeReaderArray<float> Electron_eta(reader,"Electron_eta");
   TTreeReaderArray<float> Electron_phi(reader,"Electron_phi");
   TTreeReaderArray<float> Electron_mass(reader,"Electron_mass");
   TTreeReaderArray<int> Electron_charge(reader,"Electron_charge");
   TTreeReaderArray<int> Electron_cutBased(reader,"Electron_cutBased");
   TTreeReaderArray<float> Electron_RelIso(reader,"Electron_pfRelIso03_all");
   TTreeReaderArray<int> Electron_jetIdx(reader,"Electron_jetIdx");
   TTreeReaderArray<int> Electron_pdfId(reader,"Electron_pdgId");

   TTreeReaderValue<unsigned int> nJet(reader,"nJet");
   TTreeReaderArray<float> Jet_pt(reader,"Jet_pt");
   TTreeReaderArray<float> Jet_eta(reader,"Jet_eta");
   TTreeReaderArray<float> Jet_phi(reader,"Jet_phi");
   TTreeReaderArray<float> Jet_mass(reader,"Jet_mass");
   TTreeReaderArray<int> Jet_jetId(reader,"Jet_jetId");
   TTreeReaderArray<float> Jet_btagCSVV2(reader,"Jet_btagCSVV2");
   TTreeReaderArray<float> Jet_btagDeepB(reader,"Jet_btagDeepB");
   TTreeReaderArray<unsigned char> Jet_cleanmask(reader,"Jet_cleanmask");
   TTreeReaderArray<int> Jet_nMuons(reader,"Jet_nMuons");
   TTreeReaderArray<int> Jet_nElectrons(reader,"Jet_nElectrons");

 
   //if(!CheckValue(MuonPt)) return false;


   TFile *fout = new TFile("histo.root","RECREATE");
 
   TNtuple *tup = new TNtuple("Craneen__MuMu", "Craneen__MuMu", 
		"nLeps:Lep1Pt:Lep2Pt:Lep1Eta:Lep2Eta:Lep1Phi:Lep2Phi:Lep1Charge:Lep2Charge:Lep1RelIso:Lep2RelIso:Lep1Id:Lep2Id:Lep1jetIdx:Lep2jetIdx:Lep1pdfId:Lep2pdfId:" //17
		"nJets:nMtags:1stJetPt:2ndJetPt:3rdJetPt:4thJetPt:1stbtagPt:2ndbtagPt:" //8
                "diLepMass" //1

);

   unsigned long long Nevent = 0;
   while(reader.Next()) {
      Nevent++; 
      if(Nevent % 100000 == 0) cout << "Processing the " << Nevent << "th event" << flush << "\r" << endl;
 
      float Lep1Pt = -99., Lep1Eta = -99., Lep1Phi = -99., Lep1Mass = -99., Lep1Charge = -99., Lep1RelIso = -99., Lep1Id = -99., Lep1jetIdx = -99., Lep1pdfId = -99.,
            Lep2Pt = -99., Lep2Eta = -99., Lep2Phi = -99., Lep2Mass = -99., Lep2Charge = -99., Lep2RelIso = -99., Lep2Id = -99., Lep2jetIdx = -99., Lep2pdfId = -99.;
      int nMu = 0, nEl = 0, nLeps = 0;
      float nJets = -1., nMtags = -1, Jet1Pt = -99., Jet2Pt = -99., Jet3Pt = -99., Jet4Pt = -99., Tag1Pt = -99., Tag2Pt = -99.;

      TLorentzVector lep1, lep2, diLep, jets;
      vector<TLorentzVector> jets_TLV, btags_TLV;

      nMu = *nMuon; nEl = *nElectron; nLeps = nMu+nEl;

      if(!(nMu==2 && nEl==0 && *nJet>=4)) continue;
      Lep1Pt = Muon_pt[0]; Lep1Eta = Muon_eta[0]; Lep1Phi = Muon_phi[0]; Lep1Charge = Muon_charge[0]; Lep1RelIso = Muon_RelIso[0]; Lep1Mass = Muon_mass[0]; Lep1jetIdx = Muon_jetIdx[0]; Lep1pdfId = Muon_pdfId[0];
      Lep2Pt = Muon_pt[1]; Lep2Eta = Muon_eta[1]; Lep2Phi = Muon_phi[1]; Lep2Charge = Muon_charge[1]; Lep2RelIso = Muon_RelIso[1]; Lep2Mass = Muon_mass[1]; Lep2jetIdx = Muon_jetIdx[1]; Lep2pdfId = Muon_pdfId[1];
      lep1.SetPtEtaPhiM(Lep1Pt, Lep1Eta, Lep1Phi, Lep1Mass); lep2.SetPtEtaPhiM(Lep2Pt, Lep2Eta, Lep2Phi, Lep2Mass);
      diLep = lep1 + lep2;

      if(! ( Lep1Pt>=25 && Lep2Pt>=25 && abs(Lep1Eta)<2.4 && abs(Lep2Eta)<2.4 && Lep1RelIso<0.15 && Lep2RelIso<0.15)) continue;
      if(Lep1Charge==Lep2Charge) continue;
      
      for(int nj=0; nj<*nJet; nj++){
         if(Jet_pt[nj]<30 || abs(Jet_eta[nj])>2.4 || Jet_cleanmask[nj]!=1) continue;
         jets.SetPtEtaPhiM(Jet_pt[nj], Jet_eta[nj], Jet_phi[nj], Jet_mass[nj]);
         jets_TLV.push_back(jets);
         if(Jet_btagDeepB[nj]>0.4941) btags_TLV.push_back(jets);
         //cout << Jet_cleanmask[nj]<< endl;
      }
      nJets = jets_TLV.size(); nMtags = btags_TLV.size();
      if(! (nJets>=4 && nMtags>=2) ) continue;










      float vals[26] = {(float)nLeps, Lep1Pt, Lep2Pt, Lep1Eta, Lep2Eta, Lep1Phi, Lep2Phi, Lep1Charge, Lep2Charge, Lep1RelIso, Lep2RelIso, Lep1Id, Lep2Id, Lep1jetIdx, Lep2jetIdx, Lep1pdfId, Lep2pdfId, nJets, nMtags, (float)jets_TLV[0].Pt(), (float)jets_TLV[1].Pt(), (float)jets_TLV[2].Pt(), (float)jets_TLV[3].Pt(), (float)btags_TLV[0].Pt(), (float)btags_TLV[1].Pt(), (float)diLep.M()};
      tup->Fill(vals);
   }





   fout->cd();
   tup->Write();
   fout->Close();

}


