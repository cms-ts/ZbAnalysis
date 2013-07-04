// -*- C++ -*-
//
// Package: GenbAnalyzer
// Class: GenbAnalyzer
//
/**\class GenbAnalyzer GenbAnalyzer.cc ZbAnalysis/GenbAnalyzer/src/GenbAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]

*/
//
// Original Author: Vieri Candelise
// Created: Thu Jan 10 15:57:03 CET 2013
// $Id: GenbAnalyzer.cc,v 1.28 2013/07/02 13:40:44 dellaric Exp $
//
//

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stddef.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/View.h"
#include "Math/VectorUtil.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "PhysicsTools/CandUtils/interface/CandCombiner.h"
#include "CommonTools/Utils/interface/MassRangeSelector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"

//
// class declaration
//

class GenbAnalyzer:public  edm::EDAnalyzer {

public:

  explicit GenbAnalyzer (const edm::ParameterSet &);
  ~GenbAnalyzer ();

private:

  virtual void beginJob ();
  virtual void analyze (const edm::Event &, const edm::EventSetup &);
  virtual void endJob ();

  virtual void beginRun (edm::Run const &, edm::EventSetup const &);
  virtual void endRun (edm::Run const &, edm::EventSetup const &);
  virtual void beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);
  virtual void endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);

  struct order_leptons { bool operator() (const TLorentzVector &p1, const TLorentzVector &p2) const {
      return (p1.Pt() > p2.Pt());
    }
  };
  struct order_jets { bool operator() (const pat::Jet &j1, const pat::Jet &j2) const {
      return (j1.pt() > j2.pt());
    }
  };

  void fill(TH1F* histogram, double value, double weight=1.0) {
    TAxis* axis = histogram->GetXaxis();
    Int_t nx = histogram->GetNbinsX();
    if (axis->FindBin(value) <= 0) {
      histogram->Fill(histogram->GetBinCenter(1), weight);
    } else if (axis->FindBin(value) >= nx+1) {
      histogram->Fill(histogram->GetBinCenter(nx), weight);
    } else {
      histogram->Fill(value, weight);
    }
  };

  // ----------member data ---------------------------

  std::string pileup_;
  std::string lepton_;
  std::string path_;

  edm::LumiReWeighting LumiWeights_;

  TH1F*     h_gen_weights;

  TH1F*     h_jetmultiplicity;
  TH1F*     h_jet_pt;
  TH1F*     h_ele_pt;
  TH1F*     h_muon_pt;
  TH1F*     w_pt_Z_ee;
  TH1F*     w_pt_Z_mm;
  TH1F*     w_pt_Z_ee_b;
  TH1F*     w_pt_Z_mm_b;

  TH1F*     w_jetmultiplicity;
  TH1F*     w_first_jet_pt;      // leading jet of any type
  TH1F*     w_first_jet_eta;
  TH1F*     w_first_bjet_pt;
  TH1F*     w_first_bjet_eta;
  TH1F*     w_bjetmultiplicity;
  TH1F*     w_first_jet_pt_b;   // leading jet with at least one b jet in the event
  TH1F*     w_first_jet_eta_b;
  TH1F*     w_first_ele_pt;
  TH1F*     w_second_ele_pt;
  TH1F*     w_dressed_ele_pt;
  TH1F*     w_first_muon_pt;
  TH1F*     w_mass_ee;
  TH1F*     w_mass_mm;
  TH1F*     w_mass_ee_b;  // at least one b jet in the event
  TH1F*     w_mass_mm_b;
  TH1F*     w_delta_ee_b;
  TH1F*     w_delta_mm_b;
  TH1F*     w_delta_ee;
  TH1F*     w_delta_mm;
  TH1F*     w_Ht;
  TH1F*     w_Ht_b;

};

using namespace  pat;

//
// constants, enums and typedefs
//


//
// static data member definitions
//


//
// constructors and destructor
//
GenbAnalyzer::GenbAnalyzer (const edm::ParameterSet & iConfig) {

  pileup_ = iConfig.getUntrackedParameter < std::string > ("pileup", "S7");
  lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
  path_ =   iConfig.getUntrackedParameter < std::string > ("path", "/gpfs/cms/users/candelis/work/ZbSkim/test");

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  h_gen_weights     =   fs->make < TH1F > ("h_gen_weights",      "h_gen_weights", 2, 0, 2);

  w_jetmultiplicity =   fs->make < TH1F > ("w_jetmultiplicity",  "w_jetmultiplicity;", 8, 0.5, 8.5);
  w_first_ele_pt =      fs->make < TH1F > ("w_first_ele_pt",     "w_first_ele_pt; [GeV]", 50, 0., 450.);
  w_second_ele_pt =     fs->make < TH1F > ("w_second_ele_pt",    "w_second_ele_pt; [GeV]", 50, 0., 450.);
  w_dressed_ele_pt =    fs->make < TH1F > ("w_dressed_ele_pt",   "w_dressed_ele_pt; [GeV]", 50, 0., 450.);
  w_first_muon_pt =     fs->make < TH1F > ("w_first_muon_pt",    "w_first_muon_pt; [GeV]", 50, 0., 450.);
  w_bjetmultiplicity =  fs->make < TH1F > ("w_bjetmultiplicity", "w_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  w_first_bjet_pt =     fs->make < TH1F > ("w_first_bjet_pt",    "w_first_bjet_pt; [GeV]", 50, 30., 700.);
  w_first_jet_pt =      fs->make < TH1F > ("w_first_jet_pt",     "w_first_jet_pt; [GeV]", 50, 30., 700.);
  w_first_jet_eta =     fs->make < TH1F > ("w_first_jet_eta",    "w_first_jet_eta", 16, -2.5, 2.5);
  w_first_bjet_eta =    fs->make < TH1F > ("w_first_bjet_eta",   "w_first_bjet_eta", 16, -2.5, 2.5);
  w_first_jet_pt_b =    fs->make < TH1F > ("w_first_jet_pt_b",   "w_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  w_first_jet_eta_b =   fs->make < TH1F > ("w_first_jet_eta_b",  "w_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  w_pt_Z_ee =           fs->make < TH1F > ("w_pt_Z_ee",          "w_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_mm =           fs->make < TH1F > ("w_pt_Z_mm",          "w_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_ee_b =         fs->make < TH1F > ("w_pt_Z_ee_b",        "w_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_mm_b =         fs->make < TH1F > ("w_pt_Z_mm_b",        "w_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  w_mass_ee =           fs->make < TH1F > ("w_mass_ee",          "w_mass_ee;Mass [GeV]", 80, 71., 111.);
  w_mass_mm =           fs->make < TH1F > ("w_mass_mm",          "w_mass_mm;Mass [GeV]", 80, 71., 111.);
  w_mass_ee_b =         fs->make < TH1F > ("w_mass_ee_b",        "w_mass_ee_b;Mass [GeV]", 80, 71., 111.);
  w_mass_mm_b =         fs->make < TH1F > ("w_mass_mm_b",        "w_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  w_Ht =                fs->make < TH1F > ("w_Ht",               "w_Ht [GeV]", 50, 30., 1000.);
  w_Ht_b =              fs->make < TH1F > ("w_Ht_b",             "w_Ht_b [GeV]", 50, 30., 1000.);
  w_delta_mm =          fs->make < TH1F > ("w_delta_phi_mm",     "w_delta_phi_mm", 12, 0, TMath::Pi ());
  w_delta_mm_b =        fs->make < TH1F > ("w_delta_phi_mm_b",   "w_delta_phi_mm_b", 12, 0, TMath::Pi ());
  w_delta_ee =          fs->make < TH1F > ("w_delta_phi_ee",     "w_delta_phi_ee", 12, 0, TMath::Pi ());
  w_delta_ee_b =        fs->make < TH1F > ("w_delta_phi_ee_b",   "w_delta_phi_ee_b", 12, 0, TMath::Pi ());

}

GenbAnalyzer::~GenbAnalyzer () {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event ------------
void GenbAnalyzer::analyze (const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace edm;
  using namespace std;

  edm::Handle<vector<reco::GenParticle> > genPart;
  iEvent.getByLabel ("genParticles", genPart);

  edm::Handle<vector<reco::GenJet> > gJets;
  iEvent.getByLabel(edm::InputTag("ak5GenJetsNoNu",""), gJets);
  if (! gJets.isValid()) {
    iEvent.getByLabel(edm::InputTag("selectedPatJetsPFlow","genJets"), gJets);
  }

  bool ee_event = false;
  bool mm_event = false;

  int Nj = 0;
  int Nb = 0;

  double diele_mass = 0;
  double diele_phi = 0;
  double diele_pt = 0;

  double dimuon_mass = 0;
  double dimuon_phi = 0;
  double dimuon_pt = 0;

  double MyWeight = 1;

  double Ht = 0;

  // ++++++ Pile-Up

  MyWeight = 1.0;

  Handle < vector < PileupSummaryInfo > > PupInfo;

  if (iEvent.getByLabel (edm::InputTag ("addPileupInfo"), PupInfo))  {

    float Tnpv = -1;

    for (vector < PileupSummaryInfo >::const_iterator PVI = PupInfo->begin (); PVI != PupInfo->end (); ++PVI) {
      int BX = PVI->getBunchCrossing ();
      if (BX == 0) {
        Tnpv = PVI->getTrueNumInteractions ();
        continue;
      }
    }

    MyWeight = LumiWeights_.weight (Tnpv);

  }

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;

  if (iEvent.getByLabel ("generator", genEventInfoHandle)) {

    double mcWeight = genEventInfoHandle->weight();

    h_gen_weights->Fill(0.5, 1.0);
    h_gen_weights->Fill(1.5, mcWeight);

    MyWeight = MyWeight*mcWeight;

  }

  // +++++++++ ELECTRONS

  vector < TLorentzVector > vect_ele;

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {

    if (fabs(itgen->pdgId())==11 && itgen->status()==1) { // loop over gen electrons
      TLorentzVector ele;
      ele.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());

      // Loop over photons: FSR dressing for electrons
      for (vector<reco::GenParticle>::const_iterator itgen2=genPart->begin(); itgen2!=genPart->end(); itgen2++) {
        if (fabs(itgen2->pdgId())==22 && itgen2->status()==1) { // loop over primary gen photon
       	  TLorentzVector gam;
          gam.SetPtEtaPhiM(itgen2->pt(),itgen2->eta(),itgen2->phi(),itgen2->mass());
	  double deltaPhi_eg = fabs(gam.Phi()- ele.Phi());
	  if (deltaPhi_eg > acos(-1)) deltaPhi_eg= 2*acos(-1) - deltaPhi_eg;
	  double deltaR_eg = sqrt( deltaPhi_eg*deltaPhi_eg + pow(fabs(gam.Eta()-ele.Eta()),2));
	  if (deltaR_eg< 0.1) {
	    ele += gam;
	  }
	}
      }
      if (ele.Pt()>20 && fabs(ele.Eta())<2.4) {
        ele.SetE(itgen->charge() * ele.E());
        vect_ele.push_back(ele);
      }
    }
  }

  // Sort electrons in pT
  std::sort( vect_ele.begin(), vect_ele.end(), order_leptons() );

  int iele0=0;
  int iele1=0;

  for (unsigned int i=1; i<vect_ele.size(); ++i) {
    if ((vect_ele[i].E()*vect_ele[iele0].E())<0 && iele1==0) iele1=i;
    vect_ele[i].SetE(fabs(vect_ele[i].E()));
  }
  if (vect_ele.size()!=0) vect_ele[iele0].SetE(fabs(vect_ele[iele0].E()));

  if (iele1 != 0) {
    TLorentzVector y;
    y = vect_ele[iele0] + vect_ele[iele1];
    diele_mass = y.M();
    diele_pt =   y.Pt();
    diele_phi =  y.Phi();
    if (diele_mass>71 && diele_mass<111) {
      ee_event = true;
    }
  }

  // +++++++++ MUONS

  vector < TLorentzVector > vect_muon;

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
    if (fabs(itgen->pdgId())==13 && itgen->status()==1) { // loop over gen muons
      TLorentzVector muon;
      muon.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());

      // Loop over photons: FSR dressing for muons
      for (vector<reco::GenParticle>::const_iterator itgen2=genPart->begin(); itgen2!=genPart->end(); itgen2++) {
        if (fabs(itgen2->pdgId())==22 && itgen2->status()==1) { // loop over primary gen photon
	  TLorentzVector gam;
	  gam.SetPtEtaPhiM(itgen2->pt(),itgen2->eta(),itgen2->phi(),itgen2->mass());
	  double deltaPhi_mg = fabs(gam.Phi()- muon.Phi());
	  if (deltaPhi_mg > acos(-1)) deltaPhi_mg= 2*acos(-1) - deltaPhi_mg;
	  double deltaR_mg = sqrt( deltaPhi_mg*deltaPhi_mg + pow(fabs(gam.Eta()-muon.Eta()),2));
	  if (deltaR_mg< 0.1) {
	    muon += gam;
	  }
        }
      }

      if (muon.Pt()>20 && fabs(muon.Eta())<2.4) {
        muon.SetE(itgen->charge() * muon.E());
        vect_muon.push_back(muon);
      }
    }
  }

  // Sort muons in pT
  std::sort( vect_muon.begin(), vect_muon.end(), order_leptons() );

  int imuon0=0;
  int imuon1=0;

  for (unsigned int i=1; i<vect_muon.size(); ++i) {
    if (vect_muon[i].E()*vect_muon[imuon0].E()<0 && imuon1==0) imuon1=i;
    vect_muon[i].SetE(fabs(vect_muon[i].E()));
  }
  if (vect_muon.size()!=0) vect_muon[imuon0].SetE(fabs(vect_muon[imuon0].E()));

  if (imuon1 != 0) {

    TLorentzVector y;
    y = vect_muon[imuon0] + vect_muon[imuon1];
    dimuon_mass = y.M();
    dimuon_pt = y.Pt();
    dimuon_phi = y.Phi();
    if (dimuon_mass>71 && dimuon_mass<111) {
      mm_event = true;
    }
  }

  ee_event = ee_event && (lepton_ == "electron");
  mm_event = mm_event && (lepton_ == "muon");

  // ++++++++ JETS

  vector < pat::Jet > vect_jets;

  for (vector<reco::GenJet>::const_iterator jet=gJets->begin(); jet!=gJets->end(); ++jet) {

    double jet_pt  = jet->pt ();
    double jet_eta = jet->eta();

    /*Lepton removal from jets inside sqrt Delta2+Dphi2 <0.3*/

    double deltaPhi1 = 0;
    if (ee_event) deltaPhi1 = fabs(jet->phi()- vect_ele[iele0]. Phi());
    if (mm_event) deltaPhi1 = fabs(jet->phi()- vect_muon[imuon0].Phi());

    if (deltaPhi1 > acos(-1)) deltaPhi1= 2*acos(-1) - deltaPhi1;

    double deltaR1 = 0;
    if (ee_event) deltaR1= sqrt( deltaPhi1*deltaPhi1  + pow(jet->eta()-vect_ele[iele0]. Eta(),2) );
    if (mm_event) deltaR1= sqrt( deltaPhi1*deltaPhi1  + pow(jet->eta()-vect_muon[imuon0].Eta(),2) );

    double deltaPhi2 = 0;
    if (ee_event) deltaPhi2 = fabs(jet->phi()-vect_ele[iele1]. Phi());
    if (mm_event) deltaPhi2 = fabs(jet->phi()-vect_muon[imuon1].Phi());

    if (deltaPhi2 > acos(-1)) deltaPhi2= 2*acos(-1) - deltaPhi2;

    double deltaR2 = 0;
    if (ee_event) deltaR2= sqrt( deltaPhi2*deltaPhi2  + pow(jet->eta()-vect_ele[iele1]. Eta(),2) );
    if (mm_event) deltaR2= sqrt( deltaPhi2*deltaPhi2  + pow(jet->eta()-vect_muon[imuon1].Eta(),2) );

    if (jet_pt > 30 && fabs(jet_eta) < 2.5 && (deltaR1>0.1 && deltaR2>0.1)) {

      ++Nj;

      Ht += jet_pt;

      vect_jets.push_back(*jet);

     }

  }

  // Sort jets in pT
  std::sort( vect_jets.begin(), vect_jets.end(), order_jets() );

  /*loop over gen particles, find the b*/

  vector < pat::Jet > vect_bjets;

  int nb=0;

  for (std::vector <reco::GenParticle>::const_iterator thepart =genPart->begin(); thepart != genPart->end(); thepart++) {

    if ((int) ((abs(thepart->pdgId())/100)%10) == 5 || (int) ((abs(thepart->pdgId())/1000)%10) == 5 ) {
      nb++;
      bool bdaughter = false;
      for (int i=0; i < abs(thepart->numberOfDaughters()); i++) {
        if ((int) ((abs(thepart->daughter(i)->pdgId())/100)%10) == 5 || (int) ((abs(thepart->daughter(i)->pdgId())/1000)%10) == 5) {
          bdaughter = true; // b daughter found
        }
      }
      if (!bdaughter) {
        TLorentzVector B;
        B.SetPtEtaPhiM(thepart->pt(),thepart->eta(),thepart->phi(),thepart->mass());
	int j = 0;
     	double Rmin = 9999.;
        for (unsigned int i=0; i<vect_jets.size(); ++i) {
          if (ROOT::Math::VectorUtil::DeltaR(vect_jets[i].momentum(), B) < Rmin) {
	    j = i;
	    Rmin = ROOT::Math::VectorUtil::DeltaR(vect_jets[i].momentum(), B);
	  }
        }
        if (Rmin < 0.4) {
	  Nb++;
	  vect_bjets.push_back(vect_jets[j]);
	}
      }
    }
  }

  // Sort b-jets in pT
  std::sort( vect_bjets.begin(), vect_bjets.end(), order_jets() );

  //cout<<"gen jet: "<<vect_bjets[0].pt()<<endl;

  /*
  // Get reco jet collection and print the reco b jets for check
  edm::Handle < vector < pat::Jet > > jets;
  iEvent.getByLabel ("goodJets", jets);

  for (vector < pat::Jet >::const_iterator jet = jets->begin (); jet != jets->end (); ++jet) {
    if (fabs(jet->partonFlavour ()) == 5) {
      cout<<"reco jet: "<<jet->momentum()<<endl;
    }
  }
  */

  if (ee_event && Nj > 0) {
    w_first_ele_pt->Fill (vect_ele[iele0].Pt(), MyWeight);
    w_second_ele_pt->Fill (vect_ele[iele1].Pt(), MyWeight);
    w_mass_ee->Fill(diele_mass, MyWeight);
    w_pt_Z_ee->Fill(diele_pt, MyWeight);
    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi());
    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
    w_delta_ee->Fill(delta_phi_ee, MyWeight);
  }

  if (mm_event && Nj > 0) {
    w_first_muon_pt->Fill (vect_muon[imuon0].Pt(), MyWeight);
    w_mass_mm->Fill(dimuon_mass, MyWeight);
    w_pt_Z_mm->Fill(dimuon_pt, MyWeight);
    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi());
    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
    w_delta_mm->Fill(delta_phi_mm, MyWeight);
  }

  if ((ee_event || mm_event) && Nj > 0) {
    w_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight);
    w_jetmultiplicity->Fill (Nj, MyWeight);
    w_Ht->Fill (Ht, MyWeight);
  }

  if (mm_event && Nb>0 && Nj>0) {
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight);
    w_pt_Z_mm_b->Fill (dimuon_pt, MyWeight);
    w_mass_mm_b->Fill (dimuon_mass, MyWeight);
  }
  if (ee_event && Nb>0 && Nj>0) {
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight);
    w_pt_Z_ee_b->Fill (diele_pt, MyWeight);
    w_mass_ee_b->Fill (diele_mass, MyWeight);
  }

  if ((ee_event || mm_event) && Nb>0 && Nj>0) {
    w_bjetmultiplicity->Fill (Nb, MyWeight);
    w_Ht_b->Fill (Ht, MyWeight);
  }
  if (ee_event && Nj > 0 && Nb > 0) {
    double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi());
    if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
    w_delta_ee_b->Fill(delta_phi_ee_b, MyWeight);
  }
  if (mm_event && Nj > 0 && Nb > 0) {
    double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi());
    if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
    w_delta_mm_b->Fill(delta_phi_mm_b, MyWeight);
  }

}

// ------------ method called once each job just before starting event loop ------------
void GenbAnalyzer::beginJob () {
  LumiWeights_ = edm::LumiReWeighting(path_ + "/" + "pileup_" + pileup_ + ".root", path_ + "/" + "pileup_2012.root", "pileup", "pileup");

}

// ------------ method called once each job just after ending the event loop ------------
void GenbAnalyzer::endJob () {
}

// ------------ method called when starting to processes a run ------------
void GenbAnalyzer::beginRun (edm::Run const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a run ------------
void GenbAnalyzer::endRun (edm::Run const &, edm::EventSetup const &) {
}

// ------------ method called when starting to processes a luminosity block ------------
void GenbAnalyzer::beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a luminosity block ------------
void GenbAnalyzer::endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// define this as a plug-in
DEFINE_FWK_MODULE (GenbAnalyzer);
