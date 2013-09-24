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
// $Id: GenbAnalyzer.cc,v 1.37 2013/07/22 10:27:10 dellaric Exp $
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
#include "FWCore/Framework/interface/EDProducer.h"
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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "Rivet/Projections/FastJets.hh"

//
// class declaration
//

class GenbAnalyzer:public  edm::EDProducer {

public:

  explicit GenbAnalyzer (const edm::ParameterSet &);
  ~GenbAnalyzer ();

private:

  virtual void beginJob ();
  virtual void produce (edm::Event &, const edm::EventSetup &);
  virtual void endJob ();

  virtual void beginRun (edm::Run const &, edm::EventSetup const &);
  virtual void endRun (edm::Run const &, edm::EventSetup const &);
  virtual void beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);
  virtual void endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);

  struct order_leptons { bool operator() (const TLorentzVector &p1, const TLorentzVector &p2) const {
      return (p1.Pt() > p2.Pt());
    }
  };
  struct order_jets { bool operator() (const fastjet::PseudoJet &j1, const fastjet::PseudoJet &j2) const {
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

  bool isB (reco::GenParticle gp) {
    int pid = gp.pdgId();
    if ((abs(pid)/100)%10 == 5 || (abs(pid)/1000)%10 == 5) {
      return true;
    } else {
      return false;
    }
  }

  bool hasBAncestors (reco::GenParticle gp) {
    if (gp.status()<1||gp.status()>3) return false;
    if (isB(gp)) return true;
    if (gp.numberOfMothers()==0) return false;
    bool found = false;
    for (unsigned int im=0;im<gp.numberOfMothers();im++) {
      found = found || hasBAncestors(*gp.motherRef(im));
    }
    return found;
  }

  // ----------member data ---------------------------

  std::string pileup_;
  std::string lepton_;
  std::string path_;
  bool rivet_;

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
  TH1F*     w_first_muon_pt;
  TH1F*     w_second_muon_pt;
  TH1F*     w_mass_ee;
  TH1F*     w_mass_mm;
  TH1F*     w_mass_ee_b;  // at least one b jet in the event
  TH1F*     w_mass_mm_b;
  TH1F*     w_delta_ee;
  TH1F*     w_delta_mm;
  TH1F*     w_delta_ee_b;
  TH1F*     w_delta_mm_b;
  TH1F*     w_Ht;
  TH1F*     w_Ht_b;

  TH1F*     w_single_bjet_pt;
  TH1F*     w_single_bjet_eta;
  TH1F*     w_single_pt_Z_ee_b;
  TH1F*     w_single_pt_Z_mm_b;
  TH1F*     w_single_delta_ee_b;
  TH1F*     w_single_delta_mm_b;
  TH1F*     w_single_Ht_b;

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

  rivet_  = iConfig.getUntrackedParameter < bool > ("rivet", false);

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  h_gen_weights     =   fs->make < TH1F > ("h_gen_weights",      "h_gen_weights", 2, 0, 2);

  w_jetmultiplicity =   fs->make < TH1F > ("w_jetmultiplicity",  "w_jetmultiplicity;", 8, 0.5, 8.5);
  w_first_ele_pt =      fs->make < TH1F > ("w_first_ele_pt",     "w_first_ele_pt; [GeV]", 50, 0., 450.);
  w_second_ele_pt =     fs->make < TH1F > ("w_second_ele_pt",    "w_second_ele_pt; [GeV]", 50, 0., 450.);
  w_first_muon_pt =     fs->make < TH1F > ("w_first_muon_pt",    "w_first_muon_pt; [GeV]", 50, 0., 450.);
  w_second_muon_pt =    fs->make < TH1F > ("w_second_muon_pt",   "w_second_muon_pt; [GeV]", 50, 0., 450.);
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
  w_delta_ee =          fs->make < TH1F > ("w_delta_phi_ee",     "w_delta_phi_ee", 12, 0, TMath::Pi ());
  w_delta_mm =          fs->make < TH1F > ("w_delta_phi_mm",     "w_delta_phi_mm", 12, 0, TMath::Pi ());
  w_delta_ee_b =        fs->make < TH1F > ("w_delta_phi_ee_b",   "w_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_delta_mm_b =        fs->make < TH1F > ("w_delta_phi_mm_b",   "w_delta_phi_mm_b", 12, 0, TMath::Pi ());

  w_single_bjet_pt =           fs->make < TH1F > ("w_single_bjet_pt",         "w_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  w_single_bjet_eta =          fs->make < TH1F > ("w_single_bjet_eta",        "w_single_bjet_eta", 16, -2.5, 2.5);
  w_single_pt_Z_ee_b =         fs->make < TH1F > ("w_single_pt_Z_ee_b",       "w_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  w_single_pt_Z_mm_b =         fs->make < TH1F > ("w_single_pt_Z_mm_b",       "w_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  w_single_delta_ee_b =        fs->make < TH1F > ("w_single_delta_phi_ee_b",  "w_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_single_delta_mm_b =        fs->make < TH1F > ("w_single_delta_phi_mm_b",  "w_single_delta_phi_mm_b", 12, 0, TMath::Pi ());
  w_single_Ht_b =              fs->make < TH1F > ("w_single_Ht_b",            "w_single_Ht [GeV]", 50, 30., 1000.);

  produces<std::vector<double>>("myEventWeight");

  produces<std::vector<math::XYZTLorentzVector>>("myElectrons");
  produces<std::vector<math::XYZTLorentzVector>>("myMuons");

  produces<std::vector<double>>("myPtZ");
  produces<std::vector<double>>("myPtZb");

  produces<std::vector<math::XYZTLorentzVector>>("myJets");
  produces<std::vector<math::XYZTLorentzVector>>("myJets2");
  produces<std::vector<double>>("myDeltaPhi");

  produces<std::vector<double>>("myHt");
  produces<std::vector<double>>("myHtb");

  produces<std::vector<math::XYZTLorentzVector>>("myBJets");
  produces<std::vector<math::XYZTLorentzVector>>("myBJets2");
  produces<std::vector<double>>("myBDeltaPhi");

}

GenbAnalyzer::~GenbAnalyzer () {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event ------------
void GenbAnalyzer::produce (edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace edm;
  using namespace std;

 edm::Handle<vector<reco::GenParticle> > genPart;
  iEvent.getByLabel ("genParticles", genPart);

  std::auto_ptr<std::vector<double>> myEventWeight( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myElectrons( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myMuons( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myPtZ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myPtZb( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets2( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<double>> myDeltaPhi( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myHt( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myHtb( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBJets( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBJets2( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<double>> myBDeltaPhi( new std::vector<double> );

  bool ee_event = false;
  bool mm_event = false;

  int Nj = 0;
  int Nj2 = 0;
  int Nb = 0;
  int Nb2 = 0;

  double diele_mass = 0;
  double diele_phi = 0;
  double diele_pt = 0;

  double dimuon_mass = 0;
  double dimuon_phi = 0;
  double dimuon_pt = 0;

  double MyWeight = 1;

  double Ht = 0;

  struct pt_and_particles {
    TLorentzVector p_part;
    vector<unsigned int> lepton_photon;
  };

  struct pt_and_particles ele_dres;
  struct pt_and_particles pos_dres;
  vector<unsigned int> ele_photons;
  vector<unsigned int> pos_photons;

  struct pt_and_particles mu_dres;
  struct pt_and_particles antimu_dres;
  vector<unsigned int> mu_photons;
  vector<unsigned int> antimu_photons;

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

  if (rivet_) MyWeight = 1.0;

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;

  if (iEvent.getByLabel ("generator", genEventInfoHandle)) {

    double mcWeight = genEventInfoHandle->weight();

    h_gen_weights->Fill(0.5, 1.0);
    h_gen_weights->Fill(1.5, mcWeight);

    MyWeight = MyWeight*mcWeight;

  }

  // +++++++++ ELECTRONS

  //vector < unsigned int > lepton_photon;
  vector <reco::GenParticle> part_vect;

  unsigned int index_ele=0;

  vector < TLorentzVector > vect_ele;

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
    if (fabs(itgen->pdgId())==11 && itgen->status()==1) { // loop over gen electrons
      ele_photons.clear();
      pos_photons.clear();
      TLorentzVector ele;
      ele.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());

      if (itgen->pdgId()==11) {
      	ele_photons.push_back(index_ele);
      } else {
	pos_photons.push_back(index_ele);
      }

      unsigned int index_gamma=0;
      // Loop over photons: FSR dressing for electrons
      for (vector<reco::GenParticle>::const_iterator itgen2=genPart->begin(); itgen2!=genPart->end(); itgen2++) {
        if (fabs(itgen2->pdgId())==22 && itgen2->status()==1) { // loop over primary gen photon
       	  TLorentzVector gam;
          gam.SetPtEtaPhiM(itgen2->pt(),itgen2->eta(),itgen2->phi(),itgen2->mass());
	  double deltaPhi_eg = fabs(gam.Phi()- itgen->phi());
	  if (deltaPhi_eg > acos(-1)) deltaPhi_eg= 2*acos(-1) - deltaPhi_eg;
	  double deltaR_eg = sqrt( deltaPhi_eg*deltaPhi_eg + pow(fabs(gam.Eta()-itgen->eta()),2));
	  if (deltaR_eg< 0.1) {
	    if (itgen->pdgId()==11) {
	    	ele_photons.push_back(index_gamma);
	    } else {
		pos_photons.push_back(index_gamma);
	    }
	    ele += gam;
	  }
	}
	index_gamma++;
      }
      if (ele.Pt()>20 && fabs(ele.Eta())<2.4) {
	if (itgen->pdgId()==11) { // e' un elettrone
		if (ele_dres.lepton_photon.empty()) {
			ele_dres.p_part = ele;
			ele_dres.lepton_photon = ele_photons;
        	} else {
			if (ele.Pt()>ele_dres.p_part.Pt()) {
				ele_dres.p_part = ele;
				ele_dres.lepton_photon = ele_photons;
			}
		}
      	} else { // e' un positrone
		if (pos_dres.lepton_photon.empty()) {
			pos_dres.p_part = ele;
			pos_dres.lepton_photon = pos_photons;
		} else {
			if (ele.Pt()>pos_dres.p_part.Pt()) {
				pos_dres.p_part = ele;
				pos_dres.lepton_photon = pos_photons;
			}
		}
	}
      }
    }
    index_ele++;
  }

  if ( !ele_dres.lepton_photon.empty() && !pos_dres.lepton_photon.empty() ) {
        vect_ele.push_back(ele_dres.p_part);
        vect_ele.push_back(pos_dres.p_part);
  	TLorentzVector y;
  	y = ele_dres.p_part + pos_dres.p_part;
  	diele_mass = y.M();
  	diele_pt =   y.Pt();
	diele_phi =  y.Phi();
	if (diele_mass>71 && diele_mass<111) ee_event = true;
  }

  // +++++++++ MUONS

  vector < TLorentzVector > vect_muon;
  unsigned int index_mu = 0;

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
    if (fabs(itgen->pdgId())==13 && itgen->status()==1) { // loop over gen muons
      mu_photons.clear();
      antimu_photons.clear();
      TLorentzVector muon;
      muon.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());

      if (itgen->pdgId()==13) {
        mu_photons.push_back(index_mu);
      } else {
        antimu_photons.push_back(index_mu);
      }

      // Loop over photons: FSR dressing for muons
      unsigned int index_gammamu = 0;
      for (vector<reco::GenParticle>::const_iterator itgen2=genPart->begin(); itgen2!=genPart->end(); itgen2++) {
        if (fabs(itgen2->pdgId())==22 && itgen2->status()==1) { // loop over primary gen photon
	  TLorentzVector gam;
	  gam.SetPtEtaPhiM(itgen2->pt(),itgen2->eta(),itgen2->phi(),itgen2->mass());
	  double deltaPhi_mg = fabs(gam.Phi()- itgen->phi());
	  if (deltaPhi_mg > acos(-1)) deltaPhi_mg= 2*acos(-1) - deltaPhi_mg;
	  double deltaR_mg = sqrt( deltaPhi_mg*deltaPhi_mg + pow(fabs(gam.Eta()-itgen->eta()),2));
	  if (deltaR_mg< 0.1) {
	    if (itgen->pdgId()==13) {
                mu_photons.push_back(index_gammamu);
	    } else {
                antimu_photons.push_back(index_gammamu);
	    }
	    muon += gam;
	  }
        }
	index_gammamu++;
      }
      if (muon.Pt()>20 && fabs(muon.Eta())<2.4) {
        if (itgen->pdgId()==13) { // e' un muone
                if (mu_dres.lepton_photon.empty()) {
                        mu_dres.p_part = muon;
                        mu_dres.lepton_photon = mu_photons;
                } else {
                        if (muon.Pt()>mu_dres.p_part.Pt()) {
                                mu_dres.p_part = muon;
                                mu_dres.lepton_photon = mu_photons;
                        }
                }
        } else { // e' un antimuone
                if (antimu_dres.lepton_photon.empty()) {
                        antimu_dres.p_part = muon;
                        antimu_dres.lepton_photon = antimu_photons;
                } else {
                        if (muon.Pt()>antimu_dres.p_part.Pt()) {
                                antimu_dres.p_part = muon;
                                antimu_dres.lepton_photon = antimu_photons;
                        }
                }
        }
      }
    }
    index_mu++;
  }

  if ( !mu_dres.lepton_photon.empty() && !antimu_dres.lepton_photon.empty() ) {
        vect_muon.push_back(mu_dres.p_part);
        vect_muon.push_back(antimu_dres.p_part);
        TLorentzVector y;
        y = mu_dres.p_part + antimu_dres.p_part;
        dimuon_mass = y.M();
        dimuon_pt =   y.Pt();
        dimuon_phi =  y.Phi();
        if (dimuon_mass>71 && dimuon_mass<111) mm_event = true;
  }

  vector<reco::GenParticle> part_jets;
  vector<reco::GenParticle> part_jets_st1;

  unsigned int index_nu = 0;
  unsigned index_ch = 0;

  vector<unsigned int> neutrino;
  vector<unsigned int> ch_part;

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
	if (itgen->status()==1 && (fabs(itgen->pdgId())==12 || fabs(itgen->pdgId())==14 || fabs(itgen->pdgId())==16)) {
		neutrino.push_back(index_nu);
	}
	index_nu++;
/*
	if (itgen->charge()!=0 && itgen->pt()<0.25 && itgen->status()==1) {
		ch_part.push_back(index_ch);
	}
*/
	index_ch++;
  }

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
	part_jets.push_back(*itgen);
  }

  vector<unsigned int> canc_part;

  if (ee_event) {
	canc_part = ele_dres.lepton_photon;
	for (unsigned int i =0; i < pos_dres.lepton_photon.size(); i++)
		canc_part.push_back(pos_dres.lepton_photon.at(i));
  }

  if (mm_event) {
	canc_part = mu_dres.lepton_photon;
        for (unsigned int l = 0; l < antimu_dres.lepton_photon.size(); l++)
		canc_part.push_back(antimu_dres.lepton_photon.at(l));
  }

  if (!neutrino.empty()) {
	for (unsigned int j=0; j < neutrino.size(); j++ )
		canc_part.push_back(neutrino.at(j));
  }

  if (!ch_part.empty()) {
	for (unsigned int p = 0; p < ch_part.size(); p++ )
		canc_part.push_back(ch_part.at(p));
  }

  //ordino il vettore in modo crescente
  stable_sort(canc_part.begin(), canc_part.end());

  //scorro il vettore dal basso e per ogni elemento cancello il corrispondente elemento da part_jets
  for (int i = (canc_part.size() - 1); i >= 0; i--)
	part_jets.erase(part_jets.begin() + canc_part.at(i));

  for (unsigned int j=0; j<part_jets.size(); j++) {
	if (part_jets.at(j).status()==1)
		part_jets_st1.push_back(part_jets.at(j));
  }

  // part_jets_st1 contiene ora tutte le particelle su cui fare la riclusterizzazione
  std::vector<fastjet::PseudoJet> vecs;

  for (unsigned int l = 0; l < part_jets_st1.size(); l++) {
	TLorentzVector part;
	part.SetPtEtaPhiM(part_jets_st1.at(l).pt(),part_jets_st1.at(l).eta(),part_jets_st1.at(l).phi(),part_jets_st1.at(l).mass());

	fastjet::PseudoJet pseudoJet(part_jets_st1.at(l).px(), part_jets_st1.at(l).py(), part_jets_st1.at(l).pz(), part.E());
	pseudoJet.set_user_index(l);
	vecs.push_back(pseudoJet);
  }

  // ++++++++ JETS

  vector<fastjet::PseudoJet> vect_jets;
  vector<fastjet::PseudoJet> vect_jets2;
  fastjet::ClusterSequence cseq(vecs, fastjet::JetDefinition(fastjet:: antikt_algorithm, 0.5));
  vector<fastjet::PseudoJet> jets = sorted_by_pt(cseq.inclusive_jets(0.0));
  for (unsigned int i = 0; i < jets.size(); i++) {
	double etaj = jets[i].eta();
	double ptj = jets[i].perp();

	if (fabs(etaj) < 2.5 && ptj > 30) {
		Nj++;
		vect_jets.push_back(jets[i]);
		Ht = Ht + jets[i].perp();
        }
	Nj2++;
	vect_jets2.push_back(jets[i]);
  }

  // ++++++++ BJETS

  /*loop over gen particles, find the b*/

  vector<fastjet::PseudoJet> vect_bjets;
  vector<fastjet::PseudoJet> vect_bjets2;

  bool Bjet_found = false;
  bool Bjet2_found = false;

  for(unsigned int k = 0; k < vect_jets.size() ; k++) {
    Bjet_found = false;
    vector<fastjet::PseudoJet> constituents = cseq.constituents(vect_jets[k]);

    for (unsigned int c = 0; c < constituents.size() && !Bjet_found; c++) {
      int index = constituents.at(c).user_index();
      reco::GenParticle gp = part_jets_st1.at(index);
      if (hasBAncestors(gp)) {
        Bjet_found = true;
      }
    }

    if(Bjet_found) {
      vect_bjets.push_back(vect_jets[k]);
      Nb++;
    }
  }

  for(unsigned int k = 0; k < vect_jets2.size() ; k++) {
    Bjet2_found = false;
    vector<fastjet::PseudoJet> constituents = cseq.constituents(vect_jets2[k]);

    for (unsigned int c = 0; c < constituents.size() && !Bjet2_found; c++) {
      int index = constituents.at(c).user_index();
      reco::GenParticle gp = part_jets_st1.at(index);
      if (hasBAncestors(gp)){
        Bjet2_found = true;
      }
    }

    if(Bjet2_found) {
      vect_bjets2.push_back(vect_jets2[k]);
      Nb2++;
    }
  }

  // Sort b-jets in pT
  std::sort( vect_bjets.begin(), vect_bjets.end(), order_jets() );
  std::sort( vect_bjets2.begin(), vect_bjets2.end(), order_jets() );

  ee_event = ee_event && (lepton_ == "electron");
  mm_event = mm_event && (lepton_ == "muon");

  if (ee_event && Nj > 0) {
    w_first_ele_pt->Fill (ele_dres.p_part.Pt(), MyWeight);
    w_second_ele_pt->Fill (pos_dres.p_part.Pt(), MyWeight);
    w_mass_ee->Fill(diele_mass, MyWeight);
    w_pt_Z_ee->Fill(diele_pt, MyWeight);
    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi_std());
    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
    w_delta_ee->Fill(delta_phi_ee, MyWeight);
  }

  if (mm_event && Nj > 0) {
    w_first_muon_pt->Fill (mu_dres.p_part.Pt(), MyWeight);
    w_second_muon_pt->Fill (antimu_dres.p_part.Pt(), MyWeight);
    w_mass_mm->Fill(dimuon_mass, MyWeight);
    w_pt_Z_mm->Fill(dimuon_pt, MyWeight);
    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi_std());
    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
    w_delta_mm->Fill(delta_phi_mm, MyWeight);
  }

  if ((ee_event || mm_event) && Nj > 0) {
    w_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight);
    w_jetmultiplicity->Fill (vect_jets.size(), MyWeight);
    w_Ht->Fill (Ht, MyWeight);
  }

  if (ee_event && Nb > 0 && Nj > 0) {
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight);
    w_pt_Z_ee_b->Fill (diele_pt, MyWeight);
    w_mass_ee_b->Fill (diele_mass, MyWeight);
  }

  if (mm_event && Nb > 0 && Nj > 0) {
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight);
    w_pt_Z_mm_b->Fill (dimuon_pt, MyWeight);
    w_mass_mm_b->Fill (dimuon_mass, MyWeight);
  }

  if ((ee_event || mm_event) && Nb > 0 && Nj > 0) {
    w_bjetmultiplicity->Fill (Nb, MyWeight);
    w_Ht_b->Fill (Ht, MyWeight);
  }

  if (ee_event && Nb > 0 && Nj > 0) {
    double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi_std());
    if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
    w_delta_ee_b->Fill(delta_phi_ee_b, MyWeight);
    if (Nj == 1) {
      w_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight);
      w_single_pt_Z_ee_b->Fill (diele_pt, MyWeight);
    }
  }
  if (mm_event && Nb > 0 && Nj > 0) {
    double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi_std());
    if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
    w_delta_mm_b->Fill(delta_phi_mm_b, MyWeight);
    if (Nj == 1) {
      w_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight);
      w_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight);
    }
  }

  if ((ee_event || mm_event) && Nb > 0 && Nj == 1) {
    w_single_Ht_b->Fill (Ht, MyWeight);
    w_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
  }

  // ++++++++ OUTPUT COLLECTIONS

  if ((ee_event || mm_event) && Nj2 > 0) {
     myEventWeight->push_back(MyWeight);
  }

  if (ee_event && Nj > 0) {
    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[0].Px(),vect_ele[0].Py(),vect_ele[0].Pz(),vect_ele[0].E()));
    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[1].Px(),vect_ele[1].Py(),vect_ele[1].Pz(),vect_ele[1].E()));
    myPtZ->push_back(diele_pt);
    if (Nb > 0) {
      myPtZb->push_back(diele_pt);
    }
  }

  if (mm_event && Nj > 0) {
    myMuons->push_back(math::XYZTLorentzVector(vect_muon[0].Px(),vect_muon[0].Py(),vect_muon[0].Pz(),vect_muon[0].E()));
    myMuons->push_back(math::XYZTLorentzVector(vect_muon[1].Px(),vect_muon[1].Py(),vect_muon[1].Pz(),vect_muon[1].E()));
    myPtZ->push_back(dimuon_pt);
    if (Nb > 0) {
      myPtZb->push_back(dimuon_pt);
    }
  }

  if (ee_event || mm_event) {
    for (unsigned int i=0; i<vect_jets.size(); ++i) {
      myJets->push_back(math::XYZTLorentzVector(vect_jets[i].px(),vect_jets[i].py(),vect_jets[i].pz(),vect_jets[i].e()));
    }
    for (unsigned int i=0; i<vect_bjets.size(); ++i) {
      myBJets->push_back(math::XYZTLorentzVector(vect_bjets[i].px(),vect_bjets[i].py(),vect_bjets[i].pz(),vect_bjets[i].e()));
    }
  }

  if (ee_event || mm_event) {
    for (unsigned int i=0; i<vect_jets2.size(); ++i) {
      myJets2->push_back(math::XYZTLorentzVector(vect_jets2[i].px(),vect_jets2[i].py(),vect_jets2[i].pz(),vect_jets2[i].e()));
    }
    for (unsigned int i=0; i<vect_bjets2.size(); ++i) {
      myBJets2->push_back(math::XYZTLorentzVector(vect_bjets2[i].px(),vect_bjets2[i].py(),vect_bjets2[i].pz(),vect_bjets2[i].e()));
    }
  }

  if ((ee_event || mm_event) && Nj > 0) {
    myHt->push_back(Ht);
    if (Nb > 0) {
      myHtb->push_back(Ht);
    }
  }

  if (ee_event && Nj > 0) {
    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi_std());
    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
    myDeltaPhi->push_back(delta_phi_ee);
    if (Nb > 0) {
      double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi_std());
      if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
      myBDeltaPhi->push_back(delta_phi_ee_b);
    }
  }

  if (mm_event && Nj > 0) {
    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi_std());
    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
    myDeltaPhi->push_back(delta_phi_mm);
    if (Nb > 0) {
      double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi_std());
      if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
      myBDeltaPhi->push_back(delta_phi_mm_b);
    }
  }

  iEvent.put( myEventWeight, "myEventWeight" );

  iEvent.put( myElectrons, "myElectrons" );
  iEvent.put( myMuons, "myMuons" );

  iEvent.put( myPtZ, "myPtZ" );
  iEvent.put( myPtZb, "myPtZb" );

  iEvent.put( myJets, "myJets" );
  iEvent.put( myJets2, "myJets2" );

  iEvent.put( myDeltaPhi, "myDeltaPhi" );

  iEvent.put( myHt, "myHt" );
  iEvent.put( myHtb, "myHtb" );

  iEvent.put( myBJets, "myBJets" );
  iEvent.put( myBJets2, "myBJets2" );

  iEvent.put( myBDeltaPhi, "myBDeltaPhi" );

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
