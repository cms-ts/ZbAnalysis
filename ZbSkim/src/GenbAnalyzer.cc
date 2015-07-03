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
#include "FWCore/Framework/interface/Run.h"
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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
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

  const reco::GenParticle* getBAncestors (const reco::GenParticle* gp) {
    if (gp->status()<1||gp->status()>3) return NULL;
    if (isB(*gp)) return gp;
    if (gp->numberOfMothers()==0) return NULL;
    for (unsigned int im=0; im<gp->numberOfMothers(); im++) {
      const reco::GenParticle *p = gp->motherRef(im).get();
      const reco::GenParticle* mom = getBAncestors(p);
      if (mom != NULL) return mom; 
    }
    return NULL;
  }

  // ----------member data ---------------------------

  std::string pileupMC_;
  std::string pileupDT_;
  std::string lepton_;
  std::string path_;
  double numB_;
  bool rivet_;

  edm::LumiReWeighting LumiWeights_;

  TH1F*     h_gen_weights;
  TH1F*     h_gen2_weights;

  TH1F*     h_jetmultiplicity;
  TH1F*     h_jet_pt;
  TH1F*     h_ele_pt;
  TH1F*     h_muon_pt;
  TH1F*     w_pt_Z_ee;
  TH1F*     w_y_Z_ee;
  TH1F*     w_y_Z_ee_abs;
  TH1F*     w_pt_Z_mm;
  TH1F*     w_y_Z_mm;
  TH1F*     w_y_Z_mm_abs;
  TH1F*     w_pt_Z_ee_b;
  TH1F*     w_pt_Z_mm_b;
  TH1F*     w_y_Z_ee_b;
  TH1F*     w_y_Z_ee_b_abs;
  TH1F*     w_y_Z_mm_b;
  TH1F*     w_y_Z_mm_b_abs;

  TH1F*     w_jetmultiplicity;
  TH1F*     w_first_jet_pt;      // leading jet of any type
  TH1F*     w_first_jet_eta;
  TH1F*     w_first_jet_eta_abs;
  TH1F*     w_first_bjet_pt;
  TH1F*     w_first_bjet_eta;
  TH1F*     w_first_bjet_eta_abs;
  TH1F*     w_second_bjet_pt;
  TH1F*     w_second_bjet_eta;
  TH1F*     w_second_bjet_eta_abs;
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
  TH1F*     w_mass_Zj_ee;
  TH1F*     w_mass_Zj_mm;
  TH1F*     w_mass_Zj_ee_b;
  TH1F*     w_mass_Zj_mm_b;
  TH1F*     w_delta_ee;
  TH1F*     w_delta_mm;
  TH1F*     w_delta_ee_b;
  TH1F*     w_delta_mm_b;
  TH1F*     w_Ht;
  TH1F*     w_Ht_b;
  TH1F*     w_delta_phi_2b;
  TH1F*     w_DR_bb;

  TH1F*     w_single_bjet_pt;
  TH1F*     w_single_bjet_eta;
  TH1F*     w_single_pt_Z_ee_b;
  TH1F*     w_single_pt_Z_mm_b;
  TH1F*     w_single_delta_ee_b;
  TH1F*     w_single_delta_mm_b;
  TH1F*     w_single_Ht_b;

  //Distr. Angolari
  TH1F*     w_DR_eeb_min;
  TH1F*     w_DR_mmb_min;
  TH1F*     w_DR_eeb_max;
  TH1F*     w_DR_mmb_max;
  TH1F*     w_A_eeb;  
  TH1F*     w_A_mmb;  

  //Ist. mass. inv bb e Zbb
  TH1F*     w_bb_mass;
  TH1F*     w_eebb_mass;
  TH1F*     w_mmbb_mass;

  //Ist. Phi*
  TH1F* w_Phi_star_ee;
  TH1F* w_Phi_star_mm;
  TH1F* w_Phi_star_ee_b;
  TH1F* w_Phi_star_mm_b;

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

  pileupMC_ = iConfig.getUntrackedParameter < std::string > ("pileupMC", "S10");
  pileupDT_ = iConfig.getUntrackedParameter < std::string > ("pileupDT", "");
  lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
  path_ =   iConfig.getUntrackedParameter < std::string > ("path", "/gpfs/cms/users/candelis/work/ZbSkim/test");
  numB_ =  iConfig.getUntrackedParameter <double> ("numB", 0);

  rivet_  = iConfig.getUntrackedParameter < bool > ("rivet", false);

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  TH1::SetDefaultSumw2();

  h_gen_weights     =   fs->make < TH1F > ("h_gen_weights",      "h_gen_weights", 4, 0, 4);
  h_gen2_weights    =   fs->make < TH1F > ("h_gen2_weights",     "h_gen2_weights", 4, 0, 4);

  w_jetmultiplicity =   fs->make < TH1F > ("w_jetmultiplicity",  "w_jetmultiplicity;", 8, 0.5, 8.5);
  w_first_ele_pt =      fs->make < TH1F > ("w_first_ele_pt",     "w_first_ele_pt; [GeV]", 50, 0., 450.);
  w_second_ele_pt =     fs->make < TH1F > ("w_second_ele_pt",    "w_second_ele_pt; [GeV]", 50, 0., 450.);
  w_first_muon_pt =     fs->make < TH1F > ("w_first_muon_pt",    "w_first_muon_pt; [GeV]", 50, 0., 450.);
  w_second_muon_pt =    fs->make < TH1F > ("w_second_muon_pt",   "w_second_muon_pt; [GeV]", 50, 0., 450.);
  w_bjetmultiplicity =  fs->make < TH1F > ("w_bjetmultiplicity", "w_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  w_first_bjet_pt =     fs->make < TH1F > ("w_first_bjet_pt",    "w_first_bjet_pt; [GeV]", 50, 30., 700.);
  w_first_jet_pt =      fs->make < TH1F > ("w_first_jet_pt",     "w_first_jet_pt; [GeV]", 50, 30., 700.);
  w_first_jet_eta =     fs->make < TH1F > ("w_first_jet_eta",    "w_first_jet_eta", 16, -2.5, 2.5);
  w_first_jet_eta_abs = fs->make < TH1F > ("w_first_jet_eta_abs","w_first_jet_eta_abs", 8, 0, 2.5);
  w_first_bjet_eta =    fs->make < TH1F > ("w_first_bjet_eta",   "w_first_bjet_eta", 16, -2.5, 2.5);
  w_first_bjet_eta_abs =fs->make < TH1F > ("w_first_bjet_eta_abs","w_first_bjet_eta_abs", 8, 0, 2.5);
  w_first_jet_pt_b =    fs->make < TH1F > ("w_first_jet_pt_b",   "w_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  w_first_jet_eta_b =   fs->make < TH1F > ("w_first_jet_eta_b",  "w_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  w_second_bjet_pt =    fs->make < TH1F > ("w_second_bjet_pt",   "w_second_bjet_pt; [GeV]", 50, 30., 500.);
  w_second_bjet_eta =   fs->make < TH1F > ("w_second_bjet_eta",  "w_second_bjet_eta", 16, -2.5, 2.5);
  w_second_bjet_eta_abs=fs->make < TH1F > ("w_second_bjet_eta_abs","w_second_bjet_eta_abs", 8, 0, 2.5);
  w_pt_Z_ee =           fs->make < TH1F > ("w_pt_Z_ee",          "w_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_mm =           fs->make < TH1F > ("w_pt_Z_mm",          "w_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_ee_b =         fs->make < TH1F > ("w_pt_Z_ee_b",        "w_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_mm_b =         fs->make < TH1F > ("w_pt_Z_mm_b",        "w_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  w_y_Z_ee =            fs->make < TH1F > ("w_y_Z_ee",           "w_y_Z_ee;y_Z",   20, -2.0, 2.0);
  w_y_Z_ee_abs =        fs->make < TH1F > ("w_y_Z_ee_abs",       "w_y_Z_ee_abs;abs(y_Z)",   10, 0, 2.0);
  w_y_Z_mm =            fs->make < TH1F > ("w_y_Z_mm",           "w_y_Z_mm;y_Z",   20, -2.0, 2.0);
  w_y_Z_mm_abs =        fs->make < TH1F > ("w_y_Z_mm_abs",       "w_y_Z_mm_abs;abs(y_Z)",   10, 0, 2.0);
  w_y_Z_ee_b =          fs->make < TH1F > ("w_y_Z_ee_b",         "w_y_Z_ee_b;y_Z", 20, -2.0, 2.0);
  w_y_Z_ee_b_abs =      fs->make < TH1F > ("w_y_Z_ee_b_abs",     "w_y_Z_ee_b_abs;abs(y_Z)", 10, 0, 2.0);
  w_y_Z_mm_b =          fs->make < TH1F > ("w_y_Z_mm_b",         "w_y_Z_mm_b;y_Z", 20, -2.0, 2.0);
  w_y_Z_mm_b_abs =      fs->make < TH1F > ("w_y_Z_mm_b_abs",     "w_y_Z_mm_b_abs;abs(y_Z)", 10, 0, 2.0);
  w_mass_ee =           fs->make < TH1F > ("w_mass_ee",          "w_mass_ee;Mass [GeV]", 80, 71., 111.);
  w_mass_mm =           fs->make < TH1F > ("w_mass_mm",          "w_mass_mm;Mass [GeV]", 80, 71., 111.);
  w_mass_ee_b =         fs->make < TH1F > ("w_mass_ee_b",        "w_mass_ee_b;Mass [GeV]", 80, 71., 111.);
  w_mass_mm_b =         fs->make < TH1F > ("w_mass_mm_b",        "w_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  w_mass_Zj_ee =        fs->make < TH1F > ("w_mass_Zj_ee",       "w_mass_Zj_ee", 15, 100., 330.);
  w_mass_Zj_mm =        fs->make < TH1F > ("w_mass_Zj_mm",       "w_mass_Zj_mm", 15, 100., 330.);
  w_mass_Zj_ee_b =      fs->make < TH1F > ("w_mass_Zj_ee_b",     "w_mass_Zj_ee_b", 15, 100., 330.);
  w_mass_Zj_mm_b =      fs->make < TH1F > ("w_mass_Zj_mm_b",     "w_mass_Zj_mm_b", 15, 100., 330.);
  w_Ht =                fs->make < TH1F > ("w_Ht",               "w_Ht [GeV]", 50, 30., 1000.);
  w_Ht_b =              fs->make < TH1F > ("w_Ht_b",             "w_Ht_b [GeV]", 50, 30., 1000.);
  w_delta_ee =          fs->make < TH1F > ("w_delta_phi_ee",     "w_delta_phi_ee", 12, 0, TMath::Pi ());
  w_delta_mm =          fs->make < TH1F > ("w_delta_phi_mm",     "w_delta_phi_mm", 12, 0, TMath::Pi ());
  w_delta_ee_b =        fs->make < TH1F > ("w_delta_phi_ee_b",   "w_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_delta_mm_b =        fs->make < TH1F > ("w_delta_phi_mm_b",   "w_delta_phi_mm_b", 12, 0, TMath::Pi ());
  w_delta_phi_2b =      fs->make < TH1F > ("w_delta_phi_2b",     "w_delta_phi_2b",   12, 0, TMath::Pi ());
  w_DR_bb =             fs->make < TH1F > ("w_DR_bb",            "w_DR_bb", 10, 0, 5);

  w_single_bjet_pt =           fs->make < TH1F > ("w_single_bjet_pt",         "w_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  w_single_bjet_eta =          fs->make < TH1F > ("w_single_bjet_eta",        "w_single_bjet_eta", 16, -2.5, 2.5);
  w_single_pt_Z_ee_b =         fs->make < TH1F > ("w_single_pt_Z_ee_b",       "w_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  w_single_pt_Z_mm_b =         fs->make < TH1F > ("w_single_pt_Z_mm_b",       "w_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  w_single_delta_ee_b =        fs->make < TH1F > ("w_single_delta_phi_ee_b",  "w_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_single_delta_mm_b =        fs->make < TH1F > ("w_single_delta_phi_mm_b",  "w_single_delta_phi_mm_b", 12, 0, TMath::Pi ());
  w_single_Ht_b =              fs->make < TH1F > ("w_single_Ht_b",            "w_single_Ht [GeV]", 50, 30., 1000.);
  
  //Distr. Angolari
  w_DR_eeb_min =     fs->make < TH1F > ("w_DR_eeb_min",   "w_DR_eeb_min; Delta_R",  15, 0., 4.);  
  w_DR_mmb_min =     fs->make < TH1F > ("w_DR_mmb_min",   "w_DR_mmb_min; Delta_R",  15, 0., 4.);   
  w_DR_eeb_max =     fs->make < TH1F > ("w_DR_eeb_max",   "w_DR_eeb_max; Delta_R",  15, 1.5, 5.);  
  w_DR_mmb_max =     fs->make < TH1F > ("w_DR_mmb_max",   "w_DR_mmb_max; Delta_R",  15, 1.5, 5.);    
  w_A_eeb =     fs->make < TH1F > ("w_A_eeb",   "w_A_eeb; A", 10, 0., 1.); 
  w_A_mmb =     fs->make < TH1F > ("w_A_mmb",   "w_A_mmb; A", 10, 0., 1.);  

  //Mass. Inv
  w_bb_mass =     fs->make < TH1F > ("w_bb_mass",   "w_bb_mass;Mass [GeV]", 13, 30., 400);
  w_eebb_mass =     fs->make < TH1F > ("w_eebb_mass",   "w_eebb_mass;Mass [GeV]", 12, 200., 500.);
  w_mmbb_mass =     fs->make < TH1F > ("w_mmbb_mass",   "w_mmbb_mass;Mass [GeV]", 12, 200., 500.);

  //Phi*
  w_Phi_star_ee =     fs->make < TH1F > ("w_Phi_star_ee",   "w_Phi_star_ee; Phi*", 10, 0, 1);
  w_Phi_star_mm =     fs->make < TH1F > ("w_Phi_star_mm",   "w_Phi_star_mm; Phi*", 10, 0, 1);
  w_Phi_star_ee_b =     fs->make < TH1F > ("w_Phi_star_ee_b",   "w_Phi_star_ee_b; Phi*", 10, 0, 1);
  w_Phi_star_mm_b =     fs->make < TH1F > ("w_Phi_star_mm_b",   "w_Phi_star_mm_b; Phi*", 10, 0, 1);

  produces<std::vector<double>>("myEventWeight");

  produces<std::vector<math::XYZTLorentzVector>>("myElectrons");
  produces<std::vector<math::XYZTLorentzVector>>("myMuons");

  produces<std::vector<double>>("myPtZ");
  produces<std::vector<double>>("myPtZb");
  
  produces<std::vector<double>>("myYZ");
  produces<std::vector<double>>("myYZb");

  produces<std::vector<double>>("myMassZj");
  produces<std::vector<double>>("myMassZb");

  produces<std::vector<math::XYZTLorentzVector>>("myJets");
  produces<std::vector<math::XYZTLorentzVector>>("myJets2");
  produces<std::vector<double>>("myDeltaPhi");

  produces<std::vector<double>>("myHt");
  produces<std::vector<double>>("myHtb");

  produces<std::vector<math::XYZTLorentzVector>>("myBJets");
  produces<std::vector<math::XYZTLorentzVector>>("myBJets2");
  produces<std::vector<double>>("myBDeltaPhi");

  produces<std::vector<double>>("myDRbb");
  produces<std::vector<double>>("myDeltaPhibb");

  produces<std::vector<double>>("myDRZbMin");
  produces<std::vector<double>>("myDRZbMax");

  produces<std::vector<double>>("myAZb");

  produces<std::vector<double>>("myPhiStar");
  produces<std::vector<double>>("myPhiStarb");

  produces<std::vector<double>>("mybbMass");
  produces<std::vector<double>>("mybbZMass");

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
  
  std::auto_ptr<std::vector<double>> myYZ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myYZb( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myMassZj( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myMassZb( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets2( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<double>> myDeltaPhi( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myHt( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myHtb( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBJets( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBJets2( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<double>> myBDeltaPhi( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myDRbb( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaPhibb( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myDRZbMin( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDRZbMax( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myAZb( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myPhiStar( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myPhiStarb( new std::vector<double> );

  std::auto_ptr<std::vector<double>> mybbMass( new std::vector<double> );
  std::auto_ptr<std::vector<double>> mybbZMass( new std::vector<double> );

  bool ee_event = false;
  bool mm_event = false;
  bool ist = false;

  bool b_selection = true;

  int Nj = 0;
  int Nj2 = 0;
  int Nb = 0;
  int Nb2 = 0;

  double diele_mass = 0;
  double diele_phi = 0;
  double diele_pt = 0;
  double diele_y = 0;
  double diele_eta = 0;

  double dimuon_mass = 0;
  double dimuon_phi = 0;
  double dimuon_pt = 0;
  double dimuon_y = 0;
  double dimuon_eta = 0;

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

  double lepton1_eta = -9999;
  double lepton1_phi = -9999;
  double lepton2_eta = -9999;
  double lepton2_phi = -9999;

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

    std::vector<double> mcWeights = genEventInfoHandle->weights();
    double mcWeight = genEventInfoHandle->weight();

    h_gen_weights->Fill(0.5, 1.0);
    h_gen_weights->Fill(1.5, mcWeight);

    MyWeight = MyWeight*mcWeight;

  }

  edm::Handle<LHEEventProduct> lheEventHandle;

  if (iEvent.getByLabel ("externalLHEProducer", lheEventHandle)) {

/*
    <weightgroup combine="envelope" type="scale_variation">
      <weight id="1001"> muR=0.10000E+01 muF=0.10000E+01 </weight>
      <weight id="1002"> muR=0.10000E+01 muF=0.20000E+01 </weight>
      <weight id="1003"> muR=0.10000E+01 muF=0.50000E+00 </weight>
      <weight id="1004"> muR=0.20000E+01 muF=0.10000E+01 </weight>
      <weight id="1005"> muR=0.20000E+01 muF=0.20000E+01 </weight>
      <weight id="1006"> muR=0.20000E+01 muF=0.50000E+00 </weight>
      <weight id="1007"> muR=0.50000E+00 muF=0.10000E+01 </weight>
      <weight id="1008"> muR=0.50000E+00 muF=0.20000E+01 </weight>
      <weight id="1009"> muR=0.50000E+00 muF=0.50000E+00 </weight>
    </weightgroup>
    <weightgroup combine="gaussian" type="PDF_variation">
      <weight id="2001"> pdfset=260001 </weight>
      <weight id="2002"> pdfset=260002 </weight>
      <weight id="2003"> pdfset=260003 </weight>
...
      <weight id="2098"> pdfset=260098 </weight>
      <weight id="2099"> pdfset=260099 </weight>
      <weight id="2100"> pdfset=260100 </weight>
    </weightgroup>
*/

    //int whichWeight = 1000;
    //double mcWeight2 = lheEventHandle->weights()[whichWeight].wgt/lheEventHandle->originalXWGTUP();
    double mcWeight2 = 1.0;

    h_gen_weights->Fill(2.5, mcWeight2);

    MyWeight = MyWeight*mcWeight2;

  }

  edm::Handle <std::vector<double>>  genEventInfoHandle2;

  if (iEvent.getByLabel ("GenBDWeight", genEventInfoHandle2)) {

    double mcWeight3 = genEventInfoHandle2->empty() ? 1 : (*genEventInfoHandle2)[0];

    h_gen_weights->Fill(3.5, mcWeight3);

    MyWeight = MyWeight*mcWeight3;

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

  TLorentzVector z_ee;

  if ( !ele_dres.lepton_photon.empty() && !pos_dres.lepton_photon.empty() ) {
        vect_ele.push_back(ele_dres.p_part);
        vect_ele.push_back(pos_dres.p_part);
  	z_ee = ele_dres.p_part + pos_dres.p_part;
  	diele_mass = z_ee.M();
  	diele_pt =   z_ee.Pt();
	diele_phi =  z_ee.Phi();
        diele_y = z_ee.Rapidity();
        diele_eta = z_ee.Eta();
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

  TLorentzVector z_mm;

  if ( !mu_dres.lepton_photon.empty() && !antimu_dres.lepton_photon.empty() ) {
        vect_muon.push_back(mu_dres.p_part);
        vect_muon.push_back(antimu_dres.p_part);
        z_mm = mu_dres.p_part + antimu_dres.p_part;
        dimuon_mass = z_mm.M();
        dimuon_pt =   z_mm.Pt();
        dimuon_phi =  z_mm.Phi();
        dimuon_y = z_mm.Rapidity();
        dimuon_eta = z_mm.Eta();
        if (dimuon_mass>71 && dimuon_mass<111) mm_event = true;
  }

  if (ee_event) {
    lepton1_eta = ele_dres.p_part.Eta();
    lepton2_eta = pos_dres.p_part.Eta();
    lepton1_phi = ele_dres.p_part.Phi();
    lepton2_phi = pos_dres.p_part.Phi();
   }

  if (mm_event) {
    lepton1_eta = mu_dres.p_part.Eta();
    lepton2_eta = antimu_dres.p_part.Eta();
    lepton1_phi = mu_dres.p_part.Phi();
    lepton2_phi = antimu_dres.p_part.Phi();
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
	double phij = jets[i].phi_std();
	double ptj = jets[i].perp();
	        
	if (fabs(etaj) < 2.4 && ptj > 30) {
          double delta_eta1 = lepton1_eta - etaj;
          double delta_phi1 = fabs(lepton1_phi - phij);
          if (delta_phi1 > acos(-1)) delta_phi1 = 2*acos(-1) - delta_phi1;
          double delta_eta2 = lepton2_eta - etaj;
          double delta_phi2 = fabs(lepton2_phi - phij);
          if (delta_phi2 > acos(-1)) delta_phi2 = 2*acos(-1) - delta_phi2;

          double deltaR_jl1 = sqrt(pow(delta_eta1,2) + pow(delta_phi1,2));
          double deltaR_jl2 = sqrt(pow(delta_eta2,2) + pow(delta_phi2,2));
	  
	  if (deltaR_jl1 > 0.5 && deltaR_jl2 > 0.5 && (ee_event || mm_event)) { 
	    Nj++;
            vect_jets.push_back(jets[i]);
            Ht = Ht + jets[i].perp();
          }
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

    for (unsigned int c = 0; c < constituents.size(); c++) {
      int index = constituents.at(c).user_index();
      reco::GenParticle gp = part_jets_st1.at(index);
      const reco::GenParticle* bpart = getBAncestors(&gp);
      if (bpart!=NULL) {     
        double deltaEta_Bb = vect_jets[k].eta() - bpart->eta(); 
        double deltaPhi_Bb = fabs(vect_jets[k].phi_std() - bpart->phi());
        if (deltaPhi_Bb > acos(-1)) deltaPhi_Bb = 2*acos(-1) - deltaPhi_Bb;
        double deltaR_Bb = sqrt(pow(deltaEta_Bb,2) + pow(deltaPhi_Bb,2)); 
        if (deltaR_Bb < 0.5) Bjet_found = true;
      }
    }
    if (Bjet_found) {
      vect_bjets.push_back(vect_jets[k]);
      Nb++;
    }  
  }    

  if (Nb != 1 && numB_ == 1) {
    b_selection = false;
  }

  if (Nb < 2 && numB_ == 2) {
    b_selection = false;
  }

  for(unsigned int k = 0; k < vect_jets2.size() ; k++) {
    Bjet2_found = false;
   
    vector<fastjet::PseudoJet> constituents = cseq.constituents(vect_jets2[k]);

    for (unsigned int c = 0; c < constituents.size(); c++) {
      int index = constituents.at(c).user_index();
      reco::GenParticle gp = part_jets_st1.at(index);
      const reco::GenParticle * bpart2 = getBAncestors(&gp);
      if (bpart2!=NULL) {
        double deltaEta_Bb2 = vect_jets2[k].eta() - bpart2->eta();
        double deltaPhi_Bb2 = fabs(vect_jets2[k].phi_std() - bpart2->phi());
        if (deltaPhi_Bb2 > acos(-1)) deltaPhi_Bb2 = 2*acos(-1) - deltaPhi_Bb2;
        double deltaR_Bb2 = sqrt(pow(deltaEta_Bb2,2) + pow(deltaPhi_Bb2,2));
        if (deltaR_Bb2 < 0.5) Bjet2_found = true;
      }
    }
    if (Bjet2_found) {
      vect_bjets2.push_back(vect_jets2[k]);
      Nb2++;
    }
  }

  // Sort b-jets in pT
  std::sort( vect_bjets.begin(), vect_bjets.end(), order_jets() );
  std::sort( vect_bjets2.begin(), vect_bjets2.end(), order_jets() );

  double delta_phi_2b = 0;
  double delta_eta_2b = 0;
  double DR_bb = 999;

  for (std::vector <reco::GenParticle>::const_iterator thepart = genPart->begin(); thepart != genPart->end(); thepart++) {
    if (thepart->pdgId()==23) {
      for (UInt_t i=0; i<thepart->numberOfDaughters(); i++){
        if (abs(thepart->daughter(i)->pdgId())==15 && thepart->daughter(i)->status()==3){
          ist = true;
	}
      }
    }
  }

  // ++++++++ SPECIAL WEIGHTS

#if 0
  if (!ist) {
    if (lepton_ == "electron") {
// MadGraph
      float N1b = 31313.60;
      float N2b = 2561.60;
      float w1b = 0.898136;
      float w2b = 1.185000;
// MadGraph+aMC@NLO
//      float N1b = 3.00179131500000000e+08;
//      float N2b = 2.69809477500000000e+07;
//      float w1b = 0.852622;
//      float w2b = 0.842971;
      if (Nb == 1) {
        MyWeight = MyWeight * w1b * (N1b+N2b)/(w1b*N1b+w2b*N2b);
      }
      if (Nb > 1) {
        MyWeight = MyWeight * w2b * (N1b+N2b)/(w1b*N1b+w2b*N2b);
      }
    }
    if (lepton_ == "muon") {
// MadGraph
      float N1b = 31260.50;
      float N2b = 2485.13;
      float w1b = 0.885228;
      float w2b = 1.182360;
// MadGraph+aMC@NLO
//      float N1b = 2.99066704000000000e+08;
//      float N2b = 2.80502260625000000e+07;
//      float w1b = 0.837129;
//      float w2b = 0.822508;
      if (Nb == 1) {
        MyWeight = MyWeight * w1b * (N1b+N2b)/(w1b*N1b+w2b*N2b);
      }
      if (Nb > 1) {
        MyWeight = MyWeight * w2b * (N1b+N2b)/(w1b*N1b+w2b*N2b);
      }
    }
  }
#endif

  ee_event = ee_event && (lepton_ == "electron") && !ist;
  mm_event = mm_event && (lepton_ == "muon") && !ist;

  double delta_phi_eeb = 0;
  double delta_eta_eeb = 0;
  double DR_eeb = 0;
  double DR_eeb_min = 999;
  double DR_eeb_max = -1;
  double A_eeb = 0;

  double delta_phi_mmb = 0;
  double delta_eta_mmb = 0;
  double DR_mmb = 0;
  double DR_mmb_min = 999;
  double DR_mmb_max = -1;
  double A_mmb = 0;

  double bb_mass = 0;
  double eebb_mass = 0;
  double mmbb_mass = 0;

  if (Nb > 1) {
      delta_phi_2b = fabs(vect_bjets[0].phi_std() - vect_bjets[1].phi_std());
      if (delta_phi_2b > acos(-1)) delta_phi_2b = 2 * acos(-1) - delta_phi_2b;
      delta_eta_2b = fabs(vect_bjets[0].eta() - vect_bjets[1].eta());
      DR_bb = sqrt(delta_phi_2b*delta_phi_2b + delta_eta_2b*delta_eta_2b);
      math::XYZTLorentzVector bb(vect_bjets[0].px() + vect_bjets[1].px() , vect_bjets[0].py() + vect_bjets[1].py() , vect_bjets[0].pz() + vect_bjets[1].pz() , vect_bjets[0].e() + vect_bjets[1].e());
      bb_mass = bb.mass();
      if (ee_event) {
        math::XYZTLorentzVector eebb(bb.Px() + z_ee.Px() , bb.Py() + z_ee.Py() , bb.Pz() + z_ee.Pz() , bb.E() + z_ee.E());
        eebb_mass = eebb.mass();
      }
      if (mm_event) {
        math::XYZTLorentzVector mmbb(bb.Px() + z_mm.Px() , bb.Py() + z_mm.Py() , bb.Pz() + z_mm.Pz() , bb.E() + z_mm.E());
        mmbb_mass = mmbb.mass();
      }
      for (unsigned int i=0; i<vect_bjets.size(); i++) {
        if (ee_event) {
  	  delta_phi_eeb = fabs(diele_phi - vect_bjets[i].phi_std());
          if (delta_phi_eeb > acos(-1)) delta_phi_eeb = 2 * acos(-1) - delta_phi_eeb;
      	  delta_eta_eeb = fabs(diele_eta - vect_bjets[i].eta());
          DR_eeb = TMath::Sqrt(delta_phi_eeb*delta_phi_eeb + delta_eta_eeb*delta_eta_eeb);
          if (DR_eeb <= DR_eeb_min) DR_eeb_min = DR_eeb;
          if (DR_eeb >= DR_eeb_max) DR_eeb_max = DR_eeb;
        }
        if (mm_event) {
	  delta_phi_mmb = fabs(dimuon_phi - vect_bjets[i].phi_std());
          if (delta_phi_mmb > acos(-1)) delta_phi_mmb = 2 * acos(-1) - delta_phi_mmb;
      	  delta_eta_mmb = fabs(dimuon_eta - vect_bjets[i].eta());
          DR_mmb = TMath::Sqrt(delta_phi_mmb*delta_phi_mmb + delta_eta_mmb*delta_eta_mmb);
          if (DR_mmb <= DR_mmb_min) DR_mmb_min = DR_mmb;
          if (DR_mmb >= DR_mmb_max) DR_mmb_max = DR_mmb;
        }
      }
      A_eeb = (DR_eeb_max - DR_eeb_min)/(DR_eeb_min + DR_eeb_max);
      A_mmb = (DR_mmb_max - DR_mmb_min)/(DR_mmb_min + DR_mmb_max);
  }
 
  //Phi*
  double DEta_ee = 999;
  double DPhi_ee = 999;
  double Phi_star_ee = 999;
  double DEta_mm = 999;
  double DPhi_mm = 999;
  double Phi_star_mm = 999;

  if (ee_event){
    DEta_ee = lepton1_eta - lepton2_eta;
    DPhi_ee = fabs(lepton1_phi - lepton2_phi);
    if (DPhi_ee > acos(-1)) DPhi_ee = 2 * acos(-1) - DPhi_ee;
    Phi_star_ee = tan( (acos(-1) - DPhi_ee) / 2 ) * sqrt( 1 - ( tanh( DEta_ee / 2 ) )*( tanh( DEta_ee / 2 ) ) );   
  }
  if (mm_event){
    DEta_mm = lepton1_eta - lepton2_eta; 
    DPhi_mm = fabs(lepton1_phi - lepton2_phi);
    if (DPhi_mm > acos(-1)) DPhi_mm = 2 * acos(-1) - DPhi_mm;
    Phi_star_mm = tan( (acos(-1) - DPhi_mm) / 2 ) * sqrt( 1 - ( tanh( DEta_mm / 2 ) )*( tanh( DEta_mm / 2 ) ) );   
  }

  if (ee_event && Nj > 0) {
    w_first_ele_pt->Fill (ele_dres.p_part.Pt(), MyWeight);
    w_second_ele_pt->Fill (pos_dres.p_part.Pt(), MyWeight);
    w_mass_ee->Fill(diele_mass, MyWeight);
    w_pt_Z_ee->Fill(diele_pt, MyWeight);
    w_y_Z_ee->Fill(diele_y, MyWeight);
    w_y_Z_ee_abs->Fill(fabs(diele_y), MyWeight);
    w_Phi_star_ee->Fill (Phi_star_ee, MyWeight);
    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi_std());
    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
    w_delta_ee->Fill(delta_phi_ee, MyWeight);
    math::XYZTLorentzVector zj_ee_p(vect_jets[0].px()+z_ee.Px(), vect_jets[0].py()+z_ee.Py(), vect_jets[0].pz()+z_ee.Pz(), vect_jets[0].e()+z_ee.E());
    double zj_ee_mass = zj_ee_p.mass();
    w_mass_Zj_ee->Fill (zj_ee_mass, MyWeight);
  }

  if (mm_event && Nj > 0) {
    w_first_muon_pt->Fill (mu_dres.p_part.Pt(), MyWeight);
    w_second_muon_pt->Fill (antimu_dres.p_part.Pt(), MyWeight);
    w_mass_mm->Fill(dimuon_mass, MyWeight);
    w_pt_Z_mm->Fill(dimuon_pt, MyWeight);
    w_y_Z_mm->Fill(dimuon_y, MyWeight);
    w_y_Z_mm_abs->Fill(fabs(dimuon_y), MyWeight);
    w_Phi_star_mm->Fill (Phi_star_mm, MyWeight);
    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi_std());
    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
    w_delta_mm->Fill(delta_phi_mm, MyWeight);
    math::XYZTLorentzVector zj_mm_p(vect_jets[0].px()+z_mm.Px(), vect_jets[0].py()+z_mm.Py(), vect_jets[0].pz()+z_mm.Pz(), vect_jets[0].e()+z_mm.E());
    double zj_mm_mass = zj_mm_p.mass();
    w_mass_Zj_mm->Fill (zj_mm_mass, MyWeight);
  }

  if ((ee_event || mm_event) && Nj > 0) {
    w_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight);
    w_first_jet_eta_abs->Fill (fabs(vect_jets[0].eta()), MyWeight);
    w_jetmultiplicity->Fill (vect_jets.size(), MyWeight);
    w_Ht->Fill (Ht, MyWeight);
  }

  if (ee_event && Nb > 0 && Nj > 0 && b_selection) {
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
    w_first_bjet_eta_abs->Fill (fabs(vect_bjets[0].eta()), MyWeight);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight);
    if (Nb > 1) {
      w_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight);
      w_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight);
      w_second_bjet_eta_abs->Fill (fabs(vect_bjets[1].eta()), MyWeight);
    }
    w_pt_Z_ee_b->Fill (diele_pt, MyWeight);
    w_y_Z_ee_b->Fill (diele_y, MyWeight);
    w_y_Z_ee_b_abs->Fill (fabs(diele_y), MyWeight);
    w_Phi_star_ee_b->Fill (Phi_star_ee, MyWeight);
    w_mass_ee_b->Fill (diele_mass, MyWeight);
    math::XYZTLorentzVector zb_ee_p(vect_bjets[0].px()+z_ee.Px(), vect_bjets[0].py()+z_ee.Py(), vect_bjets[0].pz()+z_ee.Pz(), vect_bjets[0].e()+z_ee.E());
    double zb_ee_mass = zb_ee_p.mass();
    w_mass_Zj_ee_b->Fill (zb_ee_mass, MyWeight);
  }

  if (mm_event && Nb > 0 && Nj > 0 && b_selection) {
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
    w_first_bjet_eta_abs->Fill (fabs(vect_bjets[0].eta()), MyWeight);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight);
    if (Nb > 1) {
      w_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight);
      w_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight);
      w_second_bjet_eta_abs->Fill (fabs(vect_bjets[1].eta()), MyWeight);
    }
    w_pt_Z_mm_b->Fill (dimuon_pt, MyWeight);
    w_y_Z_mm_b->Fill (dimuon_y, MyWeight);
    w_y_Z_mm_b_abs->Fill (fabs(dimuon_y), MyWeight);
    w_Phi_star_mm_b->Fill (Phi_star_mm, MyWeight);
    w_mass_mm_b->Fill (dimuon_mass, MyWeight);
    math::XYZTLorentzVector zb_mm_p(vect_bjets[0].px()+z_mm.Px(), vect_bjets[0].py()+z_mm.Py(), vect_bjets[0].pz()+z_mm.Pz(), vect_bjets[0].e()+z_mm.E());
    double zb_mm_mass = zb_mm_p.mass();
    w_mass_Zj_mm_b->Fill (zb_mm_mass, MyWeight);
  }

  if ((ee_event || mm_event) && Nb > 0 && Nj > 0 && b_selection) {
    w_bjetmultiplicity->Fill (Nb, MyWeight);
    w_Ht_b->Fill (Ht, MyWeight);
    if (Nb > 1) {
      w_delta_phi_2b->Fill (delta_phi_2b, MyWeight);
      w_DR_bb->Fill (DR_bb, MyWeight);
    }
  }

  if (ee_event && Nb > 0 && Nj > 0 && b_selection) {
    double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi_std());
    if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
    w_delta_ee_b->Fill(delta_phi_ee_b, MyWeight);
    if (Nj == 1) {
      w_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight);
      w_single_pt_Z_ee_b->Fill (diele_pt, MyWeight);
    }
  }
  if (mm_event && Nb > 0 && Nj > 0 && b_selection) {
    double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi_std());
    if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
    w_delta_mm_b->Fill(delta_phi_mm_b, MyWeight);
    if (Nj == 1) {
      w_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight);
      w_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight);
    }
  }

  if ((ee_event || mm_event) && Nb > 0 && Nj == 1 && b_selection) {
    w_single_Ht_b->Fill (Ht, MyWeight);
    w_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
  }
  
  if ((ee_event || mm_event) && Nj > 0 && Nb > 1 && b_selection) {
    w_bb_mass->Fill (bb_mass, MyWeight);
    if (ee_event) {
      w_DR_eeb_min->Fill (DR_eeb_min, MyWeight);
      w_DR_eeb_max->Fill (DR_eeb_max, MyWeight);
      w_A_eeb->Fill (A_eeb, MyWeight);
      w_eebb_mass->Fill (eebb_mass, MyWeight);
    } 
    if (mm_event) {
      w_DR_mmb_min->Fill (DR_mmb_min, MyWeight);
      w_DR_mmb_max->Fill (DR_mmb_max, MyWeight);
      w_A_mmb->Fill (A_mmb, MyWeight);
      w_mmbb_mass->Fill (mmbb_mass, MyWeight);
    }
  }
  

  // ++++++++ OUTPUT COLLECTIONS

  if ((ee_event || mm_event) && Nj2 > 0) {
     myEventWeight->push_back(MyWeight);
  }

  if (ee_event && Nj > 0) {
    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[0].Px(),vect_ele[0].Py(),vect_ele[0].Pz(),vect_ele[0].E()));
    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[1].Px(),vect_ele[1].Py(),vect_ele[1].Pz(),vect_ele[1].E()));
    myPtZ->push_back(diele_pt);
    myYZ->push_back(diele_y);
    myPhiStar->push_back(Phi_star_ee);
    if (Nb > 0 && b_selection) {
      myPtZb->push_back(diele_pt);
      myYZb->push_back(diele_y);
      myPhiStarb->push_back(Phi_star_ee);
      if (Nb > 1) {
        myDRZbMin->push_back(DR_eeb_min);
        myDRZbMax->push_back(DR_eeb_max);
        myAZb->push_back(A_eeb); 
        mybbZMass->push_back(eebb_mass);
      }
    }
  }

  if (mm_event && Nj > 0) {
    myMuons->push_back(math::XYZTLorentzVector(vect_muon[0].Px(),vect_muon[0].Py(),vect_muon[0].Pz(),vect_muon[0].E()));
    myMuons->push_back(math::XYZTLorentzVector(vect_muon[1].Px(),vect_muon[1].Py(),vect_muon[1].Pz(),vect_muon[1].E()));
    myPtZ->push_back(dimuon_pt);
    myYZ->push_back(dimuon_y);
    myPhiStar->push_back(Phi_star_mm);
    if (Nb > 0 && b_selection) {
      myPtZb->push_back(dimuon_pt);
      myYZb->push_back(dimuon_y);
      myPhiStarb->push_back(Phi_star_mm);
      if (Nb > 1) {
        myDRZbMin->push_back(DR_mmb_min); 
        myDRZbMax->push_back(DR_mmb_max);
        myAZb->push_back(A_mmb);
        mybbZMass->push_back(mmbb_mass);
      }
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
    if (Nb > 0 && b_selection) {
      myHtb->push_back(Ht);
      if (Nb > 1) {
        myDeltaPhibb->push_back(delta_phi_2b);
        mybbMass->push_back(bb_mass);
	myDRbb->push_back(DR_bb); 
      }
    }
  }

  if (ee_event && Nj > 0) {
    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi_std());
    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
    myDeltaPhi->push_back(delta_phi_ee);
    math::XYZTLorentzVector zj_ee_p(vect_jets[0].px()+z_ee.Px(), vect_jets[0].py()+z_ee.Py(), vect_jets[0].pz()+z_ee.Pz(), vect_jets[0].e()+z_ee.E());
    double zj_ee_mass = zj_ee_p.mass();
    myMassZj->push_back(zj_ee_mass);
    if (Nb > 0 && b_selection) {
      double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi_std());
      if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
      myBDeltaPhi->push_back(delta_phi_ee_b);
      math::XYZTLorentzVector zb_ee_p(vect_bjets[0].px()+z_ee.Px(), vect_bjets[0].py()+z_ee.Py(), vect_bjets[0].pz()+z_ee.Pz(), vect_bjets[0].e()+z_ee.E());
      double zb_ee_mass = zb_ee_p.mass();
      myMassZb->push_back(zb_ee_mass);
    }
  }

  if (mm_event && Nj > 0) {
    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi_std());
    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
    myDeltaPhi->push_back(delta_phi_mm);
    math::XYZTLorentzVector zj_mm_p(vect_jets[0].px()+z_mm.Px(), vect_jets[0].py()+z_mm.Py(), vect_jets[0].pz()+z_mm.Pz(), vect_jets[0].e()+z_mm.E());
    double zj_mm_mass = zj_mm_p.mass();
    myMassZj->push_back(zj_mm_mass);
    if (Nb > 0 && b_selection) {
      double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi_std());
      if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
      myBDeltaPhi->push_back(delta_phi_mm_b);
      math::XYZTLorentzVector zb_mm_p(vect_bjets[0].px()+z_mm.Px(), vect_bjets[0].py()+z_mm.Py(), vect_bjets[0].pz()+z_mm.Pz(), vect_bjets[0].e()+z_mm.E());
      double zb_mm_mass = zb_mm_p.mass();
      myMassZb->push_back(zb_mm_mass);
    }
  }

  iEvent.put( myEventWeight, "myEventWeight" );

  iEvent.put( myElectrons, "myElectrons" );
  iEvent.put( myMuons, "myMuons" );

  iEvent.put( myPtZ, "myPtZ" );
  iEvent.put( myPtZb, "myPtZb" );
  
  iEvent.put( myYZ,  "myYZ" );
  iEvent.put( myYZb, "myYZb" );

  iEvent.put( myMassZj, "myMassZj" );
  iEvent.put( myMassZb, "myMassZb" );

  iEvent.put( myJets, "myJets" );
  iEvent.put( myJets2, "myJets2" );

  iEvent.put( myDeltaPhi, "myDeltaPhi" );

  iEvent.put( myHt, "myHt" );
  iEvent.put( myHtb, "myHtb" );

  iEvent.put( myBJets, "myBJets" );
  iEvent.put( myBJets2, "myBJets2" );

  iEvent.put( myBDeltaPhi, "myBDeltaPhi" );

  iEvent.put( myDRbb, "myDRbb" );
  iEvent.put( myDeltaPhibb, "myDeltaPhibb" );

  iEvent.put( myDRZbMin, "myDRZbMin" );
  iEvent.put( myDRZbMax, "myDRZbMax" );
 
  iEvent.put( myAZb, "myAZb" );

  iEvent.put( myPhiStar, "myPhiStar" );
  iEvent.put( myPhiStarb, "myPhiStarb" );

  iEvent.put( mybbMass, "mybbMass" );
  iEvent.put( mybbZMass, "mybbZMass" );
}

// ------------ method called once each job just before starting event loop ------------
void GenbAnalyzer::beginJob () {

  LumiWeights_ = edm::LumiReWeighting(path_ + "/" + "pileup_" + pileupMC_ + ".root", path_ + "/" + "pileup_2012_" + pileupDT_ + ".root", "pileup", "pileup");

}

// ------------ method called once each job just after ending the event loop ------------
void GenbAnalyzer::endJob () {
}

// ------------ method called when starting to processes a run ------------
void GenbAnalyzer::beginRun (edm::Run const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a run ------------
void GenbAnalyzer::endRun (edm::Run const & iRun, edm::EventSetup const &) {

/*
  edm::Handle<LHERunInfoProduct> run;
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;

  iRun.getByLabel ("externalLHEProducer", run);
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());

  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
  }
*/

}

// ------------ method called when starting to processes a luminosity block ------------
void GenbAnalyzer::beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a luminosity block ------------
void GenbAnalyzer::endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// define this as a plug-in
DEFINE_FWK_MODULE (GenbAnalyzer);

