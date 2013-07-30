// -*- C++ -*-
//
// Package: ZbAnalyzer
// Class: ZbAnalyzer
//
/**\class ZbAnalyzer ZbAnalyzer.cc ZbAnalysis/ZbAnalyzer/src/ZbAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]

*/
//
// Original Author: Vieri Candelise
// Created: Thu Jan 10 15:57:03 CET 2013
// $Id: ZbAnalyzer.cc,v 1.111 2013/07/20 07:23:57 dellaric Exp $
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
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
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

#include "table.h"

//
// class declaration
//

class ZbAnalyzer : public edm::EDProducer {

public:

  explicit ZbAnalyzer (const edm::ParameterSet &);
  ~ZbAnalyzer ();

private:

  virtual void beginJob ();
  virtual void produce (edm::Event &, const edm::EventSetup &);
  virtual void endJob ();

  virtual void beginRun (edm::Run const &, edm::EventSetup const &);
  virtual void endRun (edm::Run const &, edm::EventSetup const &);
  virtual void beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);
  virtual void endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);

#define ECALDRIVEN 0

#if ECALDRIVEN>0
  struct order_ele { bool operator() (const pat::Electron &ele1, const pat::Electron &ele2) const {
      return (ele1.ecalDrivenMomentum().pt() > ele2.ecalDrivenMomentum().pt());
    }
  };
#endif

  double btagSF(bool isMC, double flavour, double pt, double eta) {
    if (isMC == false) return 1.0;
    if (flavour == 5 || flavour == 4) return BtSF_->Val(pt, eta);
    return LtSF_->Val(pt, eta);
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
  double par_;
  bool usePartonFlavour_;
  std::string path_;

  JetCorrectionUncertainty *jecUncDT_;
  JetCorrectionUncertainty *jecUncMC_;
  edm::LumiReWeighting LumiWeights_;

  table* ElSF_;
  table* ElSF2_;
  table* MuSF_;
  table* BtSF_;
  table* LtSF_;

  TRandom3 * gRandom_;

  /****************** LEGEND *************************

   w_plotname_b => at least a b in the event, Nb > 0
   b_plotname   => b quark fraction, isb
   c_plotname   => c quark fraction, isc

   ***************************************************/

  TH1F*     h_jetmultiplicity;

  TH1F*     ecaldriven;
  TProfile* ecaldriven2;

  TH1F*     h_pu_weights;
  TH1F*     h_tracks;
  TH1F*     w_tracks;
  TH1F*     h_recoVTX;
  TH1F*     w_recoVTX;

  TH1F*     w_jetmultiplicity;
  TH1F*     b_jetmultiplicity;
  TH1F*     c_jetmultiplicity;

  TH1F*     w_first_jet_pt;	// leading jet of any type
  TH1F*     b_first_jet_pt;
  TH1F*     c_first_jet_pt;
  TH1F*     w_first_jet_eta;
  TH1F*     b_first_jet_eta;
  TH1F*     c_first_jet_eta;
  TH1F*     w_second_jet_pt;
  TH1F*     b_second_jet_pt;
  TH1F*     c_second_jet_pt;
  TH1F*     w_second_jet_eta;
  TH1F*     b_second_jet_eta;
  TH1F*     c_second_jet_eta;
  TH1F*     w_third_jet_pt;
  TH1F*     b_third_jet_pt;
  TH1F*     c_third_jet_pt;
  TH1F*     w_third_jet_eta;
  TH1F*     b_third_jet_eta;
  TH1F*     c_third_jet_eta;

  TH1F*     w_first_jet_pt_b;	// leading jet with at least one b jet in the event
  TH1F*     b_first_jet_pt_b;
  TH1F*     c_first_jet_pt_b;
  TH1F*     w_first_jet_eta_b;
  TH1F*     b_first_jet_eta_b;
  TH1F*     c_first_jet_eta_b;
  TH1F*     w_second_jet_pt_b;
  TH1F*     b_second_jet_pt_b;
  TH1F*     c_second_jet_pt_b;
  TH1F*     w_second_jet_eta_b;
  TH1F*     b_second_jet_eta_b;
  TH1F*     c_second_jet_eta_b;
  TH1F*     w_third_jet_pt_b;
  TH1F*     b_third_jet_pt_b;
  TH1F*     c_third_jet_pt_b;
  TH1F*     w_third_jet_eta_b;
  TH1F*     b_third_jet_eta_b;
  TH1F*     c_third_jet_eta_b;

  TH1F*     w_bjetmultiplicity;
  TH1F*     b_bjetmultiplicity;
  TH1F*     c_bjetmultiplicity;

  TH1F*     w_first_bjet_pt;	// leading b jet
  TH1F*     b_first_bjet_pt;
  TH1F*     c_first_bjet_pt;
  TH1F*     w_first_bjet_eta;
  TH1F*     b_first_bjet_eta;
  TH1F*     c_first_bjet_eta;

  TH1F*     w_single_bjet_pt;	// only 1 b jet
  TH1F*     b_single_bjet_pt;
  TH1F*     c_single_bjet_pt;
  TH1F*     w_single_bjet_eta;
  TH1F*     b_single_bjet_eta;
  TH1F*     c_single_bjet_eta;

  TH1F*     w_second_bjet_pt;
  TH1F*     b_second_bjet_pt;
  TH1F*     c_second_bjet_pt;
  TH1F*     w_second_bjet_eta;
  TH1F*     b_second_bjet_eta;
  TH1F*     c_second_bjet_eta;
  TH1F*     w_third_bjet_pt;
  TH1F*     b_third_bjet_pt;
  TH1F*     c_third_bjet_pt;
  TH1F*     w_third_bjet_eta;
  TH1F*     b_third_bjet_eta;
  TH1F*     c_third_bjet_eta;

  TH1F*     w_first_ele_pt;
  TH1F*     w_first_ele_pt_b;
  TH1F*     b_first_ele_pt;
  TH1F*     c_first_ele_pt;
  TH1F*     w_second_ele_pt;
  TH1F*     b_second_ele_pt;
  TH1F*     c_second_ele_pt;
  TH1F*     w_first_muon_pt;
  TH1F*     w_first_muon_pt_b;
  TH1F*     b_first_muon_pt;
  TH1F*     c_first_muon_pt;
  TH1F*     w_second_muon_pt;
  TH1F*     b_second_muon_pt;
  TH1F*     c_second_muon_pt;
  TH1F*     w_first_ele_eta;
  TH1F*     b_first_ele_eta;
  TH1F*     c_first_ele_eta;
  TH1F*     w_second_ele_eta;
  TH1F*     b_second_ele_eta;
  TH1F*     c_second_ele_eta;
  TH1F*     w_first_muon_eta;
  TH1F*     b_first_muon_eta;
  TH1F*     c_first_muon_eta;
  TH1F*     w_second_muon_eta;
  TH1F*     b_second_muon_eta;
  TH1F*     c_second_muon_eta;

  TH1F*     w_numberOfZ;
  TH1F*     b_numberOfZ;
  TH1F*     c_numberOfZ;

  TH1F*     w_mass_ee_wide;
  TH1F*     b_mass_ee_wide;
  TH1F*     c_mass_ee_wide;
  TH1F*     w_mass_mm_wide;
  TH1F*     b_mass_mm_wide;
  TH1F*     c_mass_mm_wide;

  TH1F*     h_mass_ee;
  TH1F*     w_mass_ee;
  TH1F*     b_mass_ee;
  TH1F*     c_mass_ee;
  TH1F*     h_mass_mm;
  TH1F*     w_mass_mm;
  TH1F*     b_mass_mm;
  TH1F*     c_mass_mm;
  TH1F*     w_pt_Z_ee;
  TH1F*     b_pt_Z_ee;
  TH1F*     c_pt_Z_ee;
  TH1F*     w_pt_Z_mm;
  TH1F*     b_pt_Z_mm;
  TH1F*     c_pt_Z_mm;

  TH1F*     w_single_pt_Z_ee_b;
  TH1F*     b_single_pt_Z_ee_b;
  TH1F*     c_single_pt_Z_ee_b;
  TH1F*     w_single_pt_Z_mm_b;
  TH1F*     b_single_pt_Z_mm_b;
  TH1F*     c_single_pt_Z_mm_b;

  TH1F*     w_mass_ee_b_wide;	// at least one b jet in the event
  TH1F*     b_mass_ee_b_wide;
  TH1F*     c_mass_ee_b_wide;
  TH1F*     w_mass_mm_b_wide;
  TH1F*     b_mass_mm_b_wide;
  TH1F*     c_mass_mm_b_wide;

  TH1F*     w_mass_ee_b;	// at least one b jet in the event
  TH1F*     b_mass_ee_b;
  TH1F*     c_mass_ee_b;
  TH1F*     w_mass_mm_b;
  TH1F*     b_mass_mm_b;
  TH1F*     c_mass_mm_b;
  TH1F*     w_pt_Z_ee_b;
  TH1F*     b_pt_Z_ee_b;
  TH1F*     c_pt_Z_ee_b;
  TH1F*     w_pt_Z_mm_b;
  TH1F*     b_pt_Z_mm_b;
  TH1F*     c_pt_Z_mm_b;
  TH1F*     w_delta_ee;
  TH1F*     w_delta_ee_b;
  TH1F*     b_delta_ee_b;
  TH1F*     c_delta_ee_b;
  TH1F*     w_delta_mm;
  TH1F*     w_delta_mm_b;
  TH1F*     b_delta_mm_b;
  TH1F*     c_delta_mm_b;

  TH1F*     w_single_delta_ee_b;
  TH1F*     b_single_delta_ee_b;
  TH1F*     c_single_delta_ee_b;
  TH1F*     w_single_delta_mm_b;
  TH1F*     b_single_delta_mm_b;
  TH1F*     c_single_delta_mm_b;

  TH1F*     h_secondvtx_N;
  TH1F*     w_secondvtx_N;
  TH1F*     w_secondvtx_N_zoom;
  TH1F*     w_secondvtx_N_mass;
  TH1F*     b_secondvtx_N;
  TH1F*     b_secondvtx_N_zoom;
  TH1F*     b_secondvtx_N_mass;
  TH1F*     c_secondvtx_N;
  TH1F*     c_secondvtx_N_zoom;
  TH1F*     c_secondvtx_N_mass;

  TH1F*     w_SVTX_mass_jet;
  TH1F*     b_SVTX_mass_jet;
  TH1F*     c_SVTX_mass_jet;
  TH1F*	    w_SVTX_mass_trk;
  TH1F*     b_SVTX_mass_trk;
  TH1F*     c_SVTX_mass_trk;
  TH1F*     w_SVTX_mass;
  TH1F*     b_SVTX_mass;
  TH1F*     c_SVTX_mass;

  TH1F*     w_BJP;
  TH1F*     w_JBP;
  TH1F*     b_BJP;
  TH1F*     b_JBP;
  TH1F*     c_BJP;
  TH1F*     c_JBP;

  TH1F*     w_BJP_mass;
  TH1F*     w_JBP_mass;
  TH1F*     b_BJP_mass;
  TH1F*     b_JBP_mass;
  TH1F*     c_BJP_mass;
  TH1F*     c_JBP_mass;

  TH1F*     w_Ht;
  TH1F*     b_Ht;
  TH1F*     c_Ht;

  TH1F*     w_single_Ht_b;
  TH1F*     b_single_Ht_b;
  TH1F*     c_single_Ht_b;

  TH1F*     w_Ht_b; // at least one b jet in the event
  TH1F*     b_Ht_b;
  TH1F*     c_Ht_b;

  TH1F*     w_MET;
  TH1F*     b_MET;
  TH1F*     c_MET;
  TH1F*     w_MET_sign;
  TH1F*     b_MET_sign;
  TH1F*     c_MET_sign;

  TH1F*     w_MET_b;
  TH1F*     b_MET_b;
  TH1F*     c_MET_b;
  TH1F*     w_MET_sign_b;
  TH1F*     b_MET_sign_b;
  TH1F*     c_MET_sign_b;

  TH1F*     w_Afb;

  TH1F*     h_scaleFactor_first_ele;
  TH1F*     b_scaleFactor_first_ele;
  TH1F*     h_scaleFactor_first_muon;
  TH1F*     b_scaleFactor_first_muon;
  TH1F*     h_scaleFactor_second_ele;
  TH1F*     b_scaleFactor_second_ele;
  TH1F*     h_scaleFactor_second_muon;
  TH1F*     b_scaleFactor_second_muon;

  TH1F*     h_JEC_uncert;

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
ZbAnalyzer::ZbAnalyzer (const edm::ParameterSet & iConfig) {

  pileup_ = iConfig.getUntrackedParameter < std::string > ("pileup", "S7");
  lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
  par_ =    iConfig.getUntrackedParameter <double> ("JEC", 0);
  path_ =   iConfig.getUntrackedParameter < std::string > ("path", "/gpfs/cms/users/candelis/work/ZbSkim/test");

  usePartonFlavour_ = iConfig.getUntrackedParameter <bool> ("usePartonFlavour", false);

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  h_jetmultiplicity =   fs->make < TH1F > ("h_jetmultiplicity", "h_jetmultiplicity;N_jets", 8, 0.5, 8.5);

  ecaldriven =          fs->make < TH1F > ("ecaldriven",        "ecaldriven - pf", 100, -5, 5);
  ecaldriven2 =         fs->make < TProfile > ("ecaldriven2",   "ecaldriven - pf versus pf", 100, -2.5, 2.5, -1, 1);

  h_pu_weights =        fs->make < TH1F > ("h_pu_weights",      "h_pu_weights;PU weight", 10, 0, 10);

  h_tracks =            fs->make < TH1F > ("h_tracks",          "h_tracks;N_tracks", 100, 0, 2500);
  w_tracks = 	        fs->make < TH1F > ("w_tracks",          "w_tracks;N_tracks", 100, 0, 2500);
  h_recoVTX =           fs->make < TH1F > ("h_recoVTX",         "h_recoVTX;N_vtx", 40, 0., 40.);
  w_recoVTX =           fs->make < TH1F > ("w_recoVTX",         "w_recoVTX;N_vtx", 40, 0., 40.);

  w_jetmultiplicity =   fs->make < TH1F > ("w_jetmultiplicity", "w_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  b_jetmultiplicity =   fs->make < TH1F > ("b_jetmultiplicity", "b_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  c_jetmultiplicity =   fs->make < TH1F > ("c_jetmultiplicity", "c_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  w_first_jet_pt =      fs->make < TH1F > ("w_first_jet_pt",    "w_first_jet_pt;P_t [GeV]", 50, 30., 700.);
  b_first_jet_pt =      fs->make < TH1F > ("b_first_jet_pt",    "b_first_jet_pt;P_t [GeV]", 50, 30., 700.);
  c_first_jet_pt =      fs->make < TH1F > ("c_first_jet_pt",    "c_first_jet_pt;P_t [GeV]", 50, 30., 700.);
  w_first_jet_eta =     fs->make < TH1F > ("w_first_jet_eta",   "w_first_jet_eta;Eta", 16, -2.5, 2.5);
  b_first_jet_eta =     fs->make < TH1F > ("b_first_jet_eta",   "b_first_jet_eta;Eta", 16, -2.5, 2.5);
  c_first_jet_eta =     fs->make < TH1F > ("c_first_jet_eta",   "c_first_jet_eta;Eta", 16, -2.5, 2.5);
  w_second_jet_pt =     fs->make < TH1F > ("w_second_jet_pt",   "w_second_jet_pt;P_t [GeV]", 50, 30., 500.);
  b_second_jet_pt =     fs->make < TH1F > ("b_second_jet_pt",   "b_second_jet_pt;P_t [GeV]", 50, 30., 500.);
  c_second_jet_pt =     fs->make < TH1F > ("c_second_jet_pt",   "c_second_jet_pt;P_t [GeV]", 50, 30., 500.);
  w_second_jet_eta =    fs->make < TH1F > ("w_second_jet_eta",  "w_second_jet_eta;Eta", 16, -2.5, 2.5);
  b_second_jet_eta =    fs->make < TH1F > ("b_second_jet_eta",  "b_second_jet_eta;Eta", 16, -2.5, 2.5);
  c_second_jet_eta =    fs->make < TH1F > ("c_second_jet_eta",  "c_second_jet_eta;Eta", 16, -2.5, 2.5);
  w_third_jet_pt =      fs->make < TH1F > ("w_third_jet_pt",    "w_third_jet_pt;P_t [GeV]", 50, 30., 200.);
  b_third_jet_pt =      fs->make < TH1F > ("b_third_jet_pt",    "b_third_jet_pt;P_t [GeV]", 50, 30., 200.);
  c_third_jet_pt =      fs->make < TH1F > ("c_third_jet_pt",    "c_third_jet_pt;P_t [GeV]", 50, 30., 200.);
  w_third_jet_eta =     fs->make < TH1F > ("w_third_jet_eta",   "w_third_jet_eta;Eta", 16, -2.5, 2.5);
  b_third_jet_eta =     fs->make < TH1F > ("b_third_jet_eta",   "b_third_jet_eta;Eta", 16, -2.5, 2.5);
  c_third_jet_eta =     fs->make < TH1F > ("c_third_jet_eta",   "c_third_jet_eta;Eta", 16, -2.5, 2.5);

  w_first_jet_pt_b =    fs->make < TH1F > ("w_first_jet_pt_b",   "w_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  b_first_jet_pt_b =    fs->make < TH1F > ("b_first_jet_pt_b",   "b_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  c_first_jet_pt_b =    fs->make < TH1F > ("c_first_jet_pt_b",   "c_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  w_first_jet_eta_b =   fs->make < TH1F > ("w_first_jet_eta_b",  "w_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  b_first_jet_eta_b =   fs->make < TH1F > ("b_first_jet_eta_b",  "b_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  c_first_jet_eta_b =   fs->make < TH1F > ("c_first_jet_eta_b",  "c_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  w_second_jet_pt_b =   fs->make < TH1F > ("w_second_jet_pt_b",  "w_second_jet_pt_b;P_t [GeV]", 50, 30., 500.);
  b_second_jet_pt_b =   fs->make < TH1F > ("b_second_jet_pt_b",  "b_second_jet_pt_b;P_t [GeV]", 50, 30., 500.);
  c_second_jet_pt_b =   fs->make < TH1F > ("c_second_jet_pt_b",  "c_second_jet_pt_b;P_t [GeV]", 50, 30., 500.);
  w_second_jet_eta_b =  fs->make < TH1F > ("w_second_jet_eta_b", "w_second_jet_eta_b;Eta", 16, -2.5, 2.5);
  b_second_jet_eta_b =  fs->make < TH1F > ("b_second_jet_eta_b", "b_second_jet_eta_b;Eta", 16, -2.5, 2.5);
  c_second_jet_eta_b =  fs->make < TH1F > ("c_second_jet_eta_b", "c_second_jet_eta_b;Eta", 16, -2.5, 2.5);
  w_third_jet_pt_b =    fs->make < TH1F > ("w_third_jet_pt_b",   "w_third_jet_pt_b;P_t [GeV]", 50, 30., 200.);
  b_third_jet_pt_b =    fs->make < TH1F > ("b_third_jet_pt_b",   "b_third_jet_pt_b;P_t [GeV]", 50, 30., 200.);
  c_third_jet_pt_b =    fs->make < TH1F > ("c_third_jet_pt_b",   "c_third_jet_pt_b;P_t [GeV]", 50, 30., 200.);
  w_third_jet_eta_b =   fs->make < TH1F > ("w_third_jet_eta_b",  "w_third_jet_eta_b;Eta", 16, -2.5, 2.5);
  b_third_jet_eta_b =   fs->make < TH1F > ("b_third_jet_eta_b",  "b_third_jet_eta_b;Eta", 16, -2.5, 2.5);
  c_third_jet_eta_b =   fs->make < TH1F > ("c_third_jet_eta_b",  "c_third_jet_eta_b;Eta", 16, -2.5, 2.5);

  w_bjetmultiplicity =  fs->make < TH1F > ("w_bjetmultiplicity", "w_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  b_bjetmultiplicity =  fs->make < TH1F > ("b_bjetmultiplicity", "b_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  c_bjetmultiplicity =  fs->make < TH1F > ("c_bjetmultiplicity", "c_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  w_first_bjet_pt =     fs->make < TH1F > ("w_first_bjet_pt",    "w_first_bjet_pt;P_t [GeV]", 50, 30., 700.);
  b_first_bjet_pt =     fs->make < TH1F > ("b_first_bjet_pt",    "b_first_bjet_pt;P_t [GeV]", 50, 30., 700.);
  c_first_bjet_pt =     fs->make < TH1F > ("c_first_bjet_pt",    "c_first_bjet_pt;P_t [GeV]", 50, 30., 700.);
  w_first_bjet_eta =    fs->make < TH1F > ("w_first_bjet_eta",   "w_first_bjet_eta;Eta", 16, -2.5, 2.5);
  b_first_bjet_eta =    fs->make < TH1F > ("b_first_bjet_eta",   "b_first_bjet_eta;Eta", 16, -2.5, 2.5);
  c_first_bjet_eta =    fs->make < TH1F > ("c_first_bjet_eta",   "c_first_bjet_eta;Eta", 16, -2.5, 2.5);
  w_single_bjet_pt =    fs->make < TH1F > ("w_single_bjet_pt",    "w_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  b_single_bjet_pt =    fs->make < TH1F > ("b_single_bjet_pt",    "b_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  c_single_bjet_pt =    fs->make < TH1F > ("c_single_bjet_pt",    "c_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  w_single_bjet_eta =   fs->make < TH1F > ("w_single_bjet_eta",   "w_single_bjet_eta;Eta", 16, -2.5, 2.5);
  b_single_bjet_eta =   fs->make < TH1F > ("b_single_bjet_eta",   "b_single_bjet_eta;Eta", 16, -2.5, 2.5);
  c_single_bjet_eta =   fs->make < TH1F > ("c_single_bjet_eta",   "c_single_bjet_eta;Eta", 16, -2.5, 2.5);
  w_second_bjet_pt =    fs->make < TH1F > ("w_second_bjet_pt",   "w_second_bjet_pt;P_t [GeV]", 50, 30., 500.);
  b_second_bjet_pt =    fs->make < TH1F > ("b_second_bjet_pt",   "b_second_bjet_pt;P_t [GeV]", 50, 30., 500.);
  c_second_bjet_pt =    fs->make < TH1F > ("c_second_bjet_pt",   "c_second_bjet_pt;P_t [GeV]", 50, 30., 500.);
  w_second_bjet_eta =   fs->make < TH1F > ("w_second_bjet_eta",  "w_second_bjet_eta;Eta", 16, -2.5, 2.5);
  b_second_bjet_eta =   fs->make < TH1F > ("b_second_bjet_eta",  "b_second_bjet_eta;Eta", 16, -2.5, 2.5);
  c_second_bjet_eta =   fs->make < TH1F > ("c_second_bjet_eta",  "c_second_bjet_eta;Eta", 16, -2.5, 2.5);
  w_third_bjet_pt =     fs->make < TH1F > ("w_third_bjet_pt",    "w_third_bjet_pt;P_t [GeV]", 50, 30., 200.);
  b_third_bjet_pt =     fs->make < TH1F > ("b_third_bjet_pt",    "b_third_bjet_pt;P_t [GeV]", 50, 30., 200.);
  c_third_bjet_pt =     fs->make < TH1F > ("c_third_bjet_pt",    "c_third_bjet_pt;P_t [GeV]", 50, 30., 200.);
  w_third_bjet_eta =    fs->make < TH1F > ("w_third_bjet_eta",   "w_third_bjet_eta;Eta", 16, -2.5, 2.5);
  b_third_bjet_eta =    fs->make < TH1F > ("b_third_bjet_eta",   "b_third_bjet_eta;Eta", 16, -2.5, 2.5);
  c_third_bjet_eta =    fs->make < TH1F > ("c_third_bjet_eta",   "c_third_bjet_eta;Eta", 16, -2.5, 2.5);

  w_first_ele_pt =      fs->make < TH1F > ("w_first_ele_pt",    "w_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_first_ele_pt_b =    fs->make < TH1F > ("w_first_ele_pt_b",  "w_first_ele_pt_b;P_t [GeV]", 50, 0., 450.);
  b_first_ele_pt =      fs->make < TH1F > ("b_first_ele_pt",    "b_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  c_first_ele_pt =      fs->make < TH1F > ("c_first_ele_pt",    "c_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_second_ele_pt =     fs->make < TH1F > ("w_second_ele_pt",   "w_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  b_second_ele_pt =     fs->make < TH1F > ("b_second_ele_pt",   "b_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  c_second_ele_pt =     fs->make < TH1F > ("c_second_ele_pt",   "c_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_first_muon_pt =     fs->make < TH1F > ("w_first_muon_pt",   "w_first_muon_pt;P_t [GeV]", 50, 0., 450.);
  w_first_muon_pt_b =   fs->make < TH1F > ("w_first_muon_pt_b", "w_first_muon_pt_b [GeV]", 50, 0., 450.);
  b_first_muon_pt =     fs->make < TH1F > ("b_first_muon_pt",   "b_first_muon_pt;P_t [GeV]", 50, 0., 450.);
  c_first_muon_pt =     fs->make < TH1F > ("c_first_muon_pt",   "c_first_muon_pt;P_t [GeV]", 50, 0., 450.);
  w_second_muon_pt =    fs->make < TH1F > ("w_second_muon_pt",  "w_second_muon_pt;P_t [GeV]", 50, 0., 450.);
  b_second_muon_pt =    fs->make < TH1F > ("b_second_muon_pt",  "b_second_muon_pt;P_t [GeV]", 50, 0., 450.);
  c_second_muon_pt =    fs->make < TH1F > ("c_second_muon_pt",  "c_second_muon_pt;P_t [GeV]", 50, 0., 450.);
  w_first_ele_eta =     fs->make < TH1F > ("w_first_ele_eta",   "w_first_ele_eta;Eta", 16, -2.5, 2.5);
  b_first_ele_eta =     fs->make < TH1F > ("b_first_ele_eta",   "b_first_ele_eta;Eta", 16, -2.5, 2.5);
  c_first_ele_eta =     fs->make < TH1F > ("c_first_ele_eta",   "c_first_ele_eta;Eta", 16, -2.5, 2.5);
  w_second_ele_eta =    fs->make < TH1F > ("w_second_ele_eta",  "w_second_ele_eta;Eta", 16, -2.5, 2.5);
  b_second_ele_eta =    fs->make < TH1F > ("b_second_ele_eta",  "b_second_ele_eta;Eta", 16, -2.5, 2.5);
  c_second_ele_eta =    fs->make < TH1F > ("c_second_ele_eta",  "c_second_ele_eta;Eta", 16, -2.5, 2.5);
  w_first_muon_eta =    fs->make < TH1F > ("w_first_muon_eta",  "w_first_muon_eta;Eta", 16, -2.5, 2.5);
  b_first_muon_eta =    fs->make < TH1F > ("b_first_muon_eta",  "b_first_muon_eta;Eta", 16, -2.5, 2.5);
  c_first_muon_eta =    fs->make < TH1F > ("c_first_muon_eta",  "c_first_muon_eta;Eta", 16, -2.5, 2.5);
  w_second_muon_eta =   fs->make < TH1F > ("w_second_muon_eta", "w_second_muon_eta;Eta", 16, -2.5, 2.5);
  b_second_muon_eta =   fs->make < TH1F > ("b_second_muon_eta", "b_second_muon_eta;Eta", 16, -2.5, 2.5);
  c_second_muon_eta =   fs->make < TH1F > ("c_second_muon_eta", "c_second_muon_eta;Eta", 16, -2.5, 2.5);

  w_numberOfZ =         fs->make < TH1F > ("w_numberOfZ",       "w_numberOfZ;N_Z", 5, 0, 5);
  b_numberOfZ =         fs->make < TH1F > ("b_numberOfZ",       "b_numberOfZ;N_Z", 5, 0, 5);
  c_numberOfZ =         fs->make < TH1F > ("c_numberOfZ",       "c_numberOfZ;N_Z", 5, 0, 5);

  w_mass_ee_wide =      fs->make < TH1F > ("w_mass_ee_wide",    "w_mass_ee_wide;Mass [GeV]", 80, 60., 120.);
  b_mass_ee_wide =      fs->make < TH1F > ("b_mass_ee_wide",    "b_mass_ee_wide;Mass [GeV]", 80, 60., 120.);
  c_mass_ee_wide =      fs->make < TH1F > ("c_mass_ee_wide",    "c_mass_ee_wide;Mass [GeV]", 80, 60., 120.);
  w_mass_mm_wide =      fs->make < TH1F > ("w_mass_mm_wide",    "w_mass_mm_wide;Mass [GeV]", 80, 60., 120.);
  b_mass_mm_wide =      fs->make < TH1F > ("b_mass_mm_wide",    "b_mass_mm_wide;Mass [GeV]", 80, 60., 120.);
  c_mass_mm_wide =      fs->make < TH1F > ("c_mass_mm_wide",    "c_mass_mm_wide;Mass [GeV]", 80, 60., 120.);

  h_mass_ee =           fs->make < TH1F > ("h_mass_ee",         "h_mass_ee;Mass [GeV]", 80, 71., 111.);
  w_mass_ee = 	        fs->make < TH1F > ("w_mass_ee",         "w_mass_ee;Mass [GeV]", 80, 71., 111.);
  b_mass_ee =           fs->make < TH1F > ("b_mass_ee",         "b_mass_ee;Mass [GeV]", 80, 71., 111.);
  c_mass_ee =           fs->make < TH1F > ("c_mass_ee",         "c_mass_ee;Mass [GeV]", 80, 71., 111.);
  h_mass_mm =           fs->make < TH1F > ("h_mass_mm",         "h_mass_mm;Mass [GeV]", 80, 71., 111.);
  w_mass_mm = 	        fs->make < TH1F > ("w_mass_mm",         "w_mass_mm;Mass [GeV]", 80, 71., 111.);
  b_mass_mm =           fs->make < TH1F > ("b_mass_mm",         "b_mass_mm;Mass [GeV]", 80, 71., 111.);
  c_mass_mm =           fs->make < TH1F > ("c_mass_mm",         "c_mass_mm;Mass [GeV]", 80, 71., 111.);
  w_pt_Z_ee =           fs->make < TH1F > ("w_pt_Z_ee",         "w_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_ee =           fs->make < TH1F > ("b_pt_Z_ee",         "b_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_ee =           fs->make < TH1F > ("c_pt_Z_ee",         "c_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_mm =           fs->make < TH1F > ("w_pt_Z_mm",         "w_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_mm =           fs->make < TH1F > ("b_pt_Z_mm",         "b_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_mm =           fs->make < TH1F > ("c_pt_Z_mm",         "c_pt_Z_mm;P_t [GeV]", 40, 0., 400.);

  w_mass_ee_b_wide =    fs->make < TH1F > ("w_mass_ee_b_wide",  "w_mass_ee_b_wide;Mass [GeV]", 80, 60., 120.);
  b_mass_ee_b_wide =    fs->make < TH1F > ("b_mass_ee_b_wide",  "b_mass_ee_b_wide;Mass [GeV]", 80, 60., 120.);
  c_mass_ee_b_wide =    fs->make < TH1F > ("c_mass_ee_b_wide",  "c_mass_ee_b_wide;Mass [GeV]", 80, 60., 120.);
  w_mass_mm_b_wide =    fs->make < TH1F > ("w_mass_mm_b_wide",  "w_mass_mm_b_wide;Mass [GeV]", 80, 60., 120.);
  b_mass_mm_b_wide =    fs->make < TH1F > ("b_mass_mm_b_wide",  "b_mass_mm_b_wide;Mass [GeV]", 80, 60., 120.);
  c_mass_mm_b_wide =    fs->make < TH1F > ("c_mass_mm_b_wide",  "c_mass_mm_b_wide;Mass [GeV]", 80, 60., 120.);

  w_mass_ee_b =         fs->make < TH1F > ("w_mass_ee_b",       "w_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  b_mass_ee_b =         fs->make < TH1F > ("b_mass_ee_b",       "b_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  c_mass_ee_b =         fs->make < TH1F > ("c_mass_ee_b",       "c_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  w_mass_mm_b =         fs->make < TH1F > ("w_mass_mm_b",       "w_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  b_mass_mm_b =         fs->make < TH1F > ("b_mass_mm_b",       "b_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  c_mass_mm_b =         fs->make < TH1F > ("c_mass_mm_b",       "c_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  w_pt_Z_ee_b =         fs->make < TH1F > ("w_pt_Z_ee_b",       "w_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_ee_b =         fs->make < TH1F > ("b_pt_Z_ee_b",       "b_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_ee_b =         fs->make < TH1F > ("c_pt_Z_ee_b",       "c_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_mm_b =         fs->make < TH1F > ("w_pt_Z_mm_b",       "w_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_mm_b =         fs->make < TH1F > ("b_pt_Z_mm_b",       "b_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_mm_b =         fs->make < TH1F > ("c_pt_Z_mm_b",       "c_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  w_single_pt_Z_ee_b =  fs->make < TH1F > ("w_single_pt_Z_ee_b",       "w_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  b_single_pt_Z_ee_b =  fs->make < TH1F > ("b_single_pt_Z_ee_b",       "b_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  c_single_pt_Z_ee_b =  fs->make < TH1F > ("c_single_pt_Z_ee_b",       "c_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  w_single_pt_Z_mm_b =  fs->make < TH1F > ("w_single_pt_Z_mm_b",       "w_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  b_single_pt_Z_mm_b =  fs->make < TH1F > ("b_single_pt_Z_mm_b",       "b_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  c_single_pt_Z_mm_b =  fs->make < TH1F > ("c_single_pt_Z_mm_b",       "c_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  w_delta_ee =          fs->make < TH1F > ("w_delta_phi_ee",    "w_delta_phi_ee", 12, 0, TMath::Pi ());
  w_delta_ee =          fs->make < TH1F > ("w_delta_phi_ee",    "w_delta_phi_ee", 12, 0, TMath::Pi ());
  w_delta_ee_b =        fs->make < TH1F > ("w_delta_phi_ee_b",  "w_delta_phi_ee_b", 12, 0, TMath::Pi ());
  b_delta_ee_b =        fs->make < TH1F > ("b_delta_phi_ee_b",  "b_delta_phi_ee_b", 12, 0, TMath::Pi ());
  c_delta_ee_b =        fs->make < TH1F > ("c_delta_phi_ee_b",  "c_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_delta_mm =          fs->make < TH1F > ("w_delta_phi_mm",    "w_delta_phi_mm", 12, 0, TMath::Pi ());
  w_delta_mm_b =        fs->make < TH1F > ("w_delta_phi_mm_b",  "w_delta_phi_mm_b", 12, 0, TMath::Pi ());
  b_delta_mm_b =        fs->make < TH1F > ("b_delta_phi_mm_b",  "b_delta_phi_mm_b", 12, 0, TMath::Pi ());
  c_delta_mm_b =        fs->make < TH1F > ("c_delta_phi_mm_b",  "c_delta_phi_mm_b", 12, 0, TMath::Pi ());

  w_single_delta_ee_b =        fs->make < TH1F > ("w_single_delta_phi_ee_b",  "w_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  b_single_delta_ee_b =        fs->make < TH1F > ("b_single_delta_phi_ee_b",  "b_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  c_single_delta_ee_b =        fs->make < TH1F > ("c_single_delta_phi_ee_b",  "c_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_single_delta_mm_b =        fs->make < TH1F > ("w_single_delta_phi_mm_b",  "w_single_delta_phi_mm_b", 12, 0, TMath::Pi ());
  b_single_delta_mm_b =        fs->make < TH1F > ("b_single_delta_phi_mm_b",  "b_single_delta_phi_mm_b", 12, 0, TMath::Pi ());
  c_single_delta_mm_b =        fs->make < TH1F > ("c_single_delta_phi_mm_b",  "c_single_delta_phi_mm_b", 12, 0, TMath::Pi ());

  h_secondvtx_N =         fs->make < TH1F > ("h_secondvtx_N",     "h_secondvtx_N", 50, 0, 1);
  w_secondvtx_N =         fs->make < TH1F > ("w_secondvtx_N",     "w_secondvtx_N", 50, 0, 1);
  w_secondvtx_N_zoom =    fs->make < TH1F > ("w_secondvtx_N_zoom",  "w_secondvtx_N_zoom", 20, 0.898, 1);
  w_secondvtx_N_mass =    fs->make < TH1F > ("w_secondvtx_N_mass",  "w_secondvtx_N_mass", 20, 0.898, 1);

  b_secondvtx_N =         fs->make < TH1F > ("b_secondvtx_N",       "b_secondvtx_N", 50, 0, 1);
  b_secondvtx_N_zoom =    fs->make < TH1F > ("b_secondvtx_N_zoom",  "b_secondvtx_N_zoom", 20, 0.898, 1);
  b_secondvtx_N_mass =    fs->make < TH1F > ("b_secondvtx_N_mass",  "b_secondvtx_N_mass", 20, 0.898, 1);

  c_secondvtx_N =         fs->make < TH1F > ("c_secondvtx_N",       "c_secondvtx_N", 50, 0, 1);
  c_secondvtx_N_zoom =    fs->make < TH1F > ("c_secondvtx_N_zoom",  "c_secondvtx_N_zoom", 20, 0.898, 1);
  c_secondvtx_N_mass =    fs->make < TH1F > ("c_secondvtx_N_mass",  "c_secondvtx_N_mass", 20, 0.898, 1);

  w_SVTX_mass_jet =     fs->make < TH1F > ("w_SVTX_mass_jet",   "w_SVTX_mass_jet;Mass [GeV]", 80, 0, 6);
  b_SVTX_mass_jet =     fs->make < TH1F > ("b_SVTX_mass_jet",   "b_SVTX_mass_jet;Mass [GeV]", 80, 0, 6);
  c_SVTX_mass_jet =     fs->make < TH1F > ("c_SVTX_mass_jet",   "c_SVTX_mass_jet;Mass [GeV]", 80, 0, 6);
  w_SVTX_mass_trk =     fs->make < TH1F > ("w_SVTX_mass_trk",   "w_SVTX_mass_trk;Mass [GeV]", 160, 0, 80);
  b_SVTX_mass_trk =     fs->make < TH1F > ("b_SVTX_mass_trk",   "b_SVTX_mass_trk;Mass [GeV]", 160, 0, 80);
  c_SVTX_mass_trk =     fs->make < TH1F > ("c_SVTX_mass_trk",   "c_SVTX_mass_trk;Mass [GeV]", 160, 0, 80);
  w_SVTX_mass     =     fs->make < TH1F > ("w_SVTX_mass",       "w_SVTX_mass;Mass [GeV]", 80, 0, 6);
  b_SVTX_mass     =     fs->make < TH1F > ("b_SVTX_mass",       "b_SVTX_mass;Mass [GeV]", 80, 0, 6);
  c_SVTX_mass     =     fs->make < TH1F > ("c_SVTX_mass",       "c_SVTX_mass;Mass [GeV]", 80, 0, 6);

  w_BJP       =     fs->make < TH1F > ("w_BJP",   "w_BJP", 80, 0, 10);
  w_JBP       =     fs->make < TH1F > ("w_JBP",   "w_JBP", 50, 0, 3);
  b_BJP       =     fs->make < TH1F > ("b_BJP",   "b_BJP", 80, 0, 10);
  b_JBP       =     fs->make < TH1F > ("b_JBP",   "b_JBP", 50, 0, 3);
  c_BJP       =     fs->make < TH1F > ("c_BJP",   "c_BJP", 80, 0, 10);
  c_JBP       =     fs->make < TH1F > ("c_JBP",   "c_JBP", 50, 0, 3);

  w_BJP_mass  =     fs->make < TH1F > ("w_BJP_mass",   "w_BJP_mass", 80, 0, 10);
  w_JBP_mass  =     fs->make < TH1F > ("w_JBP_mass",   "w_JBP_mass", 50, 0, 3);
  b_BJP_mass  =     fs->make < TH1F > ("b_BJP_mass",   "b_BJP_mass", 80, 0, 10);
  b_JBP_mass  =     fs->make < TH1F > ("b_JBP_mass",   "b_JBP_mass", 50, 0, 3);
  c_BJP_mass  =     fs->make < TH1F > ("c_BJP_mass",   "c_BJP_mass", 80, 0, 10);
  c_JBP_mass  =     fs->make < TH1F > ("c_JBP_mass",   "c_JBP_mass", 50, 0, 3);

  w_Ht =                fs->make < TH1F > ("w_Ht",              "w_Ht [GeV]", 50, 30., 1000.);
  b_Ht =                fs->make < TH1F > ("b_Ht",              "b_Ht [GeV]", 50, 30., 1000.);
  c_Ht =                fs->make < TH1F > ("c_Ht",              "c_Ht [GeV]", 50, 30., 1000.);
  w_Ht_b =              fs->make < TH1F > ("w_Ht_b",            "w_Ht [GeV]", 50, 30., 1000.);
  b_Ht_b =              fs->make < TH1F > ("b_Ht_b",            "b_Ht [GeV]", 50, 30., 1000.);
  c_Ht_b =              fs->make < TH1F > ("c_Ht_b",            "c_Ht [GeV]", 50, 30., 1000.);

  w_single_Ht_b =       fs->make < TH1F > ("w_single_Ht_b",            "w_single_Ht [GeV]", 50, 30., 1000.);
  b_single_Ht_b =       fs->make < TH1F > ("b_single_Ht_b",            "b_single_Ht [GeV]", 50, 30., 1000.);
  c_single_Ht_b =       fs->make < TH1F > ("c_single_Ht_b",            "c_single_Ht [GeV]", 50, 30., 1000.);

  w_MET =               fs->make < TH1F > ("w_MET",             "w_MET;MET [GeV]", 50, 0., 250.);
  b_MET =               fs->make < TH1F > ("b_MET",             "b_MET;MET [GeV]", 50, 0., 250.);
  c_MET =               fs->make < TH1F > ("c_MET",             "c_MET;MET [GeV]", 50, 0., 250.);
  w_MET_sign = 	        fs->make < TH1F > ("w_MET_sign",        "w_MET_sign;MET significance [GeV]", 50, 0., 100.);
  b_MET_sign = 	        fs->make < TH1F > ("b_MET_sign",        "b_MET_sign;MET significance [GeV]", 50, 0., 100.);
  c_MET_sign = 	        fs->make < TH1F > ("c_MET_sign",        "c_MET_sign;MET significance [GeV]", 50, 0., 100.);

  w_MET_b =             fs->make < TH1F > ("w_MET_b",         "w_MET_b;MET [GeV]", 50, 0., 250.);
  b_MET_b =             fs->make < TH1F > ("b_MET_b",         "b_MET_b;MET [GeV]", 50, 0., 250.);
  c_MET_b =             fs->make < TH1F > ("c_MET_b",         "c_MET_b;MET [GeV]", 50, 0., 250.);
  w_MET_sign_b = 	fs->make < TH1F > ("w_MET_sign_b",    "w_MET_sign_b;MET significance [GeV]", 50, 0., 100.);
  b_MET_sign_b = 	fs->make < TH1F > ("b_MET_sign_b",    "b_MET_sign_b;MET significance [GeV]", 50, 0., 100.);
  c_MET_sign_b = 	fs->make < TH1F > ("c_MET_sign_b",    "c_MET_sign_b;MET significance [GeV]", 50, 0., 100.);

  w_Afb =               fs->make < TH1F > ("b_asymmetry",       "b_asymmetry", 10, -1, 1);

  h_JEC_uncert =        fs->make < TH1F > ("JEC uncert", "JEC uncert", 10, -0.5, 0.5);

  h_scaleFactor_first_ele =   fs->make < TH1F > ("h_scaleFactor_first_ele",   "h_scaleFactor_first_ele", 50, 0.95, 1.05);
  b_scaleFactor_first_ele =   fs->make < TH1F > ("b_scaleFactor_first_ele",   "b_scaleFactor_first_ele", 50, 0.95, 1.05);
  h_scaleFactor_first_muon =  fs->make < TH1F > ("h_scaleFactor_first_muon",  "h_scaleFactor_first_muon", 50, 0.95, 1.05);
  b_scaleFactor_first_muon =  fs->make < TH1F > ("b_scaleFactor_first_muon",  "b_scaleFactor_first_muon", 50, 0.95, 1.05);
  h_scaleFactor_second_ele =  fs->make < TH1F > ("h_scaleFactor_second_ele",  "h_scaleFactor_second_ele", 50, 0.95, 1.05);
  b_scaleFactor_second_ele =  fs->make < TH1F > ("b_scaleFactor_second_ele",  "b_scaleFactor_second_ele", 50, 0.95, 1.05);
  h_scaleFactor_second_muon = fs->make < TH1F > ("h_scaleFactor_second_muon", "h_scaleFactor_second_muon", 50, 0.95, 1.05);
  b_scaleFactor_second_muon = fs->make < TH1F > ("b_scaleFactor_second_muon", "b_scaleFactor_second_muon", 50, 0.95, 1.05);

  produces<std::vector<double>>("myEventWeight");

  produces<std::vector<math::XYZTLorentzVector>>("myElectrons");
  produces<std::vector<math::XYZTLorentzVector>>("myMuons");

  produces<std::vector<double>>("myPtZ");

  produces<std::vector<math::XYZTLorentzVector>>("myJets");
  produces<std::vector<double>>("myDeltaPhi");

  produces<std::vector<double>>("myHt");

  produces<std::vector<double>>("myBjetsWeights");

  produces<std::vector<math::XYZTLorentzVector>>("myBjets");
  produces<std::vector<double>>("myBdeltaPhi");

}

ZbAnalyzer::~ZbAnalyzer () {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event ------------
void ZbAnalyzer::produce (edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace edm;
  using namespace std;

  // Get muon collection
  edm::Handle < pat::MuonCollection > muons;
  iEvent.getByLabel ("matchedMuons", muons);

  // Get muon collection
  edm::Handle < pat::ElectronCollection > electrons;
  iEvent.getByLabel ("matchedElectrons", electrons);

  // Get jet collection
  edm::Handle < vector < pat::Jet > > jets;
  iEvent.getByLabel ("goodJets", jets);

  // Get the Z->mm collection
  edm::Handle < reco::CompositeCandidateCollection > zmm;
  iEvent.getByLabel ("zmuMatchedmuMatched", zmm);

  // Get the Z->ee collection
  edm::Handle < reco::CompositeCandidateCollection > zee;
  iEvent.getByLabel ("zeleMatchedeleMatched", zee);

  // Get tracks
  edm::Handle < vector < reco::Track > > tracks;
  iEvent.getByLabel ("generalTracks", tracks);

  // Get METs
  edm::Handle < vector < pat::MET > > mets;
  iEvent.getByLabel (edm::InputTag ("patMETsPFlow"), mets);

  edm::Handle<vector<reco::GenParticle> > genPart;
  iEvent.getByLabel ("genParticles", genPart);

  edm::Handle<vector<reco::GenJet> > gJets;
  iEvent.getByLabel(edm::InputTag("goodJets","genJets"), gJets);

  std::auto_ptr<std::vector<double>> myEventWeight( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myElectrons( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myMuons( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myPtZ( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myDeltaPhi( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myHt( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myBjetsWeights( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBjets( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myBdeltaPhi( new std::vector<double> );

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

  double Ht = 0;

  double MyWeight = 1;

  double scalFac_first_e = 1;
  double scalFac_second_e = 1;
  double scalFac_first_m = 1;
  double scalFac_second_m = 1;
  double scalFac_b = 1;

  // ++++++ Pile-Up

  bool isMC = false;
  JetCorrectionUncertainty* jecUnc;

  jecUnc = jecUncDT_;

  MyWeight = 1.0;

  Handle < vector < PileupSummaryInfo > > PupInfo;

  if (iEvent.getByLabel (edm::InputTag ("addPileupInfo"), PupInfo))  {

    isMC = true;
    jecUnc = jecUncMC_;

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

    MyWeight = MyWeight*mcWeight;

  }

  // +++++++++ ELECTRONS

  vector < pat::Electron > vect_ele;

  for (pat::ElectronCollection::const_iterator ele = electrons->begin (); ele != electrons->end (); ++ele) {

#if ECALDRIVEN>0
    if (ele->ecalDrivenMomentum().pt()>=20) {
#else
    if (ele->pt()>=20) {
#endif
      vect_ele.push_back (*ele);
    }

    ecaldriven->Fill(ele->pt() - ele->ecalDrivenMomentum().pt());
    ecaldriven2->Fill(ele->eta(), ele->pt() - ele->ecalDrivenMomentum().pt());

  }

#if ECALDRIVEN>0
  std::sort( vect_ele.begin(), vect_ele.end(), order_ele() );
#endif

  int iele0=0;
  int iele1=0;

  for (unsigned int i=1; i<vect_ele.size(); ++i) {
    if (vect_ele[i].charge()*vect_ele[iele0].charge()<0 && iele1==0) iele1=i;
  }

  if (iele1!=0) {
    math::XYZTLorentzVector z;
#if ECALDRIVEN>1
    z = vect_ele[iele0].ecalDrivenMomentum() + vect_ele[iele1].ecalDrivenMomentum();
#else
    z = vect_ele[iele0].p4() + vect_ele[iele1].p4();
#endif
    diele_mass = z.mass();
    diele_phi = z.phi();
    diele_pt = z.pt();
//
//    if (isMC) {
//      double smear = gRandom_->Gaus(0, 1.0);
//      diele_mass = 1.0023 * diele_mass + smear;
//      diele_pt = 1.0023 * diele_pt + smear;
//    }
//
    if (diele_mass>71 && diele_mass<111) ee_event = true;
  }

  // +++++++++ MUONS

  vector < pat::Muon > vect_muon;

  for (pat::MuonCollection::const_iterator muon = muons->begin (); muon != muons->end (); ++muon) {

    if (muon->pt()>=20) {
      vect_muon.push_back (*muon);
    }

  }

  int imuon0=0;
  int imuon1=0;

  for (unsigned int i=1; i<vect_muon.size(); ++i) {
    if (vect_muon[i].charge()*vect_muon[imuon0].charge()<0 && imuon1==0) imuon1=i;
  }

  if (imuon1!=0) {
    math::XYZTLorentzVector z;
    z = vect_muon[imuon0].p4() + vect_muon[imuon1].p4();
    dimuon_mass = z.mass();
    dimuon_phi = z.phi();
    dimuon_pt = z.pt();
    if (dimuon_mass>71 && dimuon_mass<111) mm_event = true;
  }

  ee_event = ee_event && (lepton_ == "electron");
  mm_event = mm_event && (lepton_ == "muon");

  if (isMC) {
    if (ee_event) {
      scalFac_first_e  =  ElSF_->Val (vect_ele[iele0].pt(), vect_ele[iele0].eta()) * ElSF2_->Val (vect_ele[iele0].pt(), vect_ele[iele0].eta());
      scalFac_second_e =  ElSF_->Val (vect_ele[iele1].pt(), vect_ele[iele1].eta()) * ElSF2_->Val (vect_ele[iele1].pt(), vect_ele[iele1].eta());
      MyWeight = MyWeight * scalFac_first_e * scalFac_second_e;
    }
    if (mm_event) {
      scalFac_first_m  = MuSF_->Val (vect_muon[imuon0].pt(), vect_muon[imuon0].eta());
      scalFac_second_m = MuSF_->Val (vect_muon[imuon1].pt(), vect_muon[imuon1].eta());
      MyWeight = MyWeight * scalFac_first_m * scalFac_second_m;
      //cout<<vect_muon[imuon0].pt()<<vect_muon[imuon0].eta()<< " mu  SF =" << scalFac_first_m <<endl;
    }
  }

  // ++++++++ VERTICES

  bool vtx_cut = true;

  edm::Handle < vector < reco::Vertex > > vertices;
  iEvent.getByLabel (edm::InputTag ("goodOfflinePrimaryVertices"), vertices);

  if (vertices->size() > 0) {
    const reco::Vertex* theVertex = &(vertices->front());
    if (theVertex->ndof() < 5) vtx_cut = false;
    if (fabs(theVertex->z()) > 24.0) vtx_cut = false;
    if (fabs(theVertex->position().rho()) > 2.0) vtx_cut = false;
  } else {
    vtx_cut = false;
  }

  int NVtx = 0;

  if (vtx_cut) {
    for (vector < reco::Vertex >::const_iterator itv = vertices->begin (); itv != vertices->end (); ++itv) {
      if (itv->ndof() < 5) continue;
      if (fabs(itv->z()) > 50.0) continue;
      if (fabs(itv->position().rho()) > 2.0) continue;
      ++NVtx;
    }
  }

  // ++++++++ LOOP OVER GEN PARTICLES

  bool isb = false;
  bool isc = false;

  if (isMC) {
    for (std::vector <reco::GenParticle>::const_iterator thepart = genPart->begin(); thepart != genPart->end(); thepart++) {
      if ((int) (abs(thepart->pdgId() / 100) % 10) == 5 || (int) (abs(thepart->pdgId() / 1000) % 10) == 5) {
        isb = true;
      }
      if ((int) (abs(thepart->pdgId() / 100) % 10 ) == 4 || (int) (abs(thepart->pdgId() / 1000) % 10) == 4) {
        isc = true;
      }
    }
  }

  // ++++++++ MET CUT

   bool met_cut = mets->empty() ? true : (*mets)[0].significance() < 30.;

  // ++++++++ JETS

  vector < pat::Jet > vect_jets;
  vector < pat::Jet > vect_bjets;

  double sumVertexMassJet = 0.;
  double sumVertexMassTrk = 0.;
  double sumVertexMass = 0.;

  for (vector < pat::Jet >::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {

    double jet_pt  = jet->pt ();
    double jet_eta = jet->eta();

    // JEC Uncertainty

    jecUnc->setJetPt(jet_pt);
    jecUnc->setJetEta(jet_eta);
    double unc = jecUnc->getUncertainty(true);
    double cor = (1.0+unc*par_);
    h_JEC_uncert->Fill (unc);
    //cout<< "JEC syst =" << unc << endl;

    jet_pt = jet_pt * cor;

    if (jet_pt > 30) {

      ++Nj;

      Ht += jet_pt;

      vect_jets.push_back (*jet);

      double discrCSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
      //cout << discrCSV << endl;

      if ((ee_event || mm_event) && vtx_cut && met_cut) {
        h_secondvtx_N->Fill (discrCSV);
        w_secondvtx_N->Fill (discrCSV, MyWeight);
	if (isb) {
	  b_secondvtx_N->Fill (discrCSV, MyWeight);
	}
	if (isc && !isb) {
	  c_secondvtx_N->Fill (discrCSV, MyWeight);
	}
      }

      if (discrCSV > 0.898 || (usePartonFlavour_ && jet->partonFlavour()==5)) {

	++Nb;
	//cout << Nb << endl;

        vect_bjets.push_back (*jet);

        scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
        if ((ee_event || mm_event) && vtx_cut && met_cut) {
	  w_secondvtx_N_zoom->Fill (discrCSV, MyWeight*scalFac_b);
	  if (isb) {
	    b_secondvtx_N_zoom->Fill (discrCSV, MyWeight*scalFac_b);
	  }
	  if (isc && !isb) {
	    c_secondvtx_N_zoom->Fill (discrCSV, MyWeight*scalFac_b);
	  }
        }

        reco::SecondaryVertexTagInfo const * svTagInfos = jet->tagInfoSecondaryVertex("secondaryVertex");

        if ((ee_event || mm_event) && vtx_cut && met_cut) {
	  if (svTagInfos && svTagInfos->nVertices() > 0) {
	    ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > sumVecJet;
	    for (reco::Vertex::trackRef_iterator track = svTagInfos->secondaryVertex(0).tracks_begin(); track != svTagInfos->secondaryVertex(0).tracks_end(); ++track) {
	      const double kPionMass = 0.13957018;
	      ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > vec;
	      vec.SetPx( (*track)->px() );
	      vec.SetPy( (*track)->py() );
	      vec.SetPz( (*track)->pz() );
	      vec.SetM (kPionMass);
	      sumVecJet += vec;
	    }
	    sumVertexMassJet += sumVecJet.M();
	  }
        }

        if ((ee_event || mm_event) && vtx_cut && met_cut) {
	  ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > sumVecTrk;
	  for (size_t itrack=0; itrack < jet->associatedTracks().size(); ++itrack) {
	    const double kPionMass = 0.13957018;
	    ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > vec;
	    vec.SetPx( jet->associatedTracks()[itrack]->px() );
	    vec.SetPy( jet->associatedTracks()[itrack]->py() );
	    vec.SetPz( jet->associatedTracks()[itrack]->pz() );
	    vec.SetM (kPionMass);
	    sumVecTrk += vec;
	  }
	  sumVertexMassTrk += sumVecTrk.M();
        }

        if ((ee_event || mm_event) && vtx_cut && met_cut) {
	  if (svTagInfos && svTagInfos->nVertices() > 0) {
    	    const reco::Vertex &vertex = svTagInfos->secondaryVertex(0);
	    reco::TrackKinematics vertexKinematics(vertex);
	    math::XYZTLorentzVector vertexSum = vertexKinematics.weightedVectorSum();
	    sumVertexMass += vertexSum.M();
	  }
        }

      }

    }

  }

  if (Nb > 0) {
    sumVertexMassJet /= Nb;
    sumVertexMassTrk /= Nb;
    sumVertexMass /= Nb;
    //cout << "VTX mass JET = " << sumVertexMassJet << endl;
    //cout << "VTX mass TRK = " << sumVertexMassTrk << endl;
    //cout << "VTX mass NEW = " << sumVertexMass << endl;
  }

  // ++++++++ MET PLOTS

  if ((ee_event || mm_event) && Nj > 0 && vtx_cut) {
    w_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight);
    w_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight);
    if (isb) {
      b_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight);
      b_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight);
    }
    if (isc && !isb) {
      c_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight);
      c_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight);
    }
    if (Nb > 0) {
      scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
      //cout << vect_bjets[0].pt() << " " << vect_bjets[0].eta() <<"   SFb = " << scalFac_b << endl;
      w_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      w_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      if (isb) {
        b_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        b_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (isc && !isb) {
        c_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        c_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ HT PLOTS

  if ((ee_event || mm_event) && Nj > 0 && vtx_cut && met_cut) {
    w_Ht->Fill (Ht, MyWeight);
    if (isb) {
      b_Ht->Fill (Ht, MyWeight);
    }
    if (isc && !isb) {
      c_Ht->Fill (Ht, MyWeight);
    }
    if (Nb > 0) {
      scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
      w_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      if (Nj == 1) w_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      if (isb) {
        b_Ht_b->Fill (Ht, MyWeight*scalFac_b);
        if (Nj == 1) b_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      }
      if (isc && !isb) {
        c_Ht_b->Fill (Ht, MyWeight*scalFac_b);
        if (Nj == 1) c_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ DIELECTRON Z PLOTS

  if (iele1!=0 && Nj > 0 && vtx_cut && met_cut) {
    w_mass_ee_wide->Fill (diele_mass, MyWeight);
    if (isb) {
      b_mass_ee_wide->Fill (diele_mass, MyWeight);
    }
    if (isc && !isb) {
      c_mass_ee_wide->Fill (diele_mass, MyWeight);
    }
    if (Nb > 0) {
      scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
      w_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      if (isb) {
        b_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      }
      if (isc && !isb) {
        c_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      }
    }
  }

  if (ee_event && vtx_cut && met_cut) {
    w_numberOfZ->Fill (zee->size(), MyWeight);
    if (isb) {
      b_numberOfZ->Fill (zee->size(), MyWeight);
    }
    if (isc && !isb) {
      c_numberOfZ->Fill (zee->size(), MyWeight);
    }
  }

  if (ee_event && Nj > 0 && vtx_cut && met_cut) {
    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi());
    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
    h_mass_ee->Fill (diele_mass);
    w_mass_ee->Fill (diele_mass, MyWeight);
    w_pt_Z_ee->Fill (diele_pt, MyWeight);
    w_delta_ee->Fill (delta_phi_ee, MyWeight);
    if (isb) {
      b_mass_ee->Fill (diele_mass, MyWeight);
      b_pt_Z_ee->Fill (diele_pt, MyWeight);
    }
    if (isc && !isb) {
      c_mass_ee->Fill (diele_mass, MyWeight);
      c_pt_Z_ee->Fill (diele_pt, MyWeight);
    }
    if (Nb > 0) {
      scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
      w_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
      w_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
      double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi());
      if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
      w_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
      if (Nj == 1) {
        w_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        w_single_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
      }
      if (isb) {
        b_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
        b_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
        b_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
	if (Nj == 1) {
	  b_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
          b_single_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
	}
      }
      if (isc && !isb) {
        c_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
        c_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
        c_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        if (Nj == 1) {
	  c_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
          c_single_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
      }
    }
  }
}
  // ++++++++ DIMUON Z PLOTS

  if (imuon1!=0 && Nj > 0 && vtx_cut && met_cut) {
    w_mass_mm_wide->Fill (dimuon_mass, MyWeight);
    if (isb) {
      b_mass_mm_wide->Fill (dimuon_mass, MyWeight);
    }
    if (isc && !isb) {
      c_mass_mm_wide->Fill (dimuon_mass, MyWeight);
    }
    if (Nb > 0) {
      scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
      w_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      if (isb) {
        b_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      }
      if (isc && !isb) {
        c_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      }
    }
  }

  if (mm_event && vtx_cut && met_cut) {
    w_numberOfZ->Fill (zmm->size(), MyWeight);
    if (isb) {
      b_numberOfZ->Fill (zmm->size(), MyWeight);
    }
    if (isc && !isb) {
      c_numberOfZ->Fill (zmm->size(), MyWeight);
    }
  }

  if (mm_event && Nj > 0 && vtx_cut && met_cut) {
    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi());
    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
    h_mass_mm->Fill (dimuon_mass);
    w_mass_mm->Fill (dimuon_mass, MyWeight);
    w_pt_Z_mm->Fill (dimuon_pt, MyWeight);
    w_delta_mm->Fill (delta_phi_mm, MyWeight);
    if (isb) {
      b_mass_mm->Fill (dimuon_mass, MyWeight);
      b_pt_Z_mm->Fill (dimuon_pt, MyWeight);
    }
    if (isc && !isb) {
      c_mass_mm->Fill (dimuon_mass, MyWeight);
      c_pt_Z_mm->Fill (dimuon_pt, MyWeight);
    }
    if (Nb > 0) {
      scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
      w_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
      w_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
      double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi());
      if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
      w_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
      if (Nj == 1) {
        w_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        w_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
      }
      if (isb) {
        b_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
        b_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        b_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        if (Nj == 1) {
          b_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
          b_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
	}
      }
      if (isc && !isb) {
        c_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
        c_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        c_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        if (Nj == 1) {
          c_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
          c_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        }
      }
    }
  }

  // ++++++++ MISC PLOTS

  if ((ee_event || mm_event) && Nj > 0 && vtx_cut && met_cut) {
    h_pu_weights->Fill (MyWeight);
    h_tracks->Fill (tracks->size());
    w_tracks->Fill (tracks->size(), MyWeight);
    h_recoVTX->Fill (NVtx);
    w_recoVTX->Fill (NVtx, MyWeight);
  }

  // ++++++++  ELECTRONS PLOTS

  if (ee_event && Nj > 0 && vtx_cut && met_cut) {
#if ECALDRIVEN>1
    w_first_ele_pt->Fill (vect_ele[iele0].ecalDrivenMomentum().pt(), MyWeight);
    w_first_ele_eta->Fill (vect_ele[iele0].ecalDrivenMomentum().eta(), MyWeight);
    w_second_ele_pt->Fill (vect_ele[iele1].ecalDrivenMomentum().pt(), MyWeight);
    w_second_ele_eta->Fill (vect_ele[iele1].ecalDrivenMomentum().eta(), MyWeight);
#else
    w_first_ele_pt->Fill (vect_ele[iele0].pt(), MyWeight);
    w_first_ele_eta->Fill (vect_ele[iele0].eta(), MyWeight);
    w_second_ele_pt->Fill (vect_ele[iele1].pt(), MyWeight);
    w_second_ele_eta->Fill (vect_ele[iele1].eta(), MyWeight);
#endif
    if (isb) {
#if ECALDRIVEN>1
      b_first_ele_pt->Fill (vect_ele[iele0].ecalDrivenMomentum().pt(), MyWeight);
      b_first_ele_eta->Fill (vect_ele[iele0].ecalDrivenMomentum().eta(), MyWeight);
      b_second_ele_pt->Fill (vect_ele[iele1].ecalDrivenMomentum().pt(), MyWeight);
      b_second_ele_eta->Fill (vect_ele[iele1].ecalDrivenMomentum().eta(), MyWeight);
#else
      b_first_ele_pt->Fill (vect_ele[iele0].pt(), MyWeight);
      b_first_ele_eta->Fill (vect_ele[iele0].eta(), MyWeight);
      b_second_ele_pt->Fill (vect_ele[iele1].pt(), MyWeight);
      b_second_ele_eta->Fill (vect_ele[iele1].eta(), MyWeight);
#endif
    }
    if (isc && !isb) {
#if ECALDRIVEN>1
      c_first_ele_pt->Fill (vect_ele[iele0].ecalDrivenMomentum().pt(), MyWeight);
      c_first_ele_eta->Fill (vect_ele[iele0].ecalDrivenMomentum().eta(), MyWeight);
      c_second_ele_pt->Fill (vect_ele[iele1].ecalDrivenMomentum().pt(), MyWeight);
      c_second_ele_eta->Fill (vect_ele[iele1].ecalDrivenMomentum().eta(), MyWeight);
#else
      c_first_ele_pt->Fill (vect_ele[iele0].pt(), MyWeight);
      c_first_ele_eta->Fill (vect_ele[iele0].eta(), MyWeight);
      c_second_ele_pt->Fill (vect_ele[iele1].pt(), MyWeight);
      c_second_ele_eta->Fill (vect_ele[iele1].eta(), MyWeight);
#endif
    }
  }

  if (ee_event && Nj > 0 && Nb > 0 && vtx_cut && met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
    w_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
#if ECALDRIVEN>1
    w_first_ele_pt_b->Fill (vect_ele[iele0].ecalDrivenMomentum().pt(), MyWeight*scalFac_b);
#else
    w_first_ele_pt_b->Fill (vect_ele[iele0].pt(), MyWeight*scalFac_b);
#endif
    if (isb) {
      b_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
    }
    if (isc && !isb) {
      c_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
    }
  }

  // ++++++++ MUONS PLOTS

  if (mm_event && Nj > 0 && vtx_cut && met_cut) {
    w_first_muon_pt->Fill (vect_muon[imuon0].pt(), MyWeight);
    w_first_muon_eta->Fill (vect_muon[imuon0].eta(), MyWeight);
    w_second_muon_pt->Fill (vect_muon[imuon1].pt(), MyWeight);
    w_second_muon_eta->Fill (vect_muon[imuon1].eta(), MyWeight);
    if (isb) {
      b_first_muon_pt->Fill (vect_muon[imuon0].pt(), MyWeight);
      b_first_muon_eta->Fill (vect_muon[imuon0].eta(), MyWeight);
      b_second_muon_pt->Fill (vect_muon[imuon1].pt(), MyWeight);
      b_second_muon_eta->Fill (vect_muon[imuon1].eta(), MyWeight);
    }
    if (isc && !isb) {
      c_first_muon_pt->Fill (vect_muon[imuon0].pt(), MyWeight);
      c_first_muon_eta->Fill (vect_muon[imuon0].eta(), MyWeight);
      c_second_muon_pt->Fill (vect_muon[imuon1].pt(), MyWeight);
      c_second_muon_eta->Fill (vect_muon[imuon1].eta(), MyWeight);
    }
  }

  if (mm_event && Nj > 0 && Nb > 0 && vtx_cut && met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
    w_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
    w_first_muon_pt_b ->Fill (vect_muon[imuon0].pt(), MyWeight*scalFac_b);
    if (isb) {
      b_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
    }
    if (isc && !isb) {
      c_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
    }
  }

  // ++++++++ SVTX MASS PLOTS

  if ((ee_event || mm_event) && Nj > 0 && Nb > 0 && vtx_cut && met_cut) {

    w_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
    w_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
    w_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    if (isb) {
      b_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      b_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      b_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
    if (isc && !isb) {
      c_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      c_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      c_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
  }

  // ++++++++ CSV PLOTS

  if ((ee_event || mm_event) && Nj > 0 && Nb > 0 && vtx_cut && met_cut && sumVertexMass > 0.0 ) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
    w_secondvtx_N_mass->Fill (vect_bjets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
    if (isb) {
      b_secondvtx_N_mass->Fill (vect_bjets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
    }
    if (isc && !isb) {
      c_secondvtx_N_mass->Fill (vect_bjets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
    }
  }

  // ++++++++ BJP/JBP PLOTS

  if ((ee_event || mm_event) && Nj > 0 && vtx_cut && met_cut && Nb > 0) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
    w_BJP->Fill (vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
    w_JBP->Fill (vect_bjets[0].bDiscriminator("jetProbabilityBJetTags"), MyWeight*scalFac_b);
    if (sumVertexMass > 0.0) {
      w_BJP_mass->Fill (vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
      w_JBP_mass->Fill (vect_bjets[0].bDiscriminator("jetProbabilityBJetTags"), MyWeight*scalFac_b);
    }
    if (isb) {
      b_BJP->Fill (vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
      b_JBP->Fill (vect_bjets[0].bDiscriminator("jetProbabilityBJetTags"), MyWeight*scalFac_b);
      if (sumVertexMass > 0.0) {
        b_BJP_mass->Fill (vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
        b_JBP_mass->Fill (vect_bjets[0].bDiscriminator("jetProbabilityBJetTags"), MyWeight*scalFac_b);
      }
    }
    if (isc && !isb) {
      c_BJP->Fill (vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
      c_JBP->Fill (vect_bjets[0].bDiscriminator("jetProbabilityBJetTags"), MyWeight*scalFac_b);
      if (sumVertexMass > 0.0) {
        c_BJP_mass->Fill (vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
        c_JBP_mass->Fill (vect_bjets[0].bDiscriminator("jetProbabilityBJetTags"), MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ JETS PLOTS

  if ((ee_event || mm_event) && Nj > 0 && vtx_cut && met_cut) {
    h_jetmultiplicity->Fill (Nj);
    w_jetmultiplicity->Fill (Nj, MyWeight);
    w_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight);
    if (isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_jetmultiplicity->Fill (Nj, MyWeight);
      b_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight);
      b_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight);
    }
    if (isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_jetmultiplicity->Fill (Nj, MyWeight);
      c_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight);
      c_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight);
    }
  }

  if ((ee_event || mm_event) && Nj > 1 && vtx_cut && met_cut) {
    w_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight);
    w_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight);
    if (isMC && fabs(vect_jets[1].partonFlavour()) == 5) {
      b_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight);
      b_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight);
    }
    if (isMC && fabs(vect_jets[1].partonFlavour()) == 4) {
      c_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight);
      c_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight);
    }
  }

  if ((ee_event || mm_event) && Nj > 2 && vtx_cut && met_cut) {
    w_third_jet_pt->Fill (vect_jets[2].pt(), MyWeight);
    w_third_jet_eta->Fill (vect_jets[2].eta(), MyWeight);
    if (isMC && fabs(vect_jets[2].partonFlavour()) == 5) {
      b_third_jet_pt->Fill (vect_jets[2].pt(), MyWeight);
      b_third_jet_eta->Fill (vect_jets[2].eta(), MyWeight);
    }
    if (isMC && fabs(vect_jets[2].partonFlavour()) == 4) {
      c_third_jet_pt->Fill (vect_jets[2].pt(), MyWeight);
      c_third_jet_eta->Fill (vect_jets[2].eta(), MyWeight);
    }
  }

  // ++++++++ B JETS PLOTS

  if ((ee_event || mm_event) && Nj > 0 && Nb > 0 && vtx_cut && met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
    w_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    if (isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      b_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      b_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      b_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      c_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      c_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      c_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
  }

  if ((ee_event || mm_event) && Nj > 1 && Nb > 1 && vtx_cut && met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta())*btagSF(isMC, vect_bjets[1].partonFlavour(), vect_bjets[1].pt(), vect_bjets[1].eta());
    w_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
    w_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    w_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
    w_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
    if (isMC && fabs(vect_jets[1].partonFlavour()) == 5) {
      b_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      b_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_bjets[1].partonFlavour()) == 5) {
      b_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
      b_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_jets[1].partonFlavour()) == 4) {
      c_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      c_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_bjets[1].partonFlavour()) == 4) {
      c_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
      c_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
    }
  }

  if ((ee_event || mm_event) && Nj > 2 && Nb > 2 && vtx_cut && met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta())*btagSF(isMC, vect_bjets[1].partonFlavour(), vect_bjets[1].pt(), vect_bjets[1].eta())*btagSF(isMC, vect_bjets[2].partonFlavour(), vect_bjets[2].pt(), vect_bjets[2].eta());
    w_third_jet_pt_b->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
    w_third_jet_eta_b->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
    w_third_bjet_pt->Fill (vect_bjets[2].pt(), MyWeight*scalFac_b);
    w_third_bjet_eta->Fill (vect_bjets[2].eta(), MyWeight*scalFac_b);
    if (isMC && fabs(vect_jets[2].partonFlavour()) == 5) {
      b_third_jet_pt_b->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
      b_third_jet_eta_b->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_bjets[2].partonFlavour()) == 5) {
      b_third_bjet_pt->Fill (vect_bjets[2].pt(), MyWeight*scalFac_b);
      b_third_bjet_eta->Fill (vect_bjets[2].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_jets[2].partonFlavour()) == 4) {
      c_third_jet_pt_b->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
      c_third_jet_eta_b->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_bjets[2].partonFlavour()) == 4) {
      c_third_bjet_pt->Fill (vect_bjets[2].pt(), MyWeight*scalFac_b);
      c_third_bjet_eta->Fill (vect_bjets[2].eta(), MyWeight*scalFac_b);
    }
  }

  // ++++++++ SINGLE BJET

  if ((ee_event || mm_event) && Nj == 1 && Nb > 0 && vtx_cut && met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
    w_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
    w_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
    w_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);

    if (isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      b_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
    if (isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      c_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
  }

  // ++++++++ EXTRA PLOTS

  int Nf = 0;
  int Nbk = 0;

  double Afb = 0;

  if ((ee_event || mm_event) && Nj > 0 && Nb > 0 && vtx_cut && met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets[0].partonFlavour(), vect_bjets[0].pt(), vect_bjets[0].eta());
    if (fabs (vect_bjets[0].eta()) > 0) Nf++;
    if (fabs (vect_bjets[0].eta()) < 0) Nbk++;
    if ((Nf+Nbk) != 0) Afb = (Nf - Nbk) / (Nf + Nbk);
    w_Afb->Fill (Afb, MyWeight*scalFac_b);
  }

  if (ee_event && Nj > 0 && vtx_cut && met_cut) {
    h_scaleFactor_first_ele->Fill (scalFac_first_e, MyWeight / (scalFac_first_e * scalFac_second_e));
    h_scaleFactor_second_ele->Fill (scalFac_second_e, MyWeight / (scalFac_first_e * scalFac_second_e));
    if (isb) {
      b_scaleFactor_first_ele->Fill (scalFac_first_e, MyWeight / (scalFac_first_e * scalFac_second_e));
      b_scaleFactor_second_ele->Fill (scalFac_second_e, MyWeight / (scalFac_first_e * scalFac_second_e));
    }
  }
  if (mm_event && Nj > 0 && vtx_cut && met_cut) {
    h_scaleFactor_first_muon->Fill (scalFac_first_m, MyWeight / (scalFac_first_m * scalFac_second_m));
    h_scaleFactor_second_muon->Fill (scalFac_second_m, MyWeight / (scalFac_first_m * scalFac_second_m));
    if (isb) {
      b_scaleFactor_first_muon->Fill (scalFac_first_m, MyWeight / (scalFac_first_m * scalFac_second_m));
      b_scaleFactor_second_muon->Fill (scalFac_second_m, MyWeight / (scalFac_first_m * scalFac_second_m));
    }
  }

  // ++++++++ OUTPUT COLLECTIONS

  if ((ee_event || mm_event) && Nj > 0 && vtx_cut && met_cut) {
    myEventWeight->push_back(MyWeight);
  }

  if (ee_event && Nj > 0 && vtx_cut && met_cut) {
    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[iele0].px(),vect_ele[iele0].py(),vect_ele[iele0].pz(),vect_ele[iele0].energy()));
    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[iele1].px(),vect_ele[iele1].py(),vect_ele[iele1].pz(),vect_ele[iele1].energy()));
    myPtZ->push_back(diele_pt);
  }

  if (mm_event && Nj > 0 && vtx_cut && met_cut) {
    myMuons->push_back(math::XYZTLorentzVector(vect_muon[imuon0].px(),vect_muon[imuon0].py(),vect_muon[imuon0].pz(),vect_muon[imuon0].energy()));
    myMuons->push_back(math::XYZTLorentzVector(vect_muon[imuon1].px(),vect_muon[imuon1].py(),vect_muon[imuon1].pz(),vect_muon[imuon1].energy()));
    myPtZ->push_back(dimuon_pt);
  }

  if ((ee_event || mm_event) && Nj > 0 && vtx_cut && met_cut) {
    for (unsigned int i=0; i<vect_jets.size(); ++i) {
      myJets->push_back(math::XYZTLorentzVector(vect_jets[i].px(),vect_jets[i].py(),vect_jets[i].pz(),vect_jets[i].energy()));
    }
    for (unsigned int i=0; i<vect_bjets.size(); ++i) {
      scalFac_b = btagSF(isMC, vect_bjets[i].partonFlavour(), vect_bjets[i].pt(), vect_bjets[i].eta());
      myBjetsWeights->push_back(scalFac_b);
      myBjets->push_back(math::XYZTLorentzVector(vect_bjets[i].px(),vect_bjets[i].py(),vect_bjets[i].pz(),vect_bjets[i].energy()));
    }
  }

  if ((ee_event || mm_event) && Nj > 0 && vtx_cut && met_cut) {
    myHt->push_back(Ht);
  }

  if (ee_event && Nj > 0 && vtx_cut && met_cut) {
    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi());
    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
    myDeltaPhi->push_back(delta_phi_ee);
    if (Nb > 0) {
      double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi());
      if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
      myBdeltaPhi->push_back(delta_phi_ee_b);
    }
  }

  if (mm_event && Nj > 0 && vtx_cut && met_cut) {
    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi());
    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
    myDeltaPhi->push_back(delta_phi_mm);
    if (Nb > 0) {
      double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi());
      if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
      myBdeltaPhi->push_back(delta_phi_mm_b);
    }
  }

  iEvent.put( myEventWeight, "myEventWeight" );

  iEvent.put( myElectrons, "myElectrons" );
  iEvent.put( myMuons, "myMuons" );

  iEvent.put( myPtZ, "myPtZ" );

  iEvent.put( myJets, "myJets" );
  iEvent.put( myDeltaPhi, "myDeltaPhi" );

  iEvent.put( myHt, "myHt" );

  iEvent.put( myBjetsWeights, "myBjetsWeights" );

  iEvent.put( myBjets, "myBjets" );
  iEvent.put( myBdeltaPhi, "myBdeltaPhi" );

}

// ------------ method called once each job just before starting event loop ------------
void ZbAnalyzer::beginJob () {
  jecUncDT_ = new JetCorrectionUncertainty(path_ + "/" + "Fall12_V7_DATA_Uncertainty_AK5PFchs.txt");
  jecUncMC_ = new JetCorrectionUncertainty(path_ + "/" + "Fall12_V7_MC_Uncertainty_AK5PFchs.txt");
  LumiWeights_ = edm::LumiReWeighting(path_ + "/" + "pileup_" + pileup_ + ".root", path_ + "/" + "pileup_2012.root", "pileup", "pileup");

  ElSF_  = new table(path_ + "/" + "ele_eff.txt");
  ElSF2_ = new table(path_ + "/" + "ele_eff2.txt");
  MuSF_  = new table(path_ + "/" + "muon_eff.txt");
  BtSF_  = new table(path_ + "/" + "btag_eff.txt");   //btagging scale factors SFb = SFc
  LtSF_  = new table(path_ + "/" + "light_eff.txt");  //light flavour scale factors

  gRandom_ = new TRandom3();
}

// ------------ method called once each job just after ending the event loop ------------
void ZbAnalyzer::endJob () {
  delete jecUncDT_;
  delete jecUncMC_;

  delete ElSF_;
  delete ElSF2_;
  delete MuSF_;
  delete BtSF_;
  delete LtSF_;

  delete gRandom_;
}

// ------------ method called when starting to processes a run ------------
void ZbAnalyzer::beginRun (edm::Run const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a run ------------
void ZbAnalyzer::endRun (edm::Run const &, edm::EventSetup const &) {
}

// ------------ method called when starting to processes a luminosity block ------------
void ZbAnalyzer::beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a luminosity block ------------
void ZbAnalyzer::endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// define this as a plug-in
DEFINE_FWK_MODULE (ZbAnalyzer);
