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
// $Id: ZbAnalyzer.cc,v 1.55 2013/05/15 07:46:08 dellaric Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/ESHandle.h>

// system include files
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stddef.h>

// user include files
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
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "PhysicsTools/CandUtils/interface/CandCombiner.h"
#include "CommonTools/Utils/interface/MassRangeSelector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"

#include "table.h"
table ElSF ("/gpfs/cms/users/candelis/work/ZbSkim/test/ele_eff.txt");
table MuSF ("/gpfs/cms/users/candelis/work/ZbSkim/test/muon_eff.txt");
table BtSF ("/gpfs/cms/users/candelis/work/ZbSkim/test/btag_eff.txt");

class TTree;

//
// class declaration
//

class ZbAnalyzer:public  edm::EDAnalyzer {

public:

  explicit ZbAnalyzer (const edm::ParameterSet &);
  ~ZbAnalyzer ();

private:

  virtual void beginJob ();
  virtual void analyze (const edm::Event &, const edm::EventSetup &);
  virtual void endJob ();

  virtual void beginRun (edm::Run const &, edm::EventSetup const &);
  virtual void endRun (edm::Run const &, edm::EventSetup const &);
  virtual void beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);
  virtual void endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);

  std::string pileup_;
  std::string lepton_;
  double par_;
  JetCorrectionUncertainty *jecUncDT_;
  JetCorrectionUncertainty *jecUncMC_;
  edm::LumiReWeighting LumiWeights_;

  // ----------member data ---------------------------


  /****************** LEGEND *************************

   w_plotname_b => at least a b in the event, Nb > 0
   b_plotname   => b quark fraction, isb
   c_plotname   => c quark fraction, isc

   ***************************************************/

  int       Nj, Nb;
  double    jet_pt;
  double    jet_eta;
  double    jet_phi;
  double    bjet_pt;
  double    bjet_eta;
  double    bjet_phi;
  double    ele_pt;
  double    ele_eta;
  double    muon_pt;
  double    muon_eta;
  double    diele_mass;
  double    diele_phi;
  double    diele_pt;
  double    dimuon_mass;
  double    dimuon_phi;
  double    dimuon_pt;
  double    delta_phi_ee;
  double    delta_phi_mm;
  double    b_leading_pt;
  double    b_leading_eta;

  double    Ht, Ht_b;
  double    discrCSV;
  double    discrBJP, discrJBP; 
  double    MyWeight;
  double    sFac;
  double    sFacErr;
  double    scalFac_first_e;
  double    scalFac_second_e;
  double    scalFac_first_m;
  double    scalFac_second_m;
  double    scalFac_b;
  double    Nf, Nbk;
  double    Afb;

  TTree*    treeZb_;

  TH1F*     h_jetmultiplicity;
  TH1F*     h_jet_pt;
  TH1F*     h_ele_pt;
  TH1F*     h_muon_pt;

  TH1F*     h_pu_weights;
  TH1F*     h_tracks;
  TH1F*     w_tracks;
  TH1F*     h_recoVTX;
  TH1F*     w_recoVTX;

  TH1F*     w_jetmultiplicity;
  TH1F*     b_jetmultiplicity;
  TH1F*     c_jetmultiplicity;
  TH1F*     w_first_jet_pt;      // leading jet of any type
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

  TH1F*     w_first_jet_pt_b;    // leading jet with at least one b jet in the event
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
  TH1F*     w_first_bjet_pt;     // leading b jet
  TH1F*     b_first_bjet_pt;
  TH1F*     c_first_bjet_pt;
  TH1F*     w_first_bjet_eta;
  TH1F*     b_first_bjet_eta;
  TH1F*     c_first_bjet_eta;
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
  TH1F*     b_first_ele_pt;
  TH1F*     c_first_ele_pt;
  TH1F*     w_second_ele_pt;
  TH1F*     b_second_ele_pt;
  TH1F*     c_second_ele_pt;
  TH1F*     w_first_muon_pt;
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

  TH1F*     w_mass_ee_b;  // at least one b jet in the event
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
  TH1F*     w_delta_ee_b;
  TH1F*     b_delta_ee_b;
  TH1F*     c_delta_ee_b;
  TH1F*     w_delta_mm_b;
  TH1F*     b_delta_mm_b;
  TH1F*     c_delta_mm_b;

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

  TH1F*     w_Ht;
  TH1F*     b_Ht;
  TH1F*     c_Ht;
  TH1F*     w_Ht_b; // at least one b jet in the event
  TH1F*     b_Ht_b;
  TH1F*     c_Ht_b;

  TH1F*     w_MET;
  TH1F*     b_MET;
  TH1F*     c_MET;
  TH1F*     w_MET_sign;

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

enum Flavour {
  ALL_JETS = 0,
  UDSG_JETS,
  C_JETS,
  B_JETS,
  NONID_JETS,
  N_JET_TYPES
};

struct Plots {
  TH1 * nVertices;
  TH1 * deltaR, * mass, * dist, * distErr, * distSig;
  TH1 * nTracks, * chi2;
} plots_[N_JET_TYPES];


//
// static data member definitions
//


//
// constructors and destructor
//
ZbAnalyzer::ZbAnalyzer (const edm::ParameterSet & iConfig) {

  pileup_ = iConfig.getUntrackedParameter < std::string > ("pileup", "S7");
  lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
  par_ =    iConfig.getUntrackedParameter <double> ("JEC",  0);

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  h_jetmultiplicity =   fs->make < TH1F > ("h_jetmultiplicity", "h_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  h_jet_pt =            fs->make < TH1F > ("h_jet_pt",          "h_jet_pt;P_t [GeV]", 20, 30, 530);
  h_ele_pt =            fs->make < TH1F > ("h_ele_pt",          "h_ele_pt;P_t [GeV]", 20, 30, 530);
  h_muon_pt =           fs->make < TH1F > ("h_muon_pt",         "h_muon_pt;P_t [GeV]", 100, 0, 250);

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
  b_first_ele_pt =      fs->make < TH1F > ("b_first_ele_pt",    "b_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  c_first_ele_pt =      fs->make < TH1F > ("c_first_ele_pt",    "c_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_second_ele_pt =     fs->make < TH1F > ("w_second_ele_pt",   "w_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  b_second_ele_pt =     fs->make < TH1F > ("b_second_ele_pt",   "b_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  c_second_ele_pt =     fs->make < TH1F > ("c_second_ele_pt",   "c_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_first_muon_pt =     fs->make < TH1F > ("w_first_muon_pt",   "w_first_muon_pt;P_t [GeV]", 50, 0., 450.);
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

  h_mass_mm =           fs->make < TH1F > ("h_mass_mm",         "h_mass_mm;Mass [GeV]", 80, 71, 111);
  w_mass_ee = 	        fs->make < TH1F > ("w_mass_ee",         "w_mass_ee;Mass [GeV]", 80, 71, 111);
  b_mass_ee =           fs->make < TH1F > ("b_mass_ee",         "b_mass_ee;Mass [GeV]", 80, 71, 111);
  c_mass_ee =           fs->make < TH1F > ("c_mass_ee",         "c_mass_ee;Mass [GeV]", 80, 71, 111);
  h_mass_ee =           fs->make < TH1F > ("h_mass_ee",         "h_mass_ee;Mass [GeV]", 80, 71, 111);
  w_mass_mm = 	        fs->make < TH1F > ("w_mass_mm",         "w_mass_mm;Mass [GeV]", 80, 71, 111);
  b_mass_mm =           fs->make < TH1F > ("b_mass_mm",         "b_mass_mm;Mass [GeV]", 80, 71, 111);
  c_mass_mm =           fs->make < TH1F > ("c_mass_mm",         "c_mass_mm;Mass [GeV]", 80, 71, 111);
  w_pt_Z_ee =           fs->make < TH1F > ("w_pt_Z_ee",         "w_pt_Z_ee;P_t [GeV]", 70, 0., 700.);
  b_pt_Z_ee =           fs->make < TH1F > ("b_pt_Z_ee",         "b_pt_Z_ee;P_t [GeV]", 70, 0., 700.);
  c_pt_Z_ee =           fs->make < TH1F > ("c_pt_Z_ee",         "c_pt_Z_ee;P_t [GeV]", 70, 0., 700.);
  w_pt_Z_mm =           fs->make < TH1F > ("w_pt_Z_mm",         "w_pt_Z_mm;P_t [GeV]", 70, 0., 700.);
  b_pt_Z_mm =           fs->make < TH1F > ("b_pt_Z_mm",         "b_pt_Z_mm;P_t [GeV]", 70, 0., 700.);
  c_pt_Z_mm =           fs->make < TH1F > ("c_pt_Z_mm",         "c_pt_Z_mm;P_t [GeV]", 70, 0., 700.);

  w_mass_ee_b =         fs->make < TH1F > ("w_mass_ee_b",       "w_mass_mm_b;Mass [GeV]", 80, 71, 111);
  b_mass_ee_b =         fs->make < TH1F > ("b_mass_ee_b",       "b_mass_mm_b;Mass [GeV]", 80, 71, 111);
  c_mass_ee_b =         fs->make < TH1F > ("c_mass_ee_b",       "c_mass_mm_b;Mass [GeV]", 80, 71, 111);
  w_mass_mm_b =         fs->make < TH1F > ("w_mass_mm_b",       "w_mass_mm_b;Mass [GeV]", 80, 71, 111);
  b_mass_mm_b =         fs->make < TH1F > ("b_mass_mm_b",       "b_mass_mm_b;Mass [GeV]", 80, 71, 111);
  c_mass_mm_b =         fs->make < TH1F > ("c_mass_mm_b",       "c_mass_mm_b;Mass [GeV]", 80, 71, 111);
  w_pt_Z_ee_b =         fs->make < TH1F > ("w_pt_Z_ee_b",       "w_pt_Z_ee_b;P_t [GeV]", 70, 0., 700.);
  b_pt_Z_ee_b =         fs->make < TH1F > ("b_pt_Z_ee_b",       "b_pt_Z_ee_b;P_t [GeV]", 70, 0., 700.);
  c_pt_Z_ee_b =         fs->make < TH1F > ("c_pt_Z_ee_b",       "c_pt_Z_ee_b;P_t [GeV]", 70, 0., 700.);
  w_pt_Z_mm_b =         fs->make < TH1F > ("w_pt_Z_mm_b",       "w_pt_Z_mm_b;P_t [GeV]", 70, 0., 700.);
  b_pt_Z_mm_b =         fs->make < TH1F > ("b_pt_Z_mm_b",       "b_pt_Z_mm_b;P_t [GeV]", 70, 0., 700.);
  c_pt_Z_mm_b =         fs->make < TH1F > ("c_pt_Z_mm_b",       "c_pt_Z_mm_b;P_t [GeV]", 70, 0., 700.);
  w_delta_ee_b =        fs->make < TH1F > ("w_delta_phi_ee_b",  "w_delta_phi_ee_b", 12, 0, TMath::Pi ());
  b_delta_ee_b =        fs->make < TH1F > ("b_delta_phi_ee_b",  "b_delta_phi_ee_b", 12, 0, TMath::Pi ());
  c_delta_ee_b =        fs->make < TH1F > ("c_delta_phi_ee_b",  "c_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_delta_mm_b =        fs->make < TH1F > ("w_delta_phi_mm_b",  "w_delta_phi_mm_b", 12, 0, TMath::Pi ());
  b_delta_mm_b =        fs->make < TH1F > ("b_delta_phi_mm_b",  "b_delta_phi_mm_b", 12, 0, TMath::Pi ());
  c_delta_mm_b =        fs->make < TH1F > ("c_delta_phi_mm_b",  "c_delta_phi_mm_b", 12, 0, TMath::Pi ());

  h_secondvtx_N =       fs->make < TH1F > ("h_secondvtx_N",     "h_secondvtx_N", 50, 0, 1);
  w_secondvtx_N =       fs->make < TH1F > ("w_secondvtx_N",     "w_secondvtx_N", 50, 0, 1);
  w_secondvtx_N_zoom =  fs->make < TH1F > ("w_secondvtx_N_zoom",  "w_secondvtx_N_zoom", 25, 0.5, 1);
  w_secondvtx_N_mass =  fs->make < TH1F > ("w_secondvtx_N_mass",  "w_secondvtx_N_mass", 50, 0, 1);
  
  b_secondvtx_N =         fs->make < TH1F > ("b_secondvtx_N",       "b_secondvtx_N", 50, 0, 1);
  b_secondvtx_N_mass =    fs->make < TH1F > ("b_secondvtx_N_mass",  "b_secondvtx_N_mass", 50, 0, 1);
  b_secondvtx_N_zoom =    fs->make < TH1F > ("b_secondvtx_N_zoom",  "b_secondvtx_N_zoom", 25, 0.5, 1);
  
  c_secondvtx_N =         fs->make < TH1F > ("c_secondvtx_N",       "c_secondvtx_N", 50, 0, 1);
  c_secondvtx_N_mass =    fs->make < TH1F > ("c_secondvtx_N_mass",  "c_secondvtx_N_mass", 50, 0, 1);
  c_secondvtx_N_zoom =    fs->make < TH1F > ("c_secondvtx_N_zoom",  "c_secondvtx_N_zoom", 25, 0.5, 1);
  

  w_SVTX_mass_jet =     fs->make < TH1F > ("w_SVTX_mass_jet",   "w_SVTX_mass_jet;Mass [GeV]", 70, 0, 7);
  b_SVTX_mass_jet =     fs->make < TH1F > ("b_SVTX_mass_jet",   "b_SVTX_mass_jet;Mass [GeV]", 70, 0, 7);
  c_SVTX_mass_jet =     fs->make < TH1F > ("c_SVTX_mass_jet",   "c_SVTX_mass_jet;Mass [GeV]", 70, 0, 7);
  w_SVTX_mass_trk =     fs->make < TH1F > ("w_SVTX_mass_trk",   "w_SVTX_mass_trk;Mass [GeV]", 160, 0, 80);
  b_SVTX_mass_trk =     fs->make < TH1F > ("b_SVTX_mass_trk",   "b_SVTX_mass_trk;Mass [GeV]", 160, 0, 80);
  c_SVTX_mass_trk =     fs->make < TH1F > ("c_SVTX_mass_trk",   "c_SVTX_mass_trk;Mass [GeV]", 160, 0, 80);
  w_SVTX_mass     =     fs->make < TH1F > ("w_SVTX_mass",       "w_SVTX_mass;Mass [GeV]", 70, 0, 7);
  b_SVTX_mass     =     fs->make < TH1F > ("b_SVTX_mass",       "b_SVTX_mass;Mass [GeV]", 70, 0, 7);
  c_SVTX_mass     =     fs->make < TH1F > ("c_SVTX_mass",       "c_SVTX_mass;Mass [GeV]", 70, 0, 7);

  w_BJP  =     fs->make < TH1F > ("w_BJP",   "w_BJP", 80, 0, 10);    
  w_JBP  =     fs->make < TH1F > ("w_JBP",   "w_JBP", 50, 0, 3);
  b_BJP  =     fs->make < TH1F > ("b_BJP",   "b_BJP", 80, 0, 10);
  b_JBP  =     fs->make < TH1F > ("b_JBP",   "b_JBP", 50, 0, 3);
  c_BJP  =     fs->make < TH1F > ("c_BJP",   "c_BJP", 80, 0, 10);
  c_JBP  =     fs->make < TH1F > ("c_JBP",   "c_JBP", 50, 0, 3);

  w_Ht =                fs->make < TH1F > ("w_Ht",              "w_Ht [GeV]", 50, 30., 1500.);
  b_Ht =                fs->make < TH1F > ("b_Ht",              "b_Ht [GeV]", 50, 30., 1500.);
  c_Ht =                fs->make < TH1F > ("c_Ht",              "c_Ht [GeV]", 50, 30., 1500.);
  w_Ht_b =              fs->make < TH1F > ("w_Ht_b",            "w_Ht [GeV]", 50, 30., 1500.);
  b_Ht_b =              fs->make < TH1F > ("b_Ht_b",            "b_Ht [GeV]", 50, 30., 1500.);
  c_Ht_b =              fs->make < TH1F > ("c_Ht_b",            "c_Ht [GeV]", 50, 30., 1500.);

  w_MET =               fs->make < TH1F > ("w_MET",             "w_MET;MET [GeV]", 50, 0, 250);
  b_MET =               fs->make < TH1F > ("b_MET",             "b_MET;MET [GeV]", 50, 0, 250);
  c_MET =               fs->make < TH1F > ("c_MET",             "c_MET;MET [GeV]", 50, 0, 250);
  w_MET_sign = 	        fs->make < TH1F > ("w_MET_sign",        "w_MET_sign", 50, 0, 50);

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

  treeZb_ = fs->make < TTree > ("ZbTree", "ZbTree");
  treeZb_->Branch ("Nj", &Nj);
  treeZb_->Branch ("jet_pt", &jet_pt);
  treeZb_->Branch ("muon_pt", &muon_pt);
  treeZb_->Branch ("dimuon_mass", &dimuon_mass);
  treeZb_->Branch ("diele_mass", &diele_mass);
  treeZb_->Branch ("weights_pu", &MyWeight);

  for (unsigned int i = 0; i < N_JET_TYPES; i++)    {
      Plots & plots = plots_[i];
      const char *flavour, *name;

      switch ((Flavour) i) {
	case ALL_JETS:
	  flavour = "all jets";
	  name = "all";
	  break;
	case UDSG_JETS:
	  flavour = "light flavour jets";
	  name = "udsg";
	  break;
	case C_JETS:
	  flavour = "charm jets";
	  name = "c";
	  break;
	case B_JETS:
	  flavour = "bottom jets";
	  name = "b";
	  break;
	default:
	  flavour = "unidentified jets";
	  name = "ni";
	  break;
      }
      plots.dist = fs->make < TH1F > (Form ("dist_%s", name), Form ("Transverse distance between PV and SV in %s", flavour), 100, 0, 2);

    }
}

ZbAnalyzer::~ZbAnalyzer () {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event ------------
void ZbAnalyzer::analyze (const edm::Event & iEvent, const edm::EventSetup & iSetup) {

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

  bool ee_event = false;
  bool mm_event = false;

  int Ntracks = 0;
  Nj = 0;
  Nbk = 0;
  jet_pt = 0;
  jet_eta = 0;
  jet_phi = 0;
  bjet_pt = 0;
  bjet_eta = 0;
  bjet_phi = 0;
  ele_pt = 0;
  ele_eta = 0;
  muon_pt = 0;
  muon_eta = 0;
  diele_mass = 0;
  dimuon_mass = 0;
  discrCSV = 0;
  b_leading_pt = 0;
  b_leading_eta = 0;
  delta_phi_ee = 0;
  delta_phi_mm = 0;
  Nf = 0;
  Nb = 0;
  Afb = 0;
  scalFac_first_e = 1;
  scalFac_second_e = 1;
  scalFac_first_m = 1;
  scalFac_second_m = 1;
  scalFac_b = 1;

  Ht = 0;
  Ht_b = 0;

  // +++++++++ ELECTRONS

  vector < double > vect_ele_pt;
  vector < double > vect_ele_eta;

  for (pat::ElectronCollection::const_iterator ele = electrons->begin (); ele != electrons->end (); ++ele) {
    ele_pt = ele->pt ();
    ele_eta = ele->eta ();

//    ele_pt = ele->ecalDrivenMomentum().pt ();
//    ele_eta = ele->ecalDrivenMomentum().eta ();

    vect_ele_pt.push_back (ele_pt);
    vect_ele_eta.push_back (ele_eta);
  }

  if (zee->size () != 0) ee_event = true;

  // +++++++++ MUONS

  vector < double > vect_muon_pt;
  vector < double > vect_muon_eta;

  for (pat::MuonCollection::const_iterator muon = muons->begin (); muon != muons->end (); ++muon) {
    muon_pt = muon->pt ();
    muon_eta = muon->eta ();
    vect_muon_pt.push_back (muon_pt);
    vect_muon_eta.push_back (muon_eta);
  }

  if (zmm->size () != 0) mm_event = true;

  if (lepton_ == "electron" && !ee_event) return;
  if (lepton_ == "muon" && !mm_event)     return;

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

    if (ee_event) {
      scalFac_first_e  =  ElSF.Val (vect_ele_pt[0], vect_ele_eta[0]);
      scalFac_second_e =  ElSF.Val (vect_ele_pt[1], vect_ele_eta[1]);
    }
    if (mm_event) {
      scalFac_first_m  = MuSF.Val (vect_muon_pt[0], vect_muon_eta[0]);
      scalFac_second_m = MuSF.Val (vect_muon_pt[1], vect_muon_eta[1]);
      //cout<<vect_muon_pt[0]<<vect_muon_eta[0]<< " mu  SF =" << scalFac_first_m <<endl;

    }

  }

  if (ee_event) MyWeight = MyWeight * scalFac_first_e * scalFac_second_e;
  if (mm_event) MyWeight = MyWeight * scalFac_first_m * scalFac_second_m;

  // ++++++++ VERTICES

  edm::Handle < vector < reco::Vertex > > vertices_h;
  iEvent.getByLabel (edm::InputTag ("goodOfflinePrimaryVertices"), vertices_h);

  // require in the event that there is at least one reconstructed vertex
  if (vertices_h->size () <= 0)    return;

  // pick the first (i.e. highest sum pt) vertex
  const reco::Vertex * theVertex = &(vertices_h->front ());

  // require that the vertex meets certain criteria
  if (theVertex->ndof () < 5)    return;
  if (fabs (theVertex->z ()) > 24.0)    return;
  if (fabs (theVertex->position ().rho ()) > 2.0)    return;

  // now, count vertices
  int NVtx = 0;
  for (vector < reco::Vertex >::const_iterator itv = vertices_h->begin (); itv != vertices_h->end (); ++itv) {
      // require that the vertex meets certain criteria
      if (itv->ndof () < 5)	continue;
      if (fabs (itv->z ()) > 50.0)	continue;
      if (fabs (itv->position ().rho ()) > 2.0)	continue;
      ++NVtx;
  }

  // ++++++++ JETS

  vector < double > vect_jets_pt;
  vector < double > vect_jets_phi;
  vector < double > vect_jets_eta;
  vector < bool >   vect_jets_isb;
  vector < bool >   vect_jets_isc;

  vector < double > vect_bjets_pt;
  vector < double > vect_bjets_phi;
  vector < double > vect_bjets_eta;
  vector < bool > vect_bjets_isb;
  vector < bool > vect_bjets_isc;

  bool isb = false;
  bool isc = false;

  double sumVertexMassJet = 0.;
  double sumVertexMassTrk = 0.;
  double sumVertexMass = 0.;

  if (ee_event || mm_event) {

    for (vector < pat::Jet >::const_iterator jet = jets->begin (); jet != jets->end (); ++jet) {

      jet_pt  = jet->pt ();
      jet_eta = jet->eta();
      jet_phi = jet->phi();

      Ht += jet_pt;

      // JEC Uncertainty

      jecUnc->setJetPt(jet_pt);
      jecUnc->setJetEta(jet_eta);
      double unc = jecUnc->getUncertainty(true);
      double cor = (1.0+unc*par_);
      h_JEC_uncert->Fill (unc);
//      cout<< "JEC syst =" << unc << endl;

      jet_pt = jet->pt () * cor;

      if (jet_pt > 30) {

        ++Nj;

        vect_jets_pt.push_back (jet_pt);
        vect_jets_phi.push_back (jet_phi);
        vect_jets_eta.push_back (jet_eta);

        if (isMC && fabs (jet->partonFlavour ()) == 5) {
          isb = true;
          vect_jets_isb.push_back (true);
        } else {
          vect_jets_isb.push_back (false);
        }
        if (isMC && fabs (jet->partonFlavour ()) == 4) {
          isc = true;
          vect_jets_isc.push_back (true);
        } else {
          vect_jets_isc.push_back (false);
        }

        discrCSV = jet->bDiscriminator ("combinedSecondaryVertexBJetTags");
	discrBJP = jet->bDiscriminator ("jetBProbabilityBJetTags");
        discrJBP = jet->bDiscriminator ("jetProbabilityBJetTags");

	reco::SecondaryVertexTagInfo const * svTagInfos = jet->tagInfoSecondaryVertex("secondaryVertex");

        h_secondvtx_N->Fill (discrCSV);
        w_secondvtx_N->Fill (discrCSV, MyWeight);
        w_secondvtx_N_zoom->Fill (discrCSV, MyWeight);

	if (isb){
	  b_secondvtx_N->Fill (discrCSV, MyWeight);
	  b_secondvtx_N_zoom->Fill (discrCSV, MyWeight);
	}
	if (isc){
	  c_secondvtx_N->Fill (discrCSV, MyWeight);
	  c_secondvtx_N_zoom->Fill (discrCSV, MyWeight);
	}
  	  
	if ( svTagInfos && svTagInfos->nVertices() > 0 ) {
	  w_secondvtx_N_mass->Fill (discrCSV, MyWeight);

 	  if (isb){
	    b_secondvtx_N_mass->Fill (discrCSV, MyWeight);
	  }
  	  if (isc){
	    c_secondvtx_N_mass->Fill (discrCSV, MyWeight);
	  }
	}

	//cout << discrCSV << endl;

        if (discrCSV > 0.89) {
        
          w_BJP->Fill (discrBJP, MyWeight); 
          w_JBP->Fill (discrJBP, MyWeight);

	  if(isb){
  	     b_BJP->Fill (discrBJP, MyWeight);
  	     b_JBP->Fill (discrJBP, MyWeight);
	  }
	  if(isc){
	     c_BJP->Fill (discrBJP, MyWeight);
	     c_JBP->Fill (discrJBP, MyWeight);
	  }

	  ++Nb;
	  //cout << Nb << endl;

          bjet_pt  = jet->pt ();
          bjet_eta = jet->eta();
          bjet_phi = jet->phi();

          vect_bjets_pt.push_back(bjet_pt);
          vect_bjets_phi.push_back(bjet_phi);
          vect_bjets_eta.push_back(bjet_eta);

          if (isMC && fabs (jet->partonFlavour ()) == 5) {
            vect_bjets_isb.push_back (true);
          } else {
            vect_bjets_isb.push_back (false);
          }
          if (isMC && fabs (jet->partonFlavour ()) == 4) {
            vect_bjets_isc.push_back (true);
          } else {
            vect_bjets_isc.push_back (false);
          }

	  if ( svTagInfos && svTagInfos->nVertices() > 0 ) {
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

	  if ( svTagInfos && svTagInfos->nVertices() > 0 ) {
    	    const reco::Vertex &vertex = svTagInfos->secondaryVertex(0);
	    reco::TrackKinematics vertexKinematics(vertex);
	    bool useTrackWeights = true;
	    math::XYZTLorentzVector vertexSum = useTrackWeights
	      ? vertexKinematics.weightedVectorSum() : vertexKinematics.vectorSum();
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
//    cout << "VTX mass JET = " << sumVertexMassJet << endl;
//    cout << "VTX mass TRK = " << sumVertexMassTrk << endl;
//    cout << "VTX mass NEW = " << sumVertexMass << endl;
  }

  // ++++++++ MET & HT PLOTS

  float met = mets->empty() ? 0 : (*mets)[0].et();

  if ((ee_event || mm_event) && Nj > 0) {
    w_MET->Fill (met, MyWeight);
    if (isb) b_MET->Fill (met, MyWeight);
    if (isc) c_MET->Fill (met, MyWeight);
    w_MET_sign->Fill (met, MyWeight);
  }

  if ((ee_event || mm_event) && Nj > 0) {
    w_Ht->Fill (Ht, MyWeight);
    if (isb) b_Ht->Fill (Ht, MyWeight);
    if (isc) c_Ht->Fill (Ht, MyWeight);
    if (Nb > 0) {
      scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0]) : 1;
      //cout << vect_bjets_pt[0] << " " << vect_bjets_eta[0] <<"   SFb = " << scalFac_b << endl;
      w_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      if (isb) b_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      if (isc) c_Ht_b->Fill (Ht, MyWeight*scalFac_b);
    }
  }

  // ++++++++ DIELECTRON Z PLOTS

  if (ee_event) {
    w_numberOfZ->Fill (zee->size(), MyWeight);
    if (isb) b_numberOfZ->Fill (zee->size(), MyWeight);
    if (isc) c_numberOfZ->Fill (zee->size(), MyWeight);
  }

  if (ee_event && Nj > 0) {
    diele_mass = (*zee)[0].mass();
    diele_phi = (*zee)[0].phi();
    diele_pt = (*zee)[0].pt();

//    ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > z;
//    z = (*electrons)[0].ecalDrivenMomentum()+(*electrons)[1].ecalDrivenMomentum();
//    diele_mass = z.mass();
//    diele_phi = z.phi();
//    diele_pt = z.pt();

    if (diele_mass != 0) {
      h_mass_ee->Fill (diele_mass);
      w_mass_ee->Fill (diele_mass, MyWeight);
      w_pt_Z_ee->Fill (diele_pt, MyWeight);
      if (isb) {
        b_mass_ee->Fill (diele_mass, MyWeight);
        b_pt_Z_ee->Fill (diele_pt, MyWeight);
      }
      if (isc) {
        c_mass_ee->Fill (diele_mass, MyWeight);
        c_pt_Z_ee->Fill (diele_pt, MyWeight);
      }
      if (Nb > 0) {
        scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0]) : 1;
        w_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
        w_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
        delta_phi_ee = fabs (diele_phi - vect_bjets_phi[0]);
        if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
        w_delta_ee_b->Fill (delta_phi_ee, MyWeight*scalFac_b);
        if (isb) {
          b_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
          b_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
          b_delta_ee_b->Fill (delta_phi_ee, MyWeight*scalFac_b);
        }
        if (isc) {
          c_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
          c_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
          c_delta_ee_b->Fill (delta_phi_ee, MyWeight*scalFac_b);
        }
      }
    }
  }

  // ++++++++ DIMUON Z PLOTS

  if (mm_event) {
    w_numberOfZ->Fill (zmm->size(), MyWeight);
    if (isb) b_numberOfZ->Fill (zmm->size(), MyWeight);
    if (isc) c_numberOfZ->Fill (zmm->size(), MyWeight);
  }

  if (mm_event && Nj > 0) {
    dimuon_mass = (*zmm)[0].mass();
    dimuon_phi = (*zmm)[0].phi();
    dimuon_pt = (*zmm)[0].pt();

    if (dimuon_mass != 0)  {
      h_mass_mm->Fill (dimuon_mass);
      w_mass_mm->Fill (dimuon_mass, MyWeight);
      w_pt_Z_mm->Fill (dimuon_pt, MyWeight);
      if (isb) {
        b_mass_mm->Fill (dimuon_mass, MyWeight);
        b_pt_Z_mm->Fill (dimuon_pt, MyWeight);
      }
      if (isc) {
        c_mass_mm->Fill (dimuon_mass, MyWeight);
        c_pt_Z_mm->Fill (dimuon_pt, MyWeight);
      }
      if (Nb > 0) {
        scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0]) : 1;
        w_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
        w_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        delta_phi_mm = fabs (dimuon_phi - vect_bjets_phi[0]);
        if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
        w_delta_mm_b->Fill (delta_phi_mm, MyWeight*scalFac_b);
        if (isb) {
          b_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
          b_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
          b_delta_mm_b->Fill (delta_phi_mm, MyWeight*scalFac_b);
        }
        if (isc) {
          c_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
          c_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
          c_delta_mm_b->Fill (delta_phi_mm, MyWeight*scalFac_b);
        }
      }
    }
  }

  // ++++++++ MISC PLOTS

  if ((ee_event || mm_event) && Nj > 0) {
    h_pu_weights->Fill (MyWeight);
    Ntracks = tracks->size();
    h_tracks->Fill (Ntracks);
    w_tracks->Fill (Ntracks, MyWeight);
    h_recoVTX->Fill (NVtx);
    w_recoVTX->Fill (NVtx, MyWeight);
  }

  // ++++++++  ELECTRONS PLOTS

  if (ee_event && Nj > 0) {
    w_first_ele_pt->Fill (vect_ele_pt[0], MyWeight);
    w_first_ele_eta->Fill (vect_ele_eta[0], MyWeight);
    w_second_ele_pt->Fill (vect_ele_pt[1], MyWeight);
    w_second_ele_eta->Fill (vect_ele_eta[1], MyWeight);
    if (isb) {
      b_first_ele_pt->Fill (vect_ele_pt[0], MyWeight);
      b_first_ele_eta->Fill (vect_ele_eta[0], MyWeight);
      b_second_ele_pt->Fill (vect_ele_pt[1], MyWeight);
      b_second_ele_eta->Fill (vect_ele_eta[1], MyWeight);
    }
    if (isc) {
      c_first_ele_pt->Fill (vect_ele_pt[0], MyWeight);
      c_first_ele_eta->Fill (vect_ele_eta[0], MyWeight);
      c_second_ele_pt->Fill (vect_ele_pt[1], MyWeight);
      c_second_ele_eta->Fill (vect_ele_eta[1], MyWeight);
    }
  }

  if (ee_event && Nj > 0 && Nb > 0) {
    scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0]) : 1;
    w_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
    if (isb) b_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
    if (isc) c_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
  }

  // ++++++++ MUONS PLOTS

  if (mm_event && Nj > 0) {
    w_first_muon_pt->Fill (vect_muon_pt[0], MyWeight);
    w_first_muon_eta->Fill (vect_muon_eta[0], MyWeight);
    w_second_muon_pt->Fill (vect_muon_pt[1], MyWeight);
    w_second_muon_eta->Fill (vect_muon_eta[1], MyWeight);
    if (isb) {
      b_first_muon_pt->Fill (vect_muon_pt[0], MyWeight);
      b_first_muon_eta->Fill (vect_muon_eta[0], MyWeight);
      b_second_muon_pt->Fill (vect_muon_pt[1], MyWeight);
      b_second_muon_eta->Fill (vect_muon_eta[1], MyWeight);
    }
    if (isc) {
      c_first_muon_pt->Fill (vect_muon_pt[0], MyWeight);
      c_first_muon_eta->Fill (vect_muon_eta[0], MyWeight);
      c_second_muon_pt->Fill (vect_muon_pt[1], MyWeight);
      c_second_muon_eta->Fill (vect_muon_eta[1], MyWeight);
    }
  }

  if (mm_event && Nj > 0 && Nb > 0) {
    scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0]) : 1;
    w_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
    if (isb) b_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
    if (isc) c_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
  }

  // ++++++++ SVTX MASS PLOTS

  if ((ee_event || mm_event) && Nj > 0 && Nb > 0) {
    scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0]) : 1;
    w_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
    w_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
    w_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    if (isb) {
      b_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      b_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      b_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
    if (isc) {
      c_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      c_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      c_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
  }

  // ++++++++ JETS PLOTS

  if ((ee_event || mm_event) && Nj > 0) {
    h_jetmultiplicity->Fill (Nj);
    w_jetmultiplicity->Fill (Nj, MyWeight);
    w_first_jet_pt->Fill (vect_jets_pt[0], MyWeight);
    w_first_jet_eta->Fill (vect_jets_eta[0], MyWeight);
    if (vect_jets_isb[0]) {
      b_jetmultiplicity->Fill (Nj, MyWeight);
      b_first_jet_pt->Fill (vect_jets_pt[0], MyWeight);
      b_first_jet_eta->Fill (vect_jets_eta[0], MyWeight);
    }
    if (vect_jets_isc[0]) {
      c_jetmultiplicity->Fill (Nj, MyWeight);
      c_first_jet_pt->Fill (vect_jets_pt[0], MyWeight);
      c_first_jet_eta->Fill (vect_jets_eta[0], MyWeight);
    }
  }

  if ((ee_event || mm_event) && Nj > 1) {
    w_second_jet_pt->Fill (vect_jets_pt[1], MyWeight);
    w_second_jet_eta->Fill (vect_jets_eta[1], MyWeight);
    if (vect_jets_isb[1]) {
      b_second_jet_pt->Fill (vect_jets_pt[1], MyWeight);
      b_second_jet_eta->Fill (vect_jets_eta[1], MyWeight);
    }
    if (vect_jets_isc[1]) {
      c_second_jet_pt->Fill (vect_jets_pt[1], MyWeight);
      c_second_jet_eta->Fill (vect_jets_eta[1], MyWeight);
    }
  }

  if ((ee_event || mm_event) && Nj > 2) {
    w_third_jet_pt->Fill (vect_jets_pt[2], MyWeight);
    w_third_jet_eta->Fill (vect_jets_eta[2], MyWeight);
    if (vect_jets_isb[2]) {
      b_third_jet_pt->Fill (vect_jets_pt[2], MyWeight);
      b_third_jet_eta->Fill (vect_jets_eta[2], MyWeight);
    }
    if (vect_jets_isb[2]) {
      c_third_jet_pt->Fill (vect_jets_pt[2], MyWeight);
      c_third_jet_eta->Fill (vect_jets_eta[2], MyWeight);
    }
  }

  // ++++++++ B JETS PLOTS

  if ((ee_event || mm_event) && Nj > 0 && Nb > 0) {
    scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0]) : 1;
    w_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
    w_first_jet_pt_b->Fill (vect_jets_pt[0], MyWeight*scalFac_b);
    w_first_jet_eta_b->Fill (vect_jets_eta[0], MyWeight*scalFac_b);
    w_first_bjet_pt->Fill (vect_bjets_pt[0], MyWeight*scalFac_b);
    w_first_bjet_eta->Fill (vect_bjets_eta[0], MyWeight*scalFac_b);
    if (vect_jets_isb[0]) {
      b_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      b_first_jet_pt_b->Fill (vect_jets_pt[0], MyWeight*scalFac_b);
      b_first_jet_eta_b->Fill (vect_jets_eta[0], MyWeight*scalFac_b);
    }
    if (vect_bjets_isb[0]) {
      b_first_bjet_pt->Fill (vect_bjets_pt[0], MyWeight*scalFac_b);
      b_first_bjet_eta->Fill (vect_bjets_eta[0], MyWeight*scalFac_b);
    }
    if (vect_jets_isc[0]) {
      c_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      c_first_jet_pt_b->Fill (vect_jets_pt[0], MyWeight*scalFac_b);
      c_first_jet_eta_b->Fill (vect_jets_eta[0], MyWeight*scalFac_b);
    }
    if (vect_bjets_isc[0]) {
      c_first_bjet_pt->Fill (vect_bjets_pt[0], MyWeight*scalFac_b);
      c_first_bjet_eta->Fill (vect_bjets_eta[0], MyWeight*scalFac_b);
    }
  }

  if ((ee_event || mm_event) && Nj > 1 && Nb > 1) {
    scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0])*BtSF.Val(vect_bjets_pt[1], vect_bjets_eta[1]) : 1;
    w_second_jet_pt_b->Fill (vect_jets_pt[1], MyWeight*scalFac_b);
    w_second_jet_eta_b->Fill (vect_jets_eta[1], MyWeight*scalFac_b);
    w_second_bjet_pt->Fill (vect_bjets_pt[1], MyWeight*scalFac_b);
    w_second_bjet_eta->Fill (vect_bjets_eta[1], MyWeight*scalFac_b);
    if (vect_jets_isb[1]) {
      b_second_jet_pt_b->Fill (vect_jets_pt[1], MyWeight*scalFac_b);
      b_second_jet_eta_b->Fill (vect_jets_eta[1], MyWeight*scalFac_b);
    }
    if (vect_bjets_isb[1]) {
      b_second_bjet_pt->Fill (vect_bjets_pt[1], MyWeight*scalFac_b);
      b_second_bjet_eta->Fill (vect_bjets_eta[1], MyWeight*scalFac_b);
    }
    if (vect_jets_isc[1]) {
      c_second_jet_pt_b->Fill (vect_jets_pt[1], MyWeight*scalFac_b);
      c_second_jet_eta_b->Fill (vect_jets_eta[1], MyWeight*scalFac_b);
    }
    if (vect_bjets_isc[1]) {
      c_second_bjet_pt->Fill (vect_bjets_pt[1], MyWeight*scalFac_b);
      c_second_bjet_eta->Fill (vect_bjets_eta[1], MyWeight*scalFac_b);
    }
  }

  if ((ee_event || mm_event) && Nj > 2 && Nb > 2) {
    scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0])*BtSF.Val(vect_bjets_pt[1], vect_bjets_eta[1])*BtSF.Val(vect_bjets_pt[2], vect_bjets_eta[2]) : 1;
    w_third_jet_pt_b->Fill (vect_jets_pt[2], MyWeight*scalFac_b);
    w_third_jet_eta_b->Fill (vect_jets_eta[2], MyWeight*scalFac_b);
    w_third_bjet_pt->Fill (vect_bjets_pt[2], MyWeight*scalFac_b);
    w_third_bjet_eta->Fill (vect_bjets_eta[2], MyWeight*scalFac_b);
    if (vect_jets_isb[2]) {
      b_third_jet_pt_b->Fill (vect_jets_pt[2], MyWeight*scalFac_b);
      b_third_jet_eta_b->Fill (vect_jets_eta[2], MyWeight*scalFac_b);
    }
    if (vect_bjets_isb[2]) {
      b_third_bjet_pt->Fill (vect_bjets_pt[2], MyWeight*scalFac_b);
      b_third_bjet_eta->Fill (vect_bjets_eta[2], MyWeight*scalFac_b);
    }
    if (vect_jets_isc[2]) {
      c_third_jet_pt_b->Fill (vect_jets_pt[2], MyWeight*scalFac_b);
      c_third_jet_eta_b->Fill (vect_jets_eta[2], MyWeight*scalFac_b);
    }
    if (vect_bjets_isc[2]) {
      c_third_bjet_pt->Fill (vect_bjets_pt[2], MyWeight*scalFac_b);
      c_third_bjet_eta->Fill (vect_bjets_eta[2], MyWeight*scalFac_b);
    }
  }

  // ++++++++ EXTRA PLOTS

  if ((ee_event || mm_event) && Nj > 0 && Nb > 0) {
    scalFac_b = isMC ? BtSF.Val(vect_bjets_pt[0], vect_bjets_eta[0]) : 1;
    if (fabs (vect_bjets_eta[0]) > 0) Nf++;
    if (fabs (vect_bjets_eta[0]) < 0) Nbk++;
    if ((Nf+Nbk) != 0) Afb = (Nf - Nbk) / (Nf + Nbk);
    w_Afb->Fill (Afb, MyWeight*scalFac_b);
  }

  if (ee_event && Nj > 0) {
    h_scaleFactor_first_ele->Fill (scalFac_first_e, MyWeight / (scalFac_first_e * scalFac_second_e));
    h_scaleFactor_second_ele->Fill (scalFac_second_e, MyWeight / (scalFac_first_e * scalFac_second_e));
    if (isb) {
      b_scaleFactor_first_ele->Fill (scalFac_first_e, MyWeight / (scalFac_first_e * scalFac_second_e));
      b_scaleFactor_second_ele->Fill (scalFac_second_e, MyWeight / (scalFac_first_e * scalFac_second_e));
    }
  }
  if (mm_event && Nj > 0) {
    h_scaleFactor_first_muon->Fill (scalFac_first_m, MyWeight / (scalFac_first_m * scalFac_second_m));
    h_scaleFactor_second_muon->Fill (scalFac_second_m, MyWeight / (scalFac_first_m * scalFac_second_m));
    if (isb) {
      b_scaleFactor_first_muon->Fill (scalFac_first_m, MyWeight / (scalFac_first_m * scalFac_second_m));
      b_scaleFactor_second_muon->Fill (scalFac_second_m, MyWeight / (scalFac_first_m * scalFac_second_m));
    }
  }

  treeZb_->Fill();

}

// ------------ method called once each job just before starting event loop ------------
void ZbAnalyzer::beginJob () {
  jecUncDT_ = new JetCorrectionUncertainty("/gpfs/cms/users/candelis/work/ZbSkim/test/Fall12_V7_DATA_Uncertainty_AK5PFchs.txt");
  jecUncMC_ = new JetCorrectionUncertainty("/gpfs/cms/users/candelis/work/ZbSkim/test/Fall12_V7_MC_Uncertainty_AK5PFchs.txt");
  LumiWeights_ = edm::LumiReWeighting("/gpfs/cms/users/candelis/work/ZbSkim/test/pileup/pileup_" + pileup_ + ".root", "/gpfs/cms/users/candelis/work/ZbSkim/test/pileup/pileup_2012.root", "pileup", "pileup");

}

// ------------ method called once each job just after ending the event loop ------------
void ZbAnalyzer::endJob () {
  delete jecUncDT_;
  delete jecUncMC_;
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
