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
// $Id: ZbAnalyzer.cc,v 1.10 2013/04/17 13:28:55 vieri Exp $
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
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "table.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
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

table  MuSF  ("/gpfs/cms/users/lalicata/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/muon_eff.txt");
table  ElSF  ("/gpfs/cms/users/lalicata/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/ele_eff.txt");
table  BtSF  ("/gpfs/cms/users/candelis/CMSSW_5_3_9/src/ZbAnalysis/ZbSkim/test/btag_eff.txt");
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
  virtual void beginLuminosityBlock (edm::LuminosityBlock const &,
			edm::EventSetup const &);
  virtual void endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);

  std::string pileup_;
  std::string lepton_;
  double par;
  JetCorrectionUncertainty *jecUncDT;
  JetCorrectionUncertainty *jecUncMC;


  // ----------member data ---------------------------

  int       Nj, Nb;
  double    jet_pt;
  double    jet_eta;
  double    jet_phi;
  double    b_pt;
  double    b_eta;
  double    ele_pt;
  double    ele_eta;
  double    muon_pt;
  double    muon_eta;
  double    diele_inv;
  double    dimuon_inv;
  double    diele_phi;
  double    dimuon_phi;
  double    Delta_phi_ee;
  double    Delta_phi_mm;
  double    b_leading_pt;
  double    b_leading_eta;

  double    Ht, Ht_b;
  double    b_mm_invmass;
  double    b_ee_invmass;
  int       NZ, NZ2;
  double    discrCSV;
  double    MyWeight;
  double    sFac;
  double    sFacErr;
  double    scalFac_first_e;
  double    scalFac_second_e;
  double    scalFac_first_mu;
  double    scalFac_second_mu;
  double    scalFac_b;
  double    pt_Z;
  double    Nf, Nbk;
  double    Afb;

  TTree *   treeZb_;

  TH1F *    h_jetmultiplicity;
  TH1F *    h_jet_pt;
  TH1F *    h_ele_pt;
  TH1F *    h_muon_pt;
  TH1F *    h_mm_inv;
  TH1F *    h_ee_inv;
  TH1F *    h_secondvtx_N;
  TH1F *    h_PUweights;
  TH1F *    h_tracks;
  TH1D *    recoVTX_;
  TH1D *    recoVTX_w;

  TH1F *    w_first_jet_pt;
  TH1F *    w_first_jet_eta;
  TH1F *    w_first_jet_pt_b;
  TH1F *    w_first_jet_eta_b;
  TH1F *    w_first_ele_pt;
  TH1F *    w_second_ele_pt;
  TH1F *    w_first_muon_pt;
  TH1F *    w_second_muon_pt;
  TH1F *    w_first_ele_eta;
  TH1F *    w_second_ele_eta;
  TH1F *    w_first_muon_eta;
  TH1F *    w_second_muon_eta;
  TH1F *    h_delta_ee;
  TH1F *    h_delta_mm;
  TH1F *    h_pt_Z_ee_b;
  TH1F *    h_pt_Z_mm_b;
  TH1F *    h_pt_Z_ee;
  TH1F *    h_pt_Z_mm;
  TH1F *    w_jetmultiplicity;
  TH1F *    sf_first_ele_pt;
  TH1F *    w_jetmultiplicity_b;
  TH1F *    w_mm_inv;
  TH1F *    w_ee_inv;
  TH1F *    w_secondvtx_N;
  TH1F *    w_bquarks;
  TH1F *    w_tracks;
  TH1F *    flavours_;
  TH1F *    w_MET;
  TH1F *    w_MET_sign;
  TH1F *    w_MET_mtx;
  TH1F *    w_delta_phi_mm;
  TH1F *    numberOfZ;
  TH1F *    w_Ht;
  TH1F *    w_Ht_b;

  TH1F *    b_mm_inv;
  TH1F *    b_ee_inv;
  TH1F *    w_Afb;
  TH1F *    h_scalFactor_first_ele;
  TH1F *    h_scalFactor_first_muon;
  TH1F *    h_scalFactor_second_ele;
  TH1F *    h_scalFactor_second_muon;

  TH1F *    h_JEC_uncert;

  TH1F* b_jetmultiplicity;   
  TH1F* b_first_jet_pt;         
  TH1F* b_first_jet_eta;     
  TH1F* b_pt_Z_ee;              
  TH1F* b_pt_Z_mm;              
  TH1F* b_invMass_ee;
  TH1F* b_invMass_mm;
 
  TH1F* SVTX_mass_jet;
  TH1F* SVTX_mass_trk;
  TH1F* SVTX_mass_jet_b;
  TH1F* SVTX_mass_trk_b;

};

using namespace  pat;

// 
// constants, enums and typedefs
// 
enum Flavour
{
  ALL_JETS = 0,
  UDSG_JETS,
  C_JETS,
  B_JETS,
  NONID_JETS,
  N_JET_TYPES
};
struct Plots
{
  TH1 *    nVertices;
  TH1 *  deltaR, *  mass, *  dist, *  distErr, *    distSig;
  TH1 *  nTracks, *    chi2;
} plots_[N_JET_TYPES];



// 
// static data member definitions
// 



// 
// constructors and destructor
// 
ZbAnalyzer::ZbAnalyzer (const edm::ParameterSet & iConfig){

  pileup_ = iConfig.getUntrackedParameter < std::string > ("pileup", "S7");
  lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
  par  =    iConfig.getUntrackedParameter <double> ("JEC",  0);

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  h_jetmultiplicity = fs->make < TH1F > ("h_jetmultiplicity", "h_jetmultiplicity", 8, 0.5, 8.5);
  h_jet_pt = fs->make < TH1F > ("h_jet_pt", "h_jet_pt", 20, 30, 530);
  h_ele_pt = fs->make < TH1F > ("h_ele_pt", "h_ele_pt", 20, 30, 530);
  h_muon_pt = fs->make < TH1F > ("h_muon_pt", "h_muon_pt", 100, 0, 250);
  h_mm_inv = fs->make < TH1F > ("h_mm_inv", "h_mm_inv", 60, 60, 120);
  h_ee_inv = fs->make < TH1F > ("h_ee_inv", "h_ee_inv", 60, 60, 120);
  h_secondvtx_N = fs->make < TH1F > ("h_secondvtx_N", "h_secondvtx_N", 50, 0, 1);
  h_PUweights = fs->make < TH1F > ("h_pu_weights", "h_pu_weights", 10, 0, 10);
  recoVTX_ = fs->make < TH1D > ("recoVTX", "No. reconstructed vertices", 40, 0., 40.);
  recoVTX_w = fs->make < TH1D > ("recoVTXw", "No. reconstructed vertices weighted", 40, 0., 40.);
  h_tracks = fs->make < TH1F > ("h_tracks", "h_tracks", 50, 0, 50);

  // b fraction before btagging histograms

  b_jetmultiplicity =  fs->make < TH1F > ("b_jetmultiplicity", "b_jetmultiplicity", 8, 0.5, 8.5);
  b_first_jet_pt =     fs->make < TH1F > ("b_first_jet_pt",    "b_first_jet_pt", 20, 30., 500.);
  b_first_jet_eta =    fs->make < TH1F > ("b_first_jet_eta",   "b_first_jet_eta", 16, -2.5, 2.5);
  b_pt_Z_ee =          fs->make < TH1F > ("b_pt_Z_ee",           "b_pt_Z_ee", 35, 0., 350.);
  b_pt_Z_mm =          fs->make < TH1F > ("b_pt_Z_mm",           "b_pt_Z_mm", 35, 0., 350.);
  b_invMass_ee =       fs->make < TH1F > ("b_invMass_ee", "b_invMass_ee", 60, 71, 112);
  b_invMass_mm =       fs->make < TH1F > ("b_invMass_mm", "b_invMass_mm", 60, 71, 112);

  // weighted histograms
  w_first_jet_pt =   fs->make < TH1F > ("w_first_jet_pt", "w_first_jet_pt;P_t [GeV]", 20, 30., 500.);
  w_first_jet_eta =   fs->make < TH1F > ("w_first_jet_eta", "w_first_jet_eta;eta", 16, -2.5, 2.5);
  w_first_ele_pt =   fs->make < TH1F > ("first_ele_pt", "first_ele_pt;P_t [GeV]", 25, 0., 250.);
  w_second_ele_pt =  fs->make < TH1F > ("second_ele_pt", "second_ele_pt;P_t [GeV]", 25, 0., 250.);
  w_first_muon_pt =  fs->make < TH1F > ("first_muon_pt", "first_muon_pt;P_t [GeV]", 25, 0., 250.);
  w_second_muon_pt = fs->make < TH1F > ("second_muon_pt", "second_muon_pt;P_t [GeV]", 25, 0., 250.);
  w_first_ele_eta =  fs->make < TH1F > ("first_ele_eta", "first_ele_eta;Eta ", 16, -2.5, 2.5);
  w_second_ele_eta = fs->make < TH1F > ("second_ele_eta", "second_ele_eta;Eta ", 16, -2.5, 2.5);
  w_first_muon_eta = fs->make < TH1F > ("first_muon_eta", "first_muon_eta;Eta ", 16, -2.5, 2.5);
  w_second_muon_eta =fs->make < TH1F > ("second_muon_eta", "second_muon_eta;Eta ", 16, -2.5, 2.5);

  w_Ht = fs->make<TH1F>("w_Ht", "w_Ht [GeV]",20,30.,500.);
  w_Ht_b = fs->make<TH1F>("w_Ht_b", "w_Ht [GeV]",20,30.,500.);

  w_jetmultiplicity =  fs->make < TH1F > ("w_jetmultiplicity", "w_jetmultiplicity", 8, 0.5, 8.5);
  sf_first_ele_pt =    fs->make < TH1F > ("sf_first_ele_pt", "sf_first_ele_pt", 100, 0., 200.);
  w_jetmultiplicity_b =    fs->make < TH1F > ("w_jetmultiplicity_b", "w_jetmultiplicity_b", 5, 0.5, 5.5);
  w_first_jet_pt_b =      fs->make < TH1F > ("w_first_jet_pt_b", "w_first_jet_pt_b", 25, 25, 250);
  w_first_jet_eta_b =     fs->make < TH1F > ("w_first_jet_eta_b", "w_first_jet_eta_b", 16, -2.5, 2.5);
  w_mm_inv = 	       fs->make < TH1F > ("w_mm_inv", "w_mm_inv", 60, 71, 112);
  w_ee_inv = 	       fs->make < TH1F > ("w_ee_inv", "w_ee_inv", 60, 71, 112);
  w_secondvtx_N =      fs->make < TH1F > ("w_secondvtx_N", "w_secondvtx_N", 50, 0, 1);
  w_bquarks = 	       fs->make < TH1F > ("bquarks", "bquarks", 50, 0, 1);
  w_tracks = 	       fs->make < TH1F > ("w_tracks", "w_tracks", 50, 0, 50); 
  flavours_ = 	       fs->make < TH1F > ("flavours", "jet flavours", 5, 0, 5);
  w_MET = 	       fs->make < TH1F > ("w_MET", "w_MET", 50, 0, 250);
  w_MET_sign = 	       fs->make < TH1F > ("w_MET_sign", "w_MET_sign", 25, 0, 25);
  h_delta_ee =         fs->make < TH1F > ("w_delta_phi_ee", "w_delta_phi_ee", 12, 0, TMath::Pi ());
  h_delta_mm =         fs->make < TH1F > ("w_delta_phi_mm", "w_delta_phi_mm", 12, 0, TMath::Pi ());
  h_pt_Z_ee_b =        fs->make < TH1F > ("Z_pt_ee_b", "Z_pt_ee_b;P_t [GeV]", 35, 0., 350.);
  h_pt_Z_mm_b =        fs->make < TH1F > ("Z_pt_mm_b", "Z_pt_mm_b;P_t [GeV]", 35, 0., 350.);
  h_pt_Z_ee =          fs->make < TH1F > ("Z_pt_ee", "Z_pt_ee;P_t [GeV]", 35, 0., 350.);
  h_pt_Z_mm =          fs->make < TH1F > ("Z_pt_mm", "Z_pt_mm;P_t [GeV]", 35, 0., 350.);
  b_mm_inv = 	       fs->make < TH1F > ("b_mm_inv", "b_mm_inv", 60, 71, 112);
  b_ee_inv =           fs->make < TH1F > ("b_ee_inv", "b_ee_inv", 60, 71, 112);
  SVTX_mass_jet =      fs->make < TH1F > ("SVTX_mass_jet", "SVTX_mass_jet", 20, 0, 10);
  SVTX_mass_trk =      fs->make < TH1F > ("SVTX_mass_trk", "SVTX_mass_trk", 20, 0, 20);
  SVTX_mass_jet_b =      fs->make < TH1F > ("SVTX_mass_jet_b", "SVTX_mass_jet_b", 20, 0, 10);
  SVTX_mass_trk_b =      fs->make < TH1F > ("SVTX_mass_trk_b", "SVTX_mass_trk_b", 20, 0, 20);

  h_scalFactor_first_ele   =    fs->make < TH1F > ("scaleFactor_first_ele",   "scaleFactor_first_ele", 90, 0.6, 1.5);
  h_scalFactor_first_muon  =    fs->make < TH1F > ("scaleFactor_first_muon",  "scaleFactor_first_muon", 90, 0.6, 1.5);
  h_scalFactor_second_ele  =    fs->make < TH1F > ("scaleFactor_second_ele",  "scaleFactor_second_ele", 90, 0.6, 1.5);
  h_scalFactor_second_muon =    fs->make < TH1F > ("scaleFactor_second_muon", "scaleFactor_second_muon", 90, 0.6, 1.5);

  h_JEC_uncert             =    fs->make < TH1F > ("JEC uncert", "JEC uncert", 10, -0.5, 0.5);

  w_Afb = fs->make < TH1F > ("b_asymmetry", "b_asymmetry", 10, -1, 1);

  numberOfZ = fs->make < TH1F > ("numberOfZ",   "numberOfZ", 5, 0, 5);

  treeZb_ = fs->make < TTree > ("ZbTree", "ZbTree");
  treeZb_->Branch ("Nj", &Nj);
  treeZb_->Branch ("jet_pt", &jet_pt);
  treeZb_->Branch ("muon_pt", &muon_pt);
  treeZb_->Branch ("dimuon_inv", &dimuon_inv);
  treeZb_->Branch ("diele_inv", &diele_inv);
  treeZb_->Branch ("weights_pu", &MyWeight);

  for (unsigned int i = 0; i < N_JET_TYPES; i++)    {
      Plots & plots = plots_[i];
      const char *flavour, *name;


      switch ((Flavour) i)
	{
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
      plots.dist = fs->make < TH1F > (Form ("dist_%s", name),Form ("Transverse distance between PV and SV in %s", flavour), 100, 0, 2);

    }
}

ZbAnalyzer::~ZbAnalyzer ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

// 
// member functions
// 

// ------------ method called for each event ------------
void
ZbAnalyzer::analyze (const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace edm;
  using namespace std;


  // Get muon collection
  edm::Handle < pat::MuonCollection > Trigmuons;
  iEvent.getByLabel ("matchedMuons", Trigmuons);

  // Get muon collection
  edm::Handle < pat::ElectronCollection > Trigelectrons;
  iEvent.getByLabel ("matchedElectrons", Trigelectrons);

  // Get jet collection
  edm::Handle < std::vector < pat::Jet > >jets;
  iEvent.getByLabel ("goodJets", jets);

  // Get the Z->mm collection
  edm::Handle < reco::CompositeCandidateCollection > zmm;
  iEvent.getByLabel ("zmuMatchedmuMatched", zmm);

  // Get the Z->ee collection
  edm::Handle < reco::CompositeCandidateCollection > zee;
  iEvent.getByLabel ("zeleMatchedeleMatched", zee);

  // Get tracks
  edm::Handle < std::vector < reco::Track > >tracks;
  iEvent.getByLabel ("generalTracks", tracks);

  // Get MET
  edm::Handle < std::vector < pat::MET > >met;
  iEvent.getByLabel (edm::InputTag ("patMETsPFlow"), met);

  bool ee_event = false;
  bool mm_event = false;
  int Ntracks = 0;
  int Nmu = 0;
  Nj = 0;
  Nbk = 0;
  jet_pt = 0;
  jet_eta = 0;
  b_pt = 0;
  b_eta = 0;
  jet_phi = 0;
  ele_pt = 0;
  ele_eta = 0;
  muon_eta = 0;
  muon_pt = 0;
  diele_inv = 0;
  dimuon_inv = 0;
  NZ = 0;
  NZ2 = 0;
  discrCSV = 0;
  b_leading_pt = 0;
  b_leading_eta = 0;
  Delta_phi_mm = 0;
  Delta_phi_ee = 0;
  b_mm_invmass = 0;
  b_ee_invmass = 0;
  Nf = 0;
  Nb = 0;
  Afb = 0;
  scalFac_first_e = 1;
  scalFac_second_e = 1;
  scalFac_first_mu = 1;
  scalFac_second_mu = 1;
  scalFac_b = 1;

  Ht = 0;
  Ht_b = 0;
  // +++++++++ ELECTRONS

  vector < double >vect_ele_pt;
  vector < double >vect_ele_eta;

  for (pat::ElectronCollection::const_iterator ele = Trigelectrons->begin (); ele != Trigelectrons->end (); ++ele){   

    	  ele_pt = ele->pt ();
    	  ele_eta = ele->eta ();

    	  vect_ele_pt.push_back (ele_pt);
    	  vect_ele_eta.push_back (ele_eta);
  }

//  if (vect_ele_pt.size () >= 2){
//      	  if (fabs (vect_ele_eta[0]) < 2.4 && fabs (vect_ele_eta[1]) < 2.4 && (fabs (vect_ele_eta[0]) < 1.442 || fabs (vect_ele_eta[0]) > 1.566) && (fabs (vect_ele_eta[1]) < 1.442 || fabs (vect_ele_eta[1]) > 1.566)){
//	    	  if (vect_ele_pt[0] > 20 && vect_ele_pt[1] > 20){
//			ee_event = true;
//		  }
//	  }
//  }

  if (zee->size () != 0) {
	  ee_event = true;
  }

  // +++++++++ MUONS

  vector < double >vect_muon_pt;
  vector < double >vect_muon_eta;

  for (pat::MuonCollection::const_iterator muon = Trigmuons->begin (); muon != Trigmuons->end (); ++muon) {

      muon_pt = muon->pt ();
      muon_eta = muon->eta ();

      vect_muon_pt.push_back (muon_pt);
      vect_muon_eta.push_back (muon_eta);
  }

//  if (vect_muon_pt.size () >= 2){
//      	  if (fabs (vect_muon_eta[0]) < 2.4 && fabs (vect_muon_eta[1]) < 2.4){
//	    	  if (vect_muon_pt[0] > 20 && vect_muon_pt[1] > 20){
//			mm_event = true;
//		  }
//	  }
//  }

  if (zmm->size () != 0){
	  mm_event = true;
  }

  if (lepton_ == "electron" && !ee_event) return;
  if (lepton_ == "muon" && !mm_event)     return;

  // ++++++ Pile-Up

  bool isMC = false;
  JetCorrectionUncertainty * jecUnc;

  jecUnc = jecUncDT;

  MyWeight = 1.0;

  Handle < std::vector < PileupSummaryInfo > >PupInfo;

  if (iEvent.getByLabel (edm::InputTag ("addPileupInfo"), PupInfo))  {
      
      isMC = true;
      jecUnc=jecUncMC;

      std::vector < PileupSummaryInfo >::const_iterator PVI;

      float Tnpv = -1;

      for (PVI = PupInfo->begin (); PVI != PupInfo->end (); ++PVI) {
    	      int BX = PVI->getBunchCrossing ();
    	      if (BX == 0) {
		      Tnpv = PVI->getTrueNumInteractions ();
		      continue;
  	      }
      }

      edm::LumiReWeighting LumiWeights_;
      LumiWeights_ = edm::LumiReWeighting("/gpfs/cms/users/candelis/work/Zb/pileup/pileup_" + pileup_ + ".root", "/gpfs/cms/users/candelis/work/Zb/pileup/pileup_2012.root", "pileup", "pileup");

      MyWeight = LumiWeights_.weight (Tnpv);


      if (ee_event){
	  scalFac_first_e  =  ElSF.Val (vect_ele_pt[0], vect_ele_eta[0]);
	  scalFac_second_e =  ElSF.Val (vect_ele_pt[1], vect_ele_eta[1]);
      }
      if (mm_event){
	  scalFac_first_mu  = MuSF.Val (vect_muon_pt[0], vect_muon_eta[0]);
	  scalFac_second_mu = MuSF.Val (vect_muon_pt[1], vect_muon_eta[1]);
      } 


  }


  if (ee_event) MyWeight = MyWeight * scalFac_first_e * scalFac_second_e;
  if (mm_event) MyWeight = MyWeight * scalFac_first_mu * scalFac_second_mu;

  if (ee_event)	  numberOfZ->Fill(zee->size ());
  if (mm_event)	  numberOfZ->Fill(zmm->size());

  // ++++++++ VERTICES

  edm::Handle < std::vector < reco::Vertex > >vertices_h;
  iEvent.getByLabel (edm::InputTag ("goodOfflinePrimaryVertices"), vertices_h);
  if (!vertices_h.isValid ()) {
      // std::cout<<"empty vertex collection!!!\n";
      // return;
  }

  // require in the event that there is at least one reconstructed vertex
  if (vertices_h->size () <= 0)    return;
  // pick the first (i.e. highest sum pt) vertex
  const reco::Vertex * theVertex = &(vertices_h->front ());
  // require that the vertex meets certain criteria
  if (theVertex->ndof () < 5)    return;
  if (fabs (theVertex->z ()) > 24.0)    return;
  if (fabs (theVertex->position ().rho ()) > 2.0)    return;

  std::vector < reco::Vertex >::const_iterator itv;
  int NVtx = 0;
  // now, count vertices
  for (itv = vertices_h->begin (); itv != vertices_h->end (); ++itv) {
      // require that the vertex meets certain criteria
      if (itv->ndof () < 5)	continue;
      if (fabs (itv->z ()) > 50.0)	continue;
      if (fabs (itv->position ().rho ()) > 2.0)	continue;
      ++NVtx;
  } 

  // +++++++++ JETS

  vector < double >vect_jet_pt;
  vector < double >vect_jet_phi;
  vector < double >vect_jet_eta;
  vector < double >vect_jet_discrCSV;
  bool isb = false;

  double sumVertexMassJet = 0.;
  double sumVertexMassTrk = 0.;

  if (ee_event || mm_event) {

      for (std::vector < pat::Jet >::const_iterator jet = jets->begin (); jet != jets->end (); ++jet) {

	  jet_pt  = jet->pt ();
	  jet_eta = jet->eta();
	  jet_phi = jet->phi();
   
	  // for events with a generated b
	  b_pt    = jet->pt ();
          b_eta   = jet->eta();

	  Ht += jet_pt;

	  // JEC Uncertainty 
	      
	  jecUnc->setJetPt(jet_pt);
	  jecUnc->setJetEta(jet_eta);
	  double unc = jecUnc->getUncertainty(true);
	  double cor = (1.0+unc*par);
	  h_JEC_uncert -> Fill(unc); 
//	  cout<< "JEC syst =" << unc << endl;

	  jet_pt  = jet->pt () * cor;

	  if (jet_pt > 30) {
	      
	      ++Nj;

	      vect_jet_pt.push_back (jet_pt);
	      vect_jet_phi.push_back (jet_phi);
	      vect_jet_eta.push_back (jet_eta);
              

	      
	      /*
	       * 
	       * // Flavour studies
	       * 
	       * Flavour flavour; // find out the jet flavour (differs
	       * between quark and anti-quark)
	       * switch(std::abs(jet->partonFlavour())) { case 1: case
	       * 2: case 3: case 21: flavour = UDSG_JETS; break; case 4:
	       * flavour = C_JETS; break; case 5: flavour = B_JETS;
	       * break; default: flavour = NONID_JETS; }
	       * 
	       * 
	       * const reco::SecondaryVertexTagInfo &svTagInfo =
	       * *jet->tagInfoSecondaryVertex(); if
	       * (svTagInfo.nVertices() < 1) continue; const reco::Vertex 
	       * &sv = svTagInfo.secondaryVertex(0);
	       * 
	       * // simply count the number of accepted jets
	       * flavours_->Fill(ALL_JETS, MyWeight);
	       * flavours_->Fill(flavour, MyWeight);
	       * 
	       * // dxy Measurement1D distance =
	       * svTagInfo.flightDistance(0, true);
	       * plots_[flavour].dist->Fill(distance.value());
	       * 
	       */
	
	      // b studies
	  
	      discrCSV = jet->bDiscriminator ("combinedSecondaryVertexBJetTags");
     
	      vect_jet_discrCSV.push_back (discrCSV);

	      if (fabs (jet->partonFlavour ()) == 5) isb = true;

              h_secondvtx_N->Fill (discrCSV);
	      w_secondvtx_N->Fill (discrCSV, MyWeight);

	      //cout<<discrCSV<<endl;

	      if (discrCSV > 0.89){
		   ++Nb;
		   //cout<<Nb<<endl;

		   reco::SecondaryVertexTagInfo const * svTagInfos = jet->tagInfoSecondaryVertex("secondaryVertex");

		   if ( svTagInfos && svTagInfos->nVertices() > 0 ) {

		     ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > sumVecJet;
		     reco::Vertex::trackRef_iterator kEndTracks = svTagInfos->secondaryVertex(0).tracks_end();
		     for (reco::Vertex::trackRef_iterator trackIter = svTagInfos->secondaryVertex(0).tracks_begin(); trackIter != kEndTracks; ++trackIter ) {
		        const double kPionMass = 0.13957018;
		        ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > vec;
		        vec.SetPx( (*trackIter)->px() );
		        vec.SetPy( (*trackIter)->py() );
		        vec.SetPz( (*trackIter)->pz() );
		        vec.SetM (kPionMass);
		        sumVecJet += vec;
		     }
		     sumVertexMassJet += sumVecJet.M();

		     ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > sumVecTrk;
		     for(size_t itrack=0; itrack < jet->associatedTracks().size(); ++itrack){
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

	      }
	  }
      }

      if (Nb) {
	sumVertexMassJet /= Nb;
	sumVertexMassTrk /= Nb;

	SVTX_mass_jet -> Fill(sumVertexMassJet);
	SVTX_mass_trk -> Fill(sumVertexMassTrk);

cout << "Nb = " << Nb << endl;
cout << "sumVertexMassJet = " << sumVertexMassJet << endl;
cout << "sumVertexMassTrk = " << sumVertexMassTrk << endl;
      }

      if (Nj != 0) {
        w_first_jet_pt->Fill (vect_jet_pt[0], MyWeight);
        w_first_jet_eta->Fill (vect_jet_eta[0], MyWeight);
	w_Ht->Fill(Ht, MyWeight);
        recoVTX_->Fill (float (NVtx));
        recoVTX_w->Fill (NVtx, MyWeight);
        if (isb){
		w_bquarks->Fill (discrCSV, MyWeight);
		b_first_jet_pt->Fill (jet_pt, MyWeight);
		b_first_jet_eta->Fill (jet_eta, MyWeight);
		SVTX_mass_trk_b->Fill(sumVertexMassTrk);
		SVTX_mass_jet_b->Fill(sumVertexMassJet);
  	}
      }
	      	  
      if(isMC){
//	      scalFac_b = BtSF.Val (vect_jet_pt[0], vect_jet_eta[0]);
	      //cout<<scalFac_b<<endl;
      }

      if (Nb != 0) {
	      
        w_jetmultiplicity_b->Fill (Nb, MyWeight*scalFac_b);
	w_Ht_b->Fill(Ht, MyWeight*scalFac_b);
        w_first_jet_pt_b->Fill (vect_jet_pt[0], MyWeight*scalFac_b);
        
	w_first_jet_eta_b->Fill (vect_jet_eta[0], MyWeight*scalFac_b);
        if (fabs (vect_jet_eta[0]) > 0) Nf++;
        if (fabs (vect_jet_eta[0]) < 0) Nbk++;
        Afb = (Nf - Nbk) / (Nf + Nbk);
        w_Afb->Fill (Afb, MyWeight*scalFac_b);
      }

  }

  // MET

  double significance = 0;
  if ((ee_event || mm_event) && Nj != 0) {
    	  for (std::vector < pat::MET >::const_iterator m = met->begin (); m != met->end (); ++m) {
		  w_MET->Fill (m->et (), MyWeight);
		  significance = m->significance();
		  w_MET_sign->Fill(significance, MyWeight);
  	  }
  }

  // TRACKS
  
  if ((ee_event || mm_event) && Nj != 0) {
      for (std::vector < reco::Track >::const_iterator t = tracks->begin (); t != tracks->end (); ++t) {
    	      Ntracks = t->numberOfValidHits ();
      }
  }

  // DIMUON Z

  if (mm_event && Nj != 0) {

    std::vector < reco::CompositeCandidate >::const_iterator Zit = zmm->begin ();

    dimuon_inv = Zit->mass ();
    dimuon_phi = Zit->phi ();
    pt_Z = Zit->pt ();

    if (dimuon_inv != 0)  {
	  
	  h_mm_inv->Fill (dimuon_inv);
    	  w_mm_inv->Fill (dimuon_inv, MyWeight);
          h_pt_Z_mm->Fill (pt_Z, MyWeight);
	  
	  if(isb){
		  b_pt_Z_mm->Fill(pt_Z, MyWeight);
	          b_invMass_mm->Fill(dimuon_inv, MyWeight);
      	  }
	  if (Nb != 0) {
		  b_mm_invmass = Zit->mass ();
		  b_mm_inv->Fill (b_mm_invmass, MyWeight*scalFac_b);
		  Delta_phi_mm = fabs (dimuon_phi - vect_jet_phi[0]);
		  
		  if (Delta_phi_mm > acos (-1)) Delta_phi_mm = 2 * acos (-1) - Delta_phi_mm;

		  h_delta_mm->Fill (Delta_phi_mm, MyWeight*scalFac_b);
		  h_pt_Z_mm_b->Fill (pt_Z, MyWeight*scalFac_b);
  	  }
    }
  }
 
  // DIELECTRON Z

  if (ee_event && Nj != 0) {

    std::vector <reco::CompositeCandidate >::const_iterator Zit2 = zee->begin ();

    diele_inv = Zit2->mass ();
    diele_phi = Zit2->phi ();
    pt_Z = Zit2->pt ();

    if (diele_inv != 0) {
    	  h_ee_inv->Fill (diele_inv);
    	  w_ee_inv->Fill (diele_inv, MyWeight);
          h_pt_Z_ee->Fill (pt_Z, MyWeight);
	  
	  if(isb){
		  b_pt_Z_ee->Fill(pt_Z, MyWeight);
	          b_invMass_ee->Fill(diele_inv, MyWeight);
      	  }

	  if (Nb != 0){
    		  b_ee_invmass = Zit2->mass ();
    		  b_ee_inv->Fill (b_ee_invmass, MyWeight*scalFac_b);
    		  Delta_phi_ee = fabs (diele_phi - vect_jet_phi[0]);
    		  if (Delta_phi_ee > acos (-1))  Delta_phi_ee = 2 * acos (-1) - Delta_phi_ee;
    		  h_delta_ee->Fill (Delta_phi_ee, MyWeight*scalFac_b);
    		  h_pt_Z_ee_b->Fill (pt_Z, MyWeight*scalFac_b);
    	  }
    }
  }

  treeZb_->Fill();

  if ((ee_event || mm_event) && Nj != 0) {
    h_PUweights->Fill (MyWeight);
  }

  if ((ee_event || mm_event) && Nj != 0) {
    h_scalFactor_first_ele->Fill (scalFac_first_e);
    h_scalFactor_first_muon->Fill (scalFac_first_mu);
    h_scalFactor_second_ele->Fill (scalFac_second_e);
    h_scalFactor_second_muon->Fill (scalFac_second_mu);
  }

  if ((ee_event || mm_event) && Nj != 0) {
      h_jetmultiplicity->Fill (Nj);
      w_jetmultiplicity->Fill (Nj, MyWeight);
  }

  if (ee_event && Nj != 0) {
  	w_first_ele_pt->Fill (vect_ele_pt[0], MyWeight);
    	sf_first_ele_pt->Fill (vect_ele_pt[0], MyWeight / (scalFac_first_e * scalFac_second_e));
  	w_second_ele_pt->Fill (vect_ele_pt[1], MyWeight);
 	w_first_ele_eta->Fill (vect_ele_eta[0], MyWeight);
	w_second_ele_eta->Fill (vect_ele_eta[1], MyWeight);
  }

  if (mm_event && Nj != 0) {
      w_first_muon_pt->Fill (vect_muon_pt[0], MyWeight);
      w_second_muon_pt->Fill (vect_muon_pt[1], MyWeight);
      w_first_muon_eta->Fill (vect_muon_eta[0], MyWeight);
      w_second_muon_eta->Fill (vect_muon_eta[1], MyWeight);
  }

  if ((ee_event || mm_event) && Nj != 0) {
    h_tracks->Fill (Ntracks);
    w_tracks->Fill (Ntracks, MyWeight);
  }

}

// ------------ method called once each job just before starting event
// loop ------------
void
ZbAnalyzer::beginJob ()
{
  jecUncDT = new JetCorrectionUncertainty("/gpfs/cms/users/candelis/work/Zb/Fall12_V7_DATA_Uncertainty_AK5PFchs.txt");
  jecUncMC = new JetCorrectionUncertainty("/gpfs/cms/users/candelis/work/Zb/Fall12_V7_MC_Uncertainty_AK5PFchs.txt");
}

// ------------ method called once each job just after ending the event
// loop ------------
void
ZbAnalyzer::endJob ()
{
  delete jecUncDT;
  delete jecUncMC;
}

// ------------ method called when starting to processes a run
// ------------
void
ZbAnalyzer::beginRun (edm::Run const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a run
// ------------
void
ZbAnalyzer::endRun (edm::Run const &, edm::EventSetup const &)
{
}

// ------------ method called when starting to processes a luminosity
// block ------------
void
ZbAnalyzer::beginLuminosityBlock (edm::LuminosityBlock const &,
				  edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a luminosity
// block ------------
void
ZbAnalyzer::endLuminosityBlock (edm::LuminosityBlock const &,
				edm::EventSetup const &)
{
}

// define this as a plug-in
DEFINE_FWK_MODULE (ZbAnalyzer);
