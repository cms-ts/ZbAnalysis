// -*- C++ -*-
//
// Package:    ZbDumper
// Class:      ZbDumper
//
/**\class ZbDumper ZbDumper.cc ZbAnalysis/ZbDumper/src/ZbDumper.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vieri Candelise
//         Created:  Thu Jul 18 15:17:47 CEST 2013
// $Id$
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
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
#include "Math/VectorUtil.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class ZbDumper : public edm::EDAnalyzer {
   public:
      explicit ZbDumper(const edm::ParameterSet&);
      ~ZbDumper();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

     std::string lepton_;
     std::string pileupDT_;
     double par_;
     double par2_;
     bool usePartonFlavour_;
     bool pcut_;
     bool useDeltaR_;
     double numB_;

     TH2F* w_first_jet_pt;
     TH2F* w_first_jet_eta;
     TH2F* w_first_jet_eta_abs;
     TH2F* w_first_bjet_pt;
     TH2F* w_first_bjet_eta;
     TH2F* w_first_bjet_eta_abs;
     TH2F* w_second_bjet_pt;
     TH2F* w_second_bjet_eta;
     TH2F* w_second_bjet_eta_abs;
     TH2F* w_pt_Z_ee;
     TH2F* w_pt_Z_mm;
     TH2F* w_pt_Z_ee_b;
     TH2F* w_pt_Z_mm_b;
     TH2F* w_y_Z_ee;
     TH2F* w_y_Z_mm;
     TH2F* w_y_Z_ee_b;
     TH2F* w_y_Z_mm_b;
     TH2F* w_y_Z_ee_abs;
     TH2F* w_y_Z_mm_abs;
     TH2F* w_y_Z_ee_b_abs;
     TH2F* w_y_Z_mm_b_abs;
     TH2F* w_mass_Zj_ee;
     TH2F* w_mass_Zj_mm;
     TH2F* w_mass_Zj_ee_b;
     TH2F* w_mass_Zj_mm_b;
     TH2F* w_Ht;
     TH2F* w_Ht_b;
     TH2F* w_delta_ee;
     TH2F* w_delta_mm;
     TH2F* w_delta_ee_b;
     TH2F* w_delta_mm_b;
     TH2F* w_delta_phi_2b;
     TH2F* w_DR_bb;
     TH2F* w_DR_eeb_min;
     TH2F* w_DR_eeb_max;
     TH2F* w_DR_mmb_min;
     TH2F* w_DR_mmb_max;
     TH2F* w_Phi_star_ee;
     TH2F* w_Phi_star_ee_b;
     TH2F* w_Phi_star_mm;
     TH2F* w_Phi_star_mm_b;
     TH2F* w_bb_mass;
     TH2F* w_eebb_mass;
     TH2F* w_mmbb_mass;
     TH2F* w_A_eeb;
     TH2F* w_A_mmb;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZbDumper::ZbDumper(const edm::ParameterSet& iConfig) {

   lepton_           = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
   pileupDT_         = iConfig.getUntrackedParameter < std::string > ("pileupDT", "");
   par_              = iConfig.getUntrackedParameter <double> ("JEC", 0);
   par2_             = iConfig.getUntrackedParameter <double> ("JER", 0);
   pcut_             = iConfig.getUntrackedParameter <bool> ("pcut", false);   
   useDeltaR_        = iConfig.getUntrackedParameter <bool> ("useDeltaR", false);
   numB_             = iConfig.getUntrackedParameter <double> ("numB", 0);
   //now do what ever initialization is needed
   edm::Service < TFileService > fs;

   w_first_jet_pt    = fs->make < TH2F > ("w_first_jet_pt",    "w_first_jet_pt;P_t [GeV]", 50, 30., 700., 50, 30., 700.);
   w_first_jet_eta   = fs->make < TH2F > ("w_first_jet_eta",   "w_first_jet_eta;Eta", 16, -2.5, 2.5,16, -2.5, 2.5);
   w_first_jet_eta_abs   = fs->make < TH2F > ("w_first_jet_eta_abs",   "w_first_jet_eta_abs;abs(Eta)", 8, 0, 2.5, 8, 0, 2.5);
   w_first_bjet_pt   = fs->make < TH2F > ("w_first_bjet_pt",   "w_first_bjet_pt;P_t [GeV]", 50, 30., 700., 50, 30., 700.);
   w_first_bjet_eta  = fs->make < TH2F > ("w_first_bjet_eta",  "w_first_bjet_eta;Eta", 16, -2.5, 2.5,16, -2.5, 2.5);
   w_first_bjet_eta_abs  = fs->make < TH2F > ("w_first_bjet_eta_abs",  "w_first_bjet_eta_abs;abs(Eta)", 8, 0, 2.5, 8, 0, 2.5);
   w_second_bjet_pt   = fs->make < TH2F > ("w_second_bjet_pt",   "w_second_bjet_pt;P_t [GeV]", 50, 30., 500., 50, 30., 500.);
   w_second_bjet_eta  = fs->make < TH2F > ("w_second_bjet_eta",  "w_second_bjet_eta;Eta", 16, -2.5, 2.5,16, -2.5, 2.5);
   w_second_bjet_eta_abs  = fs->make < TH2F > ("w_second_bjet_eta_abs",  "w_second_bjet_eta_abs;abs(Eta)", 8, 0, 2.5, 8, 0, 2.5);
   w_pt_Z_ee         = fs->make < TH2F > ("w_pt_Z_ee",         "w_pt_Z_ee;P_t [GeV]", 40, 0., 400., 40, 0., 400.);
   w_pt_Z_mm         = fs->make < TH2F > ("w_pt_Z_mm",         "w_pt_Z_mm;P_t [GeV]", 40, 0., 400., 40, 0., 400.);
   w_pt_Z_ee_b 	     = fs->make < TH2F > ("w_pt_Z_ee_b",       "w_pt_Z_ee_b;P_t [GeV]", 40, 0., 400., 40, 0., 400.);
   w_pt_Z_mm_b       = fs->make < TH2F > ("w_pt_Z_mm_b",       "w_pt_Z_mm_b;P_t [GeV]", 40, 0., 400., 40, 0., 400.);
   w_y_Z_ee          = fs->make < TH2F > ("w_y_Z_ee",          "w_y_Z_ee;y",   20, -2.0, 2.0, 20, -2.0, 2.0);
   w_y_Z_mm          = fs->make < TH2F > ("w_y_Z_mm",          "w_y_Z_mm;y",   20, -2.0, 2.0, 20, -2.0, 2.0);
   w_y_Z_ee_b 	     = fs->make < TH2F > ("w_y_Z_ee_b",        "w_y_Z_ee_b;y", 20, -2.0, 2.0, 20, -2.0, 2.0);
   w_y_Z_mm_b        = fs->make < TH2F > ("w_y_Z_mm_b",        "w_y_Z_mm_b;y", 20, -2.0, 2.0, 20, -2.0, 2.0);

   w_y_Z_ee_abs          = fs->make < TH2F > ("w_y_Z_ee_abs",          "w_y_Z_ee_abs;abs(y)",   10, 0, 2.0, 10, 0, 2.0);
   w_y_Z_mm_abs          = fs->make < TH2F > ("w_y_Z_mm_abs",          "w_y_Z_mm_abs;abs(y)",   10, 0, 2.0, 10, 0, 2.0);
   w_y_Z_ee_b_abs        = fs->make < TH2F > ("w_y_Z_ee_b_abs",        "w_y_Z_ee_b_abs;abs(y)", 10, 0, 2.0, 10, 0, 2.0);
   w_y_Z_mm_b_abs        = fs->make < TH2F > ("w_y_Z_mm_b_abs",        "w_y_Z_mm_b_abs;abs(y)", 10, 0, 2.0, 10, 0, 2.0);


   w_mass_Zj_ee      = fs->make < TH2F > ("w_mass_Zj_ee",      "w_mass_Zj_ee;P_t [GeV]", 15, 100., 330., 15, 100., 330.);
   w_mass_Zj_mm      = fs->make < TH2F > ("w_mass_Zj_mm",      "w_mass_Zj_mm;P_t [GeV]", 15, 100., 330., 15, 100., 330.);
   w_mass_Zj_ee_b    = fs->make < TH2F > ("w_mass_Zj_ee_b",    "w_mass_Zj_ee_b;P_t [GeV]", 15, 100., 330., 15, 100., 330.);
   w_mass_Zj_mm_b    = fs->make < TH2F > ("w_mass_Zj_mm_b",    "w_mass_Zj_mm_b;P_t [GeV]", 15, 100., 330., 15, 100., 330.);

   w_Ht 	     = fs->make < TH2F > ("w_Ht",              "w_Ht [GeV]", 50, 30., 1000.,50, 30., 1000.);
   w_Ht_b 	     = fs->make < TH2F > ("w_Ht_b",            "w_Ht [GeV]", 50, 30., 1000.,50, 30., 1000.);
   w_delta_ee        = fs->make < TH2F > ("w_delta_phi_ee",    "w_delta_phi_ee", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());
   w_delta_mm        = fs->make < TH2F > ("w_delta_phi_mm",    "w_delta_phi_mm", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());
   w_delta_ee_b      = fs->make < TH2F > ("w_delta_phi_ee_b",  "w_delta_phi_ee_b", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());
   w_delta_mm_b      = fs->make < TH2F > ("w_delta_phi_mm_b",  "w_delta_phi_mm_b", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());

   w_delta_phi_2b    = fs->make < TH2F > ("w_delta_phi_2b",    "w_delta_phi_2b", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());
  
   w_DR_bb           = fs->make < TH2F > ("w_DR_bb",           "w_DR_bb", 10, 0, 5, 10, 0, 5);
 
   w_DR_eeb_min      = fs->make < TH2F > ("w_DR_eeb_min",      "w_DR_eeb_min", 15, 0., 4., 15, 0., 4.);
   w_DR_eeb_max      = fs->make < TH2F > ("w_DR_eeb_max",      "w_DR_eeb_max", 15, 1.5, 5., 15, 1.5, 5.);
   w_DR_mmb_min      = fs->make < TH2F > ("w_DR_mmb_min",      "w_DR_mmb_min", 15, 0., 4., 15, 0., 4.);
   w_DR_mmb_max      = fs->make < TH2F > ("w_DR_mmb_max",      "w_DR_mmb_max", 15, 1.5, 5., 15, 1.5, 5.);

   w_Phi_star_ee     = fs->make < TH2F > ("w_Phi_star_ee",     "w_Phi_star_ee; Phi*", 10, 0., 1., 10, 0., 1.);
   w_Phi_star_ee_b   = fs->make < TH2F > ("w_Phi_star_ee_b",   "w_Phi_star_ee_b; Phi*", 10, 0., 1., 10, 0., 1.);
   w_Phi_star_mm     = fs->make < TH2F > ("w_Phi_star_mm",     "w_Phi_star_mm; Phi*", 10, 0., 1., 10, 0., 1.);
   w_Phi_star_mm_b   = fs->make < TH2F > ("w_Phi_star_mm_b",   "w_Phi_star_mm_b; Phi*", 10, 0., 1., 10, 0., 1.);

   w_bb_mass         = fs->make < TH2F > ("w_bb_mass",         "w_bb_mass;Mass [GeV]", 13, 30., 400., 13, 30., 400.);   
   
   w_eebb_mass       = fs->make < TH2F > ("w_eebb_mass",       "w_eebb_mass;Mass [GeV]", 13, 180., 500., 13, 180., 500.);
   w_mmbb_mass       = fs->make < TH2F > ("w_mmbb_mass",       "w_mmbb_mass;Mass [GeV]", 13, 180., 500., 13, 180., 500.);

   w_A_eeb           = fs->make < TH2F > ("w_A_eeb",           "w_A_eeb; A", 10, 0., 1., 10, 0., 1.);
   w_A_mmb           = fs->make < TH2F > ("w_A_mmb",           "w_A_mmb; A", 10, 0., 1., 10, 0., 1.);
}

ZbDumper::~ZbDumper() {

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void ZbDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace std;

    edm::Handle <std::vector<math::XYZTLorentzVector>> electrons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> muons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> jets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> bjets;
    edm::Handle <std::vector<double>>   ptZ;
    edm::Handle <std::vector<double>>   ptZ_b;
    edm::Handle <std::vector<double>>   yZ;
    edm::Handle <std::vector<double>>   yZ_b;
    edm::Handle <std::vector<double>>   zj_mass;
    edm::Handle <std::vector<double>>   zb_mass;
    edm::Handle <std::vector<double>>   Ht;
    edm::Handle <std::vector<double>>   Ht_b;
    edm::Handle <std::vector<double>>   delta_phi;
    edm::Handle <std::vector<double>>   bdelta_phi;
    edm::Handle <std::vector<double>>   weight;
    edm::Handle <std::vector<double>>   bweight;
    edm::Handle <std::vector<double>>   delta_phi_bb;
    edm::Handle <std::vector<double>>   DR_bb;
    edm::Handle <std::vector<double>>   DR_Zb_min;
    edm::Handle <std::vector<double>>   DR_Zb_max;
    edm::Handle <std::vector<double>>   A_Zb;
    edm::Handle <std::vector<double>>   phi_star;
    edm::Handle <std::vector<double>>   phi_star_b;
    edm::Handle <std::vector<double>>   bb_Mass;
    edm::Handle <std::vector<double>>   bbZ_Mass;    
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_electrons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_muons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_jets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_jets2;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_bjets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_bjets2;
    edm::Handle <std::vector<double>>   gen_ptZ;
    edm::Handle <std::vector<double>>   gen_ptZ_b;
    edm::Handle <std::vector<double>>   gen_yZ;
    edm::Handle <std::vector<double>>   gen_yZ_b;
    edm::Handle <std::vector<double>>   gen_zj_mass;
    edm::Handle <std::vector<double>>   gen_zb_mass;
    edm::Handle <std::vector<double>>   gen_Ht;
    edm::Handle <std::vector<double>>   gen_Ht_b;
    edm::Handle <std::vector<double>>   gen_delta_phi;
    edm::Handle <std::vector<double>>   gen_bdelta_phi;
    edm::Handle <std::vector<double>>   gen_weight;
    edm::Handle <std::vector<double>>   gen_delta_phi_bb;
    edm::Handle <std::vector<double>>   gen_DR_bb;
    edm::Handle <std::vector<double>>   gen_DR_Zb_min;
    edm::Handle <std::vector<double>>   gen_DR_Zb_max;
    edm::Handle <std::vector<double>>   gen_A_Zb;
    edm::Handle <std::vector<double>>   gen_phi_star;
    edm::Handle <std::vector<double>>   gen_phi_star_b;
    edm::Handle <std::vector<double>>   gen_bb_Mass;
    edm::Handle <std::vector<double>>   gen_bbZ_Mass;

    string postfix = "";
    string postfixgen = "";
    if (pileupDT_=="ee_pup" || pileupDT_=="mm_pup") postfix = "Pup";
    if (pileupDT_=="ee_pum" || pileupDT_=="mm_pum") postfix = "Pum";
    if (par_==+1) postfix = "Up";
    if (par_==-1) postfix = "Down";
    if (par2_==+1) postfix = "JerUp";
    if (par2_==-1) postfix = "JerDown";
    if (pcut_) postfix = "Pur";
    if (useDeltaR_) postfix = "DR";
    if (numB_==1) {
      postfix = "1b";
      postfixgen = "1b";
    }
    if (numB_==2) {
      postfix = "2b";
      postfixgen = "2b";
    }

   if (lepton_== "electron") {

     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myEventWeight"), weight);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myJets"), jets);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myPtZ"), ptZ);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myPtZb"), ptZ_b);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myYZ"), yZ);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myYZb"), yZ_b);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myMassZj"), zj_mass);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myMassZb"), zb_mass);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myHtb"), Ht_b);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myDeltaPhi"), delta_phi);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myBDeltaPhi"), bdelta_phi);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myBJets"), bjets);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myBJetsWeights"), bweight);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myDeltaPhibb"), delta_phi_bb);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myDRbb"), DR_bb);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myDRZbMin"), DR_Zb_min);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myDRZbMax"), DR_Zb_max);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myAZb"), A_Zb);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myPhiStar"), phi_star);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myPhiStarb"), phi_star_b);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"mybbMass"), bb_Mass);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"mybbZMass"), bbZ_Mass);


     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myEventWeight"), gen_weight);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myElectrons"), gen_electrons);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myMuons"), gen_muons);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myJets"), gen_jets);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myJets2"), gen_jets2);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myPtZ"), gen_ptZ);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myPtZb"), gen_ptZ_b);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myYZ"), gen_yZ);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myYZb"), gen_yZ_b);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myMassZj"), gen_zj_mass);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myMassZb"), gen_zb_mass);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myHt"), gen_Ht);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myHtb"), gen_Ht_b);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myDeltaPhi"), gen_delta_phi);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myBDeltaPhi"), gen_bdelta_phi);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myBJets"), gen_bjets);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myBJets2"), gen_bjets2);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myDeltaPhibb"), gen_delta_phi_bb);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myDRbb"), gen_DR_bb);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myDRZbMin"), gen_DR_Zb_min);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myDRZbMax"), gen_DR_Zb_max);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myAZb"), gen_A_Zb);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myPhiStar"), gen_phi_star);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"myPhiStarb"), gen_phi_star_b);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"mybbMass"), gen_bb_Mass);
     iEvent.getByLabel (edm::InputTag("demoEleGen"+postfixgen,"mybbZMass"), gen_bbZ_Mass);
   }

   if (lepton_== "muon") {

     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myEventWeight"), weight);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myJets"), jets);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myPtZ"), ptZ);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myPtZb"), ptZ_b);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myYZ"), yZ);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myYZb"), yZ_b);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myMassZj"), zj_mass);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myMassZb"), zb_mass);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myHtb"), Ht_b);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myDeltaPhi"), delta_phi);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myBDeltaPhi"), bdelta_phi);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myBJets"), bjets);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myBJetsWeights"), bweight);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myDeltaPhibb"), delta_phi_bb);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myDRbb"), DR_bb);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myDRZbMin"), DR_Zb_min);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myDRZbMax"), DR_Zb_max);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myAZb"), A_Zb);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myPhiStar"), phi_star);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myPhiStarb"), phi_star_b);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"mybbMass"), bb_Mass);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"mybbZMass"), bbZ_Mass);


     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myEventWeight"), gen_weight);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myElectrons"), gen_electrons);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myMuons"), gen_muons);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myJets"), gen_jets);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myJets2"), gen_jets2);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myPtZ"), gen_ptZ);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myPtZb"), gen_ptZ_b);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myYZ"), gen_yZ);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myYZb"), gen_yZ_b);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myMassZj"), gen_zj_mass);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myMassZb"), gen_zb_mass);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myHt"), gen_Ht);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myHtb"), gen_Ht_b);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myDeltaPhi"), gen_delta_phi);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myBDeltaPhi"), gen_bdelta_phi);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myBJets"), gen_bjets);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myBJets2"), gen_bjets2);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myDeltaPhibb"), gen_delta_phi_bb);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myDRbb"), gen_DR_bb);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myDRZbMin"), gen_DR_Zb_min);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myDRZbMax"), gen_DR_Zb_max);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myAZb"), gen_A_Zb);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myPhiStar"), gen_phi_star);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"myPhiStarb"), gen_phi_star_b);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"mybbMass"), gen_bb_Mass);
     iEvent.getByLabel (edm::InputTag("demoMuoGen"+postfixgen,"mybbZMass"), gen_bbZ_Mass);
   }

   int k=-1;
   if (jets->size()>0) {
     double R = 0.5;
     for (unsigned int i=0; i<gen_jets2->size(); ++i) {
       if (ROOT::Math::VectorUtil::DeltaR((*jets)[0], (*gen_jets2)[i]) < R) {
         k=i;
         R = ROOT::Math::VectorUtil::DeltaR((*jets)[0], (*gen_jets2)[i]);
       }
     }
   }
   int k_b=-1;
   if (bjets->size()>0) {
     double R_b = 0.5;
     for (unsigned int i=0; i<gen_bjets2->size(); ++i) {
       if (ROOT::Math::VectorUtil::DeltaR((*bjets)[0], (*gen_bjets2)[i]) < R_b) {
         k_b=i;
         R_b = ROOT::Math::VectorUtil::DeltaR((*bjets)[0], (*gen_bjets2)[i]);
       }
     }
   }
   int k_b2=-1;
   if (bjets->size()>1) {
     double R_b = 0.5;
     for (unsigned int i=0; i<gen_bjets2->size(); ++i) {
       int j = i;
       if (ROOT::Math::VectorUtil::DeltaR((*bjets)[1], (*gen_bjets2)[i]) < R_b && j!=k_b) {
         k_b2=i;
         R_b = ROOT::Math::VectorUtil::DeltaR((*bjets)[1], (*gen_bjets2)[i]);
       }
     }
   }

   if (k!=-1) {
     if ((*gen_jets2)[k].pt() < 30) k=-1;
     if (fabs((*gen_jets2)[k].eta()) > 2.4) k=-1;
   }
   if (k_b!=-1) {
     if ((*gen_bjets2)[k_b].pt() < 30) k_b=-1;
     if (fabs((*gen_bjets2)[k_b].eta()) > 2.4) k_b=-1;
   }
   if (k_b2!=-1) {
     if ((*gen_bjets2)[k_b2].pt() < 30) k_b2=-1;
     if (fabs((*gen_bjets2)[k_b2].eta()) > 2.4) k_b2=-1;
   }

   double my_weight = weight->empty() ? ( gen_weight->empty() ? -1 : (*gen_weight)[0] ) : (*weight)[0];
   double my_bweight = my_weight * ( bweight->empty() ? 1 : (*bweight)[0] );

   if (my_weight>0) {

     bool ee_event = lepton_ == "electron" && (electrons->size() != 0 || gen_electrons->size() != 0);
     bool mm_event = lepton_ == "muon" && (muons->size() != 0 || gen_muons->size() != 0);

     if (ee_event || mm_event) {
       w_first_jet_pt->Fill(jets->empty() ? -1 : (*jets)[0].pt(), k<0 ? -1 : (*gen_jets2)[k].pt(), my_weight);
       w_first_jet_eta->Fill(jets->empty() ? -3 : (*jets)[0].eta(), k<0 ? -3 : (*gen_jets2)[k].eta(), my_weight);
       w_first_jet_eta_abs->Fill(jets->empty() ? -3 : fabs((*jets)[0].eta()), k<0 ? -3 : fabs((*gen_jets2)[k].eta()), my_weight);
       w_first_bjet_pt->Fill(bjets->empty() ? -1 : (*bjets)[0].pt(), k_b<0 ? -1 : (*gen_bjets2)[k_b].pt(), my_bweight);
       w_first_bjet_eta->Fill(bjets->empty() ? -3 : (*bjets)[0].eta(), k_b<0 ? -3 : (*gen_bjets2)[k_b].eta(), my_bweight);
       w_first_bjet_eta_abs->Fill(bjets->empty() ? -3 : fabs((*bjets)[0].eta()), k_b<0 ? -3 : fabs((*gen_bjets2)[k_b].eta()), my_bweight);
       if (numB_==2) {
         w_second_bjet_pt->Fill(bjets->empty() ? -1 : (*bjets)[1].pt(), k_b2<0 ? -1 : (*gen_bjets2)[k_b2].pt(), my_bweight);
         w_second_bjet_eta->Fill(bjets->empty() ? -3 : (*bjets)[1].eta(), k_b2<0 ? -3 : (*gen_bjets2)[k_b2].eta(), my_bweight);
         w_second_bjet_eta_abs->Fill(bjets->empty() ? -3 : fabs((*bjets)[1].eta()), k_b2<0 ? -3 : fabs((*gen_bjets2)[k_b2].eta()), my_bweight);
       }
     }

     if (ee_event) {
       w_pt_Z_ee->Fill(ptZ->empty() ? -1 : (*ptZ)[0], gen_ptZ->empty() ? -1 : (*gen_ptZ)[0], my_weight);
       w_pt_Z_ee_b->Fill(ptZ_b->empty() ? -1 : (*ptZ_b)[0], gen_ptZ_b->empty() ? -1 : (*gen_ptZ_b)[0], my_bweight);
     }

     if (mm_event) {
       w_pt_Z_mm->Fill(ptZ->empty() ? -1 : (*ptZ)[0], gen_ptZ->empty() ? -1 : (*gen_ptZ)[0], my_weight);
       w_pt_Z_mm_b->Fill(ptZ_b->empty() ? -1 : (*ptZ_b)[0], gen_ptZ_b->empty() ? -1 : (*gen_ptZ_b)[0], my_bweight);
     }
     
     if (ee_event) {
       w_y_Z_ee->Fill(yZ->empty() ? -3 : (*yZ)[0], gen_yZ->empty() ? -3 : (*gen_yZ)[0], my_weight);
       w_y_Z_ee_b->Fill(yZ_b->empty() ? -3 : (*yZ_b)[0], gen_yZ_b->empty() ? -3 : (*gen_yZ_b)[0], my_bweight);
       w_y_Z_ee_abs->Fill(yZ->empty() ? -3 : fabs((*yZ)[0]), gen_yZ->empty() ? -3 : fabs((*gen_yZ)[0]), my_weight);
       w_y_Z_ee_b_abs->Fill(yZ_b->empty() ? -3 : fabs((*yZ_b)[0]), gen_yZ_b->empty() ? -3 : fabs((*gen_yZ_b)[0]), my_bweight);
     }

     if (mm_event) {
       w_y_Z_mm->Fill(yZ->empty() ? -3 : (*yZ)[0], gen_yZ->empty() ? -3 : (*gen_yZ)[0], my_weight);
       w_y_Z_mm_b->Fill(yZ_b->empty() ? -3 : (*yZ_b)[0], gen_yZ_b->empty() ? -3 : (*gen_yZ_b)[0], my_bweight);
       w_y_Z_mm_abs->Fill(yZ->empty() ? -3 : fabs((*yZ)[0]), gen_yZ->empty() ? -3 : fabs((*gen_yZ)[0]), my_weight);
       w_y_Z_mm_b_abs->Fill(yZ_b->empty() ? -3 : fabs((*yZ_b)[0]), gen_yZ_b->empty() ? -3 : fabs((*gen_yZ_b)[0]), my_bweight);
     }

     if ((ee_event || mm_event) && numB_!=1) {
       w_delta_phi_2b->Fill(delta_phi_bb->empty() ? -1 : (*delta_phi_bb)[0], gen_delta_phi_bb->empty() ? -1 : (*gen_delta_phi_bb)[0], my_bweight);       
     }
     
     if ((ee_event || mm_event) && numB_==2) {
       w_DR_bb->Fill(DR_bb->empty() ? -1 : (*DR_bb)[0], gen_DR_bb->empty() ? -1 : (*gen_DR_bb)[0], my_bweight);
     }
 
     if (ee_event && numB_!=1) {
       w_DR_eeb_min->Fill(DR_Zb_min->empty() ? -1 : (*DR_Zb_min)[0], gen_DR_Zb_min->empty() ? -1 : (*gen_DR_Zb_min)[0], my_bweight);
       w_DR_eeb_max->Fill(DR_Zb_max->empty() ? -1 : (*DR_Zb_max)[0], gen_DR_Zb_max->empty() ? -1 : (*gen_DR_Zb_max)[0], my_bweight);
     }

     if (mm_event && numB_!=1) {
       w_DR_mmb_min->Fill(DR_Zb_min->empty() ? -1 : (*DR_Zb_min)[0], gen_DR_Zb_min->empty() ? -1 : (*gen_DR_Zb_min)[0], my_bweight);
       w_DR_mmb_max->Fill(DR_Zb_max->empty() ? -1 : (*DR_Zb_max)[0], gen_DR_Zb_max->empty() ? -1 : (*gen_DR_Zb_max)[0], my_bweight);
     }

     if (ee_event && numB_!=1) {
       w_A_eeb->Fill(A_Zb->empty() ? -1 : (*A_Zb)[0], gen_A_Zb->empty() ? -1 : (*gen_A_Zb)[0], my_bweight);
     }

     if (mm_event && numB_!=1) {
       w_A_mmb->Fill(A_Zb->empty() ? -1 : (*A_Zb)[0], gen_A_Zb->empty() ? -1 : (*gen_A_Zb)[0], my_bweight);
     }

     if (ee_event) {
       w_Phi_star_ee->Fill(phi_star->empty() ? -1 : (*phi_star)[0], gen_phi_star->empty() ? -1 : (*gen_phi_star)[0], my_weight);
       w_Phi_star_ee_b->Fill(phi_star_b->empty() ? -1 : (*phi_star_b)[0], gen_phi_star_b->empty() ? -1 : (*gen_phi_star_b)[0], my_weight);
     }

     if (mm_event) {
       w_Phi_star_mm->Fill(phi_star->empty() ? -1 : (*phi_star)[0], gen_phi_star->empty() ? -1 : (*gen_phi_star)[0], my_weight);
       w_Phi_star_mm_b->Fill(phi_star_b->empty() ? -1 : (*phi_star_b)[0], gen_phi_star_b->empty() ? -1 : (*gen_phi_star_b)[0], my_weight);
     }

     if (ee_event) {
       w_mass_Zj_ee->Fill(zj_mass->empty() ? -1 : (*zj_mass)[0], gen_zj_mass->empty() ? -1 : (*gen_zj_mass)[0], my_weight);
       w_mass_Zj_ee_b->Fill(zb_mass->empty() ? -1 : (*zb_mass)[0], gen_zb_mass->empty() ? -1 : (*gen_zb_mass)[0], my_bweight);
     }

     if (mm_event) {
       w_mass_Zj_mm->Fill(zj_mass->empty() ? -1 : (*zj_mass)[0], gen_zj_mass->empty() ? -1 : (*gen_zj_mass)[0], my_weight);
       w_mass_Zj_mm_b->Fill(zb_mass->empty() ? -1 : (*zb_mass)[0], gen_zb_mass->empty() ? -1 : (*gen_zb_mass)[0], my_bweight);
     }

     if ((ee_event || mm_event) && numB_!=1) {
       w_bb_mass->Fill(bb_Mass->empty() ? -1 : (*bb_Mass)[0], gen_bb_Mass->empty() ? -1 : (*gen_bb_Mass)[0], my_bweight);
     }

     if (ee_event && numB_!=1) {
       w_eebb_mass->Fill(bbZ_Mass->empty() ? -1 : (*bbZ_Mass)[0], gen_bbZ_Mass->empty() ? -1 : (*gen_bbZ_Mass)[0], my_bweight);
     }

     if (mm_event && numB_!=1) {
       w_mmbb_mass->Fill(bbZ_Mass->empty() ? -1 : (*bbZ_Mass)[0], gen_bbZ_Mass->empty() ? -1 : (*gen_bbZ_Mass)[0], my_bweight);
     }

     if (ee_event || mm_event) {
       w_Ht->Fill(Ht->empty() ? -1 : (*Ht)[0], gen_Ht->empty() ? -1 : (*gen_Ht)[0], my_weight);
       w_Ht_b->Fill(Ht_b->empty() ? -1 : (*Ht_b)[0], gen_Ht_b->empty() ? -1 : (*gen_Ht_b)[0], my_bweight);
     }

//     if (ee_event) {
//       w_delta_ee->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi->empty() ? -1 : (*gen_delta_phi)[0], my_weight);
//       w_delta_ee_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi->empty() ? -1 : (*gen_bdelta_phi)[0], my_bweight);
//     }

     if (ee_event) {
       math::XYZTLorentzVector gen_z;
       if (gen_electrons->size() != 0) gen_z = (*gen_electrons)[0] + (*gen_electrons)[1];
       double gen_delta_phi = (gen_electrons->empty() || k<0) ? -1 : fabs(gen_z.phi() - (*gen_jets2)[k].phi());
       double gen_bdelta_phi = (gen_electrons->empty() || k_b<0) ? -1 : fabs(gen_z.phi() - (*gen_bjets2)[k_b].phi());
       if (gen_delta_phi > acos (-1)) gen_delta_phi = 2 * acos (-1) - gen_delta_phi;
       if (gen_bdelta_phi > acos (-1)) gen_bdelta_phi = 2 * acos (-1) - gen_bdelta_phi;
       w_delta_ee->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi, my_weight);
       w_delta_ee_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi, my_bweight);
     }

//     if (mm_event) {
//       w_delta_mm->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi->empty() ? -1 : (*gen_delta_phi)[0], my_weight);
//       w_delta_mm_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi->empty() ? -1 : (*gen_bdelta_phi)[0], my_bweight);
//     }

     if (mm_event) {
       math::XYZTLorentzVector gen_z;
       if (gen_muons->size() != 0) gen_z = (*gen_muons)[0] + (*gen_muons)[1];
       double gen_delta_phi = (gen_muons->empty() || k<0) ? -1 : fabs(gen_z.phi() - (*gen_jets2)[k].phi());
       double gen_bdelta_phi = (gen_muons->empty() || k_b<0) ? -1 : fabs(gen_z.phi() - (*gen_bjets2)[k_b].phi());
       if (gen_delta_phi > acos (-1)) gen_delta_phi = 2 * acos (-1) - gen_delta_phi;
       if (gen_bdelta_phi > acos (-1)) gen_bdelta_phi = 2 * acos (-1) - gen_bdelta_phi;
       w_delta_mm->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi, my_weight);
       w_delta_mm_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi, my_bweight);
     }

   }
}


// ------------ method called once each job just before starting event loop  ------------
void ZbDumper::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void ZbDumper::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void ZbDumper::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void ZbDumper::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void ZbDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void ZbDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZbDumper);
