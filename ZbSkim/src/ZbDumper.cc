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

     TH2F* w_first_jet_pt;
     TH2F* w_first_jet_eta;
     TH2F* w_first_bjet_pt;
     TH2F* w_first_bjet_eta;
     TH2F* w_pt_Z_ee;
     TH2F* w_pt_Z_mm;
     TH2F* w_pt_Z_ee_b;
     TH2F* w_pt_Z_mm_b;
     TH2F* w_mass_Zj_ee;
     TH2F* w_mass_Zj_mm;
     TH2F* w_mass_Zb_ee;
     TH2F* w_mass_Zb_mm;
     TH2F* w_Ht;
     TH2F* w_Ht_b;
     TH2F* w_delta_ee;
     TH2F* w_delta_mm;
     TH2F* w_delta_ee_b;
     TH2F* w_delta_mm_b;

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

   //now do what ever initialization is needed
   lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
   edm::Service < TFileService > fs;

   w_first_jet_pt    = fs->make < TH2F > ("w_first_jet_pt",    "w_first_jet_pt;P_t [GeV]", 50, 30., 700., 50, 30., 700.);
   w_first_jet_eta   = fs->make < TH2F > ("w_first_jet_eta",   "w_first_jet_eta;Eta", 16, -2.5, 2.5,16, -2.5, 2.5);
   w_first_bjet_pt   = fs->make < TH2F > ("w_first_bjet_pt",   "w_first_bjet_pt;P_t [GeV]", 50, 30., 700., 50, 30., 700.);
   w_first_bjet_eta  = fs->make < TH2F > ("w_first_bjet_eta",  "w_first_bjet_eta;Eta", 16, -2.5, 2.5,16, -2.5, 2.5);
   w_pt_Z_ee         = fs->make < TH2F > ("w_pt_Z_ee",         "w_pt_Z_ee;P_t [GeV]", 40, 0., 400., 40, 0., 400.);
   w_pt_Z_mm         = fs->make < TH2F > ("w_pt_Z_mm",         "w_pt_Z_mm;P_t [GeV]", 40, 0., 400., 40, 0., 400.);
   w_pt_Z_ee_b 	     = fs->make < TH2F > ("w_pt_Z_ee_b",       "w_pt_Z_ee_b;P_t [GeV]", 40, 0., 400., 40, 0., 400.);
   w_pt_Z_mm_b       = fs->make < TH2F > ("w_pt_Z_mm_b",       "w_pt_Z_mm_b;P_t [GeV]", 40, 0., 400., 40, 0., 400.);

   w_mass_Zj_ee      = fs->make < TH2F > ("w_mass_Zj_ee",      "w_mass_Zj_ee;P_t [GeV]", 60, 15., 330., 60, 15., 330.);
   w_mass_Zj_mm      = fs->make < TH2F > ("w_mass_Zj_mm",      "w_mass_Zj_mm;P_t [GeV]", 60, 15., 330., 60, 15., 330.);
   w_mass_Zb_ee      = fs->make < TH2F > ("w_mass_Zb_ee",      "w_mass_Zb_ee;P_t [GeV]", 60, 15., 330., 60, 15., 330.);
   w_mass_Zb_mm      = fs->make < TH2F > ("w_mass_Zb_mm",      "w_mass_Zb_mm;P_t [GeV]", 60, 15., 330., 60, 15., 330.);

   w_Ht 	     = fs->make < TH2F > ("w_Ht",              "w_Ht [GeV]", 50, 30., 1000.,50, 30., 1000.);
   w_Ht_b 	     = fs->make < TH2F > ("w_Ht_b",            "w_Ht [GeV]", 50, 30., 1000.,50, 30., 1000.);
   w_delta_ee        = fs->make < TH2F > ("w_delta_phi_ee",    "w_delta_phi_ee", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());
   w_delta_mm        = fs->make < TH2F > ("w_delta_phi_mm",    "w_delta_phi_mm", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());
   w_delta_ee_b      = fs->make < TH2F > ("w_delta_phi_ee_b",  "w_delta_phi_ee_b", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());
   w_delta_mm_b      = fs->make < TH2F > ("w_delta_phi_mm_b",  "w_delta_phi_mm_b", 12, 0, TMath::Pi (), 12, 0, TMath::Pi ());

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
    edm::Handle <std::vector<double>>   zj_mass;
    edm::Handle <std::vector<double>>   zb_mass;
    edm::Handle <std::vector<double>>   Ht;
    edm::Handle <std::vector<double>>   Ht_b;
    edm::Handle <std::vector<double>>   delta_phi;
    edm::Handle <std::vector<double>>   bdelta_phi;
    edm::Handle <std::vector<double>>   weight;
    edm::Handle <std::vector<double>>   bweight;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_electrons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_muons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_jets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_jets2;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_bjets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_bjets2;
    edm::Handle <std::vector<double>>   gen_ptZ;
    edm::Handle <std::vector<double>>   gen_ptZ_b;
    edm::Handle <std::vector<double>>   gen_zj_mass;
    edm::Handle <std::vector<double>>   gen_zb_mass;
    edm::Handle <std::vector<double>>   gen_Ht;
    edm::Handle <std::vector<double>>   gen_Ht_b;
    edm::Handle <std::vector<double>>   gen_delta_phi;
    edm::Handle <std::vector<double>>   gen_bdelta_phi;
    edm::Handle <std::vector<double>>   gen_weight;

   if (lepton_== "electron") {

     iEvent.getByLabel (edm::InputTag("demoEle","myEventWeight"), weight);
     iEvent.getByLabel (edm::InputTag("demoEle","myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoEle","myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoEle","myJets"), jets);
     iEvent.getByLabel (edm::InputTag("demoEle","myPtZ"), ptZ);
     iEvent.getByLabel (edm::InputTag("demoEle","myPtZb"), ptZ_b);
     iEvent.getByLabel (edm::InputTag("demoEle","myMassZj"), zj_mass);
     iEvent.getByLabel (edm::InputTag("demoEle","myMassZb"), zb_mass);
     iEvent.getByLabel (edm::InputTag("demoEle","myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoEle","myHtb"), Ht_b);
     iEvent.getByLabel (edm::InputTag("demoEle","myDeltaPhi"), delta_phi);
     iEvent.getByLabel (edm::InputTag("demoEle","myBDeltaPhi"), bdelta_phi);
     iEvent.getByLabel (edm::InputTag("demoEle","myBJets"), bjets);
     iEvent.getByLabel (edm::InputTag("demoEle","myBJetsWeights"), bweight);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myEventWeight"), gen_weight);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myElectrons"), gen_electrons);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myMuons"), gen_muons);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myJets"), gen_jets);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myJets2"), gen_jets2);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myPtZ"), gen_ptZ);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myPtZb"), gen_ptZ_b);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myMassZj"), gen_zj_mass);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myMassZb"), gen_zb_mass);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myHt"), gen_Ht);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myHtb"), gen_Ht_b);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDeltaPhi"), gen_delta_phi);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myBDeltaPhi"), gen_bdelta_phi);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myBJets"), gen_bjets);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myBJets2"), gen_bjets2);

   }

   if (lepton_== "muon") {

     iEvent.getByLabel (edm::InputTag("demoMuo","myEventWeight"), weight);
     iEvent.getByLabel (edm::InputTag("demoMuo","myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoMuo","myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoMuo","myJets"), jets);
     iEvent.getByLabel (edm::InputTag("demoMuo","myPtZ"), ptZ);
     iEvent.getByLabel (edm::InputTag("demoMuo","myPtZb"), ptZ_b);
     iEvent.getByLabel (edm::InputTag("demoMuo","myMassZj"), zj_mass);
     iEvent.getByLabel (edm::InputTag("demoMuo","myMassZb"), zb_mass);
     iEvent.getByLabel (edm::InputTag("demoMuo","myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoMuo","myHtb"), Ht_b);
     iEvent.getByLabel (edm::InputTag("demoMuo","myDeltaPhi"), delta_phi);
     iEvent.getByLabel (edm::InputTag("demoMuo","myBDeltaPhi"), bdelta_phi);
     iEvent.getByLabel (edm::InputTag("demoMuo","myBJets"), bjets);
     iEvent.getByLabel (edm::InputTag("demoMuo","myBJetsWeights"), bweight);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myEventWeight"), gen_weight);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myElectrons"), gen_electrons);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myMuons"), gen_muons);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myJets"), gen_jets);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myJets2"), gen_jets2);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myPtZ"), gen_ptZ);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myPtZb"), gen_ptZ_b);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myMassZj"), gen_zj_mass);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myMassZb"), gen_zb_mass);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myHt"), gen_Ht);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myHtb"), gen_Ht_b);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDeltaPhi"), gen_delta_phi);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myBDeltaPhi"), gen_bdelta_phi);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myBJets"), gen_bjets);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myBJets2"), gen_bjets2);

   }

   int k=-1;
   if (jets->size()>0) {
     double R = 0.1;
     for (unsigned int i=0; i<gen_jets2->size(); ++i) {
       if (ROOT::Math::VectorUtil::DeltaR((*jets)[0], (*gen_jets2)[i]) < R) {
         k=i;
         R = ROOT::Math::VectorUtil::DeltaR((*jets)[0], (*gen_jets2)[i]);
       }
     }
   }
   int k_b=-1;
   if (bjets->size()>0) {
     double R_b = 0.1;
     for (unsigned int i=0; i<gen_bjets2->size(); ++i) {
       if (ROOT::Math::VectorUtil::DeltaR((*bjets)[0], (*gen_bjets2)[i]) < R_b) {
         k_b=i;
         R_b = ROOT::Math::VectorUtil::DeltaR((*bjets)[0], (*gen_bjets2)[i]);
       }
     }
   }

   if (k!=-1) {
     if ((*gen_jets2)[k].pt() < 30) k=-1;
     if (fabs((*gen_jets2)[k].eta()) > 2.5) k=-1;
   }
   if (k_b!=-1) {
     if ((*gen_bjets2)[k_b].pt() < 30) k_b=-1;
     if (fabs((*gen_bjets2)[k_b].eta()) > 2.5) k_b=-1;
   }

   double my_weight = weight->empty() ? ( gen_weight->empty() ? -1 : (*gen_weight)[0] ) : (*weight)[0];
   double my_bweight = my_weight * ( bweight->empty() ? 1 : (*bweight)[0] );

   if (my_weight>0) {

     bool ee_event = lepton_ == "electron" && (electrons->size() != 0 || gen_electrons->size() != 0);
     bool mm_event = lepton_ == "muon" && (muons->size() != 0 || gen_muons->size() != 0);

     if (ee_event || mm_event) {
       w_first_jet_pt->Fill(jets->empty() ? -1 : (*jets)[0].pt(), k<0 ? -1 : (*gen_jets2)[k].pt(), my_weight);
       w_first_jet_eta->Fill(jets->empty() ? -3 : (*jets)[0].eta(), k<0 ? -3 : (*gen_jets2)[k].eta(), my_weight);
       w_first_bjet_pt->Fill(bjets->empty() ? -1 : (*bjets)[0].pt(), k_b<0 ? -1 : (*gen_bjets2)[k_b].pt(), my_bweight);
       w_first_bjet_eta->Fill(bjets->empty() ? -3 : (*bjets)[0].eta(), k_b<0 ? -3 : (*gen_bjets2)[k_b].eta(), my_bweight);
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
       w_mass_Zj_ee->Fill(zj_mass->empty() ? -1 : (*zj_mass)[0], gen_zj_mass->empty() ? -1 : (*gen_zj_mass)[0], my_weight);
       w_mass_Zb_ee->Fill(zb_mass->empty() ? -1 : (*zb_mass)[0], gen_zb_mass->empty() ? -1 : (*gen_zb_mass)[0], my_bweight);
     }

     if (mm_event) {
       w_mass_Zj_mm->Fill(zj_mass->empty() ? -1 : (*zj_mass)[0], gen_zj_mass->empty() ? -1 : (*gen_zj_mass)[0], my_weight);
       w_mass_Zb_mm->Fill(zb_mass->empty() ? -1 : (*zb_mass)[0], gen_zb_mass->empty() ? -1 : (*gen_zb_mass)[0], my_bweight);
     }

     if (ee_event || mm_event) {
       w_Ht->Fill(Ht->empty() ? -1 : (*Ht)[0], gen_Ht->empty() ? -1 : (*gen_Ht)[0], my_weight);
       w_Ht_b->Fill(Ht_b->empty() ? -1 : (*Ht_b)[0], gen_Ht_b->empty() ? -1 : (*gen_Ht_b)[0], my_bweight);
     }

     if (ee_event) {
       w_delta_ee->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi->empty() ? -1 : (*gen_delta_phi)[0], my_weight);
       w_delta_ee_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi->empty() ? -1 : (*gen_bdelta_phi)[0], my_bweight);
     }

     if (mm_event) {
       w_delta_mm->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi->empty() ? -1 : (*gen_delta_phi)[0], my_weight);
       w_delta_mm_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi->empty() ? -1 : (*gen_bdelta_phi)[0], my_bweight);
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
