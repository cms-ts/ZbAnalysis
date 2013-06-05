
/************** WORK IN PROGRESS **************************/


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
// $Id: GenbAnalyzer.cc,v 1.69 2013/05/24 13:48:59 vieri Exp $
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
	    if (p1.Pt() < p2.Pt()) return false;
	    return true;
      }
    };

  std::string pileup_;
  std::string lepton_;
  edm::LumiReWeighting LumiWeights_;

  // ----------member data ---------------------------

  int       Nj, Nb;
  double    jet_pt;
  double    jet_eta;
  double    jet_phi;
  double    ele_pt;
  double    ele_pt2;
  double    ele_eta;
  double    muon_pt;
  double    muon_eta;
  double    diele_mass;
  double    diele_phi;
  double    diele_pt;
  double    dimuon_mass;
  double    dimuon_phi;
  double    dimuon_pt;
  double    delta_phi_ee_b;
  double    delta_phi_mm_b;
  double    delta_phi_ee;
  double    delta_phi_mm;
  double    deltaPhi1;
  double    deltaPhi2; 
  double    MyWeight;
  double    deltaR1;
  double    deltaR2;
  double    deltaR_eg;
  double    deltaR_mg;
  double    deltaPhi_eg;
  double    deltaPhi_mg;
  double    Eg;

  TH1F*     h_jetmultiplicity;
  TH1F*     h_jet_pt;
  TH1F*     h_ele_pt;
  TH1F*     h_muon_pt;

  TH1F*     h_pu_weights;

  TH1F*     w_jetmultiplicity;
  TH1F*     w_first_jet_pt;      // leading jet of any type
  TH1F*     w_first_jet_eta;
  TH1F*     w_second_jet_pt;
  TH1F*     w_second_jet_eta;
  TH1F*     w_third_jet_pt;
  TH1F*     w_third_jet_eta;
  TH1F*     w_first_jet_pt_b;    // leading jet with at least one b jet in the event
  TH1F*     w_first_jet_eta_b;
  TH1F*     w_second_jet_pt_b;
  TH1F*     w_second_jet_eta_b;
  TH1F*     w_third_jet_pt_b;
  TH1F*     w_third_jet_eta_b;

  TH1F*     w_bjetmultiplicity;


  TH1F*     w_first_ele_pt;
  TH1F*     w_dressed_ele_pt;

  TH1F*     w_first_muon_pt;

  TH1F*     w_mass_ee_b;  // at least one b jet in the event
  TH1F*     w_mass_mm_b;
  TH1F*     w_pt_Z_ee_b;
  TH1F*     w_pt_Z_mm_b;
  TH1F*     w_delta_ee_b;
  TH1F*     w_delta_mm_b;
  TH1F*     w_delta_mm;
 
  TH2F*     dR_ele1;
  TH2F*     dR_ele2;

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

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  w_jetmultiplicity =   fs->make < TH1F > ("w_jetmultiplicity", "w_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  w_first_ele_pt =      fs->make < TH1F > ("w_first_ele_pt",    "w_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_dressed_ele_pt =    fs->make < TH1F > ("w_dressed_ele_pt",  "w_dressed_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_first_muon_pt =     fs->make < TH1F > ("w_first_muon_pt",   "w_first_muon_pt;P_t [GeV]", 50, 0., 450.);
  dR_ele1 =             fs->make < TH2F > ("dR_ele1",   "dR_ele1", 10, 0., 1, 100, 0., 100.);
  dR_ele2 =             fs->make < TH2F > ("dR_ele2",   "dR_ele2", 10, 0., 1, 100, 0., 100.);
  
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
  iEvent.getByLabel(edm::InputTag("goodJets","genJets"), gJets); 

  bool ee_event = false;
  bool mm_event    = false;
  bool isDressed   = false;
  bool isDressedMu = false;
  
  jet_pt = 0;
  jet_eta = 0;
  jet_phi = 0;
  ele_pt = 0;
  muon_pt = 0;
  deltaR_eg = 0;
  deltaR_mg = 0;
  deltaPhi_eg=0;
  deltaPhi_mg=0;

  // +++++++++ ELECTRONS

  vector < TLorentzVector > vect_ele;
  
  for(vector<reco::GenParticle>::const_iterator itgen=genPart->begin();itgen!=genPart->end();itgen++){
  
      if (fabs(itgen->pdgId())==11 && itgen->status()==1){ // loop over gen electrons
	      TLorentzVector ele;
 	      ele.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());
	      
	      // Loop over photons: FSR dressing for electrons   
	      for(vector<reco::GenParticle>::const_iterator itgen2=genPart->begin();itgen2!=genPart->end();itgen2++){
		      if (fabs(itgen2->pdgId())==22 && itgen2->status()==1){ // loop over primary gen photon
       		          TLorentzVector gam;
		          gam.SetPtEtaPhiM(itgen2->pt(),itgen2->eta(),itgen2->phi(),itgen2->mass());
			  deltaPhi_eg = fabs(gam.Phi()- itgen->phi());	
			  if (deltaPhi_eg > acos(-1)) deltaPhi_eg= 2*acos(-1) - deltaPhi_eg; 
			  deltaR_eg = sqrt( deltaPhi_eg*deltaPhi_eg + pow(fabs(gam.Eta()-itgen->eta()),2));
			  if(deltaR_eg< 0.3) {
			      ele += gam;     
			  }
		      }
	      }
	      if(ele.Pt()>20 && fabs(ele.Eta())<2.4){
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
    if (diele_mass>71 && diele_mass<111) {
    ee_event = true;
    }
  }
	
  // +++++++++ MUONS
  
  vector < TLorentzVector > vect_muon;

  for(vector<reco::GenParticle>::const_iterator itgen=genPart->begin();itgen!=genPart->end();itgen++){
	  if (fabs(itgen->pdgId())==13 && itgen->status()==1){ // loop over gen muons
		  TLorentzVector muon;
		  muon.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());
	      
		  // Loop over photons: FSR dressing for muons   
		  for(vector<reco::GenParticle>::const_iterator itgen2=genPart->begin();itgen2!=genPart->end();itgen2++){
			  if (fabs(itgen2->pdgId())==22 && itgen2->status()==1){ // loop over primary gen photon
				  TLorentzVector gam;
				  gam.SetPtEtaPhiM(itgen2->pt(),itgen2->eta(),itgen2->phi(),itgen2->mass());
				  deltaPhi_mg = fabs(gam.Phi()- itgen->phi());
				  if (deltaPhi_mg > acos(-1)) deltaPhi_mg= 2*acos(-1) - deltaPhi_mg;
				  deltaR_mg = sqrt( deltaPhi_mg*deltaPhi_mg + pow(fabs(gam.Eta()-itgen->eta()),2));
				  if(deltaR_mg< 0.3) {
					  muon += gam; 
				  }
			  }
		  }

		  if(muon.Pt()>20 && fabs(muon.Eta())<2.4){
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
    if (dimuon_mass>71 && dimuon_mass<111){
	    mm_event = true;
    }
  }

  if (lepton_ == "electron" && !ee_event) return;
  if (lepton_ == "muon" && !mm_event)     return;

  // ++++++ Pile-Up reweighting

  bool isMC = false;
  MyWeight = 1.0;

  Handle < vector < PileupSummaryInfo > > PupInfo;

  if (iEvent.getByLabel (edm::InputTag ("addPileupInfo"), PupInfo))  {

    isMC = true;

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

  // ++++++++ Vertices analysis

  edm::Handle < vector < reco::Vertex > > vertices_h;
  iEvent.getByLabel (edm::InputTag ("goodOfflinePrimaryVertices"), vertices_h);

  // require in the event that there is at least one reconstructed vertex
  if (vertices_h->size () <= 0)    return;

  // pick the first (i.e. highest sum pt) vertex
  const reco::Vertex * theVertex = &(vertices_h->front ());

  // require that the vertex meets certain criteria
  if (theVertex->ndof () < 5)    return;
  if (fabs(theVertex->z ()) > 24.0)    return;
  if (fabs(theVertex->position ().rho ()) > 2.0)    return;

  // now, count vertices
  int NVtx = 0;
  for (vector < reco::Vertex >::const_iterator itv = vertices_h->begin (); itv != vertices_h->end (); ++itv) {
    // require that the vertex meets certain criteria
    if (itv->ndof () < 5)	continue;
    if (fabs(itv->z ()) > 50.0)	continue;
    if (fabs(itv->position ().rho ()) > 2.0)	continue;
    ++NVtx;
  }

  // ++++++++ JETS

  vector < pat::Jet > vect_jets;
  vector < pat::Jet > vect_bjets;

  bool isb = false;
  bool isc = false;
  Nj=0;
  Nb=0;
  deltaPhi1=0;
  deltaPhi2=0;
  deltaR1 = 0;
  deltaR2 = 0;

  if (ee_event || mm_event) {

      for (vector<reco::GenJet>::const_iterator jet=gJets->begin();jet!=gJets->end();++jet){
	
	  isb = false;  
	  jet_pt  = jet->pt ();
 	  jet_eta = jet->eta();
	  jet_phi = jet->phi();
 
	  /*Lepton removal from jets inside sqrt Deta2+Dphi2 <0.3*/
	  if (ee_event) deltaPhi1 = fabs(jet->phi()- vect_ele[iele0]. Phi());
	  if (mm_event) deltaPhi1 = fabs(jet->phi()- vect_muon[imuon0].Phi());

	  if (deltaPhi1 > acos(-1)) deltaPhi1= 2*acos(-1) - deltaPhi1;

	  if (ee_event) deltaR1= sqrt( deltaPhi1*deltaPhi1  + pow(jet->eta()-vect_ele[iele0]. Eta(),2) );
	  if (mm_event) deltaR1= sqrt( deltaPhi1*deltaPhi1  + pow(jet->eta()-vect_muon[imuon0].Eta(),2) );

	  if (ee_event) deltaPhi2 = fabs(jet->phi()-vect_ele[iele1]. Phi());
          if (mm_event) deltaPhi2 = fabs(jet->phi()-vect_muon[imuon1].Phi());

	  if (deltaPhi2 > acos(-1)) deltaPhi2= 2*acos(-1) - deltaPhi2;

	  if (ee_event) deltaR2= sqrt( deltaPhi2*deltaPhi2  + pow(jet->eta()-vect_ele[iele1]. Eta(),2) );	   
          if (mm_event) deltaR2= sqrt( deltaPhi2*deltaPhi2  + pow(jet->eta()-vect_muon[imuon1].Eta(),2) );
	  
          if (jet_pt > 30 && fabs(jet_eta) < 2.5 && deltaR1<0.3 && deltaR2<0.3) {
	      ++Nj;
	      vect_jets.push_back (*jet);
	  }

      }
              /*loop over gen particles, find the b*/

	      int nb=0;
	      for( std::vector <reco::GenParticle>::const_iterator thepart =genPart->begin();thepart != genPart->end(); thepart++) {
		    if( (int) (abs(thepart->pdgId() / 100)%10 ) == 5 || (int) (abs(thepart->pdgId() / 1000)%10 ) == 5 ){
		      nb++;
		      bool bdaughter = false; // b candidate has no daughters
		      for(int i=0; i < abs(thepart->numberOfDaughters()); i++){
			      if( (int) (abs(thepart->daughter(i)->pdgId() / 100)%10 ) == 5 || (int) (abs( thepart->daughter(i)->pdgId() / 1000)%10 ) == 5 ) {	
				      bdaughter=true; // b daughter found
			      }
		      }
 		          if(!bdaughter){
			        TLorentzVector B;
				B.SetPtEtaPhiM(thepart->pt(),thepart->eta(),thepart->phi(),thepart->mass());
    				    int njet=0;
				    for (vector<reco::GenJet>::const_iterator jet2=gJets->begin();jet2!=gJets->end();++jet2){	   
				    njet++;
				    double Rmin(9999.);
				       	    if(ROOT::Math::VectorUtil::DeltaR(jet2->momentum(),B)<0.5){
				       	       if(ROOT::Math::VectorUtil::DeltaR(jet2->momentum(),B)<Rmin){
   						       Rmin = ROOT::Math::VectorUtil::DeltaR(jet2->momentum(),B);
   					               isb=true;
					       }
					    }
				    	    if(isb) {
						    vect_bjets.push_back(*jet2);
						    cout<<"gen jet: "<<jet2->momentum()<<endl;
					    }
				    }
			  }
		    }
	      }
  }


  // Get reco jet collection and print the reco b jets for check

  edm::Handle < vector < pat::Jet > > jets;
  iEvent.getByLabel ("goodJets", jets);

	     for (vector < pat::Jet >::const_iterator jet = jets->begin (); jet != jets->end (); ++jet) {
     		     ++Nj;
     		     if(fabs(jet->partonFlavour ()) == 5){
	     		     cout<<"reco jet: "<<jet->momentum()<<endl;
     		     }
	     }

  // ++++++++ DIELECTRON Z PLOTS

  // ++++++++ DIMUON Z PLOTS

  // ++++++++ MISC PLOTS

  if(ee_event && Nj > 0) dR_ele1->Fill(deltaR1, vect_ele[iele0].Pt());
  if(ee_event && Nj > 0) dR_ele2->Fill(deltaR2, vect_ele[iele1].Pt());

  // ++++++++  ELECTRONS PLOTS
  if (ee_event && Nj > 0) w_first_ele_pt->Fill (vect_ele[iele0].  Pt(), MyWeight);
  if (ee_event && Nj > 0) w_dressed_ele_pt->Fill (vect_ele[iele0].Pt(), MyWeight);

  // ++++++++ MUONS PLOTS
  if (mm_event && Nj > 0) w_first_muon_pt->Fill (vect_muon[imuon0].Pt(), MyWeight);

  // ++++++++ JETS PLOTS
  if ((ee_event || mm_event) && Nj > 0) w_jetmultiplicity->Fill (Nj, MyWeight);

  // ++++++++ B JETS PLOTS

}

// ------------ method called once each job just before starting event loop ------------
void GenbAnalyzer::beginJob () {
  LumiWeights_ = edm::LumiReWeighting("/gpfs/cms/users/candelis/work/ZbSkim/test/pileup/pileup_" + pileup_ + ".root", "/gpfs/cms/users/candelis/work/ZbSkim/test/pileup/pileup_2012.root", "pileup", "pileup");

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