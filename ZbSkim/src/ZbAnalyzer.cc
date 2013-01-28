// -*- C++ -*-
//
// Package:    ZbAnalyzer
// Class:      ZbAnalyzer
// 
/**\class ZbAnalyzer ZbAnalyzer.cc ZbAnalysis/ZbAnalyzer/src/ZbAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vieri Candelise
//         Created:  Thu Jan 10 15:57:03 CET 2013
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

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// system include files
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stddef.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

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

class TTree;
//
// class declaration
//

class GreaterPt{
public:
  bool operator()( const math::XYZTLorentzVector& a, const math::XYZTLorentzVector& b) {
    return a.Pt() > b.Pt();
  }
};



class ZbAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ZbAnalyzer(const edm::ParameterSet&);
      ~ZbAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
   
      int Nj;
      double jet_pt;
      double ele_pt;
      double muon_pt;
      double diele_inv;
      double dimuon_inv;

      int NZ, NZ2;
      double discrCSV;
      double MyWeight;

      TTree* treeZb_;

      TH1F* h_jetmultiplicity;
      TH1F* h_jet_pt;
      TH1F* h_ele_pt;
      TH1F* h_muon_pt;
      TH1F* h_mm_inv;
      TH1F* h_ee_inv; 
      TH1F* h_secondvtx_N;
      TH1F* h_PUweights;
      TH1D* recoVTX_;
      TH1D* recoVTX_w;

};
	using namespace pat;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZbAnalyzer::ZbAnalyzer(const edm::ParameterSet& iConfig)
{

   //now do what ever initialization is needed
   edm::Service<TFileService> fs; 

   h_jetmultiplicity = fs->make<TH1F>("h_jetmultiplicity","h_jetmultiplicity",6,0,6);
   h_jet_pt = fs->make<TH1F>("h_jet_pt","h_jet_pt",100,0,100); 
   h_ele_pt = fs->make<TH1F>("h_ele_pt","h_ele_pt",100,0,100); 
   h_muon_pt = fs->make<TH1F>("h_muon_pt","h_muon_pt",100,0,100); 
   h_mm_inv = fs->make<TH1F>("h_mm_inv","h_mm_inv",60,60,120); 
   h_ee_inv = fs->make<TH1F>("h_ee_inv","h_ee_inv",60,60,120); 
   h_secondvtx_N = fs->make<TH1F>("h_secondvtx_N","h_secondvtx_N",10,0,2); 
   h_PUweights   = fs->make<TH1F>("h_pu_weights", "h_pu_weights", 10, 0, 10);
   recoVTX_ =        fs->make<TH1D>("recoVTX","No. reconstructed vertices",40,0.,40.);
   recoVTX_w =        fs->make<TH1D>("recoVTXw","No. reconstructed vertices weighted",40,0.,40.);

   treeZb_ = fs->make<TTree>("ZbTree","ZbTree");
   treeZb_ ->Branch("Nj", &Nj);
   treeZb_ ->Branch("jet_pt", &jet_pt);
   treeZb_ ->Branch("muon_pt", &muon_pt);
   treeZb_ ->Branch("dimuon_inv", &dimuon_inv);
   treeZb_ ->Branch("diele_inv", &diele_inv);
   treeZb_ ->Branch("weights_pu", &MyWeight);

}

ZbAnalyzer::~ZbAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
ZbAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   //PU

   MyWeight = 1.0;
   Handle<std::vector< PileupSummaryInfo > >  PupInfo;

   if (iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo)) {

	   std::vector<PileupSummaryInfo>::const_iterator PVI;

	   float Tnpv = -1;
	   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
		   int BX = PVI->getBunchCrossing();
		   if(BX == 0) { 
			   Tnpv = PVI->getTrueNumInteractions();
			   continue;
		   }
	   }

	   edm::LumiReWeighting LumiWeights_;
	   LumiWeights_ = edm::LumiReWeighting("/gpfs/cms/users/candelis/work/Zb/pileup/pileup_S7.root", 
			   		       "/gpfs/cms/users/candelis/work/Zb/pileup/pileup_2012.root", "pileup", "pileup");
   
	   MyWeight = LumiWeights_.weight( Tnpv );
   }

   cout<<"weight="<<MyWeight<<endl;

   //get vertices 
    edm::Handle< std::vector<reco::Vertex> > vertices_h;
    iEvent.getByLabel(edm::InputTag ("offlinePrimaryVertices"), vertices_h);
    if (!vertices_h.isValid()) {
    	    //std::cout<<"empty vertex collection!!!\n";
    	    //return;
    }
    // require in the event that there is at least one reconstructed vertex
    if(vertices_h->size()<=0) return;
    // pick the first (i.e. highest sum pt) vertex
    const reco::Vertex* theVertex=&(vertices_h->front());
    // require that the vertex meets certain criteria
    if(theVertex->ndof()<5) return;
    if(fabs(theVertex->z())>24.0) return;
    if(fabs(theVertex->position().rho())>2.0) return;
   
     std::vector<reco::Vertex>::const_iterator itv;
     int NVtx = 0;
     // now, count vertices
     for (itv = vertices_h->begin(); itv != vertices_h->end(); ++itv) {
	     // require that the vertex meets certain criteria
	     if(itv->ndof()<5) continue;
	     if(fabs(itv->z())>50.0) continue;
	     if(fabs(itv->position().rho())>2.0) continue;
	     ++NVtx;
     }
     //recoVTX_->Fill(float(NVtx));
     //recoVTX_w ->Fill(NVtx, MyWeight);

   //get electron collection
   //edm::Handle<std::vector<reco::PFCandidate> > electrons;
   //iEvent.getByLabel("pfIsolatedElectronsPFlow", electrons);

   //get muon collection
   //edm::Handle<std::vector<reco::PFCandidate>  > muons;
   //iEvent.getByLabel("pfIsolatedMuonsPFlow", muons);
   
   //get muon collection
   edm::Handle<pat::MuonCollection> Trigmuons;
   iEvent.getByLabel("matchedMuons", Trigmuons);
   
   //get muon collection
   edm::Handle<pat::ElectronCollection> Trigelectrons;
   iEvent.getByLabel("matchedElectrons", Trigelectrons);

   //get jet collection
   //edm::Handle<std::vector<reco::PFJet>  > jets;
   //iEvent.getByLabel("ak5PFJets",jets);

   //get jet collection
   edm::Handle<std::vector<pat::Jet>  > jets;
   iEvent.getByLabel("goodJets",jets);
   
   // Get the Z->mm collection
   edm::Handle<reco::CompositeCandidateCollection> zmm;
   iEvent.getByLabel("zmuMatchedmuMatched", zmm);
   
   // Get the Z->ee collection
   edm::Handle<reco::CompositeCandidateCollection> zee;
   iEvent.getByLabel("zeleMatchedeleMatched", zee);
  
   Nj=0;
   jet_pt=0;
   ele_pt=0;
   muon_pt=0;
   diele_inv=0;
   dimuon_inv=0;
   NZ=0;
   NZ2=0;
   discrCSV=0; 

    //test on muon Vs muonmatched///////////////////////////////////

	           //cout<<"mu size="<<muons->size()<<endl;
	           cout<<"Trig ele size="<<Trigelectrons->size()<<endl;
	           cout<<"Trig muo size="<<Trigmuons->size()<<endl;
  	           //cout<<"Z="<<Zbosons->size()<<endl; 
   		   //cout<<"Z="<<zmm->size()<<endl;

   ///////////////////////////////////////////////////////
   
   // Loop over pat jets
   for(std::vector<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
	    if(jet->pt() > 30 && abs(jet->eta()) < 2.5
	       && jet->chargedEmEnergyFraction()<0.99
	       && jet->neutralHadronEnergyFraction()<0.99
	       && jet->neutralEmEnergyFraction()<0.99
	       && jet->chargedHadronEnergyFraction()>0
	       && jet->chargedMultiplicity()>0){	
	
		++Nj;
	        jet_pt = jet->pt();
	        h_jet_pt -> Fill(jet_pt);

/*
const std::vector< std::pair< std::string, float > > discrPairs = jet->getPairDiscri();
for(uint pair_i = 0;pair_i <discrPairs.size(); ++pair_i ){
cout << "discr name: "<< discrPairs[pair_i].first<<" value: "<< discrPairs[pair_i].second<<endl;
      }
*/

                discrCSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
		cout<<"CSV discrim = "<< discrCSV <<endl;
		
		//const reco::SecondaryVertexTagInfo &svTagInfo = *jet->tagInfoSecondaryVertex();
		//if (svTagInfo.nVertices() < 1) continue;
	        //const reco::Vertex &sv = svTagInfo.secondaryVertex(0);

                 
    	    }
   } 

   // Loop over triggered matched muons

   for(pat::MuonCollection::const_iterator mu=Trigmuons->begin(); mu!=Trigmuons->end(); ++mu){
	   if(mu->pt()>25 && abs(mu->eta())<2.5){
	      muon_pt = mu->pt();
              if(muon_pt!=0) h_muon_pt->Fill(muon_pt);
	   }
   }

   // Loop over triggered matched muons

   for(pat::ElectronCollection::const_iterator el=Trigelectrons->begin(); el!=Trigelectrons->end(); ++el){
           if(el->pt()>25 && abs(el->eta())<2.5){
              ele_pt = el->pt();
	      if(ele_pt!=0) h_ele_pt->Fill(ele_pt);
           }
   }

   // Loop over Z bosons candidates

      for (std::vector<reco::CompositeCandidate>::const_iterator Zit = zmm->begin() ; Zit != zmm->end(); ++Zit){

	     NZ++;
  	     cout<<"mass mm="<<Zit->mass()<<endl;
	     dimuon_inv = Zit->mass();
	     if(dimuon_inv!=0) h_mm_inv->Fill(dimuon_inv);
      }
      for (std::vector<reco::CompositeCandidate>::const_iterator Zit2 = zee->begin() ; Zit2 != zee->end(); ++Zit2){

	     NZ2++;
  	     cout<<"mass ee="<<Zit2->mass()<<endl;
	     diele_inv = Zit2->mass();
	     if(diele_inv!=0) h_ee_inv->Fill(diele_inv);
      }



   // Fill the tree

	//std::cout<<Nj<<std::endl;
	treeZb_->Fill();

	h_PUweights ->Fill(MyWeight);

	h_jetmultiplicity -> Fill(Nj);
        h_secondvtx_N ->Fill(discrCSV);
	recoVTX_w->Fill(float(NVtx)-1, MyWeight);
	recoVTX_->Fill(float(NVtx)-1);

}

// ------------ method called once each job just before starting event loop  ------------
void 
ZbAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZbAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ZbAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZbAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZbAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZbAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZbAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZbAnalyzer);
