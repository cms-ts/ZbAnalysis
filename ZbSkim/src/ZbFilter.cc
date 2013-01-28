// -*- C++ -*-
//
// Package:    ZbFilter
// Class:      ZbFilter
// 
/**\class ZbFilter ZbFilter.cc Zbanalysis/ZbFilter/src/ZbFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vieri Candelise
//         Created:  Thu Nov  1 11:32:14 CET 2012
// $Id$
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
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TLorentzVector.h"

//
// class declaration
//



class ZbFilter : public edm::EDFilter {
   public:
      explicit ZbFilter(const edm::ParameterSet&);
      ~ZbFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

//
// constants, enums and typedefs
//

//
// static data member definitions
//

  edm::Service<TFileService> fs; 

//
// constructors and destructor
//
};

ZbFilter::ZbFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


ZbFilter::~ZbFilter()
{

 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZbFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

int Nj=0;

   //get electron collection
   //edm::Handle<std::vector<reco::PFCandidate> > electrons;
   //iEvent.getByLabel("pfIsolatedElectronsPFlow", electrons);
   
   //get muon collection
   //edm::Handle<std::vector<reco::PFCandidate>  > muons;
   //iEvent.getByLabel("pfIsolatedMuonsPFlow", muons);
    
   //get jet collection
   edm::Handle<std::vector<pat::Jet>  > jets;
   iEvent.getByLabel("goodJets",jets);

   //get electron collection
   edm::Handle<pat::MuonCollection> Trigmuons;
   iEvent.getByLabel("matchedMuons", Trigmuons);

   //get muon collection
   edm::Handle<pat::ElectronCollection> Trigelectrons;
   iEvent.getByLabel("matchedElectrons", Trigelectrons);

   // Get the Z->mm collection
   edm::Handle<reco::CompositeCandidateCollection> zmm;
   iEvent.getByLabel("zmuMatchedmuMatched", zmm);

   // Get the Z->ee collection
   edm::Handle<reco::CompositeCandidateCollection> zee;
   iEvent.getByLabel("zeleMatchedeleMatched", zee);

   //std::cout<<"numero j="<<jets->size()<<std::endl;

   if ( (zmm->size()==0 && zee->size()==0) || jets->size()==0) return false;

   return true;

}
// ------------ method called once each job just before starting event loop  ------------
void 
ZbFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZbFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
ZbFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ZbFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ZbFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ZbFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZbFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZbFilter);
