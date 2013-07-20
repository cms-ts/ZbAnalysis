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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class ZbDumper : public edm::EDAnalyzer {
   public:
      explicit ZbDumper(const edm::ParameterSet&);
      ~ZbDumper();

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

     std::string lepton_;
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
ZbDumper::ZbDumper(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");


}


ZbDumper::~ZbDumper()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZbDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

    edm::Handle <std::vector<math::XYZTLorentzVector>>  electrons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> muons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> jets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> bjets;
    edm::Handle <std::vector<double>>   ptZ;
    edm::Handle <std::vector<double>>   Ht;
    edm::Handle <std::vector<double>>   weight;
    edm::Handle <std::vector<double>>   bweight;

   if(lepton_== "electron") {

     iEvent.getByLabel (edm::InputTag("demoEle","myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoEle","myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoEle","myJets"), jets);
     iEvent.getByLabel (edm::InputTag("demoEle","myPtZ"), ptZ);
     iEvent.getByLabel (edm::InputTag("demoEle","myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoEle","myBjets"), bjets);
     iEvent.getByLabel (edm::InputTag("demoEle","myBjetsWeights"), bweight);

   } else if (lepton_== "muon") {
     
     iEvent.getByLabel (edm::InputTag("demoMuo","myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoMuo","myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoMuo","myJets"), jets);
     iEvent.getByLabel (edm::InputTag("demoMuo","myPtZ"), ptZ);
     iEvent.getByLabel (edm::InputTag("demoMuo","myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoMuo","myBjets"), bjets);
     iEvent.getByLabel (edm::InputTag("demoMuo","myBjetsWeights"), bweight);
     
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myJets"), jets);
     cout<<"gen j="<<jets->size()<<endl;
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myPtZ"), ptZ);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myBjets"), bjets);

     }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ZbDumper::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZbDumper::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ZbDumper::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZbDumper::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZbDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZbDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZbDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZbDumper);
