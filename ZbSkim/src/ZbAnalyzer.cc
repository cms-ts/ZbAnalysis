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
// $Id: ZbAnalyzer.cc,v 1.1 2013/01/28 12:52:11 dellaric Exp $
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
#include <string>

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

      std::string pileup_;
      

      // ----------member data ---------------------------
   
      int Nj, Nb;
      double jet_pt;
      double jet_phi;
      double ele_pt;
      double muon_pt;
      double diele_inv;
      double dimuon_inv;
      double diele_phi;
      double dimuon_phi;
      double Delta_phi_mm;
      double b_leading_pt;
      double  b_mm_invmass;
      double  b_ee_invmass;
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
      TH1F* h_tracks;
      TH1D* recoVTX_;
      TH1D* recoVTX_w;
      
      TH1F* w_jetmultiplicity;
      TH1F* w_bmultiplicity;
      TH1F* w_jet_pt;
      TH1F* w_ele_pt;
      TH1F* w_muon_pt;
      TH1F* w_mm_inv;
      TH1F* w_ee_inv; 
      TH1F* w_secondvtx_N;
      TH1F* w_bquarks;
      TH1F* w_tracks;
      TH1F* flavours_;
      TH1F* w_MET;
      TH1F* w_delta_phi_mm;

      TH1F* b_mm_inv;
      TH1F* b_ee_inv;
      TH1F* w_bleading_pt;
};
	using namespace pat;

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
			TH1 *nVertices;
			TH1 *deltaR, *mass, *dist, *distErr, *distSig;
			TH1 *nTracks, *chi2;
		} plots_[N_JET_TYPES];


//
// static data member definitions
//

//
// constructors and destructor
//
ZbAnalyzer::ZbAnalyzer(const edm::ParameterSet& iConfig)
{

   pileup_ = iConfig.getUntrackedParameter<std::string>("pileup", "S7");

   //now do what ever initialization is needed
   edm::Service<TFileService> fs; 

   h_jetmultiplicity = fs->make<TH1F>("h_jetmultiplicity","h_jetmultiplicity",8,0.5,8.5);
   h_jet_pt = fs->make<TH1F>("h_jet_pt","h_jet_pt",20,30,530); 
   h_ele_pt = fs->make<TH1F>("h_ele_pt","h_ele_pt",20,30,530); 
   h_muon_pt = fs->make<TH1F>("h_muon_pt","h_muon_pt",100,0,250); 
   h_mm_inv = fs->make<TH1F>("h_mm_inv","h_mm_inv",60,60,120); 
   h_ee_inv = fs->make<TH1F>("h_ee_inv","h_ee_inv",60,60,120); 
   h_secondvtx_N = fs->make<TH1F>("h_secondvtx_N","h_secondvtx_N",50,0,1); 
   h_PUweights   = fs->make<TH1F>("h_pu_weights", "h_pu_weights", 10, 0, 10);
   recoVTX_ =        fs->make<TH1D>("recoVTX","No. reconstructed vertices",40,0.,40.);
   recoVTX_w =        fs->make<TH1D>("recoVTXw","No. reconstructed vertices weighted",40,0.,40.);
   h_tracks  = fs->make<TH1F>("h_tracks","h_tracks",50,0,50);
    
   //weighted histograms
   w_jetmultiplicity = fs->make<TH1F>("w_jetmultiplicity","w_jetmultiplicity",8,0.5,8.5);
   w_bmultiplicity = fs->make<TH1F>("w_bjetmultiplicity","w_bjetmultiplicity",5,0.5,5.5);
   w_jet_pt = fs->make<TH1F>("w_jet_pt","w_jet_pt", 20,30,530);
   w_bleading_pt = fs->make<TH1F>("w_bleading_pt","w_bleading_pt",20,30,530);
   w_ele_pt = fs->make<TH1F>("w_ele_pt","w_ele_pt",100,0,250); 
   w_muon_pt = fs->make<TH1F>("w_muon_pt","w_muon_pt",100,0,250); 
   w_mm_inv = fs->make<TH1F>("w_mm_inv","w_mm_inv",60,60,120); 
   w_ee_inv = fs->make<TH1F>("w_ee_inv","w_ee_inv",60,60,120); 
   w_secondvtx_N = fs->make<TH1F>("w_secondvtx_N","w_secondvtx_N",50,0,1); 
   w_bquarks = fs->make<TH1F>("bquarks","bquarks",50,0,1);
   w_tracks  = fs->make<TH1F>("w_tracks","w_tracks",50,0,50);
   flavours_ = fs->make<TH1F>("flavours", "jet flavours", 5, 0, 5);
   w_MET  = fs->make<TH1F>("w_MET","w_MET",100,0,250);
   w_delta_phi_mm  = fs->make<TH1F>("w_delta_phi_mm","w_delta_phi_mm", 12,0,TMath::Pi());

   b_mm_inv = fs->make<TH1F>("b_mm_inv","b_mm_inv",60,60,120); 
   b_ee_inv = fs->make<TH1F>("b_ee_inv","b_ee_inv",60,60,120); 

   treeZb_ = fs->make<TTree>("ZbTree","ZbTree");
   treeZb_ ->Branch("Nj", &Nj);
   treeZb_ ->Branch("jet_pt", &jet_pt);
   treeZb_ ->Branch("muon_pt", &muon_pt);
   treeZb_ ->Branch("dimuon_inv", &dimuon_inv);
   treeZb_ ->Branch("diele_inv", &diele_inv);
   treeZb_ ->Branch("weights_pu", &MyWeight);

   for(unsigned int i = 0; i < N_JET_TYPES; i++) {
	   Plots &plots = plots_[i];
	   const char *flavour, *name;


	   switch((Flavour)i) {
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
	   plots.dist = fs->make<TH1F>(Form("dist_%s", name),
			   Form("Transverse distance between PV and SV in %s", flavour), 100, 0, 2);


}
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
	   LumiWeights_ = edm::LumiReWeighting("/gpfs/cms/users/candelis/work/Zb/pileup/pileup_" + pileup_ + ".root", 
			   		       "/gpfs/cms/users/candelis/work/Zb/pileup/pileup_2012.root", "pileup", "pileup");
   
	   MyWeight = LumiWeights_.weight( Tnpv );

   }

   //cout<<"weight="<<MyWeight<<endl;

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

   //get muon collection
   edm::Handle<pat::MuonCollection> Trigmuons;
   iEvent.getByLabel("matchedMuons", Trigmuons);
   
   //get muon collection
   edm::Handle<pat::ElectronCollection> Trigelectrons;
   iEvent.getByLabel("matchedElectrons", Trigelectrons);

   //get jet collection
   edm::Handle<std::vector<pat::Jet>  > jets;
   iEvent.getByLabel("goodJets",jets);
   
   // Get the Z->mm collection
   edm::Handle<reco::CompositeCandidateCollection> zmm;
   iEvent.getByLabel("zmuMatchedmuMatched", zmm);
   
   // Get the Z->ee collection
   edm::Handle<reco::CompositeCandidateCollection> zee;
   iEvent.getByLabel("zeleMatchedeleMatched", zee);

   // Get tracks
   edm::Handle<std::vector<reco::Track> >  tracks;
   iEvent.getByLabel("generalTracks",tracks);
    
   //Get MET
   edm::Handle< std::vector<pat::MET> > met;
   iEvent.getByLabel(edm::InputTag ("patMETsPFlow"), met);

   int Ntracks=0;
   int Nmu=0;
   Nj=0;
   Nb=0;
   jet_pt=0;
   jet_phi=0;
   ele_pt=0;
   muon_pt=0;
   diele_inv=0;
   dimuon_inv=0;
   NZ=0;
   NZ2=0;
   discrCSV=0; 
   b_leading_pt=0;
   Delta_phi_mm=0;
   b_mm_invmass=0;
   b_ee_invmass=0;
   
if(zmm->size()!=0 || zee->size()!=0){
   
   // Loop over pat jets
   for(std::vector<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
	    if(jet->pt() > 30 && abs(jet->eta()) < 2.4
	       && jet->chargedEmEnergyFraction()<0.99
	       && jet->neutralHadronEnergyFraction()<0.99
	       && jet->neutralEmEnergyFraction()<0.99
	       && jet->chargedHadronEnergyFraction()>0
	       && jet->chargedMultiplicity()>0){		
	       
	        ++Nj;
	        jet_pt = jet->pt();

		if(Nj>0){
	        h_jet_pt -> Fill(jet_pt);
	        w_jet_pt -> Fill(jet_pt, MyWeight);
		}

		/*
  // Flavour studies 

		Flavour flavour;
		// find out the jet flavour (differs between quark and anti-quark)
		switch(std::abs(jet->partonFlavour())) {
			case 1:
			case 2:
			case 3:
			case 21:
				flavour = UDSG_JETS;
				break;
			case 4:
				flavour = C_JETS;
				break;
			case 5:
				flavour = B_JETS;
				break;
			default:
				flavour = NONID_JETS;
		}

			
		const reco::SecondaryVertexTagInfo &svTagInfo = *jet->tagInfoSecondaryVertex();
		if (svTagInfo.nVertices() < 1) continue;
	        const reco::Vertex &sv = svTagInfo.secondaryVertex(0); 
	
		// simply count the number of accepted jets
		flavours_->Fill(ALL_JETS, MyWeight);
		flavours_->Fill(flavour, MyWeight);

		// dxy
		Measurement1D distance = svTagInfo.flightDistance(0, true);
		plots_[flavour].dist->Fill(distance.value());
 
		*/
                
   // b studies
	
                discrCSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
                	
                if(abs(jet->partonFlavour())==5 && Nj>0){
                	w_bquarks->Fill(discrCSV, MyWeight);
                }
                
                w_secondvtx_N ->Fill(discrCSV, MyWeight);
                if(discrCSV>0.89) ++Nb;
                if(discrCSV>0.89) jet_phi = jet->phi();
		if(discrCSV>0.89 && Nb!=0) b_leading_pt = jet->pt();
		cout<<"bjetpt="<<b_leading_pt<<endl;
	    }
   } 
   
   // MET

     for(std::vector<pat::MET>::const_iterator m=met->begin(); m!=met->end(); ++m){
	   w_MET->Fill(m->et());
     }

     for(std::vector<reco::Track>::const_iterator t=tracks->begin(); t!=tracks->end(); ++t){
 	   Ntracks = t->numberOfValidHits();
     }

   // Loop over triggered matched muons

     for(pat::MuonCollection::const_iterator mu=Trigmuons->begin(); mu!=Trigmuons->end(); ++mu){
	   if(mu->pt()>25 && abs(mu->eta())<2.5){
	      Nmu++;
	      muon_pt = mu->pt();
              if(muon_pt!=0 && Nj>0) {
		      h_muon_pt->Fill(muon_pt);
		      w_muon_pt->Fill(muon_pt, MyWeight);
		      cout<<"mu N="<<Nmu<<"  "<<"pt="<<muon_pt<<endl;
	      }
      		      	      
	   }
     }

   // Loop over triggered matched muons

     for(pat::ElectronCollection::const_iterator el=Trigelectrons->begin(); el!=Trigelectrons->end(); ++el){
           if(el->pt()>25 && abs(el->eta())<2.5){
              ele_pt = el->pt();
	      if(ele_pt!=0 && Nj>0){
		      h_ele_pt->Fill(ele_pt);
		      w_ele_pt->Fill(ele_pt, MyWeight);
	      }
           }
     }

   // Loop over Z bosons candidates

      for (std::vector<reco::CompositeCandidate>::const_iterator Zit = zmm->begin() ; Zit != zmm->end(); ++Zit){

	     NZ++;
  	     //cout<<"mass mm="<<Zit->mass()<<endl;
	     dimuon_inv = Zit->mass();
	     dimuon_phi = Zit->phi();

	     if(dimuon_inv!=0 && Nj>0){
		     h_mm_inv->Fill(dimuon_inv);
		     w_mm_inv->Fill(dimuon_inv, MyWeight);
		     if(Nb!=0 && discrCSV>0.89) b_mm_invmass = Zit->mass();
		     b_mm_inv->Fill(b_mm_invmass, MyWeight);
		     cout<<"bZmass="<<b_mm_invmass<<endl;
		     
	     }
      }
      for (std::vector<reco::CompositeCandidate>::const_iterator Zit2 = zee->begin() ; Zit2 != zee->end(); ++Zit2){

	     NZ2++;
  	     cout<<"mass ee="<<Zit2->mass()<<endl;

	     diele_inv = Zit2->mass();
	     diele_phi = Zit2->phi();

	     if(diele_inv!=0 && Nj>0){
		     h_ee_inv->Fill(diele_inv);
		     w_ee_inv->Fill(diele_inv, MyWeight);
		     if(Nb!=0 && discrCSV>0.89) b_ee_invmass = Zit2->mass();
		     b_ee_inv->Fill(b_ee_invmass, MyWeight);
		     cout<<"bZeeinv="<<b_ee_invmass<<endl;

	     }
      }

	//std::cout<<Nj<<std::endl;
	treeZb_->Fill();

	h_PUweights ->Fill(MyWeight);

	if(Nj>0){
		h_jetmultiplicity -> Fill(Nj);
		w_jetmultiplicity -> Fill(Nj, MyWeight);
	}

        h_secondvtx_N ->Fill(discrCSV);
        //w_secondvtx_N ->Fill(discrCSV, MyWeight);

	recoVTX_w->Fill(float(NVtx)-1, MyWeight);
	recoVTX_->Fill(float(NVtx)-1);
	h_tracks->Fill(Ntracks);
	w_tracks->Fill(Ntracks, MyWeight);   


	// Phi studies

	if(jet_phi>0){
	       	Delta_phi_mm = fabs(dimuon_phi - jet_phi);
	}

	if((Delta_phi_mm) > acos(-1)) Delta_phi_mm = 2*acos(-1) - Delta_phi_mm;

	w_delta_phi_mm -> Fill(Delta_phi_mm, MyWeight);
	if(Nb>0) w_bmultiplicity -> Fill(Nb, MyWeight);
        w_bleading_pt-> Fill(b_leading_pt, MyWeight);


   }
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
