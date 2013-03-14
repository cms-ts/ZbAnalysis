// -*- C++ -*-
//
// Package:    ZJetsAnalyzer
// Class:      ZJetsAnalyzer
// 
/**\class ZJetsAnalyzer ZJetsAnalyzer.cc analyzer/ZJetsAnalyzer/src/ZJetsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Chiara La Licata
//         Created:  Mon Feb 11 13:52:51 CET 2013
// $Id: ZJetsAnalyzer.cc,v 1.2 2013/03/13 10:07:50 clalicat Exp $
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
#include "table.h"
table MuonEff("muon_eff.txt");
table EleEff("ele_eff.txt");

//
// class declaration
//

class ZJetsAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ZJetsAnalyzer(const edm::ParameterSet&);
      ~ZJetsAnalyzer();

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

      double ele_pt;
      double ele_eta;
      double muon_pt;
      double muon_eta;
      double jet_pt;
      double jet_eta;
      double diele_inv;
      double dimuon_inv;
     
      double MyWeight;

      double ht;

      double jet_phi;
 
      int NZ_ele, NZ_muon;     

      int Nj;

      TH1F* h_NZ_ee;
      TH1F* h_NZ_mm;

      TH1F* h_nEvent;

      TH1F* h_ele_pt;
      TH1F* h_first_ele_pt;
      TH1F* h_second_ele_pt;
      TH1F* h_first_ele_eta;
      TH1F* h_second_ele_eta;

      TH1F* h_muon_pt;
      TH1F* h_first_muon_pt;
      TH1F* h_second_muon_pt;
      TH1F* h_first_muon_eta;
      TH1F* h_second_muon_eta;

      TH1F* h_jet_pt;
      TH1F* h_jet_mult;
      TH1F* h_first_jet_pt;
      TH1F* h_second_jet_pt;
      TH1F* h_third_jet_pt;
      TH1F* h_fourth_jet_pt;
      TH1F* h_fifth_jet_pt;
      TH1F* h_first_jet_eta;
      TH1F* h_second_jet_eta;
      TH1F* h_third_jet_eta;
      TH1F* h_fourth_jet_eta;
      TH1F* h_fifth_jet_eta;
      TH1F* h_jet_pt_rescaled;
      TH1F* h_ht;
      TH1D* recoVTX_;
      TH1D* recoVTX_w;
      TH1F* h_mm_inv;
      TH1F* h_ee_inv;

      TH1F* h_PUweights;

      TH1F* h_pt_Z_ee;
      TH1F* h_pt_Z_mm;
  
      TH1F* h_delta_ee;
      TH1F* h_delta_mm;
     
      TH1F* h_size_ele;
      TH2F* h_pt_first_vs_second;
      TH1F* h_check_sort;

      
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
ZJetsAnalyzer::ZJetsAnalyzer(const edm::ParameterSet& iConfig)

{

	pileup_ = iConfig.getUntrackedParameter<std::string>("pileup", "S7");
	edm::Service<TFileService> fs;
   	//now do what ever initialization is needed

	h_first_ele_pt = fs->make<TH1F>("first_ele_pt", "first_ele_pt;P_t [GeV]",100,0.,200.);
        h_second_ele_pt = fs->make<TH1F>("second_ele_pt", "second_ele_pt;P_t [GeV]",100,0.,200.);
	h_first_ele_eta =  fs->make<TH1F>("first_ele_eta", "first_ele_eta;Eta ",16,-2.5,2.5);
	h_second_ele_eta =  fs->make<TH1F>("second_ele_eta", "second_ele_eta;Eta ",16,-2.5,2.5);


	h_first_muon_pt = fs->make<TH1F>("first_muon_pt", "first_muon_pt;P_t [GeV]",100,0.,200.);
	h_second_muon_pt = fs->make<TH1F>("second_muon_pt", "second_muon_pt;P_t [GeV]",100,0.,200.);
        h_first_muon_eta =  fs->make<TH1F>("first_muon_eta", "first_muon_eta;Eta ",16,-2.5,2.5);
        h_second_muon_eta =  fs->make<TH1F>("second_muon_eta", "second_muon_eta;Eta ",16,-2.5,2.5);

	h_jet_mult = fs->make<TH1F>("jet_multiplicity", "jet_multiplicity;Nj",8,0.5,8.5);
	h_first_jet_pt = fs->make<TH1F>("first_jet_pt", "first_jet_pt;P_t [GeV]",15,30.,330.);
        h_second_jet_pt = fs->make<TH1F>("second_jet_pt", "second_jet_pt;P_t [GeV]",10,30.,330.);
        h_third_jet_pt = fs->make<TH1F>("third_jet_pt", "third_jet_pt;P_t [GeV]",8,30.,200.);
        h_fourth_jet_pt = fs->make<TH1F>("fourth_jet_pt", "fourth_jet_pt;P_t [GeV]",7,30.,100.);
	h_fifth_jet_pt = fs->make<TH1F>("fifth_jet_pt", "fifth_jet_pt;P_t [GeV]",7,30.,100.);
	h_first_jet_eta =  fs->make<TH1F>("first_jet_eta", "first_jet_eta;Eta ",24,-2.5,2.5);
	h_second_jet_eta =  fs->make<TH1F>("second_jet_eta", "second_jet_eta;Eta ",20,-2.5,2.5);
	h_third_jet_eta =  fs->make<TH1F>("third_jet_eta", "third_jet_eta;Eta ",16,-2.5,2.5);
	h_fourth_jet_eta =  fs->make<TH1F>("fourth_jet_eta", "fourth_jet_eta;Eta ",12,-2.5,2.5);
        h_fifth_jet_eta =  fs->make<TH1F>("fifth_jet_eta", "fifth_jet_eta;Eta ",12,-2.5,2.5);
	h_ht = fs->make<TH1F>("ht", "ht;P_t [GeV]",15,30.,330.);
 
	h_mm_inv = fs->make<TH1F>("h_mm_inv","h_mm_inv",60,60,120);
	h_ee_inv = fs->make<TH1F>("h_ee_inv","h_ee_inv",60,60,120);

	h_pt_Z_ee = fs->make<TH1F>("Z_pt_ee", "Z_pt_ee;P_t [GeV]",100,0.,200.);
	h_pt_Z_mm = fs->make<TH1F>("Z_pt_mm", "Z_pt_mm;P_t [GeV]",100,0.,200.);

	h_delta_ee = fs->make<TH1F>("deltaPhi_ee", "deltaPhi_ee;DeltaPhi",12,0,TMath::Pi());
	h_delta_mm = fs->make<TH1F>("deltaPhi_mm", "deltaPhi_mm;DeltaPhi",12,0,TMath::Pi());

	h_NZ_ee = fs->make<TH1F>("NZ ele", "NZ ele",4,-0.5,3.5);
	h_NZ_mm = fs->make<TH1F>("NZ muon", "NZ muon",4,-0.5,3.5);

	h_nEvent = fs->make<TH1F>("nEvent", "nEvent",4,-1.5,2.5);
	h_PUweights   = fs->make<TH1F>("h_pu_weights", "h_pu_weights", 10, 0, 10);


}


ZJetsAnalyzer::~ZJetsAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;


//PU REWEIGHTING


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
           LumiWeights_ = edm::LumiReWeighting("/gpfs/cms/users/lalicata/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/work/pileup/pileup_" + pileup_ + ".root",
                                               "/gpfs/cms/users/lalicata/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/work/pileup/pileup_2012.root", "pileup", "pileup");



           MyWeight = LumiWeights_.weight( Tnpv );
   }

cout<<"weight="<<MyWeight<<endl;
h_PUweights ->Fill(MyWeight);

/*
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
    recoVTX_->Fill(float(NVtx));
    recoVTX_w ->Fill(NVtx, MyWeight);
*/

//electrons collection
edm::Handle<pat::ElectronCollection> Trigelectrons;
iEvent.getByLabel("matchedElectrons", Trigelectrons);

//muon collection
edm::Handle<pat::MuonCollection> Trigmuons;
iEvent.getByLabel("matchedMuons", Trigmuons);

//jet collection
edm::Handle<std::vector<pat::Jet>> jets;
iEvent.getByLabel("goodJets",jets);

//Z->mm collection
edm::Handle<reco::CompositeCandidateCollection> zmm;
iEvent.getByLabel("zmuMatchedmuMatched", zmm);

//Z->ee collection
edm::Handle<reco::CompositeCandidateCollection> zee;
iEvent.getByLabel("zeleMatchedeleMatched", zee);


bool ee_event=false;
bool mm_event=false;
ele_pt=0;
muon_pt=0;
ele_eta=0;
muon_eta=0;


//+++++++++ ELECTRONS
vector <double> vect_ele_pt;
vector <double> vect_ele_eta;

for(pat::ElectronCollection::const_iterator ele=Trigelectrons->begin(); ele!=Trigelectrons->end(); ++ele){
	ele_pt=ele->pt();
	ele_eta=ele->eta();
	
	vect_ele_pt.push_back(ele_pt);
	vect_ele_eta.push_back(ele_eta);
		
}


if(vect_ele_pt.size()!=0)
	if(fabs(vect_ele_eta[0])<2.4 && fabs(vect_ele_eta[1])<2.4 && (fabs(vect_ele_eta[0])<1.44442 || fabs(vect_ele_eta[0])>1.5660) && (fabs(vect_ele_eta[1])<1.4442 || fabs(vect_ele_eta[1])>1.5660))
		if(vect_ele_pt[0]>25 && vect_ele_pt[1]>25)
			ee_event=true;


if(ee_event){
double scalFac = EleEff.Val(vect_ele_pt[0],vect_ele_eta[0]);
cout << "EVENTO ELETTRONI" << endl;
cout << "scalFac ele= " << scalFac << endl;
cout << "pt ele = " << vect_ele_pt[0] << endl;
cout << "eta ele = " << vect_ele_eta[0] << endl;
}


//+++++++++ MUONS
vector <double> vect_muon_pt;
vector <double> vect_muon_eta;

for(pat::MuonCollection::const_iterator muon=Trigmuons->begin(); muon!=Trigmuons->end(); ++muon){
	muon_pt=muon->pt();
	muon_eta=muon->eta();

	vect_muon_pt.push_back(muon_pt);
	vect_muon_eta.push_back(muon_eta);
}

if(vect_muon_pt.size()!=0)
	if(fabs(vect_muon_eta[0])<2.4 && fabs(vect_muon_eta[1])<2.4)
                        if(vect_muon_pt[0]>25 && vect_muon_pt[1]>25)
				mm_event=true;


if(mm_event){
double scalFac = MuonEff.Val(vect_muon_pt[0],vect_muon_eta[0]);
cout << "EVENTO MUONI" << endl;
cout << "scalFac muon= " << scalFac << endl;
cout << "pt muon = " << vect_muon_pt[0] << endl;
cout << "eta muon = " << vect_muon_eta[0] << endl; 
}

//++++++++ JETS
std::vector <double> vect_jet_pt;
std::vector <double> vect_jet_eta;
std::vector <double> vect_jet_phi;

Nj=0;
ht = 0;
jet_pt=0;
jet_eta=0;
jet_phi=0;

 if((ee_event || mm_event) && (zmm->size()!=0 || zee->size()!=0)){
	for(std::vector<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
		if(jet->pt() > 30 && fabs(jet->eta()) < 2.4
		   && jet->chargedEmEnergyFraction()<0.99
               	   && jet->neutralHadronEnergyFraction()<0.99
               	   && jet->neutralEmEnergyFraction()<0.99
              	   && jet->chargedHadronEnergyFraction()>0
               	   && jet->chargedMultiplicity()>0){
			jet_pt = jet->pt();
                	jet_eta = jet->eta();
                	jet_phi = jet->phi();
			++Nj;
			ht = ht+jet_pt;
			vect_jet_pt.push_back(jet_pt);
			vect_jet_eta.push_back(jet_eta);
			vect_jet_phi.push_back(jet_phi);
		}

	}

}


int nEvent=-1;

if((zmm->size()!=0 || zee->size()!=0) && (ee_event || mm_event) && Nj!=0)
nEvent=1;

h_nEvent->Fill(nEvent);

if(ht!=0)
	h_ht->Fill(ht,MyWeight);

if(vect_jet_pt.size()!=0){
	h_jet_mult->Fill(vect_jet_pt.size(),MyWeight);
}

if(vect_jet_pt.size()>0){
	h_first_jet_pt->Fill(vect_jet_pt[0],MyWeight);
	h_first_jet_eta->Fill(vect_jet_eta[0],MyWeight);
}

if(vect_jet_pt.size()>1){
	h_second_jet_pt->Fill(vect_jet_pt[1],MyWeight);
	h_second_jet_eta->Fill(vect_jet_eta[1],MyWeight);
}

if(vect_jet_pt.size()>2){
	h_third_jet_pt->Fill(vect_jet_pt[2],MyWeight);
	h_third_jet_eta->Fill(vect_jet_eta[2],MyWeight);	
}

if(vect_jet_pt.size()>3){
	h_fourth_jet_pt->Fill(vect_jet_pt[3],MyWeight);
	h_fourth_jet_eta->Fill(vect_jet_eta[3],MyWeight);
}

if(vect_jet_pt.size()>4){
	h_fifth_jet_pt->Fill(vect_jet_pt[4],MyWeight);
	h_fifth_jet_eta->Fill(vect_jet_eta[4],MyWeight);
}

if(ee_event && vect_jet_pt.size()!=0 && zee->size()!=0){
	h_first_ele_pt->Fill(vect_ele_pt[0],MyWeight);
	h_second_ele_pt->Fill(vect_ele_pt[1],MyWeight);
	h_first_ele_eta->Fill(vect_ele_eta[0],MyWeight);
	h_second_ele_eta->Fill(vect_ele_eta[1],MyWeight);
}

if(mm_event && vect_jet_pt.size()!=0 && zmm->size()!=0){
	h_first_muon_pt->Fill(vect_muon_pt[0],MyWeight);
	h_second_muon_pt->Fill(vect_muon_pt[1],MyWeight);
	h_first_muon_eta->Fill(vect_muon_eta[0],MyWeight);
	h_second_muon_eta->Fill(vect_muon_eta[1],MyWeight);
}


//+++++++++ Z in ee
NZ_ele=0;

double pt_Z;
double phi_Z;
double Delta_phi;


if(zee->size()!=0 && vect_jet_pt.size()!=0 && ee_event){

	for (std::vector<reco::CompositeCandidate>::const_iterator Zit2 = zee->begin() ; Zit2 != zee->end(); ++Zit2){

		NZ_ele++;
		diele_inv = Zit2->mass();
		pt_Z = Zit2->pt();
		phi_Z = Zit2->phi();
		

		if(diele_inv!=0){
		        
                	Delta_phi = fabs(phi_Z - vect_jet_phi[0]);

                	if(Delta_phi > acos(-1))
                        	Delta_phi = 2*acos(-1) - Delta_phi;

			h_NZ_ee->Fill(NZ_ele,MyWeight);
			h_ee_inv->Fill(diele_inv,MyWeight);
			h_pt_Z_ee -> Fill(pt_Z,MyWeight);
			h_delta_ee -> Fill(Delta_phi,MyWeight);
		}

	}

}





//++++++++ Z in mm
NZ_muon=0;

if(zmm->size()!=0 && vect_jet_pt.size()!=0 && mm_event){

	for (std::vector<reco::CompositeCandidate>::const_iterator Zit = zmm->begin() ; Zit != zmm->end(); ++Zit){

		NZ_muon++;
		dimuon_inv = Zit->mass();
		pt_Z = Zit->pt();
	       	phi_Z = Zit->phi();	

		if(dimuon_inv!=0){

			Delta_phi = fabs(phi_Z - vect_jet_phi[0]);

                	if(Delta_phi > acos(-1))
                        	Delta_phi = 2*acos(-1) - Delta_phi;

			h_NZ_mm->Fill(NZ_muon,MyWeight);
			h_mm_inv->Fill(dimuon_inv,MyWeight);
			h_pt_Z_mm -> Fill(pt_Z,MyWeight);
                	h_delta_mm -> Fill(Delta_phi,MyWeight);
		}

	}

}


}


// ------------ method called once each job just before starting event loop  ------------
void 
ZJetsAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZJetsAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ZJetsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZJetsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZJetsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZJetsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZJetsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZJetsAnalyzer);
