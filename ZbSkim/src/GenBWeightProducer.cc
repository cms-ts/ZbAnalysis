#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>

#include "Utilities/General/interface/FileInPath.h"
struct Bhadron {
  TLorentzVector        kinematics;
  reco::SecondaryVertex vertex;
  GlobalVector          flightdir;
  int                   isSelected;
};

struct compareBCandbyPt {  
  bool operator()( const Bhadron BCand1, const Bhadron BCand2 ) const {
    return BCand1.kinematics.Pt() > BCand2.kinematics.Pt();
  }
};

class GenBWeightProducer : public edm::EDProducer {
   public:
      explicit GenBWeightProducer(const edm::ParameterSet&);
      ~GenBWeightProducer();
      float Weight;	
   private:
      virtual void beginJob();
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      void MapFiller(std::ifstream &,int);
      bool isFinal(const reco::GenParticle*,int);
      std::string GetKey(const reco::GenParticle*);
      std::map<std::string, std::map<std::string,float> > theMap;
      std::map<std::string,float> bPropMap;
      std::ifstream BRelWeights;
      std::ifstream B511DecayWeights;
      std::ifstream B521DecayWeights;
      std::ifstream B531DecayWeights;
      std::ifstream B541DecayWeights;
      std::ifstream D411DecayWeights;
      std::ifstream D421DecayWeights;
      std::ifstream D431DecayWeights;
      std::ifstream D441DecayWeights;
      //GenBDWeight GBDW;
      edm::FileInPath pprop_,p511_,p521_,p531_,p541_,p411_,p421_,p431_,p441_;
};

GenBWeightProducer::GenBWeightProducer(const edm::ParameterSet& pset):
 pprop_(pset.getParameter<edm::FileInPath>("pprop")),
 p511_(pset.getParameter<edm::FileInPath>("p511")),
 p521_(pset.getParameter<edm::FileInPath>("p521")),
 p531_(pset.getParameter<edm::FileInPath>("p531")),
 p541_(pset.getParameter<edm::FileInPath>("p541")),
 p411_(pset.getParameter<edm::FileInPath>("p411")),
 p421_(pset.getParameter<edm::FileInPath>("p421")),
 p431_(pset.getParameter<edm::FileInPath>("p431")),
 p441_(pset.getParameter<edm::FileInPath>("p441"))
{
	//std::cout<<"in const"<<std::endl;	
	//Reading the weight tables
        BRelWeights.open(pprop_.fullPath().c_str());
        B511DecayWeights.open(p511_.fullPath().c_str() );
        B521DecayWeights.open(p521_.fullPath().c_str() );
        B531DecayWeights.open(p531_.fullPath().c_str() );
        B541DecayWeights.open(p541_.fullPath().c_str() );
        D411DecayWeights.open(p411_.fullPath().c_str() );
        D421DecayWeights.open(p421_.fullPath().c_str() );
        D431DecayWeights.open(p431_.fullPath().c_str() );
        D441DecayWeights.open(p441_.fullPath().c_str() );
	
        std::string key;
	float value;
	//std::cout<<"Checking 511"<<std::endl;
	MapFiller(B511DecayWeights,511);
	//std::cout<<"Checking 521"<<std::endl;
	MapFiller(B521DecayWeights,521);
	//std::cout<<"Checking 531"<<std::endl;
	MapFiller(B531DecayWeights,531);
	//std::cout<<"Checking 541"<<std::endl;
	MapFiller(B541DecayWeights,541);
	//std::cout<<"Checking 411"<<std::endl;
        MapFiller(D411DecayWeights,411);
	//std::cout<<"Checking 421"<<std::endl;
        MapFiller(D421DecayWeights,421);
	//std::cout<<"Checking 431"<<std::endl;
        MapFiller(D431DecayWeights,431);
	//std::cout<<"Checking 441"<<std::endl;
        MapFiller(D441DecayWeights,441);
	//std::cout<<"Checking BRelWeights"<<std::endl;
	if(BRelWeights.fail())throw cms::Exception("TxtFile") << BRelWeights.get()<<" not found";
       	while(!BRelWeights.eof()){
       	        BRelWeights >> key >> value;
       	        bPropMap[key]=value;
       	}
	Weight=1;
 	BRelWeights.close();
	B511DecayWeights.close();
	B521DecayWeights.close();
	B531DecayWeights.close();
	B541DecayWeights.close();	
        D411DecayWeights.close();
        D421DecayWeights.close();
        D431DecayWeights.close();
        D441DecayWeights.close();

	//produces<GenBDWeightCollection>("GenBDWeight");
	produces<std::vector<double>>("GenBDWeight");
}


GenBWeightProducer::~GenBWeightProducer()
{}

void GenBWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	Weight=1;
	//std::auto_ptr<GenBDWeightCollection> GenBDW(new GenBDWeightCollection);
	std::auto_ptr<std::vector<double>> GenBDW( new std::vector<double> );
	//std::cout<<"in produce"<<std::endl;
 	using namespace edm;

        Handle<reco::GenParticleCollection> GPC;

        if (iEvent.getByLabel("genParticles", GPC)) {
	
	        for(reco::GenParticleCollection::const_iterator it = GPC->begin(); it != GPC->end(); ++it){
			float Id=fabs(it->pdgId());
			if(!isFinal(&*it,Id))continue;
			if(Id==411 || Id==421 || Id==431 || Id==441 || Id==511 || Id==521 || Id==531 || Id==541){
				std::ostringstream opdg;
        			opdg << Id;
				std::string spdg=opdg.str();
				std::string key=GetKey(&*it);
				float weightProd=1.0,weightDecay=1.0;
				if(Id>500)weightProd=bPropMap[spdg];
				if(theMap[spdg][key]) weightDecay=theMap[spdg][key];
				Weight*=weightDecay*weightProd;
				//std::cout<<"Found a hadron: "<<Id<<" decay="<<key<<" "<<weightProd<<" "<<weightDecay<<"==> value="<<Weight<<std::endl;
			}
        	}
        }
 
	//GBDW.theWeight(Weight);	
	GenBDW->push_back(Weight);
	//std::auto_ptr<GenBDWeightCollection> GenBDWColl(GenBDW);
	//iEvent.put(GenBDWColl,"GenBDWeight");	
	
        iEvent.put(GenBDW, "GenBDWeight");	
}	

bool GenBWeightProducer::isFinal(const reco::GenParticle* gp,int pdg)
{
        bool answer=true;
        for(unsigned int nDaught=0;nDaught<gp->numberOfDaughters();nDaught++){
                const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(gp->daughter(nDaught));
                int pidd = abs(dau->pdgId());
                if(((pidd/100)==pdg || (pidd/1000)==pdg) )answer=false;
        }
        return answer;
}


std::string GenBWeightProducer::GetKey(const reco::GenParticle* gp){
        std::vector<int> decaylist;
        for(unsigned int nDaught=0;nDaught<gp->numberOfDaughters();nDaught++){
        	int pdidc = abs(gp->daughter(nDaught)->pdgId());
        	decaylist.push_back(pdidc);
        }
        sort(decaylist.begin(), decaylist.end());
        std::string Identificator="Decay";
        for (std::vector<int>::size_type i = 0; i != decaylist.size(); ++i){
        	std::ostringstream os;
                os << decaylist[i];
                Identificator=Identificator+"_"+os.str();
	}
	return Identificator;
}

void GenBWeightProducer::MapFiller(std::ifstream &textfile,int pdg)
{
	if(textfile.fail())throw cms::Exception("TxtFile") << textfile.get() << " not found";
	//std::cout<<"in MapFiller:checking "<<textfile.get()<<std::endl;
	std::string key;
	float value;
	std::ostringstream opdg;
        opdg << pdg;
        std::string spdg=opdg.str();
	std::string line;
	while(getline(textfile, line)){
		textfile >> key >> value;
		theMap[spdg][key]=value;
		//std::cout<<"line is "<<key<<" "<<value<<std::endl;
	}
}


// ------------ method called once each job just before starting event loop ------------
void
GenBWeightProducer::beginJob()
{
 

}

// ------------ method called once each job just after ending the event loop ------------
void
GenBWeightProducer::endJob() {

}


DEFINE_FWK_MODULE(GenBWeightProducer);
