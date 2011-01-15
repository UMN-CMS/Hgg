// -*- C++ -*-
//
// Package:    ZEventProducer
// Class:      ZEventProducer
// 
/**\class ZEventProducer ZEventProducer.cc ZShape/ZEventProducer/src/ZEventProducer.cc

Description: <one line class summary>
Produces pair of electrons (GenParticleCollection) from Z decay as generated with Phytia
Implementation:
<Notes on implementation>
*/
//
// Original Author:  Giovanni FRANZONI
//         Created:  Mon Feb 18 21:18:39 CET 2008
// $Id: ZEventProducer.cc,v 1.1 2010/12/14 21:26:26 mansj Exp $
//
//


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Particle.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/GeometryVector/interface/Phi.h"

#include "CommonTools/CandUtils/interface/CandCombiner.h"

#include "DataFormats/Candidate/interface/Candidate.h"


#include "TRandom3.h"
#include <TMath.h>
#include <math.h>

//
// class decleration
//


class BestZProducer : public edm::EDProducer {
public:
  explicit BestZProducer(const edm::ParameterSet&);
  ~BestZProducer();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------

  std::vector<edm::InputTag> allProbeCandTags_;
  std::string outcol_;
  int randomSeed_;
  TRandom3 *randomNum_;
};


BestZProducer::BestZProducer(const edm::ParameterSet& iConfig):
  randomSeed_(iConfig.getUntrackedParameter<int>("randomSeed",100))
{
  using namespace reco;

  //This little guy reads in the differnt InputTag collections.
  const edm::InputTag dBlankTag("Blank");
  
  std::vector< edm::InputTag > defaultAllProbeCandTags;
  allProbeCandTags_ = iConfig.getUntrackedParameter< std::vector<edm::InputTag> >( 		  
"allProbeCandTags",defaultAllProbeCandTags);
  for( int i=0; i<(int)allProbeCandTags_.size(); ++i ) 
    std::cout << allProbeCandTags_[i] << std::endl;
  
  //Setting the random number to make things reproducable
  randomNum_ = new TRandom3(randomSeed_);


  outcol_   = iConfig.getUntrackedParameter<std::string>("outcol", "");

  //register your products
  produces<reco::CandidateCollection>(outcol_);


}


BestZProducer::~BestZProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
BestZProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm; 
  using namespace std; 
  using namespace reco;
  using namespace cms;
 
  //if (fromfile_) return;
  
  edm::Handle<reco::CandidateView> tagprobes;
      if ( !iEvent.getByLabel(allProbeCandTags_[0],tagprobes) ) {
         //std::cout << " DIDN't get the darn tag " << allProbeCandTags_[i] << std::endl;
      }

  bool goodevent = false;

  if (tagprobes.isValid()) if (tagprobes->size() > 0) goodevent = true;

  double myran = tagprobes->size()*randomNum_->Uniform(0.,1.0);
  uint imyran = (uint) myran;
  if (imyran == tagprobes->size() ) imyran -= 1;


  //reco::Candidate zee;

  auto_ptr<reco::CandidateCollection> 
    mycompcol(new reco::CandidateCollection);


  if ( goodevent )
    {
      uint tpnum = 0;
      //--End picking a random Tag and Probe combination
      reco::CandidateView::const_iterator tpItr = tagprobes->begin();
      for ( ; tpItr != tagprobes->end(); ++tpItr,++tpnum)
	  { 
	    if (tpnum != imyran ) continue;

	    //zee.addDaughter(*(tpItr->daughter(0)),"electron1");
	    //zee.addDaughter(*(tpItr->daughter(1)),"electron2");
	    //mycompcol->push_back(zee);
	    mycompcol->push_back(*(tpItr));
	    
	  }
      
    }     
     


  iEvent.put(mycompcol,outcol_);

  
  
}


// ------------ Method called once each job just before starting event loop  ------------
void BestZProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BestZProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(BestZProducer);
