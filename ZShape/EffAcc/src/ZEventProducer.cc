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
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Particle.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//
// class decleration
//


class ZEventProducer : public edm::EDProducer {
public:
  explicit ZEventProducer(const edm::ParameterSet&);
  ~ZEventProducer();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------

  edm::InputTag m_srcTag;
};


ZEventProducer::ZEventProducer(const edm::ParameterSet& iConfig) : 
  m_srcTag(iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("source")))
{
  using namespace reco;

  //register your products
  produces<GenParticleCollection>("ZEventParticles");

}


ZEventProducer::~ZEventProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ZEventProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm; 
  using namespace std; 
  using namespace reco;
  
  std::auto_ptr<GenParticleCollection> pZeeParticles(  new GenParticleCollection() );


  typedef math::PtEtaPhiMLorentzVector PolarLorentzVector;
  // electric charge type
  int q=0;
  // Lorentz vector
  typedef int Charge;
  typedef math::XYZTLorentzVector LorentzVector;
  // Lorentz vector
  typedef math::PtEtaPhiMLorentzVector PolarLorentzVector;
  // point in the space
  typedef math::XYZPoint Point;
  // point in the space
  typedef math::XYZVector Vector;
  

  Handle< HepMCProduct > EvtHandle ; 
  
  // find initial (unsmeared, unfiltered,...) HepMCProduct 
  iEvent.getByLabel(m_srcTag, EvtHandle ) ; 
  
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ; 
     
  int status;
  int pid;

  ////////////////////////
  // looking for electrons
  int myPid = 11;
  int ne=0;
  

  // looping on all the particles of the simulated event
  for(HepMC::GenEvent::particle_const_iterator mcpart = Evt->particles_begin();
      mcpart != Evt->particles_end();
      ++ mcpart ) {
    
    status = (*mcpart)->status();
    pid    = (*mcpart)->pdg_id();

    if     ((*mcpart)->pdg_id() >0) q=1;
    else  q=-1;
      
    if (status == 3 && abs(pid) == myPid ){
      
      const HepMC::GenVertex * vertex_=(*mcpart)->production_vertex();
    
      if (ne==0) {
	//evt_.z0_=::math::XYZVector(vertex_->position().x(),vertex_->position().y(),vertex_->position().z() );
      }
      
      LorentzVector momentum;
      momentum.SetPx(  (*mcpart)->momentum().x() );        momentum.SetPy(  (*mcpart)->momentum().y() );
      momentum.SetPz(  (*mcpart)->momentum().z() );        momentum.SetE(   (*mcpart)->momentum().t() );
      
      if (ne<2) {
	pZeeParticles->push_back( GenParticle(q,
					      momentum,
					      Point(vertex_->position().x(),vertex_->position().y(),vertex_->position().z() ),
					      pid,
					      status,
					      true
					      )
				  );
      }     
      ne++;
    
    }
    
  }
  // end loop on particles

  
  //
  // check than we have 2 and only 2 electrons
  if (ne!= 2)
    {
      edm::LogWarning("ZEfficiencyCalculator") << " we need two electrons (pid: " << myPid << " ) while we have: " << ne;
      return;
    }
  
  
  //
  // reorder electrons in pT
  //float pt1 = evt_.elec(0).p4_.Pt();
  //float pt2 = evt_.elec(1).p4_.Pt();
  
  std::cout << "collection size: " << pZeeParticles->size() << std::endl;
  iEvent.put(pZeeParticles,"ZEventParticles");
  
  
}


// ------------ Method called once each job just before starting event loop  ------------
void ZEventProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZEventProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZEventProducer);
