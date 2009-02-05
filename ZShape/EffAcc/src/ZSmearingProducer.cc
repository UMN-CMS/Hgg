// -*- C++ -*-
//
// Package:    ZSmearingProducer
// Class:      ZSmearingProducer
// 
/**\class ZSmearingProducer ZSmearingProducer.cc ZShape/ZSmearingProducer/src/ZSmearingProducer.cc

Description: <one line class summary>
Produces pair of electrons (GenParticleCollection) from Z decay as generated with Phytia
Implementation:
<Notes on implementation>
*/



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

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/GeometryVector/interface/Phi.h"

#include "Math/VectorUtil.h"

#include "TRandom3.h"

#include "ZShape/Base/interface/ZShapeElectron.h"

//
// class decleration
//


class ZSmearingProducer : public edm::EDProducer {
public:
  explicit ZSmearingProducer(const edm::ParameterSet&);
  ~ZSmearingProducer();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void smearElectron(math::PtEtaPhiELorentzVector &electron) ;
  double smearECAL (double eta, double et, char *det);
      
  // ----------member data ---------------------------

  edm::InputTag m_srcTag;
  edm::InputTag zElectronsTag;
  bool quiet_;
  float zElectronsCone_;
  math::XYZPoint vtx_;

};


ZSmearingProducer::ZSmearingProducer(const edm::ParameterSet& iConfig) : 
    m_srcTag(iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("source"))),
    zElectronsTag(iConfig.getUntrackedParameter<edm::InputTag>("zElectronsCollection",edm::InputTag("ZIntoElectronsEventProducer:ZEventParticles"))),
    quiet_(iConfig.getUntrackedParameter<bool>("quiet",false)),
    zElectronsCone_(iConfig.getParameter<double>("zElectronsCone"))
{
  using namespace reco;

  //register my product
  produces<GenParticleCollection>("ZEventParticles");

}


ZSmearingProducer::~ZSmearingProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ZSmearingProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm; 
  using namespace std; 
  using namespace reco;

  //The collection of new particles
  std::auto_ptr<GenParticleCollection> pZeeParticlesN(  new GenParticleCollection() );

  pZeeParticlesN->clear();

  //Get the collection of orginal particles
  Handle <GenParticleCollection> pZeeParticles;
  iEvent.getByLabel(zElectronsTag,pZeeParticles);

  const GenParticleCollection * ZeeParticles = pZeeParticles.product();
  //std::cout << " Starting the event " << std::endl;
  ////////////////////////
  // looking for electrons
  //int myPid = 11;
  int ne=0;
  int ngot[2]; ngot[0]=0; ngot[1]=0;
  
  GenParticleCollection::const_iterator me[2];
  //math::XYZTLorentzVector elecs[2];
  math::PtEtaPhiELorentzVector elecs[2];

  for(GenParticleCollection::const_iterator ZeeElectron_itr = ZeeParticles->begin();
      ZeeElectron_itr != ZeeParticles->end();
      ZeeElectron_itr++
      )
    if (abs(ZeeElectron_itr->pdgId())==11) {
      if (ne<2) {
	me[ne]=ZeeElectron_itr;
	ne++;
      }
    }

   //Now we look in a cone, we can set the cone to zero to not got any here
  for(
      GenParticleCollection::const_iterator ZeeElectron_itr = ZeeParticles->begin();
      ZeeElectron_itr != ZeeParticles->end();
      ZeeElectron_itr++
      )
    {
      double dR=ROOT::Math::VectorUtil::DeltaR(me[0]->momentum(),ZeeElectron_itr->momentum());
      if (me[0]==ZeeElectron_itr || dR<zElectronsCone_) {
	elecs[0]+=ZeeElectron_itr->p4();
	vtx_=math::XYZPoint(ZeeElectron_itr->vertex().x(),ZeeElectron_itr->vertex().y(),ZeeElectron_itr->vertex().z());
	ngot[0]++;
      } else {
	dR=ROOT::Math::VectorUtil::DeltaR(me[1]->momentum(),ZeeElectron_itr->momentum());
	if (me[1]==ZeeElectron_itr || dR<zElectronsCone_) {
	  elecs[1]+=ZeeElectron_itr->p4();
	  ngot[1]++;
	}
      }
    }
  
  //now correct the energies with the "function"
  for (int electron = 0; electron < ne; ++electron)
    {
      //std::cout << " looking at electron " << electron << std::endl;
      smearElectron(elecs[electron]);
       //std::cout << " NNew E " << elecs[electron].E() << " new eta " << elecs[electron].eta() << " new et " << elecs[electron].Et()<< " new pt " << elecs[electron].pt() << std::endl;
    }  

  //put the corrected electrons back in the event
  for(
      GenParticleCollection::const_iterator ZeeElectron_itr = ZeeParticles->begin();
      ZeeElectron_itr != ZeeParticles->end();
      ZeeElectron_itr++
      )
    {
       for (int electron = 0; electron < ne; ++electron)
         {
         if (me[electron]==ZeeElectron_itr)
           {
             //pdgId = ZeeElectron_itr->
             pZeeParticlesN->push_back( GenParticle(ZeeElectron_itr->charge(),
                                            math::XYZTLorentzVector(elecs[electron]),
                                            vtx_,
                                            ZeeElectron_itr->pdgId(),
                                            ZeeElectron_itr->status(),
                                            true
                                            )
                                );
           }
         }
    }

  iEvent.put(pZeeParticlesN,"ZEventParticles");
  
  
}

//Function that does the smearing of an electron
void ZSmearingProducer::smearElectron(math::PtEtaPhiELorentzVector &electron)
{
  //Possibly smear these as a function of detector eta rather than physics eta
  //I expect those near the edges either do not get "smeared" into another detector
  //or are cut out by the acceptance cuts
  const double EBEEboundry = 1.5;
  const double EEHFbountry = 3.0;
  //double eta = electron.eta();
  double e  = electron.E();
  double et  = electron.Et();
  //double pt = electron.Pt();

  double Efrac = 1.0;
  
  ZShapeElectron myElec;
  myElec.p4_ = math::PtEtaPhiMLorentzVector(electron);
 
  double deteta = myElec.detectorEta(vtx_);
  
  if (fabs(deteta) < EBEEboundry )
  {
    //Do EB stuff
    Efrac = smearECAL(deteta,et, "EB");
   
  }
  else if (fabs(deteta) < EEHFbountry)
  {
    //Do EE Stuff
    Efrac = smearECAL(deteta,et, "EE");
  }
  else if (fabs(deteta) < 5.0)
  {
    //Do HF Stuff
    TRandom3 randomNum(0);
    double mysig = 0.13;
    if (e < 140. ) mysig = 0.34;
    else if (e < 210. ) mysig = 0.25;
    else if (e < 420. ) mysig = 0.20;
    else if (e < 560. ) mysig = 0.18;
    else if (e < 840. ) mysig = 0.16;
    else if (e < 1120. ) mysig = 0.14;
    else if (e < 1400. ) mysig = 0.15;
    Efrac = randomNum.Gaus(1.02,mysig);
  }
  

  //double pz = electron.pz();
  electron.SetE(Efrac*e);
  //electron.SetPt(pow(-e*e*(1-Efrac*Efrac)+pt*pt,0.5));
  electron.SetPt(electron.Et());

  //std::cout << " New E " << electron.E() << " new eta " << electron.eta() << " new et " << electron.Et() << " new pt " << electron.pt() <<  " mass : " << electron.M() << std::endl;
}

//Method for smearing ECAL
double ZSmearingProducer::smearECAL(double eta, double et, char *det)
{
  double p0r = 0.;
  double p1r = 0.;
  double p2r = 0.;
  double p0l = 0.;
  double p1l = 0.;
  double p2l = 0.;

  if ( det == "EB" ) 
    {
      //First the right, then the left sides sigma
      p0r = 0.0028 + 0.380/et;
      p1r = 0.0267 - 0.617/et;
      p2r = -0.0146 + 0.597/et;

      p0l = 0.0057 + 0.321/et;
      p1l = 0.0185 - 0.410/et;
      p2l = -0.010 + 0.789/et;
    }
  else if ( det == "EE" ) 
    {
      p0r = +6.29/et;
      p1r = -4.53/et;
      p2r = +0.95/et;

      p0l = +7.98/et;
      p1l = -5.56/et;
      p2l = +1.25/et;
    }
  else return -1000.;

  double sigmaL = p0l + p1l*fabs(eta) + p2l*eta*eta;
  double sigmaR = p0r + p1r*fabs(eta) + p2r*eta*eta;
  //std::cout << " Side " << det << " sigmaL " << sigmaL << " sigmaR " << sigmaR << std::endl;

  TRandom3 randomNum(0);
  double randomnumber = randomNum.Uniform(0,1); //(between 0 and 1)
  double leftfraction = (1.0/(1.0+sigmaR/sigmaL));
  double frac = -1001.;
  if (randomnumber < leftfraction)
    {
      double myfrac = randomNum.Gaus(1.01,sigmaL);
      frac = (1.0-fabs(1.0-myfrac));
      //std::cout << " Went Left Myfac is " << myfrac << " and final Frac is " << frac << std::endl; 
    }
  else
    {
      double myfrac = randomNum.Gaus(1.01,sigmaR);
      frac = (1.0+fabs(1.0-myfrac));
      //std::cout << " Went Right Myfac is " << myfrac << " and final Frac is " << frac << std::endl; 
    }

  return frac;
}

// ------------ Method called once each job just before starting event loop  ------------
void ZSmearingProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZSmearingProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZSmearingProducer);
