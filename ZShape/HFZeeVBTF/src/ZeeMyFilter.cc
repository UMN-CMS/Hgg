// -*- C++ -*-
//
// Package:    ZeeMyFilter
// Class:      ZeeMyFilter
// 
/**\class ZeeMyFilter ZeeMyFilter.cc Work/ZeeMyFilter/src/ZeeMyFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: ZeeMyFilter.cc,v 1.21 2010/12/01 15:05:59 mansj Exp $
//
//

// system include files
#include <memory>

#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TH1.h"
#include "TH2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "ZShape/HFZeeVBTF/interface/SelectElectron.h"


//////////////////////////////////////////////////////////////////
// generically maximum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}




class ZeeMyFilter : public edm::EDFilter {
public:
  explicit ZeeMyFilter(const edm::ParameterSet&);
  ~ZeeMyFilter();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }

  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  std::string ecalID1_, ecalID2_, currentFile_;
  bool dolog_;
  std::string myName_;
  double minPtZ_,maxPtZ_;
  double minEtECAL_;
  std::vector<int> acceptedElectronIDs_;
  selectElectron w;
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
ZeeMyFilter::ZeeMyFilter(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  ecalID1_              = iConfig.getParameter<std::string>("ECALid1");
  ecalID2_              = iConfig.getParameter<std::string>("ECALid2");
  dolog_               = iConfig.getParameter<bool>("DoLog");
  acceptedElectronIDs_ = iConfig.getParameter< std::vector<int> >("acceptedElectronIDs");
  minPtZ_              = iConfig.getParameter<double>("PtZ_min");
  maxPtZ_              = iConfig.getParameter<double>("PtZ_max");
  minEtECAL_           = iConfig.getParameter<double>("minEtECAL");
  myName_              = iConfig.getParameter<std::string>("myName");

  myName_+=std::string("    ");
  
  std::vector<int>::const_iterator it;
  for(it=acceptedElectronIDs_.begin(); it!=acceptedElectronIDs_.end(); it++)
    {    w.add(*it);//      std::cout << "adding: " << (*it) << std::endl;
    }

}
  
ZeeMyFilter::~ZeeMyFilter()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
bool
ZeeMyFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  bool filterOk=false;

  edm::Handle<pat::ElectronCollection> patElectron;
  iEvent.getByLabel("patElectrons", patElectron);
  if ( ! patElectron.isValid()) {
    std::cout << "No electrons found in this event with tag " 
	      << "patElectrons" << std::endl;
    return false; // RETURN if no elecs in the event
  }

  const pat::ElectronCollection& pElecs = *(patElectron.product());
  pat::ElectronCollection::const_iterator i, ecalE1=pElecs.end(),ecalE2=pElecs.end();
  
  for (i=pElecs.begin(); i!=pElecs.end(); i++) {
    reco::SuperClusterRef scr=i->superCluster();
    if (scr.isNull()) {
      std::cout << "I'm NULL!!!\n";
      continue;
    }

    // explanation of electronID value: https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
    double eta_det=scr.get()->eta();

    // acceptance cut
    if (fabs(eta_det)>1.4442 && fabs(eta_det)<1.560) continue;

    if (i->pt()<minEtECAL_) continue;

    if (ecalE1==pElecs.end() && w.doesElePass( i->electronID(ecalID1_))) {
      ecalE1=i;
    } else if (ecalE2==pElecs.end() && w.doesElePass( i->electronID(ecalID2_))) {
      ecalE2=i;
    }
    
  }// end loop over electrons


  if (ecalE1!=pElecs.end() && ecalE2!=pElecs.end()) {
  
    // make a Z
    reco::Particle::LorentzVector Z(ecalE1->p4());
    Z+=ecalE2->p4();

    // gf: 3) loose Z mass cut and min pt
    // ==> these are the CANDIATES
    if (Z.M()>60 && Z.M()<120 && Z.pt()>minPtZ_ && Z.pt()<maxPtZ_) {

      if (dolog_) {
	std::cout << myName_  << "------------------------------------------------------" << std::endl;
	
	std::cout << myName_  << "Candidate!: event " << iEvent.id().event() << " run " << iEvent.id().run();
	
	std::cout << myName_  << "  m_ee=" << Z.M() << " Y_ee=" << Z.Rapidity() << " pT_ee=" << Z.pt() <<  std::endl;
	std::cout << myName_  << " file: " << currentFile_ << std::endl;
	
	// https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID#Electron_ID_Implementation_in_Re
	std::cout << myName_  << "ECAL1 (pt,eta,phi, eta_det) : " << ecalE1->pt() << ", " << ecalE1->eta() << ", " << ecalE1->phi() << std::endl;

	std::cout << myName_  << "  ID (60,70,80,90,95) : " 
		  << ecalE1->electronID("simpleEleId60relIso") << " " 
		  << ecalE1->electronID("simpleEleId70relIso") << " " 
		  << ecalE1->electronID("simpleEleId80relIso") << " " 
		  << ecalE1->electronID("simpleEleId90relIso") << " " 
		  << ecalE1->electronID("simpleEleId95relIso") << " " 
		  << std::endl;
	
	std::cout << myName_  << "ECAL2 (pt,eta,phi, eta_det) : " << ecalE2->pt() << ", " << ecalE2->eta() << ", " << ecalE2->phi() << std::endl;

	std::cout << myName_  << "  ID (60,70,80,90,95) : " 
		  << ecalE2->electronID("simpleEleId60relIso") << " " 
		  << ecalE2->electronID("simpleEleId70relIso") << " " 
		  << ecalE2->electronID("simpleEleId80relIso") << " " 
		  << ecalE2->electronID("simpleEleId90relIso") << " " 
		  << ecalE2->electronID("simpleEleId95relIso") << " " 
		  << std::endl;
	
	std::cout << myName_  << "======================================================" << std::endl;
      }// 4) here require electronId to have passed, for different working points
      //  5) and require HF eleID loose or tight
      //  preselections on HF isolation (0.9) and 2dvar (0.3) are already applied at clustering level.
      //  loosen as much as reasonable here

      filterOk=true;

    }// if mass is in Z window
    
  }// if the is at least on ECAL and one HF candidate


  return filterOk;
  
}


// reminder re electronID(""):
//0: fails
//1: passes electron ID only
//2: passes electron Isolation only
//3: passes electron ID and Isolation only   <<< 
//4: passes conversion rejection
//5: passes conversion rejection and ID
//6: passes conversion rejection and Isolation
//7: passes the whole selection
//https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID 

// ------------ method called once each job just before starting event loop  ------------
void 
ZeeMyFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZeeMyFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZeeMyFilter);
