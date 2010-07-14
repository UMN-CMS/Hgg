#include "TProfile.h"
#include "ZShape/Base/interface/EfficiencyStatistics.h"
#include "ZShape/ZFromData/interface/IsoIdFilter.h"
///#include "ZShape/Base/interface/ZShapeUtils.h"

#include <cmath>
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "ZShape/Base/interface/CandidateAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


#include <TF1.h>
#include <TMath.h>
#include <math.h>
#include <iostream>
#include <fstream>


//Double_t CrystalBall(Double_t *x,Double_t *par);
//
// constructors and destructor
//
IsoIdFilter::IsoIdFilter(const edm::ParameterSet& iConfig) 
{
  // 
  //now do what ever initialization is needed

  //This little guy reads in the differnt InputTag collections.
  const edm::InputTag dBlankTag("Blank");
  
  std::vector< edm::InputTag > defaultAllProbeCandTags;
  allProbeCandTags_ = iConfig.getUntrackedParameter< std::vector<edm::InputTag> >(
      "allProbeCandTags",defaultAllProbeCandTags);
   for( int i=0; i<(int)allProbeCandTags_.size(); ++i ) 
      std::cout << allProbeCandTags_[i] << std::endl;
  
}


IsoIdFilter::~IsoIdFilter()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
bool IsoIdFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace cms;
  using namespace reco;
 
  //if (fromfile_) return;

  edm::Handle<pat::ElectronCollection> patElectron;
  iEvent.getByLabel("patElectrons", patElectron);
  if ( ! patElectron.isValid()) {
    //std::cout << "No electrons found in this event with tag " 
    //        << electronCollectionTag_  << std::endl;
    return false; // RETURN if no elecs in the event
  }

  const pat::ElectronCollection& pElecs = *(patElectron.product());
  pat::ElectronCollection::const_iterator i, ecalE=pElecs.end();
  
  for (i=pElecs.begin(); i!=pElecs.end(); i++) {
    reco::SuperClusterRef scr=i->superCluster();
    if (scr.isNull()) {
      std::cout << "I'm NULL!!!\n";
      continue;
    }
    ecalE=i;
    std::cout << "  IDC (60,70,80,90) : " 
                  << ecalE->electronID("simpleEleId60cIso") << " " 
                  << ecalE->electronID("simpleEleId70cIso") << " " 
                  << ecalE->electronID("simpleEleId80cIso") << " " 
                  << ecalE->electronID("simpleEleId90cIso") << " " 
                  << std::endl;
    std::cout << "  IDR (60,70,80,90) : " 
                  << ecalE->electronID("simpleEleId60relIso") << " " 
                  << ecalE->electronID("simpleEleId70relIso") << " " 
                  << ecalE->electronID("simpleEleId80relIso") << " " 
                  << ecalE->electronID("simpleEleId90relIso") << " " 
                  << std::endl;
  }



  return true;


}




// ------------ method called once each job just before starting event loop  ------------
void 
IsoIdFilter::beginJob()
{



  
}




// ------------ method called once each job just after ending the event loop  ------------
void 
IsoIdFilter::endJob() {

 


}


// ************************************************************** //

//define this as a plug-in
//DEFINE_ANOTHER_FWK_MODULE(IsoIdFilter);
