// -*- C++ -*-
//
// Package:    EleWithMet
// Class:      EleWithMet
// 
/**\class EleWithMet EleWithMet.cc ZShape/HFZeeVBTF/src/EleWithMet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  G. Franzoni
//         Created:  Wed Aug 04 00:21:26 CDT 2010
// $Id: EleWithMet.cc,v 1.1 2010/08/06 14:14:26 franzoni Exp $
//
//


// system include files
#include <memory>

#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TH1.h"
#include "TH2.h"
#include <TMath.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "ZShape/HFZeeVBTF/interface/SelectElectron.h"

//////////////////////////////////////////////////////////////////
// generically maximum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}

class EleWithMet : public edm::EDAnalyzer {
public:
  explicit EleWithMet(const edm::ParameterSet&);
  ~EleWithMet();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  std::string ecalID_, currentFile_;
  bool dolog_;
  std::string myName_;
  std::vector<int> acceptedElectronIDs_;
  selectElectron w;
  float ETCut_;
  float METCut_;
  edm::InputTag metCollectionTag_;


  // ----------member data ---------------------------

  // gf set of histograms per each == Z definition ==
  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory td,const std::string&);
    // fill all histos of the set with the two electron candidates
    void fill(pat::ElectronCollection::const_iterator ecalE,  const reco::RecoEcalCandidate& hfE, const reco::HFEMClusterShape& hfshape, const edm::ParameterSet& myPs, const float hf_2d_cut);
    TH1* mee, *meeHFP, *meeHFM, *Yee, *Ptee;
    TH1* ec_eta, *ec_phi, *ec_pt;
    TH1* hf_eta, *hf_phi, *hf_pt;
    TH1* el_1_eta, *el_2_eta;
    TH2* el_1_eta_VS_el_2_eta, *el_1_phi_VS_el_2_phi;

    
    TH1* eleEta, *elePhi, *eleEt;
    
  };

  std::vector<std::string> HLT_Names;

  // gf set of histo for all Z definitios in a stack
  struct HistStruct {
    TH1 *nelec,*nhf;
    HistPerDef /*base,*/ mee90none, zMass;
    HistPerDef mee90loose, mee80loose, mee70loose, mee60loose;
    HistPerDef mee90tight, mee80tight, mee70tight, mee60tight;
        
    HistPerDef ele90metLoose, ele80metLoose, ele70metLoose, ele60metLoose,
      ele90metTight, ele80metTight, ele70metTight,    ele60metTight;

  } hists;

};


void EleWithMet::HistPerDef::book(TFileDirectory td, const std::string& post) {
  std::string title;

//  title=std::string("M_{ee} ")+post;
//  mee=td.make<TH1D>("mee",title.c_str(),90,40,130);  
//  title=std::string("M_{ee,HF+} ")+post;
//  meeHFP=td.make<TH1D>("mee-HFP",title.c_str(),90,40,130);  
//  title=std::string("M_{ee,HF-} ")+post;
//  meeHFM=td.make<TH1D>("mee-HFM",title.c_str(),90,40,130);  
//  title=std::string("Y_{ee} ")+post;
//  Yee=td.make<TH1D>("yee",title.c_str(),50,-4,4);  
//  title=std::string("pT_{ee} ")+post;
//  Ptee=td.make<TH1D>("ptee",title.c_str(),90,0,90);  
//
//  title=std::string("eta_{ecal} ")+post;
//  ec_eta=td.make<TH1D>("etaecal",title.c_str(),30,-3,3);  
//  title=std::string("phi_{ecal} ")+post;
//  ec_phi=td.make<TH1D>("phiecal",title.c_str(),30,-3.14159,3.14159);  
//  title=std::string("pt_{ecal} ")+post;
//  ec_pt=td.make<TH1D>("ptecal",title.c_str(),90,0,90);  
//
//  title=std::string("eta_{hf} ")+post;
//  hf_eta=td.make<TH1D>("etahf",title.c_str(),50,-5,5);  
//  title=std::string("phi_{hf} ")+post;
//  hf_phi=td.make<TH1D>("phihf",title.c_str(),30,-3.14159,3.14159);
//  title=std::string("pt_{hf} ")+post;
//  hf_pt=td.make<TH1D>("pthf",title.c_str(),90,0,90);  
//
//
//  title=std::string("eta_{el_1} ")+post;
//  el_1_eta=td.make<TH1D>("etael1",title.c_str(),50,-5,5);  
//  title=std::string("eta_{el_2} ")+post;
//  el_2_eta=td.make<TH1D>("etael2",title.c_str(),50,-5,5);  
//
//  title=std::string("eta_{el_1}_VS_eta_{el_2} ")+post;
//  el_1_eta_VS_el_2_eta=td.make<TH2D>("etael1-vs-etael2",title.c_str(),50,-5,5,50,-5,5);
//  title=std::string("phi_{el_1}_VS_phi_{el_2} ")+post;
//  el_1_phi_VS_el_2_phi=td.make<TH2D>("phiel1-vs-phiel2",title.c_str(),30,-3.14159,3.14159,30,-3.14159,3.14159);
//  





  title=std::string("eta_{el} ")+post;
  eleEta=td.make<TH1D>("eta_{el}",title.c_str(),30,-3,3);
  title=std::string("phi_{el} ")+post;
  elePhi=td.make<TH1D>("phi_{el}",title.c_str(),30,-1*TMath::Pi(),TMath::Pi());
  title=std::string("et_{el} ")+post;
  eleEt=td.make<TH1D>("et_{el}",title.c_str(),50,0,100);



}

void EleWithMet::HistPerDef::fill(pat::ElectronCollection::const_iterator ecalE,  
				 const reco::RecoEcalCandidate& hfE, const reco::HFEMClusterShape& hfshape, 
				 const edm::ParameterSet& myPs, const float hf_2d_cut) {

  //pat::ElectronCollection::const_iterator ecalE,  
  eleEta-> Fill(ecalE->eta());
  elePhi-> Fill(ecalE->phi());
  eleEt-> Fill(ecalE->pt());

}// end of fill()

static const double hf_2d_loose = 0.32;
// this value for tight is lower than CMS AN-2009/106 to accont for S/L miscalib in data
// see figure 15 therein: turn-on is sharp
static const double hf_2d_tight = 0.45;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EleWithMet::EleWithMet(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  ecalID_=iConfig.getParameter<std::string>("ECALid");
  dolog_=iConfig.getParameter<bool>("DoLog");
  acceptedElectronIDs_ = iConfig.getParameter< std::vector<int> >("acceptedElectronIDs");
  myName_=iConfig.getParameter<std::string>("myName");
  myName_+=std::string("    ");

  metCollectionTag_ = iConfig.getParameter<edm::InputTag>("metCollectionTag");

  ETCut_ = iConfig.getParameter<double>("ETCut");
  METCut_ = iConfig.getParameter<double>("METCut");



  std::vector<int>::const_iterator it;
  for(it=acceptedElectronIDs_.begin(); it!=acceptedElectronIDs_.end(); it++)
    {    w.add(*it);//      std::cout << "adding: " << (*it) << std::endl;
    }

  edm::Service<TFileService> fs;
  hists.nelec=fs->make<TH1D>("met","met:mE_{T} [GeV]",100,0,100);
  hists.nelec=fs->make<TH1D>("met eta","metEta:#eta_{mE_{T}}",100,0,100);
  hists.nelec=fs->make<TH1D>("met phi","metPhi:#phi_{mE_{T}}",100,-1*TMath::Pi(),TMath::Pi());


  hists.nelec=fs->make<TH1D>("nelec","N_Elec",10,-0.5,9.5);
  hists.nhf=fs  ->make<TH1D>("nhf","N_HF",10,-0.5,9.5);
  //hists.base.book(fs->mkdir("base"),"(base)");
  //  hists.mee90none.book(fs->mkdir("mee90none"),"(mee90none)");
  //  hists.zMass.book(fs->mkdir("zMass"),"(zMass)");
  hists.ele90metLoose.book(fs->mkdir("ele90metLoose"),"(ele90,metLoose)");
  hists.ele80metLoose.book(fs->mkdir("ele80metLoose"),"(ele80,metLoose)");
  hists.ele70metLoose.book(fs->mkdir("ele70metLoose"),"(ele70,metLoose)");  
  hists.ele60metLoose.book(fs->mkdir("ele60metLoose"),"(ele60,metLoose)");  
  hists.ele90metTight.book(fs->mkdir("ele90metTight"),"(ele90,metTight)");  
  hists.ele80metTight.book(fs->mkdir("ele80metTight"),"(ele80,metTight)");  
  hists.ele70metTight.book(fs->mkdir("ele70metTight"),"(ele70,metTight)");  
  hists.ele60metTight.book(fs->mkdir("ele60metTight"),"(ele60,metTight)");  

}
  
EleWithMet::~EleWithMet()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EleWithMet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::RecoEcalCandidateCollection> HFElectrons;
  edm::Handle<reco::SuperClusterCollection> SuperClusters;
  edm::Handle<reco::HFEMClusterShapeAssociationCollection> ClusterAssociation;
  
  iEvent.getByLabel("hfEMClusters",SuperClusters);
  iEvent.getByLabel("hfEMClusters",ClusterAssociation);
  iEvent.getByLabel("hfRecoEcalCandidate",HFElectrons);
  
  edm::Handle<pat::ElectronCollection> patElectron;
  iEvent.getByLabel("patElectrons", patElectron);
  if ( ! patElectron.isValid()) {
    std::cout << "No electrons found in this event with tag " 
	      << "patElectrons" << std::endl;
    return ; // RETURN if no elecs in the event
  }

  const pat::ElectronCollection& pElecs = *(patElectron.product());
  pat::ElectronCollection::const_iterator i, ecalE=pElecs.end();
  
  for (i=pElecs.begin(); i!=pElecs.end(); i++) {
    reco::SuperClusterRef scr=i->superCluster();
    if (scr.isNull()) {
      std::cout << "I'm NULL!!!\n";
      continue;
    }

    // gf: ECAL eta acceptance cut on supercluster. What about phi acceptance and central crack?
    // explanation of electronID value: https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
    double eta_det=scr.get()->eta();
    if (ecalE==pElecs.end() && i->electronID(ecalID_)>=(.0-0.1)   // very loose to let events pass; cuts done later
	&& (fabs(eta_det)<1.4442 || fabs(eta_det)>1.560)
	&& HFElectrons->size()>0) {
      const reco::RecoEcalCandidate& hfE=(*HFElectrons)[0];
      reco::SuperClusterRef hfclusRef=hfE.superCluster();
      const reco::HFEMClusterShapeRef hfclusShapeRef=(*ClusterAssociation).find(hfclusRef)->val;
      // const reco::HFEMClusterShape& hfshape=*hfclusShapeRef;
      // hists.base.fill(i,hfE,hfshape,robust90relIsoEleIDCutsV04_ps_);//gf: for each Z definition, N-1 control plots could be filled here too
      ;
    }
    
    // ECAL acceptance cut
    if (ecalE==pElecs.end() && w.doesElePass( i->electronID(ecalID_) )  // 3: passes electron iso/Id/convRej depending on configuration; ecalID_ defines working point
    	&& (fabs(eta_det)<1.4442 || fabs(eta_det)>1.560)) {             // ECAL acceptance cut
      ecalE=i;
    }// if
    
  }// end loop over electrons
  
  hists.nelec->Fill(pElecs.size());
  hists.nhf  ->Fill(HFElectrons->size());

  edm::Handle<pat::METCollection> patMET;
  iEvent.getByLabel(metCollectionTag_,  patMET);
  const pat::METCollection *pMet = patMET.product();
  const pat::METCollection::const_iterator met = pMet->begin();
  const pat::MET theMET = *met;


//      hists.zMass.fill(ecalE,hfE,hfshape,robust90relIsoEleIDCutsV04_ps_,0);//gf: for each Z definition, N-1 control plots could be filled here too
//      if (var2d>=hf_2d_loose) {
//	if (w.doesElePass( ecalE->electronID("simpleEleId90relIso") )) hists.mee90loose.fill(ecalE,hfE,hfshape,robust90relIsoEleIDCutsV04_ps_,hf_2d_loose);


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
EleWithMet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EleWithMet::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EleWithMet);
