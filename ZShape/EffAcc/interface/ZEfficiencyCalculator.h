// -*- C++ -*-
//
// Package:    McLevelEff
// Class:      ZEfficiencyCalculator
// 
/**\class ZEfficiencyCalculator ZEfficiencyCalculator.cc ZShape/EffAcc/src/ZEfficiencyCalculator.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Giovanni FRANZONI
//         Created:  Thu Oct  4 11:30:13 CEST 2007
// $Id$
//
//


// system include files
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FastSimulation/Particle/interface/RawParticle.h"

#include <CLHEP/Vector/LorentzVector.h>

#include "ZShape/EffAcc/interface/EffHistos.h"
#include "ZShape/EffAcc/interface/EtaAcceptance.h"
#include "ZShape/EffAcc/interface/EfficiencyStore.h"
#include "ZShape/EffAcc/interface/EfficiencyCut.h"


#include "TRandom3.h"
#include "TFile.h"


//
// class decleration
//

class ZEfficiencyCalculator : public edm::EDAnalyzer {
public:
  explicit ZEfficiencyCalculator(const edm::ParameterSet&);
  ~ZEfficiencyCalculator();


private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  std::vector<RawParticle> myElectrons_;
  TRandom3 randomNum;

  std::string outFileName_;
  bool        writeHistoConservatively_;
  TFile*      histoFile_;
  TFile *     eff1File_;
  TFile *     eff2File_;
  TFile *     eff3File_;

  TFile * effHistoFile_;
  std::string inputEffFileName_;

  std::string inputEffFileNameCut1_;
  std::string inputEffFileNameCut2_;
  std::string inputEffFileNameCut3_;

  //GF: these three will go away
  std::string cut1HistoName_;
  std::string cut2HistoName_;
  std::string cut3HistoName_;

  EfficiencyStore * eff1_;
  EfficiencyStore * eff2_;
  EfficiencyStore * eff3_;
  
  TH1F*  histoSeeding_;
  TH1F*  histoTracking_;
  TH1F*  histoEReco_;

  // objects holding histograms to fill
  EffHistos allCase_, acceptance_, aftercut1_, aftercut2_, aftercut3_;
  
  // object which applies the eta cuts
  EtaAcceptance theEtaSelector_;

  // objects for the efficiency cuts
  EfficiencyCut * seedingCut_;
  EfficiencyCut * trackingCut_;
  EfficiencyCut * eRecoCut_;
  
};
