#ifndef EFFCUT_INC
#define EFFCUT_INC

#include "TH1.h"
#include <string>
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "ZShape/Base/interface/ZShapeElectron.h"
#include "PhysicsTools/TagAndProbe/interface/EffTableLoader.h"
//
// holds efficiency info and performs the cut. Recognizes well-knonw variables

class EfficiencyCut
{
  public :

  EfficiencyCut( void );
  EfficiencyCut( TH1F * histo);

  bool passesCut( float variable ) const;
  bool passesCut( int index ) const;
  //  bool passesCut( const ::math::PtEtaPhiMLorentzVector& elec) const;
  bool passesCut( const ZShapeElectron& elec) const;
  bool passesCut( const ZShapeElectron& elec, EffTableLoader* effTable) const;
 private:
  TH1F * theClonedEffHisto_;
  //  TH1F * theTmpHisto_;
  enum CutVariable { cv_DetEta, cv_Eta, cv_Phi, cv_Pt, cv_Energy, cv_DetEta_Pt, cv_Eta_Pt };
  CutVariable theCutVariable_;
};
#endif
