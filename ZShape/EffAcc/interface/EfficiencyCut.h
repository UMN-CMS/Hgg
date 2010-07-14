#ifndef EFFCUT_INC
#define EFFCUT_INC

#include "TH1.h"
#include <string>
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "ZShape/Base/interface/ZShapeElectron.h"
#include "ZShape/Base/interface/EffTableLoader.h"
//
// holds efficiency info and performs the cut. Recognizes well-knonw variables

class TRandom;

class EfficiencyCut
{
  public :

  EfficiencyCut( TH1F * histo, EffTableLoader* indexer);

  bool passesCut( const ZShapeElectron& elec) const;
  bool passesCut( const ZShapeElectron& elec, float level ) const;
  EffTableLoader* indexer() const { return theIndexer; }
  float lastRandomLevel() const { return lastRandomLevel_; }
 private:
  bool passesCut( int index, float level=-1 ) const;
  int indexOf( const ZShapeElectron& elec) const;
  EffTableLoader* theIndexer;
  TH1F * theClonedEffHisto_;
  mutable float lastRandomLevel_;

  static TRandom* randomGenerator_;
};
#endif
