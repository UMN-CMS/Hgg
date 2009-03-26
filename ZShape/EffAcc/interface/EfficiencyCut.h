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

  EfficiencyCut( TH1F * histo, EffTableLoader* indexer);

  bool passesCut( const ZShapeElectron& elec) const;
  EffTableLoader* indexer() const { return theIndexer; }
 private:
  bool passesCut( int index ) const;
  int indexOf( const ZShapeElectron& elec) const;
  EffTableLoader* theIndexer;
  TH1F * theClonedEffHisto_;
  
};
#endif
