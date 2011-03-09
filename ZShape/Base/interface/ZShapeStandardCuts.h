#ifndef ZSHAPE_BASE_ZSHAPESTANDARDCUTS_H
#define ZSHAPE_BASE_ZSHAPESTANDARDCUTS_H 1

#include "ZShape/Base/interface/ZShapeElectron.h"
#include "ZShape/Base/interface/EtaAcceptance.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

/** \class ZShapeStandardCuts
  *  
  * $Date: 2009/04/09 19:16:45 $
  * $Revision: 1.2 $
  * \author J. Mans - Minnesota
  */
class ZShapeStandardCuts {
public:
  void acceptanceCuts(ZShapeElectron& elec) const; 
  void ptCuts(ZShapeElectron& elec) const; 
  void gen_ptCuts(reco::GenParticle& gelec,ZShapeElectron& elec) const; 
  void dummyCuts(ZShapeElectron& elec) const; 
private:
  EtaAcceptance theEtaSelector_;
};

#endif
