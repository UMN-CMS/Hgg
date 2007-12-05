#ifndef ZSHAPE_BASE_ZSHAPESTANDARDCUTS_H
#define ZSHAPE_BASE_ZSHAPESTANDARDCUTS_H 1

#include "ZShape/Base/interface/ZShapeElectron.h"
#include "ZShape/Base/interface/EtaAcceptance.h"

/** \class ZShapeStandardCuts
  *  
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class ZShapeStandardCuts {
public:
  void acceptanceCuts(ZShapeElectron& elec) const; 
  void ptCuts(ZShapeElectron& elec) const; 
private:
  EtaAcceptance theEtaSelector_;
};

#endif
