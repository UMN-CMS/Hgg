#ifndef ZSHAPE_BASE_ZSHAPEEVENT_H
#define ZSHAPE_BASE_ZSHAPEEVENT_H 1

#include "DataFormats/Math/interface/Vector3D.h"
#include "ZShape/Base/interface/ZShapeElectron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <map>
#include <string>
#include <ostream>

/** \class ZShapeEvent
  *  
  * $Date: 2008/12/10 14:15:52 $
  * $Revision: 1.3 $
  * \author J. Mans - Minnesota
  */
class ZShapeEvent {
public:
  ZShapeEvent();
  // primary vertex
  ::math::XYZPoint vtx_;
  int n_elec;
  int n_TLelec;

  ZShapeElectron& elec(int i) { return ((i==0)?(*e1_):(*e2_)); }
  const ZShapeElectron& elec(int i) const { return ((i==0)?(*e1_):(*e2_)); }
  reco::GenParticle& elecTreeLevel(int i) { return ((i==0)?(*eTL1_):(*eTL2_)); }
  const reco::GenParticle& elecTreeLevel(int i) const { return ((i==0)?(*eTL1_):(*eTL2_)); }

  virtual ~ZShapeEvent();
  virtual void clear();
  virtual void dump(std::ostream& s) const;
  virtual void dump() const;
private:
  virtual ZShapeElectron* allocateElectron();
  virtual reco::GenParticle* allocateTreeLevelElectron();
  ZShapeElectron* e1_;
  ZShapeElectron* e2_;
  reco::GenParticle * eTL1_;
  reco::GenParticle * eTL2_;
};

#endif
