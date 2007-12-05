#ifndef ZSHAPE_BASE_ZSHAPEEVENT_H
#define ZSHAPE_BASE_ZSHAPEEVENT_H 1

#include "DataFormats/Math/interface/Vector3D.h"
#include "ZShape/Base/interface/ZShapeElectron.h"
#include <map>
#include <string>
#include <ostream>

/** \class ZShapeEvent
  *  
  * $Date: 2007/11/06 18:18:24 $
  * $Revision: 1.1 $
  * \author J. Mans - Minnesota
  */
class ZShapeEvent {
public:
  ZShapeEvent();
  // primary vertex
  ::math::XYZVector z0_;
  int n_elec;

  ZShapeElectron& elec(int i) { return ((i==0)?(*e1_):(*e2_)); }
  const ZShapeElectron& elec(int i) const { return ((i==0)?(*e1_):(*e2_)); }

  virtual ~ZShapeEvent();
  virtual void clear();
  virtual void dump(std::ostream& s) const;
  virtual void dump() const;
private:
  virtual ZShapeElectron* allocateElectron();
  ZShapeElectron* e1_;
  ZShapeElectron* e2_;
};

#endif
