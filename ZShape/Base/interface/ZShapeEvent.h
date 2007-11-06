#ifndef ZSHAPE_BASE_ZSHAPEEVENT_H
#define ZSHAPE_BASE_ZSHAPEEVENT_H 1

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <map>
#include <string>
#include <ostream>

/** \class ZShapeEvent
  *  
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class ZShapeEvent {
public:
  // primary vertex
  ::math::XYZVector z0_;
  int n_elec;

  // information about each electron
  struct ElectronStruct {
    bool passed(const std::string& cutName) const; // false if not present!
    void cutResult(const std::string& cutName, bool didpass);
    // pass/fail for each possible cut
    std::map<std::string, bool> cutMap_;
    // four-vector of electron
    math::PtEtaPhiMLorentzVector p4_;
  } elec_[2];

  virtual ~ZShapeEvent() { }
  virtual void clear();
  virtual void dump(std::ostream& s) const;
  virtual void dump() const;

};

#endif
