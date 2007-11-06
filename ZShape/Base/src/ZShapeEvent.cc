#include "ZShape/Base/interface/ZShapeEvent.h"
#include <iostream>

void ZShapeEvent::dump() const {
  dump(std::cout);
}
void ZShapeEvent::dump(std::ostream& s) const {
  s << "Vertex: " << z0_.x() << ", " << z0_.y() << ", " << z0_.z() << std::endl;
  for (int i=0; i<2; i++) {
    s << "Electron " << i+1 << ": (eta,phi,pt) " << elec_[i].p4_.eta() << ", " << elec_[i].p4_.phi() << ", " << elec_[i].p4_.Pt() << std::endl;
    for (std::map<std::string,bool>::const_iterator j=elec_[i].cutMap_.begin(); j!=elec_[i].cutMap_.end(); ++j) {
      s << "  \"" << j->first << "\" : " << ((j->second)?("pass"):("fail")) << std::endl;
    }
  }
}

void ZShapeEvent::clear() {
  elec_[0].cutMap_.clear();
  elec_[1].cutMap_.clear();
  n_elec=0;
}

bool ZShapeEvent::ElectronStruct::passed(const std::string& cutName) const { // false if not present!
  std::map<std::string,bool>::const_iterator i=cutMap_.find(cutName);
  return ((i==cutMap_.end())?(false):(i->second));
}

void ZShapeEvent::ElectronStruct::cutResult(const std::string& cutName, bool didpass) {
  cutMap_[cutName]=didpass;
}
