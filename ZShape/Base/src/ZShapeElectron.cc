#include "ZShape/Base/interface/ZShapeElectron.h"
#include <iostream>

void ZShapeElectron::dump(std::ostream& s) const {
  s << "(eta,phi,pt) " << p4_.eta() << ", " << p4_.phi() << ", " << p4_.Pt() << std::endl;
  for (std::map<std::string,bool>::const_iterator j=cutMap_.begin(); j!=cutMap_.end(); ++j) {
    s << "  \"" << j->first << "\" : " << ((j->second)?("pass"):("fail")) << std::endl;
  }
}

void ZShapeElectron::clear() {
  cutMap_.clear();
}

bool ZShapeElectron::passed(const std::string& cutName) const { // false if not present!
  std::map<std::string,bool>::const_iterator i=cutMap_.find(cutName);
  return ((i==cutMap_.end())?(false):(i->second));
}

void ZShapeElectron::cutResult(const std::string& cutName, bool didpass) {
  cutMap_[cutName]=didpass;
}
