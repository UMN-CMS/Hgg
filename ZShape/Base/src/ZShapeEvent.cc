#include "ZShape/Base/interface/ZShapeEvent.h"
#include <iostream>

ZShapeElectron* ZShapeEvent::allocateElectron() {
  return new ZShapeElectron(); // vanilla
}

void ZShapeEvent::dump() const {
  dump(std::cout);
}
void ZShapeEvent::dump(std::ostream& s) const {
  s << "Vertex: " << z0_.x() << ", " << z0_.y() << ", " << z0_.z() << std::endl;
  for (int i=0; i<2; i++) {
    s << "Electron " << i+1;
    elec(i).dump(s);
  }
}

void ZShapeEvent::clear() {
  if (e1_==0) e1_=allocateElectron();
  if (e2_==0) e2_=allocateElectron();
  elec(0).clear();
  elec(1).clear();
  n_elec=0;
}

ZShapeEvent::ZShapeEvent() : e1_(0),e2_(0) {
}

ZShapeEvent::~ZShapeEvent() {
  if (e1_!=0) delete e1_;
  if (e2_!=0) delete e2_;
}

