#ifndef EFFHS_INC
#define EFFHS_INC

#include <TH1.h>
#include "FastSimulation/Particle/interface/RawParticle.h"

class EffHistos {

 public:
  void Book();
  void Fill(const RawParticle& e1, const RawParticle& e2);

 private:
  TH1 *mZ_,*YZ_, *ptZ_;
  TH1 *e1eta_,*e2eta_,*e1phi_,*e2phi_,*e1pt_,*e2pt_;

}; 

#endif
