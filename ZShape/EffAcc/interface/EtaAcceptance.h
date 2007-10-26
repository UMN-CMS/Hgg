#ifndef ETAACC_INC
#define ETAACC_INC

//
// Simple class to cut on eta acceptance
//


#include "FastSimulation/Particle/interface/RawParticle.h"

class EtaAcceptance {

 public:
  EtaAcceptance( void );
  bool isInAcceptance( const RawParticle& p);

 private:

}; 

#endif
