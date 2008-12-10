#ifndef ETAACC_INC
#define ETAACC_INC

//
// Simple class to cut on eta acceptance
//

#include "DataFormats/Math/interface/LorentzVector.h"

class EtaAcceptance {

 public:
  EtaAcceptance( void );

  bool isInAcceptance(const math::PtEtaPhiMLorentzVector& p) const;
  enum Zone { zone_EB, zone_EE, zone_EE_tracker, zone_EE_notracker, zone_ECAL, zone_HF, zone_ANY };
  bool isInAcceptance(const math::PtEtaPhiMLorentzVector& p, Zone zone) const;

 private:

}; 

#endif
