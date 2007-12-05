#include "ZShape/Base/interface/EtaAcceptance.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


// for reference

EtaAcceptance::EtaAcceptance () {
}


bool EtaAcceptance::isInAcceptance(const math::PtEtaPhiMLorentzVector& p, Zone z) const {
  float eta = fabsf( p.eta() );

  switch (z) {
  case (zone_EB) :
    return (eta<1.460) &&  // EB
      (eta>0.018) && // eta=0 crack
      (eta<0.423 || eta>0.461) && // module 1-2
      (eta<0.770 || eta>0.806) && // module 2-3
      (eta<1.127 || eta>1.163); // module 3-4
  case (zone_EE) :
    return (eta>1.558) && (eta<2.850);
  case (zone_HF) :
    return (eta>3.1) && (eta<4.60);
  case (zone_ECAL) :
    return isInAcceptance(p,zone_EB) || isInAcceptance(p,zone_EE);
  case (zone_ANY) :
    return isInAcceptance(p,zone_EB) || isInAcceptance(p,zone_EE) || isInAcceptance(p,zone_HF);
  default: 
    edm::LogWarning("ZShape") << "Unknown acceptance request: " << z;
    return false;
  }
}

bool EtaAcceptance::isInAcceptance(const math::PtEtaPhiMLorentzVector& p) const {
  return isInAcceptance(p,zone_ANY);
}
