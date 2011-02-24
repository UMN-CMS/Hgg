#include "ZShape/Base/interface/EtaAcceptance.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


// for reference

EtaAcceptance::EtaAcceptance () {
}


bool EtaAcceptance::isInAcceptance(double p, Zone z) const {
  float eta = fabsf(p);

  switch (z) {
  case (zone_EB) :
    return (eta<1.4442); //&&  // EB
    //(eta>0.018) && // eta=0 crack
    //(eta<0.423 || eta>0.461) && // module 1-2
    // (eta<0.770 || eta>0.806) && // module 2-3
    // (eta<1.127 || eta>1.163); // module 3-4
  case (zone_EE) :
    return (eta>1.566) && (eta<2.850);
  case (zone_EE_tracker) :
    return (eta>1.566) && (eta<2.5);
  case (zone_EE_notracker) :
    return (eta>=2.5) && (eta<2.850);
  case (zone_EE_muon) : // joint acceptance with muons
    return (eta>1.566) && (eta<2.10); 
  case (zone_MUON) : // only acceptance with muons
    return (eta<2.10); 
  case (zone_HF) :
    return (eta>3.10) && (eta<4.60);
  case (zone_ECAL) :
    return isInAcceptance(p,zone_EB) || isInAcceptance(p,zone_EE);
  case (zone_ECAL_muon) :
    return isInAcceptance(p,zone_EB) || isInAcceptance(p,zone_EE_muon);
  case (zone_ANY) :
    return isInAcceptance(p,zone_EB) || isInAcceptance(p,zone_EE) || isInAcceptance(p,zone_HF);
  case (zone_ALL) :
    return true;
  default: 
    edm::LogWarning("ZShape") << "Unknown acceptance request: " << z;
    return false;
  }
}
/*
bool EtaAcceptance::isInAcceptance(const ZShapeElectron& p) const {
  return isInAcceptance(p.detEta_,zone_ANY);
}
*/
