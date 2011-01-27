#ifndef ETAACC_INC
#define ETAACC_INC

//
// Simple class to cut on eta acceptance
//

#include "DataFormats/Math/interface/LorentzVector.h"
#include "ZShape/Base/interface/ZShapeElectron.h"

class EtaAcceptance {

 public:
  EtaAcceptance( void );

  //  bool isInAcceptance(const math::PtEtaPhiMLorentzVector& p) const;
  enum Zone { zone_EB, zone_EE, zone_EE_tracker, zone_EE_notracker, zone_EE_muon, zone_ECAL, zone_ECAL_muon, zone_MUON, zone_HF, zone_ANY, zone_ALL };
  bool isInAcceptance(double deteta, Zone zone) const;
  bool isInAcceptance(const ZShapeElectron& zs, Zone zone) const { return isInAcceptance(zs.detEta_,zone); }

 private:

}; 

#endif
