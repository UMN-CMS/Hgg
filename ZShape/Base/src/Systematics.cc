#include "ZShape/Base/interface/Systematics.h"
#include "ZShape/Base/interface/EtaAcceptance.h"

namespace zshape {
  

  EnergyScale::EnergyScale(bool positive) : upshift_(positive) {
  }
   
  void EnergyScale::rescale(ZShapeElectron& electron) {
    EtaAcceptance accept;
    double scale=1.0;

    if (accept.isInAcceptance(electron,EtaAcceptance::zone_EB)) {
      if (upshift_) scale=1.01;
      else scale=0.99;
    } else if (accept.isInAcceptance(electron,EtaAcceptance::zone_EE)) {
      if (upshift_) scale=1.03;
      else scale=0.97;     
    } else if (accept.isInAcceptance(electron,EtaAcceptance::zone_HF)) {
      if (upshift_) scale=1.00; // for now!
      else scale=1.00; // for now!
    }

    electron.p4_= math::PtEtaPhiMLorentzVector(electron.p4_.Pt()*scale,
					       electron.p4_.eta(),
					       electron.p4_.phi(),
					       electron.p4_.M()*scale);

  }

}
