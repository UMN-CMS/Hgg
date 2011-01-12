#include "ZShape/Base/interface/Systematics.h"
#include "ZShape/Base/interface/EtaAcceptance.h"

namespace zshape {
  

  EnergyScale::EnergyScale(int ecal, int HF) : shiftECAL_(ecal), shiftHF_(HF) {
  }
   
  void EnergyScale::rescale(ZShapeElectron& electron) {
    EtaAcceptance accept;
    double scale=1.0;

    if (accept.isInAcceptance(electron,EtaAcceptance::zone_EB)) {
      if (shiftECAL_>0) scale=1.01;
      else if (shiftECAL_<0) scale=0.99;
    } else if (accept.isInAcceptance(electron,EtaAcceptance::zone_EE)) {
      if (shiftECAL_>0) scale=1.03;
      else if (shiftECAL_<0) scale=0.97;     
    } else if (accept.isInAcceptance(electron,EtaAcceptance::zone_HF)) {
      if (shiftHF_>0) scale=1.10; 
      else if (shiftHF_<0) scale=0.90; 
    }

    electron.p4_= math::PtEtaPhiMLorentzVector(electron.p4_.Pt()*scale,
					       electron.p4_.eta(),
					       electron.p4_.phi(),
					       electron.p4_.M()*scale);

  }

}
