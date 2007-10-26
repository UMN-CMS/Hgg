#include "ZShape/EffAcc/interface/EtaAcceptance.h"

typedef math::XYZTLorentzVector XYZTLorentzVector;

// for reference
// $ROOTSYS/include/Math/GenVector/LorentzVector.h

EtaAcceptance::EtaAcceptance () {
}

bool EtaAcceptance::isInAcceptance(const RawParticle& p) {

  XYZTLorentzVector p1 = p.momentum() ;

  float eta = fabsf( p.eta() );
 
  //
  // acceptance selection on eta of particle
  if    (eta < 0.018 || 
	 (eta>0.423 && eta<0.461) ||
	 (eta>0.770 && eta<0.806) ||
	 (eta>1.127 && eta<1.163) ||
	 (eta>1.460 && eta<1.558) ||
	 (eta>2.850 && eta<3.100) || // to be verified for details
	 (eta>4.60) // to be verified for details
	 )
    {
      std::cout << "eta: " << eta << ", acceptance failed" << std::endl;
      return false;
    }
  else
    {
      //edm::LogWarning("EtaAcceptance") << " - not passed ";
      std::cout << "Eta: " << eta << " acceptance passed" << std::endl;
      return true;
    }    
}
