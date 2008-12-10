#ifndef EFFEXHS_INC
#define EFFEXHS_INC

#include <TH1.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

// holds histograms for Z and e+/- quantities
// filled to carachterize candidates' pupulation at every step along selection path

class EffExtraHistos {

 public:
  void Book(TFileDirectory& tfd);
  void Fill(const ::math::PtEtaPhiMLorentzVector& e1, const ::math::PtEtaPhiMLorentzVector& e2,const ::math::PtEtaPhiMLorentzVector& em1, const ::math::PtEtaPhiMLorentzVector& em2);

 private:
  TH1 *DmZ_,*DYZ_, *DptZ_;
  TH1 *De1eta_,*De2eta_,*De1phi_,*De2phi_,*De1pt_,*De2pt_;
  TH1 *MCmZ_,*MCYZ_, *MCptZ_;
  TH1 *MCe1eta_,*MCe2eta_,*MCe1phi_,*MCe2phi_,*MCe1pt_,*MCe2pt_;
}; 

#endif