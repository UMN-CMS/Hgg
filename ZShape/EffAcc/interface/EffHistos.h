#ifndef EFFHS_INC
#define EFFHS_INC

#include <TH1.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

// holds histograms for Z and e+/- quantities
// filled to carachterize candidates' pupulation at every step along selection path

class EffHistos {

 public:
  void Book(TFileDirectory& tfd);
  void Fill(const ::math::PtEtaPhiMLorentzVector& e1, const ::math::PtEtaPhiMLorentzVector& e2);

 private:
  TH1 *mZ_,*YZ_, *ptZ_;
  TH1 *e1eta_,*e2eta_,*e1phi_,*e2phi_,*e1pt_,*e2pt_;

}; 

#endif
