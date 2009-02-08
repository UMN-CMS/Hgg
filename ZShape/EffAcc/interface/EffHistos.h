#ifndef EFFHS_INC
#define EFFHS_INC

#include <TH1.h>
#include <TH2.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

// holds histograms for Z and e+/- quantities
// filled to carachterize candidates' pupulation at every step along selection path

class EffHistos {

 public:
  void Book(TFileDirectory& tfd);
  void Fill(const ::math::PtEtaPhiMLorentzVector& e1, const ::math::PtEtaPhiMLorentzVector& e2, 
            const ::math::PtEtaPhiMLorentzVector& eTL1, const ::math::PtEtaPhiMLorentzVector& eTL2); 
  void WrapUp(void);

 private:
  TH1 *mZ_,*YZ_, *ptZ_;
  TH1 *e1eta_,*e2eta_,*e1phi_,*e2phi_,*e1pt_,*e2pt_;
  TH2 *mZ_Y_,*mZ_pt_,*pt_Y_,*e1eta_YZ_,*e2eta_YZ_,*e1eta_ptZ_,*e2eta_ptZ_,*e1eta_e2eta_,*YZTL_YZ_,*YZTL_YZ_matrix_;

}; 

#endif
