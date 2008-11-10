#include "ZShape/ZFromData/interface/EffExtraHistos.h"

typedef math::XYZTLorentzVector XYZTLorentzVector;

// for reference
// $ROOTSYS/include/Math/GenVector/LorentzVector.h

void EffExtraHistos::Book(TFileDirectory& tdf) {

  // add here all extremes, tidily
  float pi       = 3.141593;
  float maxPt    = 400;
  float maxY     = 5.5;
  float minZmass = 50;
  float maxZmass = 140;

  DmZ_  = tdf.make<TH1F>("Z0_mass_Delta","Z0_mass_Delta;Delta m_{Z0} (GeV/c^{2})", 50, -60, 60);  
  DYZ_  = tdf.make<TH1F>("Z0_Y_Delta","Z0_Y_Delta;Y_{Z0}", int((maxY*2)*4), -maxY, maxY);  
  DptZ_ = tdf.make<TH1F>("Z0_Pt_Delta","Z0_Pt_Delta;p_{T,Z0}", 100, -50., 50.); 
  De1eta_ = tdf.make<TH1F>("e1_eta_Delta","e1_eta_Delta;#eta_{e1}", 100, -5, 5);  
  De2eta_ = tdf.make<TH1F>("e2_eta_Delta","e2_eta_Delta;#eta_{e2}", 100, -5, 5);
  De1phi_ = tdf.make<TH1F>("e1_phi_Delta","e1_phi_Delta;#phi_{e1}", 100, -pi, pi);  
  De2phi_ = tdf.make<TH1F>("e2_phi_Delta","e2_phi_Delta;#phi_{e2}", 100,  -pi, pi);  
  De1pt_  = tdf.make<TH1F>("e1_P_t_Delta","e1_P_t_Delta;p_{T,e1}", 200, -50., 50.);  
  De2pt_  = tdf.make<TH1F>("e2_P_t_Delta","e2_P_t_Delta;p_{T,e2}", 200, -50., 50.);

  MCmZ_  = tdf.make<TH1F>("Z0_mass_MC","Z0_mass_MC;m_{Z0} (GeV/c^{2})", 60, 35., 150.);  
  MCYZ_  = tdf.make<TH1F>("Z0_Y_MC","Z0_Y_MC;Y_{Z0}", int((maxY*2)*4), -maxY, maxY);  
  MCptZ_ = tdf.make<TH1F>("Z0_Pt_MC","Z0_Pt_MC;p_{T,Z0}", 100, 0., maxPt); 
  MCe1eta_ = tdf.make<TH1F>("e1_eta_MC","e1_eta_MC;#eta_{e1}", 100, -5, 5);  
  MCe2eta_ = tdf.make<TH1F>("e2_eta_MC","e2_eta_MC;#eta_{e2}", 100, -5, 5);
  MCe1phi_ = tdf.make<TH1F>("e1_phi_MC","e1_phi_MC;#phi_{e1}", 100, -pi, pi);  
  MCe2phi_ = tdf.make<TH1F>("e2_phi_MC","e2_phi_MC;#phi_{e2}", 100,  -pi, pi);  
  MCe1pt_  = tdf.make<TH1F>("e1_P_t_MC","e1_P_t_MC;p_{T,e1}", 200, 0., maxPt);  
  MCe2pt_  = tdf.make<TH1F>("e2_P_t_MC","e2_P_t_MC;p_{T,e2}", 200, 0., maxPt);
}

void EffExtraHistos::Fill(const ::math::PtEtaPhiMLorentzVector& e1, const ::math::PtEtaPhiMLorentzVector& e2, const ::math::PtEtaPhiMLorentzVector& em1, const ::math::PtEtaPhiMLorentzVector& em2) {

  XYZTLorentzVector p1(e1);
  XYZTLorentzVector p2(e2);

  XYZTLorentzVector pZ = p1 + p2 ;

  XYZTLorentzVector pm1(em1);
  XYZTLorentzVector pm2(em2);

  XYZTLorentzVector pmZ = pm1 + pm2 ;
  //  debug
//  std::cout << "this is e1:"  << p1 << std::endl;
//  std::cout << "this is e2:"  << p2 << std::endl;
//  std::cout << "this is Z:"  << pZ << std::endl;
//  float rapidity = 0.5 * log( ( pZ.E()+pZ.Pz() ) / ( pZ.E()-pZ.Pz() ) );  
//  std::cout << "rapidity calculated: " << rapidity  << " while from method: " << pZ.Rapidity() << std::endl;
  
  float zMass = sqrt ( pZ.Dot(pZ) );
  float zY    = pZ.Rapidity();
  float zPt   = pZ.Pt();

  float e1eta = p1.Eta();
  float e1phi = p1.Phi();
  float e1Pt  = p1.Pt();
  
  float e2eta = p2.Eta();
  float e2phi = p2.Phi();
  float e2Pt  = p2.Pt();

  float zmMass = sqrt ( pmZ.Dot(pmZ) );
  float zmY    = pmZ.Rapidity();
  float zmPt   = pmZ.Pt();

  float em1eta = pm1.Eta();
  float em1phi = pm1.Phi();
  float em1Pt  = pm1.Pt();
  
  float em2eta = pm2.Eta();
  float em2phi = pm2.Phi();
  float em2Pt  = pm2.Pt();

  DmZ_  -> Fill(zMass-zmMass);
  DYZ_  -> Fill(zY-zmY);
  DptZ_ -> Fill(zPt-zmPt);

  De1eta_ -> Fill(e1eta-em1eta);
  De1pt_  -> Fill(e1Pt-em1Pt);
  De1phi_ -> Fill(e1phi-em1phi);

  De2eta_ -> Fill(e2eta-em2eta);
  De2pt_  -> Fill(e2Pt-em2Pt);
  De2phi_ -> Fill(e2phi-em2phi);

  MCmZ_  -> Fill(zmMass);
  MCYZ_  -> Fill(zmY);
  MCptZ_ -> Fill(zmPt);

  MCe1eta_ -> Fill(em1eta);
  MCe1pt_  -> Fill(em1Pt);
  MCe1phi_ -> Fill(em1phi);

  MCe2eta_ -> Fill(em2eta);
  MCe2pt_  -> Fill(em2Pt);
  MCe2phi_ -> Fill(em2phi);


}
