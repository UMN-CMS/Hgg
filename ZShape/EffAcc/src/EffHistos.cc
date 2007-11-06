#include "ZShape/EffAcc/interface/EffHistos.h"

typedef math::XYZTLorentzVector XYZTLorentzVector;

// for reference
// $ROOTSYS/include/Math/GenVector/LorentzVector.h

void EffHistos::Book() {

  // add here all extremes, tidily
  float pi       = 3.141593;
  float maxPt    = 600;
  float maxY     = 10;
  float maxZmass = 200;

  mZ_  = new TH1F("Z0_mass","Z0_mass", 100, 0, maxZmass);
  mZ_ ->GetXaxis()->SetTitle("m_{Z0} (GeV/c^{2}");

  YZ_  = new TH1F("Z0_Y","Z0_Y", 100, -maxY, maxY);
  YZ_ ->GetXaxis()->SetTitle("Y_{Z0}");

  ptZ_ = new TH1F("Z0_Pt","Z0_Pt", 100, 0, maxPt);
  ptZ_ ->GetXaxis()->SetTitle("p_{T,Z0}");
  


  e1eta_ = new TH1F("e1_eta","e1_eta", 100, -10, 10);
  e1eta_ ->GetXaxis()->SetTitle("#eta_{e1}");

  e2eta_ = new TH1F("e2_eta","e2_eta", 100, -10, 10);
  e2eta_ ->GetXaxis()->SetTitle("#eta_{e2}");

  e1phi_ = new TH1F("e1_phi","e1_phi", 100, -pi, pi);
  e1phi_ ->GetXaxis()->SetTitle("#phi_{e1}");

  e2phi_ = new TH1F("e2_phi","e2_phi", 100,  -pi, pi);
  e2phi_ ->GetXaxis()->SetTitle("#phi_{e2}");

  e1pt_  = new TH1F("e1_P_t","e1_P_t", 200, 0, maxPt);
  e1pt_ ->GetXaxis()->SetTitle("p_{T,e1}");

  e2pt_  = new TH1F("e2_P_t","e2_P_t", 200, 0, maxPt);
  e2pt_ ->GetXaxis()->SetTitle("p_{T,e2}");

}

void EffHistos::Fill(const ::math::PtEtaPhiMLorentzVector& e1, const ::math::PtEtaPhiMLorentzVector& e2) {

  XYZTLorentzVector p1(e1);
  XYZTLorentzVector p2(e2);

  XYZTLorentzVector pZ = p1 + p2 ;

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

  mZ_  -> Fill(zMass);
  YZ_  -> Fill(zY);
  ptZ_ -> Fill(zPt);

  e1eta_ -> Fill(e1eta);
  e1pt_  -> Fill(e1Pt);
  e1phi_ -> Fill(e1phi);

  e2eta_ -> Fill(e2eta);
  e2pt_  -> Fill(e2Pt);
  e2phi_ -> Fill(e2phi);

}
