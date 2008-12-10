#include "ZShape/EffAcc/interface/EffHistos.h"

typedef math::XYZTLorentzVector XYZTLorentzVector;

// for reference
// $ROOTSYS/include/Math/GenVector/LorentzVector.h

void EffHistos::Book(TFileDirectory& tdf) {

  // add here all extremes, tidily
  float pi       = 3.141593;
  float maxPt    = 400;
  float maxY     = 5.5;
  float minZmass = 50;
  float maxZmass = 140;
  const int yBinsPerUnitFine = 10;
  const int yBinsPerUnitCoarse = 4;

  mZ_  = tdf.make<TH1F>("Z0_mass","Z0_mass", 50, minZmass, maxZmass);
  mZ_ ->GetXaxis()->SetTitle("m_{Z0} (GeV/c^{2}");

  YZ_  = tdf.make<TH1F>("Z0_Y","Z0_Y", int((maxY*2)*yBinsPerUnitFine), -maxY, maxY);
  YZ_ ->GetXaxis()->SetTitle("Y_{Z0}");

  ptZ_ = tdf.make<TH1F>("Z0_Pt","Z0_Pt", 100, 0, maxPt);
  ptZ_ ->GetXaxis()->SetTitle("p_{T,Z0}");

  mZ_Y_ = tdf.make<TH2F>("Z0_Y_v_mass","Z0_Y_v_mass", int((maxY*2)*yBinsPerUnitCoarse), -maxY, maxY, 50, minZmass, maxZmass);
  mZ_Y_->GetXaxis()->SetTitle("Y_{Z0}");
  mZ_Y_->GetYaxis()->SetTitle("m_{Z0}");

  mZ_pt_ = tdf.make<TH2F>("Z0_pt_v_mass","Z0_pt_v_mass", 50,0, maxPt, 50, minZmass, maxZmass);
  mZ_pt_->GetXaxis()->SetTitle("p_{T,Z0}");
  mZ_pt_->GetYaxis()->SetTitle("m_{Z0}");
  
  e1eta_ = tdf.make<TH1F>("e1_eta","e1_eta", 100, -5, 5);
  e1eta_ ->GetXaxis()->SetTitle("#eta_{e1}");

  e2eta_ = tdf.make<TH1F>("e2_eta","e2_eta", 100, -5, 5);
  e2eta_ ->GetXaxis()->SetTitle("#eta_{e2}");

  e1phi_ = tdf.make<TH1F>("e1_phi","e1_phi", 100, -pi, pi);
  e1phi_ ->GetXaxis()->SetTitle("#phi_{e1}");

  e2phi_ = tdf.make<TH1F>("e2_phi","e2_phi", 100,  -pi, pi);
  e2phi_ ->GetXaxis()->SetTitle("#phi_{e2}");

  e1pt_  = tdf.make<TH1F>("e1_P_t","e1_P_t", 200, 0, maxPt*2);
  e1pt_ ->GetXaxis()->SetTitle("p_{T,e1}");

  e2pt_  = tdf.make<TH1F>("e2_P_t","e2_P_t", 200, 0, maxPt*2);
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

  mZ_Y_->Fill(zY,zMass);
  mZ_pt_->Fill(zPt,zMass);

  e1eta_ -> Fill(e1eta);
  e1pt_  -> Fill(e1Pt);
  e1phi_ -> Fill(e1phi);

  e2eta_ -> Fill(e2eta);
  e2pt_  -> Fill(e2Pt);
  e2phi_ -> Fill(e2phi);

}
