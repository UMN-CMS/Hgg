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

  mZ_  = tdf.make<TH1F>("Z0_mass","Z0_mass;m_{Z0} (GeV/c^{2})", 50, minZmass, maxZmass);
  YZ_  = tdf.make<TH1F>("Z0_Y","Z0_Y;Y_{Z0}", int((maxY*2)*yBinsPerUnitFine), -maxY, maxY);
  ptZ_ = tdf.make<TH1F>("Z0_Pt","Z0_Pt;p_{T,Z0}", 100, 0, maxPt);

  e1eta_ = tdf.make<TH1F>("e1_eta","e1_eta;#eta_{e1}", 100, -5, 5);
  e2eta_ = tdf.make<TH1F>("e2_eta","e2_eta;#eta_{e2}", 100, -5, 5);
  e1phi_ = tdf.make<TH1F>("e1_phi","e1_phi;#phi_{e1}", 100, -pi, pi);
  e2phi_ = tdf.make<TH1F>("e2_phi","e2_phi;#phi_{e2}", 100,  -pi, pi);
  e1pt_  = tdf.make<TH1F>("e1_P_t","e1_P_t;p_{T,e1}", 200, 0, maxPt*2);
  e2pt_  = tdf.make<TH1F>("e2_P_t","e2_P_t;p_{T,e2}", 200, 0, maxPt*2);

  //Now make the 2-D histograms

  mZ_Y_ = tdf.make<TH2F>("Z0_Y_v_mass","Z0_Y_v_mass;Y_{Z0};m_{Z0}", int((maxY*2)*yBinsPerUnitCoarse), -maxY, maxY, 50, minZmass, maxZmass);
  pt_Y_ = tdf.make<TH2F>("Z0_Y_v_pt","Z0_Y_v_pt;Y_{Z0};p_{T,Z0}", int((maxY*2)*yBinsPerUnitCoarse), -maxY, maxY, 50, 0, maxPt);
  mZ_pt_ = tdf.make<TH2F>("Z0_pt_v_mass","Z0_pt_v_mass;p_{T,Z0};m_{Z0}", 50,0, maxPt, 50, minZmass, maxZmass);

  e1eta_YZ_    = tdf.make<TH2F>("e1_eta_vs_Z0Y","e1_eta_vs_Z0Y;#eta_{e1};Y_{Z0}", 50, -5, 5,int((maxY*2)*yBinsPerUnitCoarse), -maxY, maxY );
  e2eta_YZ_    = tdf.make<TH2F>("e2_eta_vs_Z0Y","e2_eta_vs_Z0Y;#eta_{e2};Y_{Z0}", 50, -5, 5,int((maxY*2)*yBinsPerUnitCoarse), -maxY, maxY );
  e1eta_ptZ_   = tdf.make<TH2F>("e1_eta_vs_Ptz","e1_eta_vs_PtZ;#eta_{e1};p_{T,Z0}", 50, -5, 5,50, 0, maxPt );
  e2eta_ptZ_   = tdf.make<TH2F>("e2_eta_vs_Ptz","e2_eta_vs_PtZ;#eta_{e2};p_{T,Z0}", 50, -5, 5,50, 0, maxPt );
  e1eta_e2eta_ = tdf.make<TH2F>("e1_eta_vs_e2_eta","e1_eta_vs_e2_eta;#eta_{e1};#eta_{e2}", 50, -5, 5,50, -5, 5 );

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
  pt_Y_->Fill(zY,zPt);

  e1eta_ -> Fill(e1eta);
  e1pt_  -> Fill(e1Pt);
  e1phi_ -> Fill(e1phi);

  e2eta_ -> Fill(e2eta);
  e2pt_  -> Fill(e2Pt);
  e2phi_ -> Fill(e2phi);
  
  //Now Fill the 2-D Histograms
  mZ_Y_->Fill(zY,zMass);
  mZ_pt_->Fill(zPt,zMass);
  pt_Y_->Fill(zY,zPt);
  
  e1eta_YZ_  -> Fill(e1eta,zY);
  e2eta_YZ_  -> Fill(e2eta,zY);
  e1eta_ptZ_ -> Fill(e1eta,zPt);
  e2eta_ptZ_ -> Fill(e2eta,zPt);
  
  e1eta_e2eta_ -> Fill(e1eta,e2eta);

}
