#include "ZShape/EffAcc/interface/EffHistos.h"
#include "ZShape/Base/interface/ZShapeBinning.h"

typedef math::XYZTLorentzVector XYZTLorentzVector;

// for reference
// $ROOTSYS/include/Math/GenVector/LorentzVector.h

void EffHistos::Book(TFileDirectory& tdf) {

  booked_=true;

  // add here all extremes, tidily
  const float pi       = 3.141593;
  const float maxPt    = 400;
  const float minZmass = 50;
  const float maxZmass = 150;
  const int   binsZmass  = 50;

  float detetaBins[104] = {-5., -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4., -3.9, -3.8, -3.7, -3.6,  //14
                       -3.5, -3.4, -3.3,  -3.2, -3.1, -3., //6
                       -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.89, -1.78, -1.67,  //13
                       -1.56, -1.4442, -1.35, -1.257, -1.163, -1.127,  -1.02, -0.913, -0.806, -0.770,  //10
                       -0.667, -.564, -0.461, -0.423, -0.32, -0.22, -0.12, -0.018, 0.018, 0.12, 0.22, 0.32,  0.423, 0.461, 0.564,0.667, //16
                       0.770, 0.806, 0.913, 1.02, 1.127, 1.163,  1.257, 1.35, 1.4442, //9
                       1.56, 1.67, 1.78, 1.89, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, //14
                       3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, //10
                       4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5. }; //11
 
  mZ_  = tdf.make<TH1F>("Z0_mass","Z0_mass;m_{Z} (GeV/c^{2})", binsZmass, minZmass, maxZmass);
  YZ_  = tdf.make<TH1F>("Z0_Y","Z0_Y;Y_{Z0}", zshape::y_bins, -zshape::y_max, zshape::y_max);
  ptZ_  = tdf.make<TH1F>("Z0_Pt","Z0_Pt;p_{T,Z0}", zshape::pt_bins, zshape::pt_binning);
  ptZmon_ = tdf.make<TH1F>("Z0_PtMon","Z0_PtMon;p_{T,Z0}", 200, 0, maxPt);

  mZTL_  = tdf.make<TH1F>("Z0_massTL","Z0_M Tree Level;m_{Z0} (GeV/c^{2})", 50, minZmass, maxZmass);
  YZTL_  = tdf.make<TH1F>("Z0_YTL","Z0_Y Tree Level;Y_{Z0}", zshape::y_bins, -zshape::y_max, zshape::y_max);
  ptZTLmon_ = tdf.make<TH1F>("Z0_PtTLMon","Z0_Pt Tree Level Mon;p_{T,Z0}", 200, 0, maxPt);

  e1eta_ = tdf.make<TH1F>("e1_eta","e1_eta;detector #eta_{e1}", 50, -5, 5);
  e2eta_ = tdf.make<TH1F>("e2_eta","e2_eta;detector #eta_{e2}", 50, -5, 5);
  eeta_ = tdf.make<TH1F>("e_eta","e_eta;detector #eta_{e}", 50, -5, 5);
  e1phi_ = tdf.make<TH1F>("e1_phi","e1_phi;#phi_{e1}", 50, -pi, pi);
  e2phi_ = tdf.make<TH1F>("e2_phi","e2_phi;#phi_{e2}", 50,  -pi, pi);
  ephi_ = tdf.make<TH1F>("e_phi","e_phi;#phi_{e}", 50,  -pi, pi);
  e1pt_  = tdf.make<TH1F>("e1_P_t","e1_P_t;p_{T,e1}", 200, 0, maxPt);
  e2pt_  = tdf.make<TH1F>("e2_P_t","e2_P_t;p_{T,e2}", 200, 0, maxPt);
 
  //Now make the 2-D histograms

  mZ_Y_ = tdf.make<TH2F>("Z0_Y_v_mass","Z0_Y_v_mass;Y_{Z0};m_{Z0}", zshape::y_bins, -zshape::y_max, zshape::y_max, binsZmass, minZmass, maxZmass);
  pt_Y_ = tdf.make<TH2F>("Z0_Y_v_pt","Z0_Y_v_pt;Y_{Z0};p_{T,Z0}", zshape::y_bins, -zshape::y_max, zshape::y_max, zshape::pt_bins, zshape::pt_binning);
  mZ_pt_ = tdf.make<TH2F>("Z0_pt_v_mass","Z0_pt_v_mass;m_{Z0};p_{T,Z0}",binsZmass, minZmass, maxZmass, zshape::pt_bins,zshape::pt_binning);

  e1eta_YZ_    = tdf.make<TH2F>("e1_eta_vs_Z0Y","e1_eta_vs_Z0Y;detector #eta_{e1};Y_{Z0}", 50, -5, 5,zshape::y_bins, -zshape::y_max, zshape::y_max );
  e2eta_YZ_    = tdf.make<TH2F>("e2_eta_vs_Z0Y","e2_eta_vs_Z0Y;detector #eta_{e2};Y_{Z0}", 50, -5, 5,zshape::y_bins, -zshape::y_max, zshape::y_max );
  e1eta_ptZ_   = tdf.make<TH2F>("e1_eta_vs_Ptz","e1_eta_vs_PtZ;detector #eta_{e1};p_{T,Z0}", 50, -5, 5, zshape::pt_bins, zshape::pt_binning);
  e2eta_ptZ_   = tdf.make<TH2F>("e2_eta_vs_Ptz","e2_eta_vs_PtZ;detector #eta_{e2};p_{T,Z0}", 50, -5, 5,zshape::pt_bins, zshape::pt_binning);
  e1eta_e2eta_ = tdf.make<TH2F>("e1_eta_vs_e2_eta","e1_eta_vs_e2_eta;detector #eta_{e1};detector #eta_{e2}", 50, -5, 5,50, -5, 5 );
 
  YZTL_YZ_       = tdf.make<TH2F>("YZTL_vs_YZ","YZTL_vs_YZ;Y_{Z};Y_{ZtreeLevel}",zshape::y_bins,-zshape::y_max,zshape::y_max,zshape::y_bins,-zshape::y_max,zshape::y_max); 
  YZTL_YZ_matrix_= tdf.make<TH2F>("YZTL_vs_YZ_matrix","YZTL_vs_YZ_matrix;Y_{Z};Y_{ZtreeLevel}",zshape::y_bins,-zshape::y_max,zshape::y_max,zshape::y_bins,-zshape::y_max,zshape::y_max); 

} 

void EffHistos::Fill(const  ZShapeElectron& e1, const  ZShapeElectron& e2, 
		     const ::math::PtEtaPhiMLorentzVector& eTL1, const ::math::PtEtaPhiMLorentzVector& eTL2, double wgt, bool doMC) { 

  if (!booked_) {
    std::cerr << "Attempt to fill without booking!\n";
    return;
  }

  
  XYZTLorentzVector p1(e1.p4_);
  XYZTLorentzVector p2(e2.p4_);

  XYZTLorentzVector pTL1(eTL1); 
  XYZTLorentzVector pTL2(eTL2); 

  XYZTLorentzVector pZ = p1 + p2 ;
  XYZTLorentzVector pTLZ = pTL1 + pTL2 ; 

  //  debug
  //  std::cout << "this is e1:"  << p1 << std::endl;
  //  std::cout << "this is e2:"  << p2 << std::endl;
  //  std::cout << "this is Z:"  << pZ << std::endl;
  //  float rapidity = 0.5 * log( ( pZ.E()+pZ.Pz() ) / ( pZ.E()-pZ.Pz() ) );  
  //  std::cout << "rapidity calculated: " << rapidity  << " while from method: " << pZ.Rapidity() << std::endl;
  //   rapidity = 0.5 * log( ( pTLZ.E()+pTLZ.Pz() ) / ( pTLZ.E()-pTLZ.Pz() ) );   
  //   std::cout << "rapidity calculated: " << rapidity  << " while from method: " << pTLZ.Rapidity() << std::endl;

  float zMass = sqrt ( pZ.Dot(pZ) );
  float zY    = pZ.Rapidity();
  float zPt   = pZ.Pt();

  float e1eta = e1.detEta_;
  float e1phi = p1.Phi();
  float e1Pt  = p1.Pt();
  
  float e2eta = e2.detEta_;
  float e2phi = p2.Phi();
  float e2Pt  = p2.Pt();

  float zTLMass = sqrt ( pTLZ.Dot(pTLZ) );
  float zTLY    = pTLZ.Rapidity();
  float zTLPt   = pTLZ.Pt();
  

  mZ_  -> Fill(zMass,wgt);
  YZ_  -> Fill(zY,wgt);
  ptZ_ -> Fill(zPt,wgt);
  ptZmon_ -> Fill(zPt,wgt);

  if (doMC) 
    {
      mZTL_  -> Fill(zTLMass,wgt);
      YZTL_  -> Fill(zTLY,wgt);
      ptZTLmon_ -> Fill(zTLPt,wgt);
    }
  e1eta_ -> Fill(e1eta,wgt);
  e1pt_  -> Fill(e1Pt,wgt);
  e1phi_ -> Fill(e1phi,wgt);

  e2eta_ -> Fill(e2eta,wgt);
  e2pt_  -> Fill(e2Pt,wgt);
  e2phi_ -> Fill(e2phi,wgt);

  eeta_ -> Fill(e1eta,wgt);
  eeta_ -> Fill(e2eta,wgt);
  ephi_ -> Fill(e1phi,wgt);
  ephi_ -> Fill(e2phi,wgt);

  
  //Now Fill the 2-D Histograms
  mZ_Y_->Fill(zY,zMass,wgt);
  mZ_pt_->Fill(zPt,zMass,wgt);
  pt_Y_->Fill(zY,zPt,wgt);
  
  e1eta_YZ_  -> Fill(e1eta,zY,wgt);
  e2eta_YZ_  -> Fill(e2eta,zY,wgt);
  e1eta_ptZ_ -> Fill(e1eta,zPt,wgt);
  e2eta_ptZ_ -> Fill(e2eta,zPt,wgt);
  
  e1eta_e2eta_ -> Fill(e1eta,e2eta,wgt);

  // tree level rapidity VS fast reco-ed rapidity 
  if (doMC) YZTL_YZ_ -> Fill(zY,zTLY,wgt); 
}


void EffHistos::WrapUp(){ 

  if( YZTL_YZ_==0 &&  YZTL_YZ_matrix_==0) return; 

  if (!booked_) {
    std::cerr << "Attempt to finish without booking!\n";
    return;
  }
 
  int numBinX =  YZTL_YZ_ ->GetNbinsX(); 
  int numBinY =  YZTL_YZ_ ->GetNbinsY(); 
   
  // filling migration matricex with by normalixing migration histograms 
  float integral; 
  for (int binY=1; binY<=numBinY; binY++) 
    { 
      integral=0; 
      for (int binX=1; binX<=numBinX; binX++) 
	{   integral += YZTL_YZ_->GetBinContent(binX,binY); } 
     
      if (integral==0) continue; 
      for (int binX=1; binX<=numBinX; binX++) 
	{   YZTL_YZ_matrix_->SetBinContent(binX,binY,   YZTL_YZ_->GetBinContent(binX,binY) / integral); } 
 
    } 
  
  // debug 
  //    std::cout << "\n num entries HIS: " << numBinX << " " << numBinY << " " << YZTL_YZ_->GetEntries() << std::endl; 
  //    std::cout << "\n num entries MATR: " << YZTL_YZ_matrix_->GetNbinsX() << " " << YZTL_YZ_matrix_->GetNbinsY()  
  //        << " " << YZTL_YZ_matrix_->GetEntries() << std::endl; 
  
} 
