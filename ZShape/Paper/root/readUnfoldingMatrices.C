#include <iostream> 
//#include <stdlib.h>
#include "TH1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TStyle.h"

// In this example, you pick up an unfolding matrix (made with makeUnfoldingMatrices.C) from a file (unfoldingMatrix_theOutPut.root)
// and apply to a rapidity distribution you get from an effAcc file (effAccSource.root)

// to execute this example
//      wget http://homepages.spa.umn.edu/~franzoni/reps/11/jan12/effAccSource.root
//      (getting this file is only needed to run the example)
//      pb-d-128-141-34-30:2011-01-12-migMatrixReleaseFRESH franzoni$ root -l 
//      root [0] .L readUnfoldingMatrices.C                                                  
//      root [1] readUnfoldingMatrices("unfoldingMatrix_theOutPut.root","effAccSource.root")



TH1* unfold(const TH1* source, const char* unfoldingMatrixIF)  {
  // load the matrix...

  TFile theunfoldingMatrixInputFile(unfoldingMatrixIF,"read");

  TMatrixD * theUnfoldingMatrix = (TMatrixD*)theunfoldingMatrixInputFile.Get("unsmearMatrices/unfoldingMatrixTotal");

  unfold(source, theUnfoldingMatrix);
}


TH1* unfold(const TH1* source, TMatrixD * theUnfoldingMatrix)  {
  Double_t arrayRapiditySmeared[100]; // this needs be hardcoded...
  for(int bin=1; bin<=100; bin++)
    {
      arrayRapiditySmeared[bin-1]=source->GetBinContent(bin);
      if (0) std::cout << source->GetBinContent(bin) << std::endl;
    }

  TVectorD vectorRapiditySmeared; vectorRapiditySmeared.Use(100,arrayRapiditySmeared);

  // this is the multiplication for clusure test
  TVectorD vectorRapidityUNSmeared = (*theUnfoldingMatrix) * vectorRapiditySmeared;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  char name[1024];
  sprintf(name,"%s_unfolded",source->GetName());
  TH1 * h_RapidityUNSmeared = (TH1*)(source->Clone(name));
  h_RapidityUNSmeared->SetDirectory(0);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // turning the unsmeared vector into a histogram
  for(int bin=1; bin<=100; bin++)
  {
   h_RapidityUNSmeared->SetBinContent(bin,vectorRapidityUNSmeared[bin-1]);
   if(fabs(vectorRapiditySmeared[bin-1]) > 1.e-9) h_RapidityUNSmeared->SetBinError(bin,vectorRapidityUNSmeared[bin-1]/vectorRapiditySmeared[bin-1]*source->GetBinError(bin));
   else                                           h_RapidityUNSmeared->SetBinError(bin,source->GetBinError(bin));
  }
  return h_RapidityUNSmeared;
}




int readUnfoldingMatrices(std::string unfoldingMatrixFileInputFile, std::string effAccFileForTest) 

{

    gStyle->SetOptStat(110011);
    gStyle->SetPaintTextFormat(".2f");

    std::cout << "\n\n\tReading test histogram  from: " << effAccFileForTest          << std::endl;
    std::cout << "\tReading unfolding matrices in:  "   << unfoldingMatrixFileInputFile <<"\n\n"<< std::endl;

    TFile *theEffAccInpuntFile = new TFile(effAccFileForTest.c_str(),"read");

    TFile *theunfoldingMatrixInputFile = new TFile(unfoldingMatrixFileInputFile.c_str(),"read");
    theunfoldingMatrixInputFile->cd();
    std::cout << "The file: " << theunfoldingMatrixInputFile->GetName() << " has been opened"<< std::endl;
    
    // this is the matrix you want to use for to unfold
    TMatrixD * theUnfoldingMatrix = (TMatrixD*)theunfoldingMatrixInputFile->Get("unsmearMatrices/unfoldingMatrixTotal");

    /*
    //////////GFfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    // now pull up a test ECAL-ECAL histogram 
    std::string prefix("mcEff/");
    std::string definition("ECAL80-ECAL95");
    std::string cut("C08-m(70,110)");
    std::string plot("Z0_Y");                  // bin migration histogram
    std::string plotName = prefix + definition;
    plotName = plotName + std::string("/");    plotName = plotName + cut;
    plotName = plotName + std::string("/");    plotName = plotName + plot;
    std::cout << "\n\n\tUnsing ECAL-ECAL plot : "<< plotName << endl;
    TH1F* h_smearedRapidityECALECAL = (TH1F*) theEffAccInpuntFile->Get(plotName.c_str());
    std::cout << "\tname h_smearedRapidityECALECAL: " << h_smearedRapidityECALECAL->GetName() << std::endl;
    plot=std::string("Z0_YTL");                  // tree level histogram
    plotName = prefix + definition + std::string("/") +  cut + std::string("/") + plot;
    TH1F* h_treeLevelRapidityECALECAL = (TH1F*) theEffAccInpuntFile->Get(plotName.c_str());

    std::cout <<  plotName << "    " << h_treeLevelRapidityECALECAL << std::endl;

    // now pull up a test ECAL-HF histogram 
    definition=std::string("ECAL80-HF");
    plot=std::string("Z0_Y"); // bin migration histogram
    plotName = prefix + definition;
    plotName = plotName + std::string("/");    plotName = plotName + cut;
    plotName = plotName + std::string("/");    plotName = plotName + plot;
    std::cout << "\tUsing ECAL-HF plot:     "<< plotName << endl;
    TH1F* h_smearedRapidityECALHF = (TH1F*) theEffAccInpuntFile->Get(plotName.c_str());
    std::cout << "\tname h_smearedRapidityECALHF: " <<  h_smearedRapidityECALHF->GetName() <<  std::endl;
    plot=std::string("Z0_YTL");                  // tree level histogram
    plotName = prefix + definition + std::string("/") +  cut + std::string("/") + plot;
    TH1F* h_treeLevelRapidityECALHF = (TH1F*) theEffAccInpuntFile->Get(plotName.c_str());

   // building total migration histogram, summing the  ECAL-ECAL and the ECAL-HF histograms
   TH1F * h_fullTest_smeared  = new TH1F("h_fullTest_smeared","h_fullTest_smeared", 
		h_smearedRapidityECALHF->GetNbinsX(),
		h_smearedRapidityECALHF->GetXaxis()->GetXmin(),
		h_smearedRapidityECALHF->GetXaxis()->GetXmax());
	
   h_fullTest_smeared->Add(h_smearedRapidityECALHF,h_smearedRapidityECALECAL);
   h_fullTest_smeared->GetXaxis()->SetTitle("Y_{RECO} for test");


   // building total tree level histogram, summing the  ECAL-ECAL and the ECAL-HF histograms
   TH1F * h_fullTest_treeLevel  = new TH1F("h_fullTest_TL","h_fullTest_TL", 
		h_smearedRapidityECALHF->GetNbinsX(),
		h_smearedRapidityECALHF->GetXaxis()->GetXmin(),
		h_smearedRapidityECALHF->GetXaxis()->GetXmax());
   
   h_fullTest_treeLevel->Add(h_treeLevelRapidityECALHF,h_treeLevelRapidityECALECAL);
   h_fullTest_treeLevel->GetXaxis()->SetTitle("Y_{RECO} for test");

   //////////GFfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    */

    //////////////////////////////////////////////////////
    // ... but to be consistent with the definitions 
    // of AN-10-367 one needs to use the YZTL_vs_YZ before Eff and Acc
    // are applied ==> move to the "All" folder
    // get hold of ECAL-HF 2d plot 
    std::string prefix("mcEff/");
    std::string definition=std::string("All");
    std::string cut("");
    std::string plot=std::string("Z0_Y"); // smeared Y before Eff&Acc
    std::string plotName = prefix + definition;
    plotName = plotName + std::string("/");    plotName = plotName + cut;
    plotName = plotName + std::string("/");    plotName = plotName + plot;
    std::cout << "\tUsing directly the total smeared Y plot:     "<< plotName << endl;
    TH2F * h_fullTest_smeared = (TH2F*) theEffAccInpuntFile->Get(plotName.c_str());

    plot=std::string("Z0_YTL"); // tree-level Y before Eff&Acc
    plotName = prefix + definition;
    plotName = plotName + std::string("/");    plotName = plotName + cut;
    plotName = plotName + std::string("/");    plotName = plotName + plot;
    std::cout << "\tComparing agains the total tree-Level Y plot:     "<< plotName << endl;
    TH2F * h_fullTest_treeLevel = (TH2F*) theEffAccInpuntFile->Get(plotName.c_str());



   //Double_t arrayRapiditySmeared[110]; // this needs be hardcoded...
  Double_t arrayRapiditySmeared[100]; // this needs be hardcoded...
  for(int bin=1; bin<=100; bin++)
  {
   arrayRapiditySmeared[bin-1]=h_fullTest_smeared->GetBinContent(bin);
  }

  TVectorD vectorRapiditySmeared; vectorRapiditySmeared.Use(100,arrayRapiditySmeared);

  // this is the multiplication for clusure test
  TVectorD vectorRapidityUNSmeared = (*theUnfoldingMatrix) * vectorRapiditySmeared;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1D * h_RapidityUNSmeared = new TH1D("Rapidity Unsmeared", "Rapidity Unsmeared", 
		h_fullTest_treeLevel->GetNbinsX(),
		h_fullTest_treeLevel->GetXaxis()->GetXmin(),
		h_fullTest_treeLevel->GetXaxis()->GetXmax());
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // turning the unsmeared vector into a histogram
  for(int bin=1; bin<=100; bin++)
  {
   h_RapidityUNSmeared->SetBinContent(bin,vectorRapidityUNSmeared[bin-1]);
   //std::cout << vectorRapidityUNSmeared[bin-1]<<std::endl;
  }

  gStyle->SetOptStat(0);
  TCanvas * theBaseCanvas = new TCanvas("check of closure","check of closure",200,200,1000,800);  
  theBaseCanvas->cd();
  h_RapidityUNSmeared->Draw();

  h_fullTest_smeared->SetLineColor(kGreen);
  h_fullTest_smeared->Draw("same");

  h_fullTest_treeLevel->SetLineColor(kRed);
  h_fullTest_treeLevel->SetLineStyle(8);
  h_fullTest_treeLevel->Draw("same");
  
  TLegend* leg = new TLegend(0.7,0.7,0.95,0.9);
  leg->SetHeader("closure: isn't it?");
  leg->AddEntry(h_fullTest_treeLevel,"tree level","lp");
  leg->AddEntry(h_RapidityUNSmeared,"unsmeared: does it overlap?","lp");
  leg->AddEntry(h_fullTest_smeared,"smeared","lp");
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.02);
  leg->Draw();  
  

  std::cout << "\n\n\n\t If you see tree plots and two are overlayed (gree and black), you've applied unfolding matrix succesfully" << std::endl;

  return 0;
}
