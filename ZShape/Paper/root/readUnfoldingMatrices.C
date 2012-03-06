//#include <iostream.h> 
//#include <stdlib.h>


// In this example, you pick up an unfolding matrix (made with makeUnfoldingMatrices.C) from a file (unfoldingMatrix_theOutPut.root)
// and apply to a rapidity distribution you get from an effAcc file (effAccSource.root)

// to execute this example
//      wget http://homepages.spa.umn.edu/~franzoni/reps/11/jan12/effAccSource.root
//      (getting this file is only needed to run the example)
//      pb-d-128-141-34-30:2011-01-12-migMatrixReleaseFRESH franzoni$ root -l 
//      root [0] .L readUnfoldingMatrices.C                                                  
//      root [1] readUnfoldingMatrices("unfoldingMatrix_theOutPut.root","effAccSource.root")




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
    TMatrixD * theUnfoldingMatrix = theunfoldingMatrixInputFile->Get("unsmearMatrices/unfoldingMatrixTotal");
    

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


    // now pull up a test ECAL-HF histogram 
    definition=std::string("ECAL80-HF");
    std::string cut("C08-m(70,110)");
    plot=std::string("Z0_Y"); // bin migration histogram
    plotName = prefix + definition;
    plotName = plotName + std::string("/");    plotName = plotName + cut;
    plotName = plotName + std::string("/");    plotName = plotName + plot;
    std::cout << "\tUsing ECAL-HF plot:     "<< plotName << endl;
    TH1F* h_smearedRapidityECALHF = (TH1F*) theEffAccInpuntFile->Get(plotName.c_str());
    std::cout << "\tname h_smearedRapidityECALHF: " <<  h_smearedRapidityECALHF->GetName() <<  std::endl;

   // building total migration histogram, summing the  ECAL-ECAL and the ECAL-HF histograms
   TH1F * h_fullTest_smeared  = new TH1F("h_fullTest_smeared","h_fullTest_smeared", 
		h_smearedRapidityECALHF->GetNbinsX(),
		h_smearedRapidityECALHF->GetXaxis()->GetXmin(),
		h_smearedRapidityECALHF->GetXaxis()->GetXmax());
	
   h_fullTest_smeared->Add(h_smearedRapidityECALHF,h_smearedRapidityECALECAL);
   h_fullTest_smeared->GetXaxis()->SetTitle("Y_{RECO} for test");



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
		h_smearedRapidityECALHF->GetNbinsX(),
		h_smearedRapidityECALHF->GetXaxis()->GetXmin(),
		h_smearedRapidityECALHF->GetXaxis()->GetXmax());
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // turning the unsmeared vector into a histogram
  for(int bin=1; bin<=100; bin++)
  {
   h_RapidityUNSmeared->SetBinContent(bin,vectorRapidityUNSmeared[bin-1]);
   //std::cout << vectorRapidityUNSmeared[bin-1]<<std::endl;
  }

  h_RapidityUNSmeared->Draw();
  h_fullTest_smeared->SetLineColor(kRed);
  h_fullTest_smeared->Draw("same");

  std::cout << "\n\n\n\t If you see two plots overlayed, you've applied unfolding matrix succesfully" << std::endl;

  return 0;
}
