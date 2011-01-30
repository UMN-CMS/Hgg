//#include <iostream.h> 
//#include <stdlib.h>

// In this example

#include <TDatime.h>

#include "tdrstyle.C"
#include "zrapidityStandard.C"

#include "readUnfoldingMatrices.C"

// note: this number is the number of parameters variation + 1 for the BASE  
#define numVariations 17
  
//int systematicsUnfolding(std::string effAccFileInputFile, std::string unfoldingMatrixOutPutFile) 
int systematicsUnfolding( ) 
{

  setTDRStyle();
  zrap_colors();

  gSystem->Load("libMatrix");
  //gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat(".2f");


   std::vector<std::string> allParamVariatios;


   allParamVariatios.push_back(std::string("MESSthisUp")); 
   //      allParamVariatios.push_back(std::string("BASE.zshape_effacc_smear")); //fail  // ********* the base needs to be the first file *********
      
//      
   allParamVariatios.push_back(std::string("a-EB-Vdown"));//ok

   allParamVariatios.push_back(std::string("a-EB-Vup"));//ok
   allParamVariatios.push_back(std::string("a-EE-Vdown"));//ok
   allParamVariatios.push_back(std::string("a-EE-Vup"));//ok
      

//        allParamVariatios.push_back(std::string("a-HF-Vdown"));//fail
   allParamVariatios.push_back(std::string("penco"));//ok


   //   allParamVariatios.push_back(std::string("a-HF-Vup"));//fail
       allParamVariatios.push_back(std::string("prenco"));//

     allParamVariatios.push_back(std::string("BASE")); //ok  // ********* the base needs to be the first file *********

      
         allParamVariatios.push_back(std::string("alpha-EB-Vdown"));//ok
         allParamVariatios.push_back(std::string("alpha-EB-Vup"));//ok
         allParamVariatios.push_back(std::string("alpha-EE-Vdown"));//ok
         allParamVariatios.push_back(std::string("alpha-EE-Vup"));//ok
         allParamVariatios.push_back(std::string("c-EB-Vdown")); // ok
         allParamVariatios.push_back(std::string("c-EB-Vup"));//ok
   
   
         //allParamVariatios.push_back(std::string("c-HF-Vdown")); // fail
   	         allParamVariatios.push_back(std::string("pinco")); // ok so-so
      
      //   allParamVariatios.push_back(std::string("c-HF-Vup")); // fail
             allParamVariatios.push_back(std::string("panco")); // ok
	     allParamVariatios.push_back(std::string("punco")); // ok
      
      std::string theDirectoryPath("effAcc-V00-05-03/");

   std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\t" 
	     << allParamVariatios.size() << " parameter variations definitions being considered "  << std::endl;


   if (numVariations != allParamVariatios.size()){
     std::cout << "\n\n\t ** Mismatch between allParamVariatios.size() : " << allParamVariatios.size() 
	       << " and numVariations: " << numVariations 
	       << " to be fixed! ==> \n bailing out.\n\n" << std::endl;
       return 0;
   }


   // open all the necessary files
   TFile * theVariationFiles[numVariations+1];
   for(int v=0; v<numVariations; v++ ){

     std::string theFileName = theDirectoryPath + allParamVariatios.at(v) + std::string(".py.hadded.root");
     theVariationFiles[v] = new TFile( theFileName.c_str(),"read");
     //std::cout << "\t opening string: " << theFileName << std::endl;

     if( theVariationFiles[v]->IsOpen() && (! theVariationFiles[v]->IsZombie()) )       {
       std::cout << "\t " << (v+1) << ".  \t" 	   //<< " for string: " <<  theFileName << "\t"
		 << theVariationFiles[v]->GetName()  	 << "\t\t is open." << std::endl;    }
     
     else{                                std::cout << "\t problems opening: " << theFileName << std::endl; }

   } // loop on parameters variation





    // get hold of ECAL-ECAL 2d plot 
    std::string prefix("mcEff/");
    std::string definition("ECAL80-ECAL95");   
    std::string cut("C08-m(70,110)");
    std::string plot("YZTL_vs_YZ");            // bin migration histogram
    std::string plotNameECALECAL = prefix + definition;
    plotNameECALECAL = plotNameECALECAL + std::string("/");    plotNameECALECAL = plotNameECALECAL + cut;
    plotNameECALECAL = plotNameECALECAL + std::string("/");    plotNameECALECAL = plotNameECALECAL + plot;
    std::cout << "\n\n\tUnsing ECAL-ECAL the 2d migration plot : "<< plotNameECALECAL << endl;

    // get hold of ECAL-HF 2d plot 
    definition=std::string("ECAL80-HF");
    std::string cut("C08-m(70,110)");
    plot=std::string("YZTL_vs_YZ"); // bin migration histogram
    plotNameECALHF = prefix + definition;
    plotNameECALHF = plotNameECALHF + std::string("/");    plotNameECALHF = plotNameECALHF + cut;
    plotNameECALHF = plotNameECALHF + std::string("/");    plotNameECALHF = plotNameECALHF + plot;
    std::cout << "\tUsing ECAL-HF 2d migration plot :   "<< plotNameECALHF << endl;
    

    // std::cout << "is file valid: " << theVariationFiles[0] << std::endl;
    // theVariationFiles[0]->ls(); 
    TH2F * histoMigrationEcalEcal[numVariations];
    TH2F * histoMigrationEcalHF[numVariations];
    
    // get pointers to all migration histograms
    for(int v=0; v<numVariations; v++ ){
      // theVariationFiles[v]->ls();   // verify that access is possible
      histoMigrationEcalEcal[v] = (TH2F*) theVariationFiles[v]->Get(plotNameECALECAL.c_str());
      histoMigrationEcalHF[v]   = (TH2F*) theVariationFiles[v]->Get(plotNameECALHF.c_str());
      //std::cout << "\n\t" << histoMigrationEcalEcal[v]->GetName() << "\t and \t" <<  histoMigrationEcalHF[v]->GetName() 
      //<< " have been opened from the file:  " << theVariationFiles[v]->GetName() << std::endl;
    }
    
    // get bins and axis limits from the ECAL-ECAL base
    // implicitely assuming (W/O checking!)
    // that for all variations these are the same
    int binsXMatrix = histoMigrationEcalEcal[0]->GetNbinsX();
    int binsYMatrix = histoMigrationEcalEcal[0]->GetNbinsY();
    float xmin      = histoMigrationEcalEcal[0]->GetXaxis()->GetXmin();
    float xmax      = histoMigrationEcalEcal[0]->GetXaxis()->GetXmax();
    float yMin      = histoMigrationEcalEcal[0]->GetYaxis()->GetXmin();
    float yMax      = histoMigrationEcalEcal[0]->GetYaxis()->GetXmax();

    TH2F *     histoMigrationTotal[numVariations];    
    TH2F *     matrixMigrationTotal[numVariations];    
    double     migrationArray[10000];    //array reordering matrixMigrationTotal, to construct TMatrix  
    TMatrixD * migrationMatrix[numVariations];
    TMatrixD * unsmearMatrix[numVariations];

    //////////////////////////////////////////////////////////////////////////
    for(int v=0; v<numVariations; v++ ){
      // adding the total migration matrix
      std::string histoMigrationTotalName = std::string("migTotal") + allParamVariatios.at(v); 
      histoMigrationTotal[v] = new TH2F(histoMigrationTotalName.c_str(), 
					histoMigrationTotalName.c_str(),
					binsXMatrix,xmin,xmax,binsYMatrix,yMin,yMax);
      
      histoMigrationTotal[v]->Add(histoMigrationEcalEcal[v],histoMigrationEcalHF[v]);

      //      histoMigrationTotal[v]->Draw(); //GF

      if(1) std::cout << "\n\n\tECAL-ECAL : " << histoMigrationEcalEcal[v]->Integral()
		<< "\tECAL-HF : "       << histoMigrationEcalHF[v]->Integral()
		<< "\tECAL-total : "         << histoMigrationTotal[v]->Integral()
		<< std::endl;

      // now manipulating to get migration matrix from migration histogram
      std::string matrixMigrationName = std::string("matrixTotal") + allParamVariatios.at(v);
      matrixMigrationTotal[v] = new TH2F(matrixMigrationName.c_str(), 
					 matrixMigrationName.c_str(),
					 binsXMatrix,xmin,xmax,binsYMatrix,yMin,yMax); 
      if(0) std::cout << "\tmatrixMigrationTotal has name:  " << matrixMigrationTotal[v]->GetName() << std::endl;

      float integral=0;
      for(int binY=1; binY<=binsYMatrix; binY++)
	{
	  integral=0;
	  for(int binX=1; binX<=binsXMatrix; binX++)
	    { 	   integral+=histoMigrationTotal[v]->GetBinContent(binX,binY);        }// binX
	  
	  if (integral==0) continue;
	  // if(1) std::cout << "\t\tintegral of column is: " << integral << std::endl;

	  // histomatrix has rows all normalized to 1
	  // x-> reco; y-> true; => sum(y=const) has to be =1 , for probabilty conservation
	  // note that in modern version of EffxAcc code there's an _matrix histogram which is normalize..
	  for(int binX=1; binX<=binsXMatrix; binX++)
	    { 	 matrixMigrationTotal[v]->SetBinContent(binX,binY, histoMigrationTotal[v]->GetBinContent(binX,binY)/integral );        }// binX

	}// binY
      
      //matrixMigrationTotal[v]->Draw("colz"); //GF

      // resetting migrationArray
      for(int t=0; t<10000; t++) {migrationArray[t]=0;}
      int counter(0);

      // now reordering the content of the migration matrix into migrationArray
      for (int u=0; u<matrixMigrationTotal[v]->GetNbinsX(); u++){//rows
	for (int z=0; z<matrixMigrationTotal[v]->GetNbinsX(); z++){//col
	  
	  migrationArray[counter]=matrixMigrationTotal[v]->GetBinContent(u+1,z+1);
  
	  // for matrix to be invertible, no diagonal elements can be 0
	  // for now set, in all rows which are empty, the diagonal element to 1
	  if(z==u &&  (matrixMigrationTotal[v]->GetBinContent(u+1,z+1)==0) ) migrationArray[counter]=1;

	  if(0 &&  migrationArray[counter]!=0 )  std::cout << u << " " << z << " " << counter << " " << migrationArray[counter] <<std::endl;  
	  // till here verified that matrix is is the same as makeUnfoldingMatrices - GF
	  counter++;
	}//col
      }//rows      

      std::string matrixMigrationName = std::string("migration") + allParamVariatios.at(v);
      migrationMatrix[v] = new TMatrixD(0, (matrixMigrationTotal[v]->GetNbinsX()-1), 0, (matrixMigrationTotal[v]->GetNbinsX()-1),migrationArray,matrixMigrationName.c_str());
      std::string inverseMigrationName = std::string("inverse") + allParamVariatios.at(v);
      unsmearMatrix[v] =  (TMatrixD*) migrationMatrix[v]->Clone( inverseMigrationName.c_str() );

      //if( (*unsmearMatrix[v]) == (*migrationMatrix[v] ) ) std::cout << "\tthey're the same" << std::endl;
      //else                                                std::cout << "\tthey're NOT the same" << std::endl;

      Double_t determinant(0);
      // do the actual matrix inversion!
      unsmearMatrix[v]->Invert(&determinant);

      TMatrixD isThisUnity( (*migrationMatrix[v])  ,TMatrixD::kMult, (*unsmearMatrix[v]) );

      //if( (*unsmearMatrix[v]) == (*migrationMatrix[v] ) ) std::cout << "\tthey're the same" << std::endl;
      //else                                                std::cout << "\tthey're NO LONGER the same" << std::endl;
      
      // create reference to the diagonal and set the diagonal to 0
      TMatrixDDiag diag2(isThisUnity); diag2=0;
      // in order to get the largest off-diagonal term
      const Double_t U1_max_offdiag = (isThisUnity.Abs()).Max();
      cout << "\n\n\n\t  Verification that inversion of matrix " << v << " was ok: maximum off-diagonal A * A^-1 =  " << U1_max_offdiag << " is it SMALL?"<< endl;
      cout << "\t  Determinant          = " << determinant << "\n\n"<< std::endl;

    }// end loop over parameters variations
    

    


    TH1D * h_RapidityTreeLevel[numVariations];
    TH1D * h_RapiditySmeared[numVariations];
    TH1D * h_RapidityUnfolded[numVariations];
    TCanvas* theControlCanvas[numVariations];

    for(int v=0; v<numVariations; v++ ){
      // build verification histograms
      h_RapidityTreeLevel[v] = histoMigrationTotal[v]->ProjectionY("Rapidity TreeLevel aa",0,-1,"");
      std::string treeLevelHistoName = std::string("treeLevelHisto-") + allParamVariatios.at(v);
      h_RapidityTreeLevel[v]->SetTitle(treeLevelHistoName.c_str());
      
      h_RapiditySmeared[v]   = histoMigrationTotal[v]->ProjectionX("Rapidity Smeared aa",0,-1,"");
      std::string smearedHistoName = std::string("smearedHisto") + allParamVariatios.at(v);
      h_RapiditySmeared[v]->SetTitle(smearedHistoName.c_str());

      /////////////////////////////////
      TMatrixD isThisUnityBis( (*migrationMatrix[v])  ,TMatrixD::kMult, (*unsmearMatrix[v]) );
      
      // create reference to the diagonal and set the diagonal to 0
      TMatrixDDiag diag2Bis(isThisUnityBis); diag2Bis=0;
      // in order to get the largest off-diagonal term
      const Double_t U1_max_offdiagBis = (isThisUnityBis.Abs()).Max();
      //cout << "\n\n\n\t  Verification **BIS** that inversion of matrix " << v << " was ok: maximum off-diagonal A * A^-1 =  " << U1_max_offdiagBis << " is it SMALL?"<< endl;
      //cout << "\t  Determinant          = " << determinant << "\n\n"<< std::endl;
      /////////////////////////////////
      // this is to verify that the unfolding matrix is the inverse of the migration matrix
      


      if(1){
	Double_t arrayRapiditySmeared[100]; // this needs be hardcoded...
	for(int bin=1; bin<=100; bin++)
	  {
	    arrayRapiditySmeared[bin-1]=h_RapiditySmeared[v]->GetBinContent(bin);
	  }
	
	TVectorD vectorRapiditySmeared; vectorRapiditySmeared.Use(100,arrayRapiditySmeared);
	// this is the multiplication for clusure test
	TVectorD vectorRapidityUNSmeared = (*unsmearMatrix[v]) * vectorRapiditySmeared;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::string nameUnfolded = std::string("Rapidity Unsmeared") +  allParamVariatios.at(v) ;
	h_RapidityUnfolded[v] = new TH1D(nameUnfolded.c_str(), nameUnfolded.c_str(), 
					 h_RapiditySmeared[v]->GetNbinsX(),
					 h_RapiditySmeared[v]->GetXaxis()->GetXmin(),
					 h_RapiditySmeared[v]->GetXaxis()->GetXmax());
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// turning the unsmeared vector into a histogram
	for(int bin=1; bin<=100; bin++)
	  {
	    h_RapidityUnfolded[v]->SetBinContent(bin,vectorRapidityUNSmeared[bin-1]);
	    //std::cout << vectorRapidityUNSmeared[bin-1]<<std::endl;
	  }

      } // end of artificial if for debugging!
      // this is to verify that the multiplication given the inverse matrix is the same as within unfold()



      //h_RapiditySmeared[v]->SetLineColor(1);
      //h_RapiditySmeared[v]->Draw();
      if(0) h_RapidityUnfolded[v] = (TH1D *) unfold(h_RapiditySmeared[v], unsmearMatrix[v]);

      if (1) std::cout << "h_RapiditySmeared integral: " << h_RapiditySmeared[v]->Integral()
		       << "\t h_RapidityTreeLevel integral: " << h_RapidityTreeLevel[v]->Integral() << "\t " <<  h_RapidityTreeLevel[v]->GetMaximum()
		       << "\t h_RapidityUnfolded integral: " <<  h_RapidityUnfolded[v]->Integral() << "\t " <<  h_RapidityUnfolded[v]->GetMaximum()
		       << std::endl;

      std::string theCanvasName = std::string("canvas-") + allParamVariatios.at(v);
      theControlCanvas [v] = new TCanvas(theCanvasName.c_str(),theCanvasName.c_str(),25*v, 25*v, 1050+25*v,750+25*v);
      theControlCanvas [v]->cd();


      h_RapidityTreeLevel[v]->Draw(""); // this (from projection of a 2d) is verified to be the same as getting histos directly

      h_RapidityUnfolded[v]->SetLineColor(kBlue);
      h_RapidityUnfolded[v]->SetLineColor(kBlue);
      h_RapidityUnfolded[v]->SetMarkerStyle(kFullDotMedium);
      h_RapidityUnfolded[v]->SetMarkerColor(kBlue);
      h_RapidityUnfolded[v]->SetLineStyle(8);
      h_RapidityUnfolded[v]->Draw("same");

      //h_RapiditySmeared[v]->Draw("");

     if(0) { // this is to verify h_RapidityTreeLevel
    definition=std::string("ECAL80-ECAL95");
    std::string cut("C08-m(70,110)");
    plot=std::string("Z0_YTL"); // bin migration histogram
    plotNameECALHF = prefix + definition;
    plotNameECALHF = plotNameECALHF + std::string("/");    plotNameECALHF = plotNameECALHF + cut;
    plotNameECALHF = plotNameECALHF + std::string("/");    plotNameECALHF = plotNameECALHF + plot;
    TH1F * theTreeLevelDirectECAL = theVariationFiles[v]->Get(plotNameECALHF.c_str()); 

    definition=std::string("ECAL80-HF");
    std::string cut("C08-m(70,110)");
    plot=std::string("Z0_YTL"); // bin migration histogram
    plotNameECALHF = prefix + definition;
    plotNameECALHF = plotNameECALHF + std::string("/");    plotNameECALHF = plotNameECALHF + cut;
    plotNameECALHF = plotNameECALHF + std::string("/");    plotNameECALHF = plotNameECALHF + plot;
    TH1F * theTreeLevelDirectHF = theVariationFiles[v]->Get(plotNameECALHF.c_str()); 

    TH1F * theTreeLevelDirectTotal = theTreeLevelDirectECAL->Clone("pinco");
    theTreeLevelDirectTotal->Add(theTreeLevelDirectHF,theTreeLevelDirectECAL);
    
    theTreeLevelDirectTotal->SetLineColor(kRed);
    theTreeLevelDirectTotal->SetLineStyle(8);
    theTreeLevelDirectTotal->Draw("same");
     }

    }

   /*
    TFile *theEffAccInpuntFile = new TFile(effAccFileInputFile.c_str());
    //theEffAccInpuntFile->ls();

    TFile *theunfoldingMatrixOutpuFile = new TFile(unfoldingMatrixOutPutFile.c_str(),"recreate");
    theunfoldingMatrixOutpuFile ->cd();

  // directory where one saves the original histograms (from EffAcc calculator) the unfolding matrices have been obtained from
  TDirectory *originalHistograms = theunfoldingMatrixOutpuFile->mkdir("originalHistograms");
  originalHistograms->cd();

    gStyle->SetOptStat(110011);



    // keep the option of using the whole statistcs as available in MC
    TH2F * histoMigrationEcalEcal = histoMigrationEcalEcalhighStat->Clone();
    TH2F * histoMigrationEcalHF = histoMigrationEcalHFhighStat->Clone();

   // building total migration histogram, summing the  ECAL-ECAL and the ECAL-HF histograms
   TH2F * histoMigration  = new TH2F("histoMigration","histoMigration",
	binsXMatrixEcalEcal,xminEcalEcal,xmaxEcalEcal,
	binsYMatrixEcalEcal,yMinEcalEcal,yMaxEcalEcal
	);
	
   histoMigration->Add(histoMigrationEcalEcalhighStat, histoMigrationEcalHFhighStat);
   histoMigration->GetXaxis()->SetTitle("Y_{ZtreeLevel}");
   histoMigration->GetYaxis()->SetTitle("Y_{ZtreeLevel}");

   TH2F * histoMigrationTotal_toBeSaved = histoMigration->Clone();
   histoMigrationTotal_toBeSaved->Write("originalHistoMigrationTotal");


   // build verification histograms
   TH1D * h_RapidityTreeLevel = histoMigration->ProjectionY("Rapidity TreeLevel",0,-1,"");
   histoMigration->SetTitle("RapidityTreeLevel");
   TH1D * h_RapidityTreeLevel_toSave =  h_RapidityTreeLevel->Clone();  h_RapidityTreeLevel_toSave->Write("original1dTreeLevelRapidity");

   TH1D * h_RapiditySmeared   = histoMigration->ProjectionX("Rapidity Smeared",0,-1,"");
   h_RapiditySmeared->SetTitle("RapiditySmeared");
   h_RapiditySmeared->SetLineColor(2);
   TH1D * h_RapiditySmeared_toSave =  h_RapiditySmeared->Clone();  h_RapiditySmeared_toSave->Write("original1dSmearedRapidity");

   TH1D * h_effect = new TH1D("Y migration","Y migration",binsXMatrixEcalEcal,xminEcalEcal,xmaxEcalEcal); 
   h_effect->Add(h_RapiditySmeared, h_RapidityTreeLevel, 1, -1);
   h_effect->Divide(h_RapidityTreeLevel);
   h_effect->SetTitle("(smeared-true)/true");
   h_effect->SetStats(0);


   TH2F * histoMatrix  = new TH2F("histoMatrix","Migration matrix",
	binsXMatrixEcalEcal,xminEcalEcal,xmaxEcalEcal,
	binsYMatrixEcalEcal,yMinEcalEcal,yMaxEcalEcal
	);
   histoMatrix->GetXaxis()->SetTitle("Y_{ZtreeLevel}");
   histoMatrix->GetYaxis()->SetTitle("Y_{ZtreeLevel}");

    float integral=0;
    for(int binY=1; binY<=binsYMatrixEcalEcal; binY++)
    {
       integral=0;

       for(int binX=1; binX<=binsXMatrixEcalEcal; binX++)
       {
        integral+=histoMigration->GetBinContent(binX,binY);
       }// binX

       if (integral==0) continue;

       // histomatrix has rows all normalized to 1
       // x-> reco; y-> true; => sum(y=const) has to be =1 , for probabilty conservation
       // note that in modern version of EffxAcc code there's an _matrix histogram which is normalize..
       for(int binX=1; binX<=binsXMatrixEcalEcal; binX++)
       {
        histoMatrix->SetBinContent(binX,binY, 1.*histoMigration->GetBinContent(binX,binY)/integral );
       }// binX

    }// binY
   

   int counter(0);
   //   double migration[1936];//array reordering histoMatrix, to construct TMatrix
   // double migration[12100];    //array reordering histoMatrix, to construct TMatrix
   double migration[10000];    //array reordering histoMatrix, to construct TMatrix
                               // this has to be hard-coded
      for (int u=0; u<histoMatrix->GetNbinsX(); u++){//rows
     for (int v=0; v<histoMatrix->GetNbinsX(); v++){//col

       migration[counter]=histoMatrix->GetBinContent(u+1,v+1);
       if(0)  std::cout << u << " " << v << " " << migration[counter] <<std::endl;

       // for matrix to be invertible, no diagonal elements can be 0
       // for now set, in all rows which are empty, the diagonal element to 1
       if(u==v && histoMatrix->GetBinContent(u+1,v+1)==0) migration[counter]=1;
       counter++;

       }	
    }

  //http://root.cern.ch/root/html/tutorials/matrix/invertMatrix.C.html
  //this is the smearing matrix: it multiplies the unsmeared vector to give the smeared vector
  gSystem->Load("libMatrix");
  TMatrixD MigrationMatrix(0, histoMatrix->GetNbinsX()-1, 0, histoMatrix->GetNbinsX()-1, migration,"");

  // directory where one saves the matrices
  TDirectory *unsmearMatrices = theunfoldingMatrixOutpuFile->mkdir("unsmearMatrices");
  unsmearMatrices->cd();
  TMatrixD inverseMigrationMatrix = MigrationMatrix;

  Double_t determinant;
  inverseMigrationMatrix.Invert(&determinant);

  // this is the name I assign to the unfolding matrix for all Z candidates: ECAL-ECAL and ECAL-HF
  inverseMigrationMatrix.Write("unfoldingMatrixTotal");


  //checking whether the inversion is sound: multiply matrix and its inverse, anche check if unity
  //   TMatrixD isThisUnity = inverseMigrationMatrix * MigrationMatrix;
  TMatrixD isThisUnity(inverseMigrationMatrix,TMatrixD::kMult,MigrationMatrix);
  // create reference to the diagonal and set the diagonal to 0
  TMatrixDDiag diag1(isThisUnity); diag1 = 0.0;
  const Double_t U1_max_offdiag = (isThisUnity.Abs()).Max();
  cout << "\n\n\n  Verification that inversion was ok: maximum off-diagonal A * A^-1 = " << U1_max_offdiag << " is it SMALL?\n\n "<< endl;
  //cout << "  Determinant          = " << determinant <<endl;

  //now apply inverse matrix to smeared histogram, and check if it is the same as unsmeared
  //http://root.cern.ch/root/html/tutorials/matrix/solveLinear.C.html
  //h_RapidityTreeLevel
  //h_RapiditySmeared

  //Double_t arrayRapiditySmeared[110]; // this needs be hardcoded...
  Double_t arrayRapiditySmeared[100]; // this needs be hardcoded...
  //  for(int bin=1; bin<=44; bin++)  // was hard-coded, but this is not necessary
  for(int bin=1; bin<=binsXMatrixEcalEcal; bin++)
  {
   arrayRapiditySmeared[bin-1]=h_RapiditySmeared->GetBinContent(bin);
  }

  TVectorD vectorRapiditySmeared; vectorRapiditySmeared.Use(binsXMatrixEcalEcal,arrayRapiditySmeared);

  // this is the multiplication for clusure test
  TVectorD vectorRapidityUNSmeared = inverseMigrationMatrix * vectorRapiditySmeared;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1D * h_RapidityUNSmeared = new TH1D("Rapidity Unsmeared", "Rapidity Unsmeared", binsXMatrixEcalEcal, xminEcalEcal, xmaxEcalEcal);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // turning the unsmeared vector into a histogram
  for(int bin=1; bin<=binsXMatrixEcalEcal; bin++)
  {
   h_RapidityUNSmeared->SetBinContent(bin,vectorRapidityUNSmeared[bin-1]);
   //std::cout << vectorRapidityUNSmeared[bin-1]<<std::endl;
  }










   */



   return 0;
}
