//#include <iostream.h> 
//#include <stdlib.h>

// In this example, you pick up migration matrix from a file from the standard EffAcc (effAccSource.root),
// produce an unsmearing matrix and save it to a file (for usage with readUnfoldingMatrices.C)

// to execute this example
//     root -l
//root [0] .L makeUnfoldingMatrices.C 
//root [1] makeUnfoldingMatrices("effAccSource.root","unfoldingMatrix_theOutPut.root")



int makeUnfoldingMatrices(std::string effAccFileInputFile, std::string unfoldingMatrixOutPutFile) 
{

    gStyle->SetPaintTextFormat(".2f");

   std::cout << "\n\n\tGetting migration matrices from: " << effAccFileInputFile          << std::endl;
   std::cout << "\tCreating unfolding matrices in:  "     << unfoldingMatrixOutPutFile <<"\n\n"<< std::endl;

    TFile *theEffAccInpuntFile = new TFile(effAccFileInputFile.c_str());
    //theEffAccInpuntFile->ls();

    TFile *theunfoldingMatrixOutpuFile = new TFile(unfoldingMatrixOutPutFile.c_str(),"recreate");
    theunfoldingMatrixOutpuFile ->cd();

    gStyle->SetOptStat(110011);

    // now pick up the migration matrix
    // open file
    // TFile * _file1 = new TFile("smearingMtx.root ");

    // get hold of ECAL-ECAL 2d plot 
    std::string prefix("mcEff/");
    std::string definition("ECAL80-ECAL95");   
    std::string cut("C08-m(70,110)");
    std::string plot("YZTL_vs_YZ");            // bin migration histogram
    std::string plotName = prefix + definition;
    plotName = plotName + std::string("/");    plotName = plotName + cut;
    plotName = plotName + std::string("/");    plotName = plotName + plot;
    std::cout << "\n\n\tUnsing ECAL-ECAL plot : "<< plotName << endl;

    TH2F* histoMigrationEcalEcalhighStat = (TH2F*) theEffAccInpuntFile->Get(plotName.c_str());
    int normEcalEcal = histoMigrationEcalEcalhighStat -> Integral();
    histoMigrationEcalEcalhighStat ->GetXaxis()->SetTitle("Y_{ZtreeLevel}");
    histoMigrationEcalEcalhighStat ->GetYaxis()->SetTitle("Y_{ZtreeLevel}");


    // get hold of ECAL-HF 2d plot 
    definition=std::string("ECAL80-HF");
    std::string cut("C08-m(70,110)");
    plot=std::string("YZTL_vs_YZ"); // bin migration histogram
    plotName = prefix + definition;
    plotName = plotName + std::string("/");    plotName = plotName + cut;
    plotName = plotName + std::string("/");    plotName = plotName + plot;
    std::cout << "\tUsing ECAL-HF plot:     "<< plotName << endl;

    TH2F* histoMigrationEcalHFhighStat = (TH2F*) theEffAccInpuntFile->Get(plotName.c_str());
    int normEcalHF = histoMigrationEcalHFhighStat -> Integral();
    histoMigrationEcalHFhighStat ->GetXaxis()->SetTitle("Y_{ZtreeLevel}");
    histoMigrationEcalHFhighStat ->GetYaxis()->SetTitle("Y_{ZtreeLevel}");
   

    ////////////////////////////////////////////////////////////////////////////////
    // evaluating statistical error of migration matrix for ECAL-ECAL
    TH2F * histoMigrationEcalEcalError = new TH2F("histoMigrationEcalEcalError", "histoMigrationEcalEcalError", 
	histoMigrationEcalEcalhighStat->GetNbinsX(), histoMigrationEcalEcalhighStat->GetXaxis()->GetXmin(), histoMigrationEcalEcalhighStat->GetXaxis()->GetXmax(),
	histoMigrationEcalEcalhighStat->GetNbinsY(), histoMigrationEcalEcalhighStat->GetYaxis()->GetXmin(), histoMigrationEcalEcalhighStat->GetYaxis()->GetXmax());
	
   for(int xi=0; xi<histoMigrationEcalEcalhighStat->GetNbinsX(); xi++){
   for(int yi=0; yi<histoMigrationEcalEcalhighStat->GetNbinsY(); yi++){
	double value = histoMigrationEcalEcalhighStat->GetBinContent(xi,yi);
        if (fabs(value)>0.000001) histoMigrationEcalEcalError->SetBinContent(xi,yi,1/sqrt(value));
  } }


    ////////////////////////////////////////////////////////////////////////////////
    // evaluating statistical error of migration matrix for ECAL-HF
    TH2F * histoMigrationEcalHFError = new TH2F("histoMigrationEcalHFError", "histoMigrationEcalHFError", 
	histoMigrationEcalHFhighStat->GetNbinsX(), histoMigrationEcalHFhighStat->GetXaxis()->GetXmin(), histoMigrationEcalHFhighStat->GetXaxis()->GetXmax(),
	histoMigrationEcalHFhighStat->GetNbinsY(), histoMigrationEcalHFhighStat->GetYaxis()->GetXmin(), histoMigrationEcalHFhighStat->GetYaxis()->GetXmax());
	

   for(int xi=0; xi<histoMigrationEcalHFhighStat->GetNbinsX(); xi++){
   for(int yi=0; yi<histoMigrationEcalHFhighStat->GetNbinsY(); yi++){
	double value = histoMigrationEcalHFhighStat->GetBinContent(xi,yi);
        if (fabs(value)>0.000001) histoMigrationEcalHFError->SetBinContent(xi,yi,1/sqrt(value));
  } }


   /////////////////////////////////////////////////////////////////////////////////////////
   // checking that sizes and bins are consistent between the ECAL-ECAL and the ECAL-HF case

   int binsXMatrixEcalEcal = histoMigrationEcalEcalhighStat->GetNbinsX();
   int binsYMatrixEcalEcal = histoMigrationEcalEcalhighStat->GetNbinsY();
   float xminEcalEcal = histoMigrationEcalEcalhighStat->GetXaxis()->GetXmin();
   float xmaxEcalEcal = histoMigrationEcalEcalhighStat->GetXaxis()->GetXmax();
   float yMinEcalEcal = histoMigrationEcalEcalhighStat->GetYaxis()->GetXmin();
   float yMaxEcalEcal = histoMigrationEcalEcalhighStat->GetYaxis()->GetXmax();
     
   int binsXMatrixEcalHF = histoMigrationEcalHFhighStat->GetNbinsX();
   int binsYMatrixEcalHF = histoMigrationEcalHFhighStat->GetNbinsY();
   float xminEcalHF = histoMigrationEcalHFhighStat->GetXaxis()->GetXmin();
   float xmaxEcalHF = histoMigrationEcalHFhighStat->GetXaxis()->GetXmax();
   float yMinEcalHF = histoMigrationEcalHFhighStat->GetYaxis()->GetXmin();
   float yMaxEcalHF = histoMigrationEcalHFhighStat->GetYaxis()->GetXmax();
     
     if (
	binsXMatrixEcalEcal ==binsYMatrixEcalEcal && yMinEcalEcal==xminEcalEcal && xmaxEcalEcal==yMaxEcalEcal &&
	binsXMatrixEcalHF   ==binsYMatrixEcalHF &&   yMinEcalHF==xminEcalHF && xmaxEcalHF==yMaxEcalHF &&
	binsXMatrixEcalEcal ==binsXMatrixEcalHF &&   yMinEcalEcal==yMinEcalHF && xmaxEcalEcal==xmaxEcalHF)
       {std::cout<< "\nbin migration histogram are square (" << binsXMatrixEcalEcal << "X" << binsYMatrixEcalEcal << ") and the same for EcalEcal/EcalHF, as it should\n" << std::endl; 	}
     else
       {std::cout<< "\nbin migration histogram are not square or not the same for EcalEcal/EcalHF, returning.\n" << std::endl;
	 return;}



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


   // build verification histograms
   TH1D * h_RapidityTreeLevel = histoMigration->ProjectionY("Rapidity TreeLevel",0,-1,"");
   histoMigration->SetTitle("RapidityTreeLevel");
   TH1D * h_RapiditySmeared   = histoMigration->ProjectionX("Rapidity Smeared",0,-1,"");
   h_RapiditySmeared->SetTitle("RapiditySmeared");
   h_RapiditySmeared->SetLineColor(2);




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









    TCanvas * c3 = new TCanvas("c3","c3",1050,750);
    c3->Divide(2);
    histoMigrationEcalEcalhighStat->GetXaxis()->SetTitle("Y_{Z,reco}");
    histoMigrationEcalEcalhighStat->GetYaxis()->SetTitle("Y_{Z,true}");
    histoMigrationEcalEcalhighStat->SetTitle("ECAL-ECAL: migration");
    c3->cd(1);    histoMigrationEcalEcalhighStat->Draw("colz");
    histoMigrationEcalHFhighStat->GetXaxis()->SetTitle("Y_{Z,reco}");
    histoMigrationEcalHFhighStat->GetYaxis()->SetTitle("Y_{Z,true}");
    histoMigrationEcalHFhighStat->SetTitle("ECAL-HF: migration");
    c3->cd(2);    histoMigrationEcalHFhighStat->Draw("colz");
    c3->Print( (definition + std::string("--") + cut + std::string(".png") ).c_str() );


    TCanvas * c3bos = new TCanvas("c3bos","c3bos",1040,740);
    c3bos->Divide(2);
    c3bos->cd(1);
    histoMigrationEcalEcalError->SetStats(0);
    histoMigrationEcalEcalError->SetTitle("ECAL-ECAL: stat error");
    histoMigrationEcalEcalError->Draw("colz");
    histoMigrationEcalEcalError->GetXaxis()->SetTitle("Y_{Z,reco}");
    histoMigrationEcalEcalError->GetYaxis()->SetTitle("Y_{Z,true}");
    c3bos->cd(2);
    histoMigrationEcalHFError->GetXaxis()->SetTitle("Y_{Z,reco}");
    histoMigrationEcalHFError->GetYaxis()->SetTitle("Y_{Z,true}");
    histoMigrationEcalHFError->SetStats(0);
    histoMigrationEcalHFError->SetTitle("ECAL-HF: stat error");
    histoMigrationEcalHFError->Draw("colz");
   c3bos->Print( (std::string("migration-stat-error") + std::string(".png") ).c_str() );
   // http://ultrahigh.org/2007/08/20/making-pretty-root-color-palettes/

   TCanvas * c4 = new TCanvas("c4","c4",1050,750);
   c4->cd(1);  c4->SetLogz();     
   histoMigration->SetTitle("Migration histogram");
   histoMigration->SetTitle("Migration Total");
   histoMigration->Draw("colz");
   c4->Print( (std::string("migration-total") + std::string(".png") ).c_str() );


   TCanvas * c4q = new TCanvas("c4q","c4q",1050,750);
   c4q->cd();
   h_effect->SetLineColor(2);
   h_effect->Draw();

   TCanvas * c5 = new TCanvas("c5","c5",1050,750);
   c5->cd(1);  c5->SetLogz();   
   histoMatrix->GetXaxis()->SetTitle("Y_{Z,reco}");
   histoMatrix->GetYaxis()->SetTitle("Y_{Z,true}");
   histoMatrix->SetStats(0);
   histoMatrix->Draw("colz");
   //histoMatrix->SetTitle("Migration matrix (normalized)");
   c5->Print("migrationMatrix.png");


   TCanvas * c7pre = new TCanvas("c7pre","c7pre",1050,750);
   c7pre->cd();
   h_RapidityTreeLevel->SetLineWidth(2);//this is black for now
   h_RapidityTreeLevel->GetYaxis()->SetTitle("events");
   h_RapidityTreeLevel->SetTitle(0);
   h_RapidityTreeLevel->Draw("");//this is black for now
  
   //h_RapidityUNSmeared->SetLineStyle(4);
   //h_RapidityUNSmeared->SetLineWidth(2);//this is black for now
   //h_RapidityUNSmeared->SetLineColor(4);
   h_RapidityUNSmeared->SetMarkerColor(4);
   h_RapidityUNSmeared->SetMarkerStyle(20);
   h_RapidityUNSmeared->SetMarkerSize(0.7);
   h_RapidityUNSmeared->SetLineWidth(1);
   h_RapidityUNSmeared->Draw("sameAP");

   TH1F * h_RapiditySmearedDummy = new TH1F("a","b",10,0,10);
   h_RapiditySmearedDummy->SetMarkerColor(3);
   h_RapiditySmearedDummy->SetMarkerStyle(23);
   h_RapiditySmearedDummy->SetMarkerSize(1.5);
   TH1F * h_RapidityUNSmearedDummy = new TH1F("c","d",10,0,10);
   h_RapidityUNSmearedDummy->SetMarkerColor(4);
   h_RapidityUNSmearedDummy->SetMarkerStyle(20);
   h_RapidityUNSmearedDummy->SetMarkerSize(1.5);

   h_RapiditySmeared->SetMarkerColor(2);
   h_RapiditySmeared->SetMarkerStyle(23);
   h_RapiditySmeared->SetMarkerSize(0.75);
   h_RapiditySmeared->Draw("sameAP");

   TLegend* leg = new TLegend(0.1,0.7,0.30,0.9);
   leg->SetHeader("Z Rapidity");
   leg->AddEntry(h_RapidityTreeLevel,"Y_{Z} treeLevel","l");
   leg->AddEntry(h_RapiditySmeared,"Y_{Z} smeared","p");
   leg->AddEntry(h_RapidityUNSmeared,"Y_{Z} unsmeared","p");
   leg->Draw();

   ///////////////////////////////////////////////////////
   c7pre->Print("rapdity-Tree-Smear-Unsmear-Closure.png");
   ///////////////////////////////////////////////////////

   TH1D * h_effectUnsmearing = new TH1D("Y unsmearing","Y unsmearing",binsXMatrixEcalEcal,xminEcalEcal,xmaxEcalEcal); 
   h_effectUnsmearing->Add(h_RapidityUNSmeared, h_RapidityTreeLevel, 1, -1);
   h_effectUnsmearing->Divide(h_RapidityTreeLevel);
   
   return 0;
}
