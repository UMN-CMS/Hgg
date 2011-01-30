//#include <iostream.h> 
//#include <stdlib.h>

// In this example, you pick up migration matrix from a file from the standard EffAcc (effAccSource.root),
// produce an unsmearing matrix and save it to a file (for usage with readUnfoldingMatrices.C)
// this macro works as well for Rapidity and transverse momentum, you have to choose when you launch it
// note if you want to _read_ the code: for how this code developed (written initially for rapidity)
// somewhere the 'variable' which is to be unsmeared is referred to as 'rapidity',
// while the macro consistently works for both Pt and Y

// to execute this example
//     root -l
//root [0] .L makeUnfoldingMatrices.C 
// if you want to produce a RAPIDITY unfolding matrix:
//root [1] makeUnfoldingMatrices("effAcc-V00-05-05/BASE.py.hadded.root","pincopincopFriday.root","Y");
// if you want to produce a TRANSVERSE MOMENTUM unfolding matrix:
//root [1] makeUnfoldingMatrices("effAcc-V00-05-05/BASE.py.hadded.root","pincopincopFriday.root","pt");

#include <TDatime.h>

#include "tdrstyle.C"
#include "zrapidityStandard.C"
int makeUnfoldingMatrices(std::string effAccFileInputFile, std::string unfoldingMatrixOutPutFile, std::string variable) 
{
  bool isRapidity(true);
  if (variable==std::string("Y")){
	std::cout << "\n\tyou have chosen Y" << std::endl;
  }  
  else if (variable==std::string("pt")){
	std::cout << "\n\tyou have chosen pt" << std::endl;
	isRapidity=false;
  }
  else{
	std::cout << "not done a valid choice of variable. Bailing out" << std::endl;
	return 0;
  }

  setTDRStyle();
  zrap_colors();
  
  //gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat(".2f");
  gStyle->SetOptStat(0);
  
  std::cout << "\n\n\tGetting migration matrices from: " << effAccFileInputFile          << std::endl;
  std::cout << "\tCreating unfolding matrices in:  "     << unfoldingMatrixOutPutFile <<"\n\n"<< std::endl;
  
  TFile *theEffAccInpuntFile = new TFile(effAccFileInputFile.c_str());
  //theEffAccInpuntFile->ls();
  
  TFile *theunfoldingMatrixOutpuFile = new TFile(unfoldingMatrixOutPutFile.c_str(),"recreate");
  theunfoldingMatrixOutpuFile ->cd();

  // directory where one saves the original histograms (from EffAcc calculator) the unfolding matrices have been obtained from
  TDirectory *originalHistograms = theunfoldingMatrixOutpuFile->mkdir("originalHistograms");
  originalHistograms->cd();

    gStyle->SetOptStat(110011);

    //////////////////////////////////////////////////////
    // to be consistent with the definitions 
    // of AN-10-367/AN-11-xxx one needs to use the YZTL_vs_YZ/PtTL_vs_PtZ before Eff and Acc
    // are applied ==> move to the "All" folder

    std::string prefix("mcEff/");
    std::string definition=std::string("All");
    std::string cut("");
    if(isRapidity) plot=std::string("YZTL_vs_YZ"); // bin migration histogram for rapidity
    else           plot=std::string("PtTL_vs_PtZ"); // bin migration histogram for Pt
    std::string plotName = prefix + definition;
    plotName = plotName + std::string("/");    plotName = plotName + cut;
    plotName = plotName + std::string("/");    plotName = plotName + plot;
    std::cout << "\tUsing directly the total plot:     "<< plotName << endl;
    TH2F * histoMigration  = (TH2F*) theEffAccInpuntFile->Get(plotName.c_str());

   int binsXMatrixEcalEcal = histoMigration->GetNbinsX();
   int binsYMatrixEcalEcal = histoMigration->GetNbinsY();
   float xminEcalEcal = histoMigration->GetXaxis()->GetXmin();
   float xmaxEcalEcal = histoMigration->GetXaxis()->GetXmax();
   float yMinEcalEcal = histoMigration->GetYaxis()->GetXmin();
   float yMaxEcalEcal = histoMigration->GetYaxis()->GetXmax();

   if(isRapidity){
   histoMigration->GetXaxis()->SetTitle("Y_{reco}");
   histoMigration->GetYaxis()->SetTitle("Y_{ZtreeLevel}");
   } else
   {
   histoMigration->GetXaxis()->SetTitle("q_{T,reco} [GeV/c]");
   histoMigration->GetYaxis()->SetTitle("q_{T,ZtreeLevel} [GeV/c]");
   histoMigration->GetXaxis()->SetTitleOffset(1.3);
   histoMigration->GetYaxis()->SetTitleOffset(1.4);
   }

   TH2F * histoMigrationTotal_toBeSaved = histoMigration->Clone();
   histoMigrationTotal_toBeSaved->Write("originalHistoMigrationTotal");

   // build verification histograms
   std::string histoMigrationTitle("");
   std::string theEffectTitle("");
   if(isRapidity){
	histoMigrationTitle = std::string("Rapidity TreeLevel");
        theEffectTitle = std::string("Y migration");
   } else
   {
	histoMigrationTitle = std::string("Pt TreeLevel");
        theEffectTitle = std::string("Pt migration");
   }
   TH1D * h_RapidityTreeLevel = histoMigration->ProjectionY(histoMigrationTitle.c_str(),0,-1,"");
   if(isRapidity) histoMigration->SetTitle("RapidityTreeLevel");
   else           histoMigration->SetTitle("PtTreeLevel");

   TH1D * h_RapidityTreeLevel_toSave =  h_RapidityTreeLevel->Clone();  
   if(isRapidity) h_RapidityTreeLevel_toSave->Write("original1dTreeLevelRapidity");
   else           h_RapidityTreeLevel_toSave->Write("original1dTreeLevelPt");

   TH1D * h_RapiditySmeared   = histoMigration->ProjectionX("Rapidity Smeared",0,-1,"");
   if(isRapidity) h_RapiditySmeared->SetTitle("RapiditySmeared");
   else           h_RapiditySmeared->SetTitle("PtSmeared");
   h_RapiditySmeared->SetLineColor(2);

   TH1D * h_RapiditySmeared_toSave =  h_RapiditySmeared->Clone();
   if(isRapidity) h_RapiditySmeared_toSave->Write("original1dSmearedRapidity");
   else           h_RapiditySmeared_toSave->Write("original1dSmearedPt");

   //TH1D * h_effect = new TH1D("Y migration","Y migration",binsXMatrixEcalEcal,xminEcalEcal,xmaxEcalEcal); 
   TH1D * h_effect = h_RapiditySmeared->Clone(theEffectTitle.c_str());
   h_effect->Reset();
   h_effect->Add(h_RapiditySmeared, h_RapidityTreeLevel, 1, -1);
   h_effect->Divide(h_RapidityTreeLevel);
   h_effect->SetTitle("(smeared-true)/true");
   h_effect->SetStats(0);


   TH2F * histoMatrix  = histoMigration->Clone("Migration matrix");
   histoMatrix->Reset();
   if(isRapidity) {
   histoMatrix->GetXaxis()->SetTitle("Y_{ZtreeLevel}");
   histoMatrix->GetYaxis()->SetTitle("Y_{ZtreeLevel}");
   } else
   {
   histoMatrix->GetXaxis()->SetTitle("q_{T,reco} [GeV/c]");
   histoMatrix->GetYaxis()->SetTitle("q_{T,true} [GeV/c]");
   histoMatrix->GetXaxis()->SetTitleOffset(1.3);
   histoMatrix->GetYaxis()->SetTitleOffset(1.4);
   }

    float integral=0;
    for(int binY=1; binY<=histoMigration->GetNbinsY(); binY++)
    {
       integral=0;

       for(int binX=1; binX<=histoMigration->GetNbinsX(); binX++)
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
   // double migration[1936];//array reordering histoMatrix, to construct TMatrix
   double migration[10000];    //array reordering histoMatrix, to construct TMatrix
                               // this has to be hard-coded
     for (int u=0; u<histoMatrix->GetNbinsX(); u++){//rows
     for (int v=0; v<histoMatrix->GetNbinsX(); v++){//col

       migration[counter]=histoMatrix->GetBinContent(u+1,v+1);

       // for matrix to be invertible, no diagonal elements can be 0
       // for now set, in all rows which are empty, the diagonal element to 1
       if(u==v && histoMatrix->GetBinContent(u+1,v+1)==0) migration[counter]=1;
       if(0 &&  migration[counter]!=0)  std::cout << u << " " << v << " " << migration[counter] <<std::endl;

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
  
  std::string tableFileName; 
  if(isRapidity) tableFileName=std::string("migration_table_Y.tex");
  else           tableFileName=std::string("migration_table_Pt.tex");
  //FILE* texTable=fopen("migration_table.tex","wt");
  FILE* texTable=fopen(tableFileName.c_str(),"wt");
  fprintf(texTable,"\\begin{tabular}{|cc||ccccc|} \\hline \n");
  if(isRapidity) fprintf(texTable,"$min(Y^\\mathrm{true}_i)$ & $max(Y^\\mathrm{true}_i)$ & $Y^\\mathrm{reco}_{i-2}$ & $Y^\\mathrm{reco}_{i-1}$ & $Y^\\mathrm{reco}_i$ & $Y^\\mathrm{reco}_{i+1}$ & $Y^\\mathrm{reco}_{i+2}$ \\\\ \\hline \n");
  else fprintf(texTable,"$min(q^\\mathrm{true}_{T,i})$ & $max(q^\\mathrm{true}_{T,i})$ & $q^\\mathrm{reco}_{T,i-2}$ & $q^\\mathrm{reco}_{T,i-1}$ & $q^\\mathrm{reco}_i$ & $q^\\mathrm{reco}_{T,i+1}$ & $q^\\mathrm{reco}_{T,i+2}$ \\\\ \\hline \n");
  for (int i=0; i<histoMigration->GetNbinsX(); i++) {
    double lv=histoMigration->GetXaxis()->GetBinLowEdge(i+1);
    double lh=histoMigration->GetXaxis()->GetBinUpEdge(i+1);

    if ((isRapidity) && ( fabs(lv)>3.8 || fabs(lh)>3.8)) continue;

    fprintf(texTable,"%5.2f & %5.2f & ",lv,lh);

    for (int j=-2; j<=2; j++) {
      int ireco=i+j;
      if (ireco<0 || ireco>=histoMigration->GetNbinsX()) {
	fprintf(texTable," -- ");
	if (j!=2) fprintf(texTable," & ");
	else fprintf(texTable,"\\\\\n");
	continue;
      } else {
	fprintf(texTable,"%.2e",MigrationMatrix(ireco,i));
	if (j!=2) fprintf(texTable," & ");
	else fprintf(texTable,"\\\\\n");
      }


    }
     
  }
  fprintf(texTable,"\\hline\n\\end{tabular}\n");  fclose(texTable);

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
  //TH1D * h_RapidityUNSmeared = new TH1D("Rapidity Unsmeared", "Rapidity Unsmeared", binsXMatrixEcalEcal, xminEcalEcal, xmaxEcalEcal);
  TH1D * h_RapidityUNSmeared = h_RapiditySmeared->Clone(theEffectTitle.c_str()); 
  h_RapidityUNSmeared->Reset();
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // turning the unsmeared vector into a histogram
  for(int bin=1; bin<=binsXMatrixEcalEcal; bin++)
  {
   h_RapidityUNSmeared->SetBinContent(bin,vectorRapidityUNSmeared[bin-1]);
   //std::cout << vectorRapidityUNSmeared[bin-1]<<std::endl;
  }







  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(1){
    
    // http://ultrahigh.org/2007/08/20/making-pretty-root-color-palettes/

   TCanvas * c4 = new TCanvas("Migration histogram","Migration histogram",1050,750);
   c4->cd(1);  c4->SetLogz();     
   if(!isRapidity) {
   c4->SetLogx();   c4->SetLogy();
   }
   histoMigration->SetTitle("Migration histogram");
   histoMigration->SetTitle("Migration Total");
   histoMigration->GetXaxis()->CenterTitle();
   histoMigration->GetYaxis()->CenterTitle();
   histoMigration->Draw("colz");
   c4->Print( (std::string("migration-total") + std::string(".png") ).c_str() );


   TCanvas * c4q = new TCanvas("c4q","c4q",1050,750);
   c4q->cd();
   h_effect->SetLineColor(2);
   h_effect->Draw();

   //new way of getting the time 
   const char* time_;
   TDatime *mytime = new TDatime();
   time_ = mytime->AsString();
   TText *plabel = new TText();
   plabel-> SetNDC();
   plabel -> SetTextFont(1);
   plabel -> SetTextColor(1);
   plabel -> SetTextSize(0.04);
   plabel -> SetTextAlign(22);
   plabel -> SetTextAngle(0);
   TText *tlabel = new TText();
   tlabel-> SetNDC();
   tlabel -> SetTextFont(1);
   tlabel -> SetTextColor(1);
   tlabel -> SetTextSize(0.02);
   tlabel -> SetTextAlign(22);
   tlabel -> SetTextAngle(0);
   std::cout << "\n\n The local time is: " << time_ << "\n\n" << std::endl;

   TCanvas * c5 = new TCanvas("migration matrix","migration matrix",1050,750);
   c5->SetTopMargin(0.02);
   c5->SetRightMargin(0.14);
   c5->SetLeftMargin(0.11);   
   c5->SetBottomMargin(0.125);
   c5->cd(1);  c5->SetLogz();   
   TH2F * dummyForRange = new TH2F("dummyForRange","dummyForRange",2,0,600,2,0,600);
   dummyForRange->GetXaxis()->CenterTitle();
   dummyForRange->GetYaxis()->CenterTitle(true);
   dummyForRange->GetYaxis()->SetTitleOffset(0.8);
   dummyForRange->SetStats(0);
   if(isRapidity) {
   histoMatrix->GetXaxis()->SetTitle("Y_{Z,reco}");
   histoMatrix->GetYaxis()->SetTitle("Y_{Z,true}");
   histoMatrix->GetXaxis()->SetRangeUser(-3.8,3.8);
   histoMatrix->GetYaxis()->SetRangeUser(-3.8,3.8);
   } else {
   histoMatrix->GetXaxis()->SetTitle("q_{T,reco} [GeV/c]");
   histoMatrix->GetYaxis()->SetTitle("q_{T,true} [GeV/c]");
   histoMatrix->GetXaxis()->SetTitleOffset(1.3);
   histoMatrix->GetYaxis()->SetTitleOffset(1.4);
   histoMatrix->GetXaxis()->SetRangeUser(0.0001,600);
   histoMatrix->GetYaxis()->SetRangeUser(0.0001,600);
   // dummyForRange->Draw();
   c5->SetLogx();   c5->SetLogy();
   }
   histoMatrix->GetXaxis()->CenterTitle();
   histoMatrix->GetYaxis()->CenterTitle(true);
   histoMatrix->GetYaxis()->SetTitleOffset(0.8);
   histoMatrix->SetMaximum(1);
   histoMatrix->SetStats(0);
   //histoMatrix->Draw("colz same");
   histoMatrix->Draw("colz");
   plabel -> DrawText(0.30, 0.93, "CMS 2010 PRELIMINARY");
   tlabel -> DrawText(0.30, 0.90, Form("%s",time_));
   histoMatrix->SetTitle("Migration matrix (normalized)");
   c5->Print("migrationMatrix.png");
   c5->Print("migrationMatrix.eps");
   c5->Print("migrationMatrix.pdf");


   TCanvas * c7pre = new TCanvas("closure test","closure test",1050,750);
   if(!isRapidity) {
   c7pre->SetLogx();   c7pre->SetLogy();
   }
   c7pre->cd();
   h_RapidityTreeLevel->SetLineWidth(2);//this is black for now
   h_RapidityTreeLevel->GetYaxis()->SetTitle("events");
   h_RapidityTreeLevel->SetTitle(0);
   h_RapidityTreeLevel->Draw("");//this is black for now
  
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
   leg->AddEntry(h_RapiditySmeared,"Y_{Z} smeared","p");
   leg->AddEntry(h_RapidityTreeLevel,"Y_{Z} treeLevel","l");
   leg->AddEntry(h_RapidityUNSmeared,"Y_{Z} unsmeared (overlap?)","p");
   leg->Draw();

   ///////////////////////////////////////////////////////
   c7pre->Print("rapdity-Tree-Smear-Unsmear-Closure.png");
   ///////////////////////////////////////////////////////

  }// end of switch on/off plotting of many canvases
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   TH1D * h_effectUnsmearing       = new TH1D("Y unsmearing","Y unsmearing",binsXMatrixEcalEcal,xminEcalEcal,xmaxEcalEcal); 
   TH1D * h_effectUnsmearingByHand = new TH1D("Y unsmearing BH","Y unsmearing BH",binsXMatrixEcalEcal,xminEcalEcal,xmaxEcalEcal); 
   h_RapidityUNSmeared->Sumw2();
   h_RapidityTreeLevel->Sumw2();
   h_effectUnsmearing->Add(h_RapidityUNSmeared, h_RapidityTreeLevel, 1, -1);
   h_effectUnsmearing->Divide(h_RapidityTreeLevel);
   h_effectUnsmearing->SetTitle(0);
   for(int a=1; a<=binsXMatrixEcalEcal; a++){
     if(h_RapidityTreeLevel->GetBinContent(a)>1){
     h_effectUnsmearingByHand->SetBinContent(a,(h_RapidityUNSmeared->GetBinContent(a)-h_RapidityTreeLevel->GetBinContent(a))/h_RapidityTreeLevel->GetBinContent(a));
     }
   }
   if(1){
     TCanvas * c8 = new TCanvas("closure test: (UNsmeared-true)/true","closure test: (UNsmeared-true)/true",1050,750);
     c8->cd();
     //h_effectUnsmearing->Draw();
     // h_effectUnsmearingByHand->Draw("e hpx HIST");
     h_effectUnsmearingByHand->Draw("hpx HIST");
     int i;
   }
   c8->Print("UNsmeared-true-over-true-Closure.png");
   std::cout << "\tFile: "<< effAccFileInputFile << "   : ";
   std::cin >> i;
   //theunfoldingMatrixOutpuFile->Write();

   return 0;
}
