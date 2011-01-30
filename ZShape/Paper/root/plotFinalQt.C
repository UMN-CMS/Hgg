#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH2.h"
#include "readUnfoldingMatrices.C"


static const int    pt_bins=19;
static const double pt_binning[]={ 0.0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };
static const double pt_binning_vis[]={ 0.7, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };

static const int BINCOUNT=pt_bins;


struct DataSeries {

  DataSeries() {
    clear();
  }
  void clear() {
    for (int i=0; i<BINCOUNT; i++) {
      xave[i]=(pt_binning[i]+pt_binning[i+1])/2;
      xwidth[i]=(pt_binning[i+1]-pt_binning[i])/2;
      y[i]=0;
      ey[i]=0;
    }
  }
  void clearContent() {
    for (int i=0; i<BINCOUNT; i++) {
      y[i]=0;
      ey[i]=0;
    }
  }
  DataSeries(const DataSeries& ds) {
    for (int i=0; i<BINCOUNT; i++) {
      xave[i]=(pt_binning[i]+pt_binning[i+1])/2;
      xwidth[i]=(pt_binning[i+1]-pt_binning[i])/2;
      y[i]=ds.y[i];
      ey[i]=ds.ey[i];
    }
  }
  DataSeries(const char* fname, int iy=0);
  DataSeries(int lumi, const char* fnameformat, int iy=0);

  TGraphErrors* makeTGE() {
    return new TGraphErrors(BINCOUNT,xave,y,xwidth,ey);
  }
  TGraph* makeTG() {
    return new TGraph(BINCOUNT,xave,y);
  }

  TH1* makeTH1(const char* name) {
    TH1F* hist=new TH1F(name,name,BINCOUNT,pt_binning_vis);
    for (int i=1; i<=BINCOUNT; i++) {
      hist->SetBinContent(i,y[i-1]);
      hist->SetBinError(i,ey[i-1]);
      //      printf("%s %d %f %f \n",name,i+startFrom,y[startFrom+i-1],ey[startFrom+i-1]);
    }
    hist->SetDirectory(0);
    return hist;
  }

  void create(const char* file, int iy);

  void normalizeToWidth() {
    for (int i=0; i<BINCOUNT; i++) {
      y[i]/=(xwidth[i]*2);
      ey[i]/=(xwidth[i]*2);
    }
  }

  double xave[BINCOUNT];
  double xwidth[BINCOUNT];
  double y[BINCOUNT];
  double ey[BINCOUNT];
};  


#include "tdrstyle.C"
#include "zrapidityStandard.C"

DataSeries::DataSeries(int lumi, const char* file, int iy) {
  char fullname[1024];
  sprintf(fullname,file,lumi);
  create(fullname,iy);
}
DataSeries::DataSeries(const char* file, int iy) {
  create(file,iy);
}
void DataSeries::create(const char* file, int iy) {
  clear();

  FILE* f=fopen(file,"r");
  if (f==0) {
    fprintf(stderr,"Cannot open '%s'\n",file);
    return;
  }
  char line[4096];
  while (!feof(f)) {
    line[0]=0;
    fgets(line,4095,f);
    if (strlen(line)<5 || line[0]=='#') continue;
    int i;
    float a,b,q[20];
    int ate, atetotal=0;
    int found=sscanf(line," %d %f %f %n",&i,&a,&b,&ate);
    if (found!=3) continue;

    atetotal+=ate;
    for (int iq=0; iq<=iy && found!=0; iq++) {
      found=sscanf(line+atetotal," %f %n ",&q[iq],&ate);
      //      printf("%d %d %d %f\n",i,iq,ate,q[iq]);
      atetotal+=ate;
    }
    if (found!=1) continue;
    if (i<1 || i>BINCOUNT) continue;
    y[i-1]=q[iy];
  }
  fclose(f);
}

void plotFinalQt(TFile* mctruth) {
  int lumi=30;
  setTDRStyle();

  TH1* truth=(TH1*)(mctruth->Get("mcEff/ECAL80-ECAL95/Acceptance/Z0_Pt_masscut")->Clone("truth"));
  truth->SetDirectory(0);
  
  DataSeries data_ee("ECAL80-ECAL95_Pt_RAW.csv");
  DataSeries data_all;
  DataSeries data_bkgd(data_all);
  
  DataSeries effAcc("eff_pt.txt"), ea_statErr("effStatPtError.csv");

  const int firsti=1;
  const int lasti=18;
  const int npoint=lasti-firsti+1;

  double tt=0, td=0;
 
  for (int i=0; i<BINCOUNT; i++) {
    data_bkgd.y[i]=data_ee.y[i];
  }


  DataSeries dataStatError(data_all);

  DataSeries backgroundAll("backgroundBkgFinalPtVals.csv",0);
  DataSeries backgroundAllUnc("backgroundBkgFinalPtVals.csv",1);

  DataSeries backgroundUncFrac(backgroundAllUnc);

  for (int i=0; i<BINCOUNT; i++) {  
    data_all.y[i]=data_bkgd.y[i]-backgroundAll.y[i];    
    data_all.ey[i]=data_bkgd.ey[i];
  }    

  DataSeries corrData(data_all);
  DataSeries totalSyst(data_all);

  for (int i=0; i<BINCOUNT; i++) {  
    if (effAcc.y[i]==0) {
      corrData.y[i]=0;
      continue;
    }
    //    printf("%2d %f %f \n",i,corrData.y[i],effAcc.y[i]);
    corrData.y[i]/=effAcc.y[i];
  }
  // unsmear


  for (int i=0; i<BINCOUNT; i++) {  
    if   (effAcc.y[i]==0) corrData.y[i]=0;
    else                  corrData.ey[i]=sqrt(data_all.y[i]+backgroundAll.y[i])/effAcc.y[i];
  }
  DataSeries corrDataClone(corrData);

  bool doUnsmearing(false); // use this to activate unsmearing + associated errors or leave it completely off
  if (doUnsmearing) {
    TH1* data_corr_smeared=corrData.makeTH1("data_corr_smeared");
    TH1* data_corr_unsmeared=unfold(data_corr_smeared,"../root/unfoldingMatrix_theOutPut.root");
    TFile*    theunfoldingMatrixInputFile = new TFile("../root/unfoldingMatrix_theOutPut.root","read");
    TMatrixD* theUnfoldingMatrix          = (TMatrixD*)theunfoldingMatrixInputFile->Get("unsmearMatrices/unfoldingMatrixTotal"); // indices [0.. 99] X [0.. 99]
    double errorCumul=0;
    for (int i=0; i<BINCOUNT; i++) {  
      corrData.y[i]=data_corr_unsmeared->GetBinContent(i+1);
      for (int s=0; s<BINCOUNT; s++) {
	errorCumul+= pow( corrDataClone.ey[s] * (*theUnfoldingMatrix)(s,i) , 2);
      }
      if (doUnsmearing) corrData.ey[i]=sqrt(errorCumul);
      errorCumul=0;
    }
  }
  DataSeries corrDataSyst(corrData);
  DataSeries corrDataBkgd(corrData);

  //  DataSeries pdfPos("pdfErrsPos.txt");
  //  DataSeries pdfNeg("pdfErrsNeg.txt");
  //DataSeries pdfTotal(data_all);
  //DataSeries pdfFrac(data_all);

  DataSeries energyScale("energyScale_pt.csv");

  DataSeries effSystGSF("effSystGSFPtError.csv");
  DataSeries effSystWP80("effSystWP80PtError.csv");
  DataSeries effSystWP95("effSystWP95PtError.csv");
  DataSeries effSystHLT("effSystHLTPtError.csv");

  DataSeries effSystErr;

  for (int i=firsti-1; i<lasti; i++) {
    effSystErr.y[i]=sqrt(pow(effSystGSF.y[i],2)+
			 pow(effSystWP80.y[i],2)+
			 pow(effSystWP95.y[i],2)+
			 pow(effSystHLT.y[i],2));
  }  

  double sumtotal=0;
  
  for (int i=firsti-1; i<lasti; i++) {
    if (effAcc.y[i]==0) {
      corrData.y[i]=0;
      printf("hih\n");
      continue;
    }
    // what is this for??
    //corrData.ey[i]=sqrt(data_all.y[i])/effAcc.y[i];
    //    corrDataSyst.xave[i]+=0.02; // uncomment to offset error bars

    backgroundUncFrac.y[i]=backgroundAllUnc.y[i]/std::max(1.0,data_bkgd.y[i]);
    // GF
    //corrDataBkgd.ey[i]=sqrt(data_all.y[i]+backgroundAll.y[i])/effAcc.y[i];
    
    dataStatError.y[i]= corrDataBkgd.ey[i]/std::max(1.0,corrDataBkgd.y[i]);

    //    corrDataSyst.y[i]/=effAcc.y[i];

    /*
    pdfTotal.y[i]=(pdfPos.y[i]+(-pdfNeg.y[i]))/2*corrDataBkgd.y[i];
    pdfFrac.y[i]=(pdfPos.y[i]+(-pdfNeg.y[i]))/2;
    */
    //  printf("%d %f %f %f %f\n",i,corrData.ey[i],corrData.y[i]*ea_statErr.y[i]);
    corrDataSyst.ey[i]=sqrt(pow(corrDataBkgd.ey[i],2)+
			    pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			    pow(backgroundAllUnc.y[i]/effAcc.y[i],2)+
			    pow(energyScale.y[i]*corrData.y[i],2)+
			    pow(effSystErr.y[i]*corrData.y[i],2)+
			    0
			    //			    pow(pdfTotal.y[i],2)
			    );

    totalSyst.ey[i]=sqrt( pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			  pow(backgroundAllUnc.y[i]/effAcc.y[i],2)+
			  pow(energyScale.y[i]*corrData.y[i],2)+
			  pow(effSystErr.y[i]*corrData.y[i],2)+
			  0
			  //			  pow(pdfTotal.y[i],2)
			  );

    totalSyst.y[i]=totalSyst.ey[i]/std::max(1.0,corrDataBkgd.y[i]);

    sumtotal+=corrDataSyst.y[i];
  }
  char fnamework[1024];

  
 
  data_all.normalizeToWidth();
  data_bkgd.normalizeToWidth();
  corrData.normalizeToWidth();
  corrDataSyst.normalizeToWidth();
  corrDataBkgd.normalizeToWidth();


  TH1* truth_vis=new TH1F("tvis","tvis",BINCOUNT,pt_binning_vis);

  for (int i=0; i<BINCOUNT; i++) {
    truth_vis->SetBinContent(i+1,truth->GetBinContent(i+1)/(2*corrData.xwidth[i]));

    tt+=truth_vis->GetBinContent(i+1);
    td+=corrData.y[i];
  }

  truth_vis->Scale(td/tt*1.0);

  sprintf(fnamework,"ZQt_final_table-%d.tex",lumi);

  FILE* ftable=fopen(fnamework,"w");

  fprintf(ftable,"\\begin{tabular}{|c|c||c|c|c|}\n\\hline\n");
  fprintf(ftable,"$|q_{T,min}|$ & $|q_{T,max}|$ & Measurement & Statistical Error & Systematic Error \\\\ \\hline\n");
  for (int j=firsti-1; j<lasti; j++) {
    fprintf(ftable,"%7.2f & %7.2f & ",
	    corrData.xave[j]-corrData.xwidth[j],
	    corrData.xave[j]+corrData.xwidth[j]);
    if (corrData.y[j]/td>0.01)
      fprintf(ftable,"%7.4f & %7.4f & %7.4f \\\\ \n",
	      corrData.y[j]/td,
	      corrData.ey[j]/td,
	      sqrt(pow(corrDataSyst.ey[j]/td,2)-pow(corrData.ey[j]/td,2))
	      );
    else 
      fprintf(ftable,"%7.2e & %7.2e & %7.2e \\\\ \n",
	      corrData.y[j]/td,
	      corrData.ey[j]/td,
	      sqrt(pow(corrDataSyst.ey[j]/td,2)-pow(corrData.ey[j]/td,2))
	      );
  }
  
  fprintf(ftable,"\\hline\n\\end{tabular}");
  fclose(ftable);
 



  
  TGraph* rawd=data_all.makeTG();
  TGraph* rawdb=data_bkgd.makeTG();
  TGraph* corrd=corrData.makeTGE();
  TGraph* corrdsys=corrDataSyst.makeTGE();
  TGraph* corrdbkgd=corrDataBkgd.makeTGE();
  // this is special
  TGraph* corrdsysb=new TGraphErrors(npoint,corrDataBkgd.xave+firsti-1,corrDataBkgd.y+firsti-1,corrDataSyst.xwidth+firsti-1,corrDataSyst.ey+firsti-1);

  double ybinage[21];
  for (int i=0; i<21; i++) {
    ybinage[i]=0.1+900/20.0*i;
  }

  TH2* dummy=new TH2F("dummy","dummy",18,pt_binning_vis,20,ybinage);

  dummy->GetYaxis()->SetTitle("Events/GeV");
  dummy->GetXaxis()->SetTitle("q_{T,Z}");
  
  TCanvas* c1=new TCanvas("finalZQt","finalZQt",800,600);
  c1->SetLogx();
  c1->SetLogy();
  dummy->Draw();
  dummy->SetDirectory(0);

  rawd->SetMarkerStyle(24);
  rawd->Draw("PSAME");
  rawdb->SetMarkerStyle(25);
  rawdb->Draw("PSAME");

  //  corrd->SetMarkerStyle(20);
  //  corrd->Draw("PSAME");

  // corrdsys->SetMarkerStyle(0);
  // corrdsys->Draw("PSAME");

  corrdbkgd->SetMarkerStyle(20);
  corrdbkgd->Draw("PSAME");

  corrdsysb->SetMarkerStyle(0);
  corrdsysb->Draw("PSAME");

  truth_vis->Draw("SAMEHIST");

  TLegend* tl=new TLegend(0.4,0.16,0.72,0.34);
  tl->AddEntry(corrd,"Corrected Data","P");
  tl->AddEntry(rawdb,"Raw Data","P");
  tl->AddEntry(rawd,"Background-Subtracted Data","P");
  tl->AddEntry(truth_vis,"POWHEG+CT10 Prediction","l");
  tl->Draw();

  zrap_Prelim(0.85,0.90,0.85,0.85);
  zrap_Lumi(0.85,0.86,lumi);

  sprintf(fnamework,"ZQt_final-%d.eps",lumi);
  c1->Print(fnamework);
  sprintf(fnamework,"ZQt_final-%d.png",lumi);
  c1->Print(fnamework);

  TCanvas* c3=new TCanvas("finalZQtErrors","finalZQtErrors",800,600);
  TH1* dummy3=new TH1F("dummy3","dummy3",10,0.7,600);
  dummy3->SetDirectory(0);
  c3->SetLogy();
  c3->SetLogx();
  dummy3->SetMaximum(1.0);
  dummy3->SetMinimum(0.001);
  dummy3->GetYaxis()->SetTitle("Fractional Error");
  dummy3->GetXaxis()->SetTitle("q_{T,Z}");

  
  //  TH1* binmiggr=binMigPos.makeTH1("bin migration");
  TH1* effacc_stat_gr=ea_statErr.makeTH1("eff stat");
  TH1* effacc_syst_gr=effSystErr.makeTH1("eff syst");
  TH1* gr_dataStatError=dataStatError.makeTH1("data stats err");
  TH1* gr_bkgd_unc=backgroundUncFrac.makeTH1("bkgd unc");
  //  TH1* gr_pdf=pdfFrac.makeTH1("pdf_err");
  TH1* gr_energyscale=energyScale.makeTH1("es_unc");
  TH1* gr_systotal=totalSyst.makeTH1("total_syst");
  dummy3->Draw();

  //  binmiggr->SetLineWidth(2);
  //  binmiggr->SetLineColor(kBlue);

  effacc_stat_gr->SetLineWidth(2);
  effacc_stat_gr->SetLineColor(kGreen+2);

  effacc_syst_gr->SetLineWidth(2);
  effacc_syst_gr->SetLineColor(kGreen-5);

  gr_dataStatError->SetLineWidth(2);
  gr_dataStatError->SetLineColor(kBlack);
  gr_dataStatError->SetLineStyle(2);

  gr_bkgd_unc->SetLineWidth(2);
  gr_bkgd_unc->SetLineColor(kRed-4);

  gr_energyscale->SetLineWidth(2);
  gr_energyscale->SetLineColor(kBlue);

  //  gr_pdf->SetLineWidth(2);
  //  gr_pdf->SetLineColor(kViolet-5);

  gr_systotal->SetLineWidth(2);
  

  //  binmiggr->Draw("HSAME");
  effacc_stat_gr->Draw("HSAME");
  effacc_syst_gr->Draw("HSAME");
  gr_dataStatError->Draw("HSAME");
  gr_bkgd_unc->Draw("HSAME");
  gr_energyscale->Draw("HSAME");
  //  gr_pdf->Draw("HSAME");
  gr_systotal->Draw("HISTSAME");

   tl=new TLegend(0.35,0.70,0.75,0.94,"Error Source Contributions");
   //  tl->AddEntry(binmiggr,"Bin Migration","L");
  tl->AddEntry(gr_dataStatError,"Data Statistics","L");
  tl->AddEntry(gr_systotal,"Total Systematic Error","L");
  tl->AddEntry(gr_bkgd_unc,"Background Uncertainty","L");
  tl->AddEntry(effacc_stat_gr,"Efficiency Statistics","L");
  //  tl->AddEntry(gr_pdf,"PDF (#epsilon #times A) Uncertainty","L");
  tl->AddEntry(gr_energyscale,"Energy Scale Uncertainty","L");
  tl->AddEntry(effacc_syst_gr,"Efficiency Systematics","L");
  tl->Draw();
  
  zrap_Prelim(0.3,0.24,0.3,0.2);
  zrap_Lumi(0.30,0.19,lumi);

  sprintf(fnamework,"ZQt_final_errors-%d.eps",lumi);
  c3->Print(fnamework);
  sprintf(fnamework,"ZQt_final_errors-%d.png",lumi);
  c3->Print(fnamework);

  sprintf(fnamework,"ZRapidity_final_fold_errtable-%d.tex",lumi);
  ftable=fopen(fnamework,"w");

  /*  
  {

    if (i==0) {
      fprintf(ftable,"\\begin{tabular}{|c|c||c|c|c|c|}\\hline\n");
      fprintf(ftable,"            &             & Efficiency &   Bin     & Energy Scale  & Background \\\\ \n");
      fprintf(ftable,"$|Y_{min}|$ & $|Y_{max}|$ & Statistics & Migration & Uncertainty & Estimation \\\\ \\hline \n");
    }
    
      fprintf(ftable," %7.2f & %7.2f & %.4f & %.4f & %.4f & %.4f \\\\ \n",
             corrDataFold.xave[j]-corrDataFold.xwidth[j],
             corrDataFold.xave[j]+corrDataFold.xwidth[j],
             staterr,systerr,energyerr,bkgderr);
        
  }
  fprintf(ftable,"\\hline\n\\end{tabular}\n");
  fclose(ftable);
  */

  //  printf(" Two models: %f %f \n",chi2_0,chi2_1);

}
