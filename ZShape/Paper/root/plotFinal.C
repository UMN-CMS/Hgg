#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH2.h"
#include "readUnfoldingMatrices.C"

static const int BINCOUNT=100;
static const double xlow=-5.0;
static const double xhigh=5.0;

struct DataSeries {

  DataSeries() {
    clear();
  }
  void clear() {
    for (int i=0; i<BINCOUNT; i++) {
      xave[i]=xlow+(xhigh-xlow)/BINCOUNT*(i+0.5);
      xwidth[i]=(xhigh-xlow)/BINCOUNT;
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
      xave[i]=xlow+(xhigh-xlow)/BINCOUNT*(i+0.5);
      xwidth[i]=(xhigh-xlow)/BINCOUNT/2;
      y[i]=ds.y[i];
      ey[i]=ds.ey[i];
    }
  }
  DataSeries(const char* fname, int iy=0);
  DataSeries(int lumi, const char* fnameformat, int iy=0);

  TGraphErrors* makeTGE(int bins, int startFrom) {
    return new TGraphErrors(bins,xave+startFrom,y+startFrom,xwidth+startFrom,ey+startFrom);
  }
  TGraph* makeTG(int bins, int startFrom) {
    return new TGraph(bins,xave+startFrom,y+startFrom);
  }

  TH1* makeTH1(const char* name, int bins, int startFrom) {
    TH1F* hist=new TH1F(name,name,bins,xave[startFrom]-xwidth[startFrom],xave[startFrom+bins-1]+xwidth[startFrom+bins-1]);
    for (int i=1; i<=bins; i++) {
      hist->SetBinContent(i,y[startFrom+i-1]);
      hist->SetBinError(i,ey[startFrom+i-1]);
      //      printf("%s %d %f %f \n",name,i+startFrom,y[startFrom+i-1],ey[startFrom+i-1]);
    }
    hist->SetDirectory(0);
    return hist;
  }

  void create(const char* file, int iy);

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

void plotFinal(TFile* mctruth) {
  int lumi=30;
  setTDRStyle();

  TH1* truth=(TH1*)(mctruth->Get("mcEff/All/Z0_Y_masscut")->Clone("truth"));
  truth->SetDirectory(0);
  
  DataSeries data_ee("ECAL80-ECAL95_RAW.csv");
  DataSeries data_ef("ECAL80-HF_RAW.csv"), data_all;
  DataSeries data_bkgd(data_all);
  
  DataSeries effAcc("effAcc_combined.csv"), ea_statErr("effStatError.csv");

  const int firsti=15;
  const int lasti=84;
  const int npoint=lasti-firsti+1;

  double tt=0, td=0;
 
  for (int i=0; i<BINCOUNT; i++) {
    data_bkgd.y[i]=data_ef.y[i]+data_ee.y[i];
  }


  DataSeries dataStatError(data_all);

  DataSeries backgroundEE("backgroundBkgFinalFixVals.txt",0);
  DataSeries backgroundHF("backgroundBkgFinalHFFixVals.txt",0);

  DataSeries backgroundEEUnc("backgroundBkgFinalFixVals.txt",1);
  DataSeries backgroundHFUnc("backgroundBkgFinalHFFixVals.txt",1);

  DataSeries backgroundAll(backgroundEE), backgroundAllUnc(backgroundEEUnc),
    backgroundUncFrac(backgroundEEUnc);


  for (int i=0; i<BINCOUNT; i++) {  
    backgroundAll.y[i]=backgroundEE.y[i]+backgroundHF.y[i];
    backgroundAllUnc.y[i]=sqrt(pow(backgroundEEUnc.y[i],2)+pow(backgroundHFUnc.y[i],2));

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
    corrData.y[i]/=effAcc.y[i];
  }
  // unsmear

  printf("A\n");
  TH1* data_corr_smeared=corrData.makeTH1("data_corr_smeared",100,0);
  TH1* data_corr_unsmeared=unfold(data_corr_smeared,"../root/unfoldingMatrix_theOutPut.root");
  for (int i=0; i<BINCOUNT; i++) {  
    //      corrData.y[i]=data_corr_unsmeared->GetBinContent(i+1);
    //    data_all.ey[i]=sqrt(data_all.y[i]);
  }
  DataSeries corrDataSyst(corrData);
  DataSeries corrDataBkgd(corrData);

  //  DataSeries pdfPos("pdfErrsPos.txt");
  //  DataSeries pdfNeg("pdfErrsNeg.txt");
  //DataSeries pdfTotal(data_all);
  //DataSeries pdfFrac(data_all);

  DataSeries energyScale("energyScaleError.csv");

  DataSeries effSystGSF("effSystGSFError.csv");
  DataSeries effSystWP80("effSystWP80Error.csv");
  DataSeries effSystWP95("effSystWP95Error.csv");
  DataSeries effSystHFEID("effSystHFEIDError.csv");
  DataSeries effSystHLT("effSystHLTError.csv");

  DataSeries effSystErr;

  for (int i=firsti-1; i<lasti; i++) {
    effSystErr.y[i]=sqrt(pow(effSystGSF.y[i],2)+
			 pow(effSystWP80.y[i],2)+
			 pow(effSystWP95.y[i],2)+
			 pow(effSystHFEID.y[i],2)+
			 pow(effSystHLT.y[i],2));
  }  

  double sumtotal=0;
  
  for (int i=firsti-1; i<lasti; i++) {
    if (effAcc.y[i]==0) {
      corrData.y[i]=0;
      continue;
    }

    corrData.ey[i]=sqrt(data_all.y[i])/effAcc.y[i];
    //    corrDataSyst.xave[i]+=0.02; // uncomment to offset error bars


    tt+=truth->GetBinContent(i+1);
    td+=corrData.y[i];

    backgroundUncFrac.y[i]=backgroundAllUnc.y[i]/std::max(1.0,data_bkgd.y[i]);

    corrDataBkgd.ey[i]=sqrt(data_all.y[i]+backgroundAll.y[i])/effAcc.y[i];

    dataStatError.y[i]= corrDataBkgd.ey[i]/std::max(1.0,corrDataBkgd.y[i]);

    //    corrDataSyst.y[i]/=effAcc.y[i];

    /*
    pdfTotal.y[i]=(pdfPos.y[i]+(-pdfNeg.y[i]))/2*corrDataBkgd.y[i];
    pdfFrac.y[i]=(pdfPos.y[i]+(-pdfNeg.y[i]))/2;
    */
    printf("%d %f %f %f %f\n",i,corrData.ey[i],corrData.y[i]*ea_statErr.y[i]);
    corrDataSyst.ey[i]=sqrt(pow(corrDataBkgd.ey[i],2)+
			    pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			    pow(backgroundAllUnc.y[i]/effAcc.y[i],2)+
			    pow(energyScale.y[i]*data_all.y[i],2)+
			    pow(effSystErr.y[i]*data_all.y[i],2)+
			    0
			    //			    pow(pdfTotal.y[i],2)
			    );

    totalSyst.ey[i]=sqrt( pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			  pow(backgroundAllUnc.y[i]/effAcc.y[i],2)+
			  pow(energyScale.y[i]*data_all.y[i],2)+
			  pow(effSystErr.y[i]*data_all.y[i],2)+
			  0
			  //			  pow(pdfTotal.y[i],2)
			  );

    totalSyst.y[i]=totalSyst.ey[i]/std::max(1.0,corrDataBkgd.y[i]);

    sumtotal+=corrDataSyst.y[i];
  }
  char fnamework[1024];

  
  DataSeries corrDataFold(corrData),corrDataSystFold(corrData);

  sprintf(fnamework,"ZRapidity_final_fold_table-%d.tex",lumi);

  FILE* ftable=fopen(fnamework,"w");

  fprintf(ftable,"\\begin{tabular}{|c|c||c|c|c|}\n\\hline\n");
  fprintf(ftable,"$|Y_{min}|$ & $|Y_{max}|$ & Measurement & Statistical Error & Systematic Error \\\\ \\hline\n");
  
  const int izero=51;
  const double binwidth=0.1;
  //  double binwidth=1.0;
  for (int i=0; i<=(lasti-firsti)/2; i++) {
    int j=izero+i-1;
    int jadd=izero-i-2;
    corrDataFold.y[j]=(corrDataBkgd.y[j]+corrDataBkgd.y[jadd])/std::max(1.0,sumtotal);


    //    printf("%d %d (%f) %d (%f) \n",i,j,corrDataBkgd.y[j],jadd,corrDataBkgd.y[jadd]);

    double ea_Ave=(effAcc.y[j]+effAcc.y[jadd])/2;

    //    double errmig=(fabs(binMig.y[j])*corrData.y[j]+fabs(binMig.y[jadd])*corrData.y[jadd])/sumtotal;
    //    double pdferr=(pdfFrac.y[j]+pdfFrac.y[jadd])/2*corrDataFold.y[j]/sumtotal;
    double bkgderr=sqrt(pow(backgroundAllUnc.y[j],2)+pow(backgroundAllUnc.y[jadd],2))/ea_Ave/sumtotal;
    double energyerr=(energyScale.y[j]+energyScale.y[jadd])/2/ea_Ave/sumtotal;


    corrDataFold.ey[j]=sqrt(data_all.y[j]+data_all.y[jadd]+backgroundAll.y[j]+backgroundAll.y[jadd])/ea_Ave/sumtotal;

    corrDataSystFold.y[j]=corrDataFold.y[j];
    corrDataSystFold.ey[j]=sqrt(pow(corrDataFold.ey[j],2)+
				pow(corrDataFold.y[j]*(ea_statErr.y[j]+ea_statErr.y[jadd])/2,2)+
				pow(corrDataFold.y[j]*(effSystErr.y[j]+effSystErr.y[jadd])/2,2)+
				//				pow(pdferr,2)+
				pow(energyerr,2)+
				pow(bkgderr,2)
				);

    corrDataFold.y[j]/=binwidth;
    corrDataFold.ey[j]/=binwidth;
    corrDataSystFold.y[j]/=binwidth;
    corrDataSystFold.ey[j]/=binwidth;


    fprintf(ftable,"%7.2f & %7.2f & %7.4f & %7.4f & %7.4f \\\\ \n",
	    corrDataFold.xave[j]-corrDataFold.xwidth[j],
	    corrDataFold.xave[j]+corrDataFold.xwidth[j],
	    corrDataFold.y[j],
	    corrDataFold.ey[j],
	    sqrt(pow(corrDataSystFold.ey[j],2)-pow(corrDataFold.ey[j],2))
	    );
  }

  fprintf(ftable,"\\hline\n\\end{tabular}");
  fclose(ftable);
 
  truth->Scale(td/tt*1.0);
 
  TGraph* rawd=data_all.makeTG(npoint,firsti-1);
  TGraph* rawdb=data_bkgd.makeTG(npoint,firsti-1);
  TGraph* corrd=corrData.makeTGE(npoint,firsti-1);
  TGraph* corrdsys=corrDataSyst.makeTGE(npoint,firsti-1);
  TGraph* corrdbkgd=corrDataBkgd.makeTGE(npoint,firsti-1);
  // this is special
  TGraph* corrdsysb=new TGraphErrors(npoint,corrDataBkgd.xave+firsti-1,corrDataBkgd.y+firsti-1,corrDataSyst.xwidth+firsti-1,corrDataSyst.ey+firsti-1);

  TH2* dummy=new TH2F("dummy","dummy",20,-3.5,3.5,30,0,600);

  dummy->GetYaxis()->SetTitle("Events/0.1 Units of Rapidity");
  dummy->GetXaxis()->SetTitle("Y_{Z}");
  
  TCanvas* c1=new TCanvas("finalZRap","finalZRap",800,600);
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

  truth->Draw("SAMEHIST");

  TLegend* tl=new TLegend(0.4,0.16,0.72,0.34);
  tl->AddEntry(corrd,"Corrected Data","P");
  tl->AddEntry(rawdb,"Raw Data","P");
  tl->AddEntry(rawd,"Background-Subtracted Data","P");
  tl->AddEntry(truth,"POWHEG+CT10 Prediction","l");
  tl->Draw();

  zrap_Prelim(0.85,0.90,0.85,0.85);
  zrap_Lumi(0.85,0.86,lumi);

  sprintf(fnamework,"ZRapidity_final-%d.eps",lumi);
  c1->Print(fnamework);
  sprintf(fnamework,"ZRapidity_final-%d.png",lumi);
  c1->Print(fnamework);


  TCanvas* c2=new TCanvas("finalZRapFold","finalZRapFold",800,600);
  TH2* dummy2=new TH2F("dummy2","dummy2",10,0.0,3.5,30,0,0.4);
  dummy2->SetDirectory(0);

  TH1* truth2=(TH1*)(truth->Clone("truth2"));
  TH1* truthalt=(TH1*)(truth->Clone("truthalt"));
  truth2->SetDirectory(0);
  truthalt->SetDirectory(0);

  double chi2_0=0;
  double chi2_1=0;
  double fudge=1.00;

  for (int i=0; i<=(lasti-firsti)/2+2; i++) {

    int j=izero+i-1;
    int jadd=izero-i-2;

    truth2->SetBinContent(j+1,(truth2->GetBinContent(j+1)+truth2->GetBinContent(jadd+1))/sumtotal/binwidth*fudge);
    /*
    double factor=1+(pdfSens25.y[i-1]+pdfSens25.y[BINCOUNT-i])/2.0;
    truthalt->SetBinContent(i,truth2->GetBinContent(i)*factor);
    if (corrDataSystFold.ey[i-1]<=0) continue;
    printf("AltUniv: %d %f %f\n",i,factor,truthalt->GetBinContent(i));
    chi2_0+=pow(truth2->GetBinContent(i)-corrDataFold.y[i-1],2)/pow(corrDataSystFold.ey[i-1],2);
    chi2_1+=pow(truthalt->GetBinContent(i)-corrDataFold.y[i-1],2)/pow(corrDataSystFold.ey[i-1],2);
    */
  }
  

  c2->cd();
  dummy2->GetYaxis()->SetTitle("1/#sigma d#sigma/dY");
  dummy2->GetXaxis()->SetTitle("|Y_{Z}|");
  
  dummy2->Draw();

  TGraph* corrdfold=new TGraphErrors(npoint/2,corrDataFold.xave+izero-1,corrDataFold.y+izero-1,corrDataFold.xwidth+izero-1,corrDataFold.ey+izero-1);
  corrdfold->Draw("PSAME");
  TGraph* corrdfoldsys=new TGraphErrors(npoint/2,corrDataFold.xave+izero-1,corrDataFold.y+izero-1,corrDataSystFold.xwidth+izero-1,corrDataSystFold.ey+izero-1);
  corrdfold->Draw("PSAME");
  corrdfoldsys->SetMarkerStyle(0);
  corrdfoldsys->Draw("PSAME");

  truth2->Draw("SAMEHIST");
  truthalt->SetLineColor(kRed);
  truthalt->Draw("SAMEHIST");

  zrap_Prelim(0.3,0.24,0.3,0.2);
  zrap_Lumi(0.30,0.19,lumi);
  

  sprintf(fnamework,"ZRapidity_final_fold-%d.eps",lumi);
  c2->Print(fnamework);
  sprintf(fnamework,"ZRapidity_final_fold-%d.png",lumi);
  c2->Print(fnamework);

  TCanvas* c3=new TCanvas("finalZRapErrors","finalZRapErrors",800,600);
  TH1* dummy3=new TH1F("dummy3","dummy3",10,-3.5,3.5);
  dummy3->SetDirectory(0);
  c3->SetLogy();
  dummy3->SetMaximum(1.0);
  dummy3->SetMinimum(0.001);
  dummy3->GetYaxis()->SetTitle("Fractional Error");
  dummy3->GetXaxis()->SetTitle("Y_{Z}");


  
  //  TH1* binmiggr=binMigPos.makeTH1("bin migration",npoint,firsti-1);
  TH1* effacc_stat_gr=ea_statErr.makeTH1("eff stat",npoint,firsti-1);
  TH1* effacc_syst_gr=effSystErr.makeTH1("eff syst",npoint,firsti-1);
  TH1* gr_dataStatError=dataStatError.makeTH1("data stats err",npoint,firsti-1);
  TH1* gr_bkgd_unc=backgroundUncFrac.makeTH1("bkgd unc",npoint,firsti-1);
  //  TH1* gr_pdf=pdfFrac.makeTH1("pdf_err",npoint,firsti-1);
  TH1* gr_energyscale=energyScale.makeTH1("es_unc",npoint,firsti-1);
  TH1* gr_systotal=totalSyst.makeTH1("total_syst",npoint,firsti-1);
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

  sprintf(fnamework,"ZRapidity_final_errors-%d.eps",lumi);
  c3->Print(fnamework);
  sprintf(fnamework,"ZRapidity_final_errors-%d.png",lumi);
  c3->Print(fnamework);

sprintf(fnamework,"ZRapidity_final_fold_errtable-%d.tex",lumi);
  ftable=fopen(fnamework,"w");

  for (int i=0; i<=(lasti-firsti)/2; i++) {
    int j=izero+i-1;
    int jadd=izero-i-2;

    double ea_Ave=(effAcc.y[j]+effAcc.y[jadd])/2;

    //    double errmig=(fabs(binMig.y[j])*corrData.y[j]+fabs(binMig.y[jadd])*corrData.y[jadd])/sumtotal;
    //    double pdferr=(pdfFrac.y[j]+pdfFrac.y[jadd])/2;//*corrDataFold.y[j]*binwidth/sumtotal;
    double bkgderr=sqrt(pow(backgroundAllUnc.y[j],2)+pow(backgroundAllUnc.y[jadd],2))/ea_Ave/sumtotal;
    double staterr=(ea_statErr.y[j]+ea_statErr.y[jadd])/2;
    double systerr=(effSystErr.y[j]+effSystErr.y[jadd])/2;
    double energyerr=(energyScale.y[j]+energyScale.y[jadd])/2;

    //    errmig/=binwidth*corrDataFold.y[j];
    //    pdferr/=binwidth*corrDataFold.y[j];
    bkgderr/=binwidth*corrDataFold.y[j];
    //    staterr/=binwidth*corrDataFold.y[j];
    //    systerr/=binwidth*corrDataFold.y[j];


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


  //  printf(" Two models: %f %f \n",chi2_0,chi2_1);

}
