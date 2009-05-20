#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH2.h"


struct DataSeries {
  DataSeries() {
    clear();
  }
  void clear() {
    for (int i=0; i<44; i++) {
      xave[i]=0;
      xwidth[i]=0;
      y[i]=0;
      ey[i]=0;
    }
  }
  DataSeries(const DataSeries& ds) {
    for (int i=0; i<44; i++) {
      xave[i]=ds.xave[i];
      xwidth[i]=ds.xwidth[i];
      y[i]=ds.y[i];
      ey[i]=ds.ey[i];
    }
  }
  DataSeries(const char* fname, int iy=0);

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
    }
    hist->SetDirectory(0);
    return hist;
  }

  double xave[44];
  double xwidth[44];
  double y[44];
  double ey[44];
};  


#include "tdrstyle.C"
#include "zrapidityStandard.C"

DataSeries::DataSeries(const char* file, int iy) {
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
      printf("%d %d %d %f\n",i,iq,ate,q[iq]);
      atetotal+=ate;
    }
    if (found!=1) continue;
    if (i<1 || i>44) continue;
    xave[i-1]=(b+a)/2;
    xwidth[i-1]=(b-a)/2;
    y[i-1]=q[iy];
  }
  fclose(f);
}

void plotFinal() {
  setTDRStyle();
  TFile mctruth("histo_40M.root");

  TH1* truth=(TH1*)(mctruth.Get("mcEff/All/Z0_Y")->Clone("truth"));
  truth->SetDirectory(0);
  
  DataSeries data_ee("Tight-ECAL-Loose-ECAL.txt"), data_ef("Tight-ECAL-HF.txt"), data_all;
  
  DataSeries effAcc("effacc.txt"), ea_statErr("staterror_100pb.txt");
  DataSeries binMig("BinMigrationEffect.txt");
  DataSeries binMigPos(binMig);

  const int firsti=8;
  const int lasti=37;
  const int npoint=lasti-firsti+1;

  double tt=0, td=0;
 
  for (int i=0; i<44; i++) {
    data_all.xave[i]=data_ef.xave[i];
    data_all.xwidth[i]=data_ef.xwidth[i];
    data_all.y[i]=data_ef.y[i]+data_ee.y[i];
    binMigPos.y[i]=fabs(binMig.y[i]);
  }


  DataSeries dataStatError(data_all);

  DataSeries backgroundExp("background.txt",0);
  DataSeries backgroundDelta("background.txt",0);
  DataSeries backgroundUncFrac(data_all);

  DataSeries pdfPos("pdfErrsPos.txt");
  DataSeries pdfNeg("pdfErrsNeg.txt");
  DataSeries pdfTotal(data_all);
  DataSeries pdfFrac(data_all);

  TRandom workRand(50);


  backgroundExp.y[8-1]=workRand.Poisson(backgroundExp.y[10-1]);
  backgroundExp.y[9-1]=workRand.Poisson(backgroundExp.y[10-1]);
  backgroundExp.y[36-1]=workRand.Poisson(backgroundExp.y[35-1]);
  backgroundExp.y[37-1]=workRand.Poisson(backgroundExp.y[35-1]); 

  for (int i=firsti-1; i<lasti; i++) {
    backgroundDelta.y[i]=workRand.Poisson(backgroundDelta.y[i])-backgroundDelta.y[i];
  }

  DataSeries corrData(data_all);
  DataSeries corrDataSyst(data_all);
  DataSeries corrDataBkgd(data_all);
  DataSeries data_bkgd(data_all);
  DataSeries totalSyst(data_all);

  double sumtotal=0;
  
  for (int i=firsti-1; i<lasti; i++) {
    if (effAcc.y[i]==0) {
      corrData.y[i]=0;
      continue;
    }
    corrData.y[i]/=effAcc.y[i];
    corrData.ey[i]=sqrt(data_all.y[i])/effAcc.y[i];
    //    corrDataSyst.xave[i]+=0.02; // uncomment to offset error bars


    tt+=truth->GetBinContent(i+1);
    td+=corrData.y[i];

    data_bkgd.y[i]+=backgroundExp.y[i];

    backgroundUncFrac.y[i]=sqrt(backgroundExp.y[i])/data_bkgd.y[i];

    corrDataBkgd.y[i]+=backgroundDelta.y[i];
    corrDataBkgd.y[i]/=effAcc.y[i];
    corrDataBkgd.ey[i]=sqrt(data_all.y[i]+backgroundExp.y[i])/effAcc.y[i];

    dataStatError.y[i]= corrDataBkgd.ey[i]/corrDataBkgd.y[i];

    corrDataSyst.y[i]/=effAcc.y[i];

    pdfTotal.y[i]=(pdfPos.y[i]+(-pdfNeg.y[i]))/2*corrDataBkgd.y[i];
    pdfFrac.y[i]=(pdfPos.y[i]+(-pdfNeg.y[i]))/2;

    double errmig=fabs(binMig.y[i])*corrDataBkgd.y[i];

    printf("%d %f %f %f %f\n",i,corrData.ey[i],corrData.y[i]*ea_statErr.y[i],errmig,pdfTotal.y[i]);
    corrDataSyst.ey[i]=sqrt(pow(corrDataBkgd.ey[i],2)+
			    pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			    fabs(backgroundExp.y[i])/effAcc.y[i]+
			    pow(errmig,2)+
			    pow(pdfTotal.y[i],2)
			    );

    totalSyst.ey[i]=sqrt( pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			  pow(errmig,2)+
			  fabs(backgroundExp.y[i])/effAcc.y[i]+
			  pow(pdfTotal.y[i],2)
			  );

    totalSyst.y[i]=totalSyst.ey[i]/corrDataBkgd.y[i];

    sumtotal+=corrDataSyst.y[i];
  }


  DataSeries corrDataFold(corrData),corrDataSystFold(corrData);

  FILE* ftable=fopen("ZRapidity_final_fold_table.tex","w");

  fprintf(ftable,"\\begin{tabular}{|c|c||c|c|c|}\n\\hline\n");
  fprintf(ftable,"$|Y_{min}|$ & $|Y_{max}|$ & Measurement & Statistical Error & Systematic Error \\\\ \\hline\n");
  
  const int izero=23;
  for (int i=0; i<=(lasti-firsti)/2; i++) {
    int j=izero+i-1;
    int jadd=izero-i-2;
    corrDataFold.y[j]=(corrDataBkgd.y[j]+corrDataBkgd.y[jadd])/sumtotal;


    //    printf("%d %d (%f) %d (%f) \n",i,j,corrDataBkgd.y[j],jadd,corrDataBkgd.y[jadd]);

    double ea_Ave=(effAcc.y[j]+effAcc.y[jadd])/2;

    double errmig=(fabs(binMig.y[j])*corrData.y[j]+fabs(binMig.y[jadd])*corrData.y[jadd])/sumtotal;
    double pdferr=(pdfFrac.y[j]+pdfFrac.y[jadd])/2*corrDataFold.y[j]/sumtotal;
    double bkgderr=sqrt(backgroundExp.y[j]+backgroundExp.y[jadd])/sumtotal;

    corrDataFold.ey[j]=sqrt(data_all.y[j]+data_all.y[jadd]+backgroundExp.y[j]+backgroundExp.y[jadd])/ea_Ave/sumtotal;

    corrDataSystFold.y[j]=corrDataFold.y[j];
    corrDataSystFold.ey[j]=sqrt(pow(corrDataFold.ey[j],2)+
				pow(corrDataFold.y[j]*(ea_statErr.y[j]+ea_statErr.y[jadd])/2,2)+
				pow(errmig,2)+
				pow(pdferr,2)+
				pow(bkgderr,2)
				);

    fprintf(ftable,"%7.2f & %7.2f & %7.5f & %7.5f & %7.5f \\\\ \n",
	    corrDataFold.xave[j]-corrDataFold.xwidth[j],
	    corrDataFold.xave[j]+corrDataFold.xwidth[j],
	    corrDataFold.y[j],
	    corrDataFold.ey[j],
	    sqrt(pow(corrDataSystFold.ey[j],2)-pow(corrDataFold.ey[j],2))
	    );
  }

  fprintf(ftable,"\\hline\n\\end{tabular}");
  fclose(ftable);

  truth->Scale(td/tt);

  TGraph* rawd=data_all.makeTG(npoint,firsti-1);
  TGraph* rawdb=data_bkgd.makeTG(npoint,firsti-1);
  TGraph* corrd=corrData.makeTGE(npoint,firsti-1);
  TGraph* corrdsys=corrDataSyst.makeTGE(npoint,firsti-1);
  TGraph* corrdbkgd=corrDataBkgd.makeTGE(npoint,firsti-1);
  // this is special
  TGraph* corrdsysb=new TGraphErrors(npoint,corrDataBkgd.xave+firsti-1,corrDataBkgd.y+firsti-1,corrDataSyst.xwidth+firsti-1,corrDataSyst.ey+firsti-1);

  TH2* dummy=new TH2F("dummy","dummy",20,-3.75,3.75,30,0,5500);

  dummy->GetYaxis()->SetTitle("Events/0.25 Units of Rapidity");
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
  tl->AddEntry(truth,"CTEQ6.1 Prediction","l");
  tl->Draw();

  zrap_Prelim(0.85,0.89,0.85,0.85);

  c1->Print("ZRapidity_final.eps");
  c1->Print("ZRapidity_final.png");

  TCanvas* c2=new TCanvas("finalZRapFold","finalZRapFold",800,600);
  TH2* dummy2=new TH2F("dummy2","dummy2",10,0.0,3.75,30,0,0.1);
  dummy2->SetDirectory(0);

  TH1* truth2=(TH1*)(truth->Clone("truth2"));
  truth2->SetDirectory(0);
  for (int i=izero; i<=lasti; i++) 
    truth2->SetBinContent(i,(truth2->GetBinContent(i)+truth2->GetBinContent(44-i+1))/sumtotal);

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

  zrap_Prelim(0.3,0.24,0.3,0.2);

  c2->Print("ZRapidity_final_fold.eps");
  c2->Print("ZRapidity_final_fold.png");


  TCanvas* c3=new TCanvas("finalZRapErrors","finalZRapErrors",800,600);
  TH1* dummy3=new TH1F("dummy3","dummy3",10,-3.76,3.75);
  dummy3->SetDirectory(0);
  c3->SetLogy();
  dummy3->SetMaximum(1.0);
  dummy3->SetMinimum(0.001);
  dummy3->GetYaxis()->SetTitle("Fractional Error");
  dummy3->GetXaxis()->SetTitle("Y_{Z}");


  
  TH1* binmiggr=binMigPos.makeTH1("bin migration",npoint,firsti-1);
  TH1* effacc_stat_gr=ea_statErr.makeTH1("eff stat",npoint,firsti-1);
  TH1* gr_dataStatError=dataStatError.makeTH1("data stats err",npoint,firsti-1);
  TH1* gr_bkgd_unc=backgroundUncFrac.makeTH1("bkgd unc",npoint,firsti-1);
  TH1* gr_pdf=pdfFrac.makeTH1("pdf_err",npoint,firsti-1);
  TH1* gr_systotal=totalSyst.makeTH1("total_syst",npoint,firsti-1);
  dummy3->Draw();

  binmiggr->SetLineWidth(2);
  binmiggr->SetLineColor(kBlue);

  effacc_stat_gr->SetLineWidth(2);
  effacc_stat_gr->SetLineColor(kGreen+2);

  gr_dataStatError->SetLineWidth(2);
  gr_dataStatError->SetLineColor(kBlack);
  gr_dataStatError->SetLineStyle(2);

  gr_bkgd_unc->SetLineWidth(2);
  gr_bkgd_unc->SetLineColor(kRed-4);

  gr_pdf->SetLineWidth(2);
  gr_pdf->SetLineColor(kViolet-5);

  gr_systotal->SetLineWidth(2);
  

  binmiggr->Draw("HSAME");
  effacc_stat_gr->Draw("HSAME");
  gr_dataStatError->Draw("HSAME");
  gr_bkgd_unc->Draw("HSAME");
  gr_pdf->Draw("HSAME");
  gr_systotal->Draw("HISTSAME");

  TLegend* tl=new TLegend(0.35,0.70,0.75,0.94,"Error Source Contributions");
  tl->AddEntry(binmiggr,"Bin Migration","L");
  tl->AddEntry(effacc_stat_gr,"Efficiency Statistics","L");
  tl->AddEntry(gr_pdf,"PDF (#epsilon #times A) Uncertainty","L");
  tl->AddEntry(gr_bkgd_unc,"Background Uncertainty","L");
  tl->AddEntry(gr_systotal,"Total Systematic Error","L");
  tl->AddEntry(gr_dataStatError,"Data Statistics","L");
  tl->Draw();
  
  zrap_Prelim(0.3,0.24,0.3,0.2);

  c3->Print("ZRapidity_final_errors.eps");
  c3->Print("ZRapidity_final_errors.png");

}
