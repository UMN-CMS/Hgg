struct DataSeries {
  DataSeries() {
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

  double xave[44];
  double xwidth[44];
  double y[44];
  double ey[44];
};  

#include "tdrstyle.C"
#include "zrapidityStandard.C"

void loadData(const char* file, DataSeries& ds) {
  setTDRStyle();
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
    float a,b,c;
    if (sscanf(line," %d %f %f %f ",&i,&a,&b,&c)!=4) continue;
    if (i<1 || i>44) continue;
    ds.xave[i-1]=(b+a)/2;
    ds.xwidth[i-1]=(b-a);
    ds.y[i-1]=c;
  }
  fclose(f);
}

void plotFinal() {
  TFile mctruth("histo_40M.root");

  TH1* truth=mctruth.Get("mcEff/All/Z0_Y")->Clone("truth");
  truth->SetDirectory(0);
  
  DataSeries data_ee, data_ef, data_all;
  
  loadData("Tight-ECAL-HF.txt",data_ef);
  loadData("Tight-ECAL-Loose-ECAL.txt",data_ee);

  DataSeries effAcc, ea_statErr;

  loadData("effacc.txt",effAcc);
  loadData("staterror_100pb.txt",ea_statErr);

  const int firsti=8;
  const int lasti=37;
  const int npoint=lasti-firsti+1;

  double tt=0, td=0;
 
  for (int i=0; i<44; i++) {
    data_all.xave[i]=data_ef.xave[i];
    data_all.xwidth[i]=data_ef.xwidth[i];
    data_all.y[i]=data_ef.y[i]+data_ee.y[i];
  }

  DataSeries corrData(data_all);
  DataSeries corrDataSyst(data_all);

  for (int i=firsti-1; i<lasti; i++) {
    corrData.y[i]/=effAcc.y[i];
    corrData.xwidth[i]/=2;
    corrData.ey[i]=sqrt(data_all.y[i])/effAcc.y[i];
    //    corrDataSyst.xave[i]+=0.02; // uncomment to offset error bars
    corrDataSyst.y[i]/=effAcc.y[i];
    corrDataSyst.xwidth[i]/=2;
    printf("%d %f %f\n",i,corrData.ey[i],corrData.y[i]*ea_statErr.y[i]);
    corrDataSyst.ey[i]=sqrt(pow(corrData.ey[i],2)+
			    pow(corrData.y[i]*ea_statErr.y[i],2)
			    );

    tt+=truth->GetBinContent(i+1);
    td+=corrData.y[i];
  }

  DataSeries corrDataFold(corrData),corrDataSystFold(corrData);

  const int izero=23;
  for (int i=0; i<=(lasti-firsti)/2; i++) {
    int j=izero+i-1;
    int jadd=izero-i-2;
    corrDataFold.y[j]=corrData.y[j]+corrData.y[jadd];

    printf("%d %d (%f) %d (%f) \n",i,j,corrData.y[j],jadd,corrData.y[jadd]);

    double ea_Ave=(effAcc.y[j]+effAcc.y[jadd])/2;

    corrDataFold.ey[j]=sqrt(data_all.y[j]+data_all.y[jadd])/ea_Ave;

    corrDataSystFold.y[j]=corrDataFold.y[j];
    corrDataSystFold.ey[j]=sqrt(pow(corrDataFold.ey[j],2)+
				pow(corrDataFold.y[j]*(ea_statErr.y[j]+ea_statErr.y[jadd])/2,2)
				
				);
  }

  truth->Scale(td/tt);

  TGraph* rawd=new TGraph(npoint,data_all.xave+firsti-1,data_all.y+firsti-1);
  TGraph* corrd=new TGraphErrors(npoint,corrData.xave+firsti-1,corrData.y+firsti-1,corrData.xwidth+firsti-1,corrData.ey+firsti-1);
  TGraph* corrdsys=new TGraphErrors(npoint,corrDataSyst.xave+firsti-1,corrDataSyst.y+firsti-1,corrDataSyst.xwidth+firsti-1,corrDataSyst.ey+firsti-1);

  TH2* dummy=new TH2F("dummy","dummy",20,-3.75,3.75,30,0,5500);

  dummy->GetYaxis()->SetTitle("Events/0.25 Units of Rapidity");
  dummy->GetXaxis()->SetTitle("Y_{Z}");
  
  TCanvas* c1=new TCanvas("finalZRap","finalZRap",800,600);
  dummy->Draw();
  dummy->SetDirectory(0);

  rawd->SetMarkerStyle(24);
  rawd->Draw("PSAME");

  corrd->SetMarkerStyle(20);
  corrd->Draw("PSAME");

  corrdsys->SetMarkerStyle(0);
  corrdsys->Draw("PSAME");

  truth->Draw("SAMEHIST");

  TLegend* tl=new TLegend(0.4,0.16,0.72,0.3);
  tl->AddEntry(corrd,"Corrected Data","P");
  tl->AddEntry(rawd,"Raw Data","P");
  tl->AddEntry(truth,"CTEQ6.1 Prediction","l");
  tl->Draw();

  zrap_Prelim(0.85,0.89,0.85,0.85);

  c1->Print("ZRapidity_final.eps");
  c1->Print("ZRapidity_final.png");

  TCanvas* c2=new TCanvas("finalZRapFold","finalZRapFold",800,600);
  TH2* dummy2=new TH2F("dummy2","dummy2",10,0.0,3.75,30,0,10000);
  dummy2->SetDirectory(0);

  TH1* truth2=truth->Clone("truth2");
  truth2->SetDirectory(0);
  for (int i=izero; i<=lasti; i++) 
    truth2->SetBinContent(i,truth2->GetBinContent(i)+truth2->GetBinContent(44-i+1));

  c2->cd();
  dummy2->GetYaxis()->SetTitle("Events/0.25 Units of Rapidity");
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

}
