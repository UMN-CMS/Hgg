#include "tdrstyle.C"
#include "zrapidityStandard.C"

TCanvas* c3;
int nplot=1;

double getChi2(int ecal, const char* pcode, double value, int sf, TH1* data) {
  char name[1024];
  const char* dcode, *item;
  
  if (ecal==0) {
    dcode="EB";
    item="mcEff/Golden-EB-EB/C07-HLT-EtaDet/Z0_mass";
  } else if (ecal==1) {
    dcode="EE";
    item="mcEff/Tight-ECAL-Loose-EE/C07-HLT-EtaDet/Z0_mass";
  } else if (ecal==2) {
    dcode="HF";
    item="mcEff/Tight-ECAL-HF/C07-HLT-EtaDet/Z0_mass";
  }

  int v=int(value);
  int sub=(int(value*pow(10,sf)+0.5))%(int(pow(10,sf)));

  //  if (v==0) 
    //    sprintf(name,"k_ecalScan%s_%s__%0*d.root",dcode,pcode,sf,sub);
  //  else 
    sprintf(name,"k3_ecalScan%s_%s_%d_%0*d.root",dcode,pcode,v,sf,sub);

  TFile* f=TFile::Open(name);
  if (f==0) return 0;
  TH1* hist=(f->Get(item))->Clone(name);
  hist->SetDirectory(0);

  double chi2=0;

  int ilow=11; //11
  int ihigh=28; //30

  hist->Scale(data->Integral(ilow,ihigh)/hist->Integral(ilow,ihigh));

  int ndof=0;
  for (int i=ilow;i<=ihigh;i++){
    double s2exp = (hist->GetBinContent(i));
    double s2obs= (data->GetBinContent(i));   
      
    if (s2exp==0){
	chi2+=0;
    }else{
      chi2+=((s2obs-s2exp)*(s2obs-s2exp))/((sqrt(s2exp))*(sqrt(s2exp)));
      ndof++;
    }
  }
  //  c3->cd();
  hist->Draw("SAMEHIST");

  delete f;
  return chi2/ndof;
}


void smearFit(TFile* mc1,int ecal,int param, const char* ofname=0){
  setTDRStyle();

  if(ecal==0){
    TH1F *da = mc1->Get("ZFromData/Golden-EB-EB/C07-HLT-EtaDet/Z0_mass")->Clone("blah");
  }else if(ecal==1){
    TH1F *da = mc1->Get("ZFromData/Tight-ECAL-Loose-EE/C07-HLT-EtaDet/Z0_mass")->Clone("blah");
  }else if(ecal==2){
    TH1F *da = mc1->Get("ZFromData/Tight-ECAL-HF/C07-HLT-EtaDet/Z0_mass")->Clone("blah");
  }

  TCanvas *c1 =new TCanvas("c1","c1",600,600);
  c3 =new TCanvas("c3","c3",600,600);

  da->SetMarkerStyle(24);
  da->Draw("E0");
  
  int i;
  int j;
  int n=7;
  double chi2[100]={0,0,0,0,0,0,0};  
  
  TF1* myF=new TF1("myF","[0]*pow(x-[1],2)+[2]");
    myF->SetParameter(1,1.40);
    myF->SetParameter(0,10.0);
    myF->SetParameter(2,1.0);
    
  if((ecal==0)&&(param==0)){
    //    Double_t xc[100]={0.014, 0.016, 0.018, 0.020,0.022,0.024,0.026,0.030,0.035,0.040,0.045, -1000};
    Double_t xc[100]={0.020,0.025,0.030,0.035,0.040,0.045,0.050,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"p0",xc[i],3,da);
    }
    n=i;
    
    c1->cd();
    myF->SetParameter(1,0.040);
    myF->SetParLimits(1,xc[0],xc[n-1]);
   
    gr = new TGraph(n,xc,chi2);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);
    gr->GetXaxis()->SetTitle("EB a");
  }else if((ecal==0)&&(param==1)){
    Double_t xc[100]={-.036, -.106,-.176,-.246, -.316, -.386,-.456};
    c1->cd();
    gr = new TGraph(n,xc,chi2); 
    gr->Fit(myF,"","",-.386,-.066);
    gr->GetXaxis()->SetTitle("EB p1");
  }else if((ecal==0)&&(param==2)){
    Double_t xc[100]={0.493, 0.593, 0.693,0.793,.893, 0.993, 1.093};
    c1->cd();
    gr = new TGraph(n,xc,chi2);
    gr->Fit(myF,"","",0.443,1.093);
    gr->GetXaxis()->SetTitle("EB p2");
  }else if((ecal==0)&&(param==3)){
    //    Double_t xc[100]={.002, .005, .008,.014, .017, .020, -1000};
    Double_t xc[100]={ 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0045, 0.0050, 0.0055, 0.0060, 0.0065, 0.0070, 0.0075, 0.0080, 0.0090, 0.0095, 0.0100,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"c",xc[i],4,da);
    }
    n=i;

    c1->cd();
    gr = new TGraph(n,xc,chi2);
    myF->SetParameter(1,.005);

    myF->SetParLimits(1,xc[0],xc[n-1]);

    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",0.001,0.02);
    gr->GetXaxis()->SetTitle("EB c");
  }else if((ecal==0)&&(param==4)){
    Double_t xc[100]={1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"alpha",xc[i],1,da);
    }
    n=i;
    

    c1->cd();
    gr = new TGraph(n,xc,chi2);
    myF->SetParameter(1,1.9);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    for (int pass=0; pass<2; pass++)
      gr->Fit(myF,"","",xc[pass],xc[n-3]);
    gr->GetXaxis()->SetTitle("EB #alpha");


  }else if((ecal==0)&&(param==5)){
    Double_t xc[100]={2.0,2.2,2.4,2.8,3.0,3.2,3.4,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"n",xc[i],1,da);
    }
    n=i;


    c1->cd();
    gr = new TGraph(n,xc,chi2);
    myF->SetParameter(1,2.40);
    myF->SetParameter(0,10.0);
    myF->SetParameter(2,1.0);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    gr->Fit(myF,"","",1.5,2.6);
    gr->GetXaxis()->SetTitle("EB n");

  }else if((ecal==0)&&(param==6)){
    Double_t xc[100]={0.992,0.993,0.994,0.995,0.996,0.997,0.998,0.999,1.000,-1000};
    //    Double_t xc[100]={0.996,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"mean",xc[i],3,da);
    }
    n=i;

    c1->cd();
    gr = new TGraph(n,xc,chi2);
    myF->SetParLimits(0,0,10000000);
    myF->SetParameter(1,1.00);
    myF->SetParameter(0,1.0);
    myF->SetParameter(2,2.2);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",0.99,1.002);
    gr->GetXaxis()->SetTitle("EB #bar{x}");
    gr->GetXaxis()->SetLabelSize(0.04);
  }else if((ecal==1)&&(param==0)){
    Double_t xc[100]={0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,-1000};
    //Double_t xc[100]={0.8,-1000};
    
    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"p0",xc[i],1,da);
    }
    n=i;
    
    c1->cd();
    myF->SetParameter(1,1.0);
    gr = new TGraph(n,xc,chi2);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);
    gr->GetXaxis()->SetTitle("EE a");
  }else if((ecal==1)&&(param==1)){
    Double_t xc[100]={0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,-1000};
    //Double_t xc[100]={0.8,-1000};
    
    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"p1",xc[i],1,da);
    }
    n=i;
    
    c1->cd();
    myF->SetParameter(1,0.8);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    gr = new TGraph(n,xc,chi2);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);


    gr->GetXaxis()->SetTitle("EE p1");
  }else if((ecal==1)&&(param==2)){
    Double_t xc[100]={0.05,0.10,0.15,0.20,0.25,0.3,-1000};
    //Double_t xc[100]={0.8,-1000};
    
    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"p2",xc[i],2,da);
    }
    n=i;
    
    c1->cd();
    myF->SetParameter(1,0.8);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    gr = new TGraph(n,xc,chi2);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);

    gr->GetXaxis()->SetTitle("EE p2");
  }else if((ecal==1)&&(param==4)){
    Double_t xc[100]={1.1,1.2, 1.3, 1.4, 1.5,1.6,1.7,1.8,1.9,-1000};
    gr = new TGraph(n,xc,chi2);
    
    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"alpha",xc[i],1,da);
    }
    n=i;

    myF->SetParameter(1,1.3);
    myF->SetParameter(0,1.0);
    myF->SetParameter(2,2.99);
    
    c1->cd();
    myF->SetParLimits(1,xc[0],xc[n-1]);
    myF->SetParLimits(0,0,100000);
    gr = new TGraph(n,xc,chi2);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);
    gr->GetXaxis()->SetTitle("EE #alpha");
  }else if((ecal==1)&&(param==5)){
    Double_t xc[100]={2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"n",xc[i],1,da);
    }
    n=i;
    
    c1->cd();
    gr = new TGraph(n,xc,chi2);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);
    gr->GetXaxis()->SetTitle("EE n");
  } else if((ecal==1)&&(param==6)){
    //    Double_t xc[100]={1.0091,1.010,1.011,1.012,1.013,1.014,1.015,1.016,1.017,-1000};
    Double_t xc[100]={0.985,0.987,0.991,0.992,0.993,0.994,0.995,0.996,0.997,0.998,0.999,1.000,1.001,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"mean",xc[i],3,da);
    }
    n=i;
    
    c1->cd();
    gr = new TGraph(n,xc,chi2);

    myF->SetParameter(1,1.013);
    myF->SetParameter(0,2.0);
    myF->SetParameter(2,2.99);
    myF->SetParLimits(1,xc[0],xc[n-1]);

    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);
    gr->GetXaxis()->SetTitle("EE #bar{x}");
    gr->GetXaxis()->SetLabelSize(0.035);
  }else if((ecal==2)&&(param==2)){
    Double_t xc[100]={0.97,0.98,0.99,1.0,1.01,1.02,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"mean",xc[i],2,da);
    }
    n=i;
    
    c1->cd();
    myF->SetParameter(1,0.99);
    gr = new TGraph(n,xc,chi2);
    myF->SetParLimits(1,xc[0],xc[n-1]);

    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);
    gr->GetXaxis()->SetTitle("HF mean");
    gr->GetXaxis()->SetLabelSize(0.035);
  }else if((ecal==2)&&(param==0)){
    Double_t xc[100]={1.90,1.95,2.00,2.05,2.10,2.15,2.20,2.25,2.30,2.35,2.40,2.45,2.50,2.55,2.60,2.65,2.70,2.75,2.80,2.85,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"stocastic",xc[i],2,da);
    }
    n=i;
    
    c1->cd();
    myF->SetParameter(1,0.055);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    gr = new TGraph(n,xc,chi2);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-1]);
    gr->GetXaxis()->SetTitle("HF stocastic");
  }else if((ecal==2)&&(param==1)){
    Double_t xc[100]={0.035,0.040,0.045,0.050,0.055,0.060,0.065,0.070,-1000};

    for (i=0; xc[i]>-999; i++) {
      chi2[i]=getChi2(ecal,"constant",xc[i],3,da);
    }
    n=i;
    
    c1->cd();
    myF->SetParameter(1,0.99);
    myF->SetParameter(0,2.0);
    myF->SetParameter(2,2.99);
    myF->SetParLimits(1,xc[0],xc[n-1]);
    gr = new TGraph(n,xc,chi2);
    for (int pass=0; pass<3; pass++)
      gr->Fit(myF,"","",xc[0],xc[n-2]);
    gr->GetXaxis()->SetTitle("HF constant");
    gr->GetXaxis()->SetLabelSize(0.035);
  } else{
    std::cout<<"EB=0, EE=1:p0-2=0-2,c=3,alpha=4,n=5, EE has no c! "<<std::endl;
    return;
  }
  
  
  
  gr->GetYaxis()->SetTitle("#chi^{2}/dof");
  gr->SetTitle(0);
  
  double pram,err;
  pram= myF->GetParameter(1);
  err=myF->GetParError(1);
  
  char text[50];
  sprintf(text,"Fit = a(x-[%.3f #pm %.1e])^{2}+b",pram,err);
  std::cout<<pram<<" "<<err<<" "<<text<<std::endl; 
  tr = new TLatex(.3,.81,text);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gr->Draw("A*");
  tr->SetNDC();
  tr->SetTextFont(42);
  tr->SetTextSize(0.04);
  tr->Draw("Same A*");
  c1->RaiseWindow();

  zrap_Prelim(0.4,0.9);

  if (ofname!=0) c1->Print(ofname);

}
