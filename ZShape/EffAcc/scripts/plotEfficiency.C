/*
  Must run ".x EfficiencyStatistics.cc+" before running
 */
void determineText(const char* name, char* title, char* xaxisTitle, char*& fname);
void plotEfficiency(const char* filename) {

  gROOT->SetStyle("Plain");

  FILE* f=fopen(filename,"r");
  char buffer[1024],title[40],xaxisTitle[40];
  char* ofilename=0;

  title[0]=0;
  xaxisTitle[0]=0;

  float x[50],y[50],exl[50],exh[50],eyl[50],eyh[50];
  int nitems=0;

  float bmin,bmax,val,systP,systM,denom;
  int bin;

  while (!feof(f)) {
    buffer[0]=0;
    fgets(buffer,1024,f);
    if (strstr(buffer,"name")!=0) {
      determineText(strstr(buffer,"name")+5,title,xaxisTitle,ofilename);
    }
    int found=sscanf(buffer, "%d %f %f %f %f %f %f %f", &bin,
		     &bmin, &bmax,&val,&systP,&systM,&denom);
    if (found==7 && denom>0) {
      x[nitems]=(bmin+bmax)/2.0;
      exl[nitems]=x[nitems]-bmin;
      exh[nitems]=-x[nitems]+bmax;
      y[nitems]=val;
      EfficiencyStatistics es(double(val),int(denom+0.5));
      eyl[nitems]=val-es.sigma(-1);
      eyh[nitems]=-val+es.sigma(1);
      //      printf("%f %f %f %f\n",x[nitems],y[nitems],eyl[nitems],eyh[nitems]);
      nitems++;
    } 
  }
  fclose(f);

  printf("%d\n",nitems);

  
  TCanvas* c1=new TCanvas("c1","c1",600,600);

  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.17);

  TGraphAsymmErrors* theGraph=new TGraphAsymmErrors(nitems,x,y,exl,exh,eyl,eyh);
  theGraph->SetMaximum(1.005);
  theGraph->SetMarkerStyle(24);
  theGraph->SetTitle(title);
  theGraph->GetXaxis()->SetTitle(xaxisTitle);
  theGraph->GetYaxis()->SetTitle("Single-Electron Efficiency");
  theGraph->GetYaxis()->SetTitleOffset(1.9);
  theGraph->GetYaxis()->CenterTitle(true);
  theGraph->GetXaxis()->CenterTitle(true);
  gStyle->SetTitleX(0.48);
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleY(0.98);
  theGraph->Draw("EAP");
  theGraph->Draw("E");

  if (ofilename!=0) {
    c1->Print(ofilename);
  }

}

void determineText(const char* name, char* title, char* xaxisTitle, char*& ofname) {
  if (strstr(name,"Iso")!=0) strcpy(title,"Isolation");
  if (strstr(name,"ElectronId")!=0) strcpy(title,"Electron Identification");
  if (strstr(name,"GsfTrack")!=0) strcpy(title,"GSF Track Matching");
  if (strstr(name,"EtaDet")!=0) strcpy(xaxisTitle,"#eta_{detector}");

  if (strstr(name,"Iso")!=0 && strstr(name,"EtaDet")!=0) {
    ofname="Eff-Iso-EtaDet.eps";
  }
  if (strstr(name,"ElectronId")!=0 && strstr(name,"EtaDet")!=0) {
    ofname="Eff-ElectronId-EtaDet.eps";
  }
  if (strstr(name,"GsfTrack")!=0 && strstr(name,"EtaDet")!=0) {
    ofname="Eff-GsfTrack-EtaDet.eps";
  }

}
