#include "readStandardFile.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"

int colorForId(int i) {
  const int colors[] = {kRed, kBlue, kGreen+2,
			kViolet, kGray
  };
  return colors[i];
}


void compareQT(bool limiNormal, const char* models,const char* outgoing=0) {
  setTDRStyle();

  const char* knownModels[] = { 
    "D6T","Pythia Tune D6T","mcFile_D6T_TL.txt",
    "P0","Pythia P0","mcFile_P0_TL.txt",
    "ProPT0","Pythia ProPT0","mcFile_ProPT0_TL.txt",
    "ProQ20","Pythia ProQ20","mcFile_ProQ20_TL.txt",
    "Z2","Pythia Tune Z2","mcFile_Z2_TL.txt",
    "R0.50","Resbos #alpha=0.50","ptShape_Resbos_a0_50.txt",
    "R0.70","Resbos #alpha=0.70","ptShape_Resbos_a0_70.txt",
    "R0.75","Resbos #alpha=0.75","ptShape_Resbos_a0_75.txt",
    "R0.80","Resbos #alpha=0.80","ptShape_Resbos_a0_80.txt",
    "R0.85","Resbos #alpha=0.85","ptShape_Resbos_a0_85.txt",
    "R0.90","Resbos #alpha=0.90","ptShape_Resbos_a0_90.txt",
    "R0.95","Resbos #alpha=0.95","ptShape_Resbos_a0_95.txt",
    "R1.00","Resbos #alpha=1.00","ptShape_Resbos_a1_0.txt",
    "R1.25","Resbos #alpha=1.25","ptShape_Resbos_a1_25.txt",
    "R1.50","Resbos #alpha=1.50","ptShape_Resbos_a1_5.txt",
    "R1.75","Resbos #alpha=1.75","ptShape_Resbos_a1_75.txt",
    "R2.00","Resbos #alpha=2.00","ptShape_Resbos_a2_0.txt",
    "R2.25","Resbos #alpha=2.25","ptShape_Resbos_a2_25.txt",
    "R2.50","Resbos #alpha=2.50","ptShape_Resbos_a2_5.txt",
    "R2.75","Resbos #alpha=2.75","ptShape_Resbos_a2_75.txt",
    "R3.00","Resbos #alpha=3.00","ptShape_Resbos_a3_0.txt",
    "R3.25","Resbos #alpha=3.25","ptShape_Resbos_a3_25.txt",
    "R3.50","Resbos #alpha=3.50","ptShape_Resbos_a3_5.txt",
    "R3.75","Resbos #alpha=3.75","ptShape_Resbos_a3_75.txt",
    "R4.00","Resbos #alpha=4.00","ptShape_Resbos_a4_0.txt",
    0,0,0
    };
  int nhists=0;
  TH1* hists[10],*histsDelta[10];
  int imodels[10];
  char names[10][10];
  
  // find the models
  int nmatch=sscanf(models,"%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,]",
		    names[0],
		    names[1],
		    names[2],
		    names[3],
		    names[4],
		    names[5],
  		    names[6],
		    names[7],
		    names[8],
		    names[9]);
  for (int i=0; i<nmatch; i++) {
    int jmodel=0;
    for (jmodel=0; knownModels[jmodel*3]!=0; jmodel++)
      if (!strcmp(knownModels[jmodel*3],names[i])) break;
    if (knownModels[jmodel*3]!=0) {
      imodels[nhists]=jmodel;
      nhists++;
    } else {
      printf("Didn't know: '%s'\n",names[i]);
    }
  }


  

  TCanvas* c1=new TCanvas("c1","c1",800,800);

  TH1* data=readStandardFilePt("data","../data/dsdqt_combined.csv");
  data->SetMaximum(0.08);
  data->SetMinimum(0.00);
  data->GetXaxis()->SetRangeUser(0.0,29.0);
  data->GetXaxis()->SetTitle("q_{T} [GeV/c]");
  data->GetYaxis()->SetTitle("1/#sigma d#sigma/dq_{T} [(GeV/c)^{-1}]");
  data->GetXaxis()->CenterTitle();
  data->Draw("E");

  char* Leg="";
  if (limiNormal) Leg="MC Normalized over 0 < q_{T} < 30";

  TLegend* tl=new TLegend(0.4,0.65,0.9,0.9,Leg);
  TLegend* tl2=new TLegend(0.4,0.75,0.9,0.9,Leg);
  //  tl->AddEntry(data,"Combined e + #mu","p");

  
  for (int i=0; i<nhists; i++) {
    char fullfname[1024];
    sprintf(fullfname,"../theory/%s",knownModels[imodels[i]*3+2]);
    printf("%s\n",fullfname);
    hists[i]=readStandardFilePt(knownModels[imodels[i]*3],fullfname);

    if (strstr(fullfname,"Resbos")) {
      double sum=hists[i]->Integral();
      for (int j=1; j<18; j++) 
	hists[i]->SetBinContent(j,hists[i]->GetBinContent(j)/hists[i]->GetBinWidth(j)/sum);
    }
    sprintf(fullfname,"%s-delta",knownModels[imodels[i]*3+2]);
    histsDelta[i]=(TH1*)(hists[i]->Clone(fullfname));
    
    for (int j=1; j<18; j++) {
      histsDelta[i]->SetBinContent(j,data->GetBinContent(j)/hists[i]->GetBinContent(j)-1);
      histsDelta[i]->SetBinError(j,data->GetBinError(j)/hists[i]->GetBinContent(j));
    }
    
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(colorForId(i));

    histsDelta[i]->SetLineWidth(2);
    histsDelta[i]->SetLineColor(colorForId(i));

    tl->AddEntry(hists[i],knownModels[imodels[i]*3+1],"L");
    tl2->AddEntry(hists[i],knownModels[imodels[i]*3+1],"L");

    if (limiNormal) 
      hists[i]->Scale(data->Integral(1,9,"width")/hists[i]->Integral(1,9,"width"));

    
    hists[i]->Draw("HIST SAME");
  }
  tl->Draw();

  zrap_Prelim(0.4,0.25,0.4,0.17);
  zrap_Lumi(0.40,0.21,36);

  TCanvas* c2=new TCanvas("c2","c2",800,400);
  TH2* dummy=new TH2F("dummy","dummy",20,0.0,30.0,30,-0.35,0.35);
  dummy->Draw();
  dummy->GetXaxis()->SetTitle("q_{T} [GeV/c]");
  dummy->GetYaxis()->SetTitleOffset(1.0);
  dummy->GetYaxis()->SetTitle("Data/Theory");

  double xv[11*2],yv[11*2];

  for (int i=1; i<=9; i++) {
    xv[i]=histsDelta[0]->GetXaxis()->GetBinCenter(i);
    yv[i]=histsDelta[0]->GetBinError(i);
  }
  xv[0]=0; yv[0]=yv[1];
  xv[10]=30; yv[10]=yv[9];
  for (int i=1; i<=9; i++) {
    xv[21-i]=histsDelta[0]->GetXaxis()->GetBinCenter(i);
    yv[21-i]=-histsDelta[0]->GetBinError(i);
  }
  xv[11]=30; yv[11]=yv[12];
  xv[21]=0; yv[21]=yv[20];

  TGraph* deltaBand=new TGraph(11*2,xv,yv);
  deltaBand->SetFillColor(kYellow);
  deltaBand->Draw("F SAME");


  for (int i=0; i<nhists; i++) {
    histsDelta[i]->Draw("SAME HIST");
  }

  tl2->SetNColumns(2);

  tl2->Draw();

  zrap_Prelim(0.4,0.25,0.4,0.17);
  zrap_Lumi(0.40,0.21,36);

  char outx[1024];
  if (outgoing!=0) {
    sprintf(outx,"%s.png",outgoing);
    //    c1->Print(outx);
    sprintf(outx,"%s.eps",outgoing);
    c1->Print(outx);
    sprintf(outx,"%s_delta.png",outgoing);
    //c2->Print(outx);
    sprintf(outx,"%s_delta.eps",outgoing);
    c2->Print(outx);
  }

}
