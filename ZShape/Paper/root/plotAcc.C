#include "tdrstyle.C"

int colorFor(int i) {
  switch (i) {
  case (0): return 1;
  case (1): return 2;
  case (2): return 3;
  case (3): return 4;
  case (4): return 6;
  case (5): return 7;
  case (6): return 8;
  default: return i+1;
  };

}

void plotAcc(TFile* f, const char* var="Z0_Y") {
  setTDRStyle();
  TDirectory* basedir=f->Get("mcEff");

  TH1* base=((TDirectory*)basedir->Get("All"))->Get(var);

  //  TDirectory* zdir=basedir->Get("Tight-ECAL-Loose-ECAL");
  // TDirectory* zdir2=basedir->Get("Tight-ECAL-HF");

  TCanvas* c1=new TCanvas("c1","c1",700,600);
  
  TLegend* tl;
  tl=new TLegend(0.35,0.15,0.70,0.35,"Pythia 6.422");

  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);

  TH1* final=0;

  int n=0;
  TH1* stage1=((TDirectory*)basedir->Get("All"))->Get("ECAL-ECAL");
  TH1* stage2=((TDirectory*)basedir->Get("All"))->Get("ECAL-HF");
	       
  TH1* work=stage1->Clone();
  work->Add(stage2);
  TH1* work2=base->Clone();

  double w2f=0;
  int w2fn=0;
  for (int iii=-3; iii<=4; iii++) {
    w2f+=work2->GetBinContent(work2->GetNbinsX()/2+iii);
    w2fn++;
  }
  w2f/=w2fn;

  work2->Scale(1.0/w2f);
  work->Divide(work,base);
  work->SetLineColor(colorFor(n));
  work->SetLineWidth(2);
  work->GetYaxis()->SetTitle("Geometrical Acceptance");
  work->GetYaxis()->SetTitleOffset(1.2);
  work->SetTitle(0);
  work->SetStats(false);
  work->SetMaximum(1.02);
  work->SetMinimum(0.0);

  if (strstr(var,"Y")!=0) {
    work->GetXaxis()->SetRangeUser(-4.5,4.5);
  } else {
    work->GetXaxis()->SetRangeUser(0.0,120.0);
  }
  work->GetXaxis()->CenterTitle(true);
  work->GetYaxis()->CenterTitle(true);
  work->Draw("HIST");

  work2->SetLineWidth(2);
  work2->SetLineStyle(2);
  work2->SetLineColor(4);
  work2->Draw("SAMEHIST");

  tl->AddEntry(work,"Acceptance");
  tl->AddEntry(work2,"Normalized d#sigma/dY");

  tl->SetFillColor(0);
  tl->Draw("SAME");

  TText* tt=new TText(0.8,0.92,"Preliminary");
  tt->SetNDC();
  tt->SetTextSize(0.03);
  tt->Draw();
  TDatime now;
  TText* tt2=new TText(0.8,0.90,now.AsString());
  tt2->SetNDC();
  tt2->SetTextSize(0.013);
  tt2->Draw();

  char stuff[1024];
  sprintf(stuff,"acc_%s.png",var);
  c1->Print(stuff);
  sprintf(stuff,"acc_%s.eps",var);
  c1->Print(stuff);

  /*
  sprintf(stuff,"effacc_%s.txt",zdef);
  FILE* stats=fopen(stuff,"wt");

  fprintf(stats," Y_min Y_max Eff*Acc\n");

  for (int i=1; i<=final->GetNbinsX(); i++) {
    fprintf(stats," %5.2f %5.2f %.5f\n",final->GetBinLowEdge(i),final->GetBinWidth(i)+final->GetBinLowEdge(i),final->GetBinContent(i)); 
  }
  fclose(stats);
  */
}
