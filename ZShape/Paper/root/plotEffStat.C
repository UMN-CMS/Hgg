#include "tdrstyle.C"

void plotEffStat(TFile* f, TFile* of=0, const char* var="Z0_Y") {
  setTDRStyle();

  TDirectory* based=(TDirectory*)f->Get("mcEff");
  
  TH1* base=((TDirectory*)based->Get("All"))->Get(var);

  TList* contents=based->GetListOfKeys();
  TDirectory* defcont=0;
  TString ctitle;
  for (int j=0; j<contents->GetSize(); j++) {
    TString q(contents->At(j)->GetName());
    printf("%s\n",q.Data());
    if (q.Index("EffStats")!=-1) {
      defcont=based->Get(q.Data());
      ctitle=q;
      ctitle.Remove(0,9);
    }
  }



  TCanvas* c1=new TCanvas("c1","c1",700,600);
  TLegend* tl=new TLegend(0.5,0.8,0.95,0.95);


  TDirectory* ecec=defcont->Get("Tight-ECAL-Loose-ECAL");
  TDirectory* echf=defcont->Get("Tight-ECAL-HF");

  if (ecec->GetListOfKeys()->GetSize()!=echf->GetListOfKeys()->GetSize()) {
    printf("Different numbers of Trials in ECAL/ECAL and HF/ECAL\n");
    return;
  }

  int ntrials=ecec->GetListOfKeys()->GetSize();

  // pass 1, we get the average number of entries
  double aveEntries=0;

  char dname[120];
  for (int i=0; i<ntrials; i++) {
    sprintf(dname,"Trial%d",i+1);
    TH1* h_ecec=((TDirectory*)ecec->Get(dname))->Get(var);
    TH1* h_echf=((TDirectory*)echf->Get(dname))->Get(var);

    aveEntries+=h_ecec->GetEntries()+h_echf->GetEntries();
  }

  aveEntries/=ntrials;
  printf("%f\n",aveEntries);

  TH1* fe=0,*s=0,*ss=0;

  for (int i=1; i<=ntrials; i++) {
    sprintf(dname,"Trial%d",i);
    TH1* h_ecec=(TH1*)((TDirectory*)ecec->Get(dname))->Get(var);
    TH1* h_echf=(TH1*)((TDirectory*)echf->Get(dname))->Get(var);

    if (i==1) {
      fe=new TH1F(ctitle.Data(),"FractionalError",h_ecec->GetNbinsX(),h_ecec->GetXaxis()->GetXmin(),h_ecec->GetXaxis()->GetXmax());
      fe->GetXaxis()->SetTitle(h_ecec->GetXaxis()->GetTitle());
      if (of!=0)       fe->SetDirectory(of);
      s=new TH1F("Sum","Sum",h_ecec->GetNbinsX(),h_ecec->GetXaxis()->GetXmin(),h_ecec->GetXaxis()->GetXmax());      
      ss=new TH1F("SumSq","SumSq",h_ecec->GetNbinsX(),h_ecec->GetXaxis()->GetXmin(),h_ecec->GetXaxis()->GetXmax());
    }

    TH1* work=h_ecec->Clone();
    work->Add(h_echf);
    work->Scale(aveEntries/work->GetEntries());

    work->Divide(work,base);    
    for (int k=1; k<=work->GetNbinsX(); k++) {
      s->Fill(work->GetBinCenter(k),work->GetBinContent(k));
      ss->Fill(work->GetBinCenter(k),work->GetBinContent(k)*work->GetBinContent(k));
    }
  }
  for (int k=1; k<=s->GetNbinsX(); k++) {
    double ave=s->GetBinContent(k)/ntrials;
    if (ave>0) {
      double rms=ss->GetBinContent(k)/ntrials;
      //      printf("bin %d %f %f %f \n",k,ave,rms,ave*ave);
      rms=sqrt(fabs(rms-ave*ave));
      //      fe->Fill(p->GetBinCenter(k),p->GetBinError(k)/p->GetBinContent(k));
      fe->Fill(s->GetBinCenter(k),rms/ave);
    }
  }
  fe->Draw();
  if (of!=0) fe->Write();
}
