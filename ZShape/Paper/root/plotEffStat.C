#include "tdrstyle.C"
#include "zrapidityStandard.C"
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


  TDirectory* ecec=defcont->Get("ECAL80-ECAL95-MUO");
  TDirectory* echf=defcont->Get("ECAL80-HF");
  TDirectory* ecnt=defcont->Get("ECAL80-NTTight");


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
    //aveEntries+=h_ecec->GetEntries();

  }

  aveEntries/=ntrials;
  printf("%f\n",aveEntries);

  TH1* fe=0,*s=0,*ss=0;

  for (int i=1; i<=ntrials; i++) {
    sprintf(dname,"Trial%d",i);
    TH1* h_ecec=(TH1*)((TDirectory*)ecec->Get(dname))->Get(var);
    TH1* h_echf=(TH1*)((TDirectory*)echf->Get(dname))->Get(var);
    TH1* h_ecnt=(TH1*)((TDirectory*)ecnt->Get(dname))->Get(var);

    if (i==1) {
      //fe=new TH1F(ctitle.Data(),"FractionalError",h_ecec->GetNbinsX(),h_ecec->GetXaxis()->GetXmin(),h_ecec->GetXaxis()->GetXmax());
      fe=new TH1F("MYH","FractionalError",h_ecec->GetNbinsX(),h_ecec->GetXaxis()->GetXmin(),h_ecec->GetXaxis()->GetXmax());

      fe->GetXaxis()->SetTitle(h_ecec->GetXaxis()->GetTitle());
      if (of!=0)       fe->SetDirectory(of);
      s=new TH1F("Sum","Sum",h_ecec->GetNbinsX(),h_ecec->GetXaxis()->GetXmin(),h_ecec->GetXaxis()->GetXmax());      
      ss=new TH1F("SumSq","SumSq",h_ecec->GetNbinsX(),h_ecec->GetXaxis()->GetXmin(),h_ecec->GetXaxis()->GetXmax());
    }

    TH1* work=h_ecec->Clone();
    work->Add(h_echf);
    work->Add(h_ecnt);

    //work->Scale(aveEntries/work->GetEntries());
    work->Scale(aveEntries/h_ecec->GetEntries());

    
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


void plotEffSyst(TFile* f, TFile *fhf,  TFile* of=0, const char* var="Z0_Pt_masscut", const char* varpm="Plus") {
  setTDRStyle();

  TDirectory* based=(TDirectory*)f->Get("mcEff");
  
  TH1* base=((TDirectory*)based->Get("All"))->Get(var);
  
  std::string mys = "C07-HLT-GSF";

  TList* contents=based->GetListOfKeys();
  TDirectory* defcont=0;
  TString ctitle;
  for (int j=0; j<contents->GetSize(); j++) {
    TString q(contents->At(j)->GetName());
    if (q.Index(varpm)!=-1 && q.Index("ECAL95")!=-1) {
      defcont=based->Get(q.Data());
      printf("%s\n",q.Data());
      ctitle=q;
      ctitle.Remove(0,9);
    }
  }

  TDirectory* basedhf=(TDirectory*)fhf->Get("mcEff");
  
  TH1* basehf=((TDirectory*)basedhf->Get("All"))->Get(var);

  
  TList* contentshf=basedhf->GetListOfKeys();
  TDirectory* defconthf=0;
  TString ctitlehf;
  for (int j=0; j<contentshf->GetSize(); j++) {
    TString q(contentshf->At(j)->GetName());
    if (q.Index(varpm)!=-1 && q.Index("HF")!=-1) {
      defconthf=basedhf->Get(q.Data());
      printf("%s\n",q.Data());
      ctitlehf=q;
      ctitlehf.Remove(0,9);
      if (q.Index("HLT")!=-1) {
	 mys = "C07-HLT-GSF_Sys";
	 if (q.Index("Plus")!=-1) mys += varpm;
	 else mys += "Min";
      }
    }
  }

  TCanvas* c1=new TCanvas("c1","c1",700,600);
  TLegend* tl=new TLegend(0.5,0.8,0.95,0.95);


  TDirectory* ecec=based->Get("ECAL80-ECAL95-MUO/C07-HLT-GSF");
  TDirectory* echf=based->Get("ECAL80-HF/C07-HLT-GSF");
  TDirectory* ecnt=based->Get("ECAL80-NTTight/C07-HLT-GSF");

  TH1* h_ececb=zpt_rebinForPlot((TH1*)ecec->Get(var));
  TH1* h_echfb=zpt_rebinForPlot((TH1*)echf->Get(var)); 
  std::cout << " my string is " << mys << " " << Form("%s",mys.c_str()) << std::endl;
  TH1* h_ecec= defcont   ? ((TDirectory*)defcont->Get(Form("%s",mys.c_str())))->Get(var) : h_ececb ;
  TH1* h_echf =defconthf ?((TDirectory*)defconthf->Get(Form("%s",mys.c_str())))->Get(var): h_echfb ;
  
  TH1* feb =  (TH1*) h_ececb->Clone("feb");
  TH1* fec =  (TH1*) h_ecec->Clone("fec");
  
  //feb->Add(h_echfb);
  //fec->Add(h_echf);

  //Scale the fec histogram since this is a SHAPE measurement... who cares about the absolute effects...
  fec->Scale((double (feb->Integral()))/ ( double (fec->Integral())));
  
  fe = (TH1*) feb->Clone("MYH");
  //fe->SetTitle("FractionalError");
  
  fe->Add(fec,-1.);
  fe->Divide(base);
  
  fe->Draw();
  
  TH1 *fea = fe->Clone("MYHA");
  for (int k=1; k<=fea->GetNbinsX(); k++) {
      fea->SetBinContent(k,fabs(fea->GetBinContent(k)));
    }
  
  if (of!=0) { of->cd(); fe->Write(); fea->Write();}
  

}


void makefinal(TFile* fp, TFile *fm,  TFile* of=0) {
  setTDRStyle();
  
  TH1 *hp = (TH1*) fp->Get("MYHA");
  TH1 *hm = (TH1*) fm->Get("MYHA");
  
  TH1 *h = hp->Clone("MYNA");
  for (int k=1; k<=h->GetNbinsX(); k++) {
    std::cout << hp->GetBinContent(k) <<  " " << hm->GetBinContent(k) << std::endl;
    h->SetBinContent(k,hp->GetBinContent(k)*.5 +hm->GetBinContent(k)*.5);
    std::cout << " NEW " << h->GetBinContent(k) << std::endl;
    //h->SetBinEntries(k,1);
  }
  
  if (of!=0) { of->cd(); h->Write();}

  
}

