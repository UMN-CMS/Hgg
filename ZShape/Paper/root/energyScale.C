TH1* extractDistribution(const char* name, TFile* f, TFile* effAccFile=0) {
  
  if (effAccFile==0) effAccFile=f;

  const char* h1name="mcEff/ECAL80-ECAL95/C08-m(70,110)/Z0_Y";
  const char* h2name="mcEff/ECAL80-HF/C08-m(70,110)/Z0_Y";

  TH1* base=effAccFile->Get("mcEff/All/Z0_YTL");
  TH1* n1=effAccFile->Get(h1name)->Clone("aa");
  TH1* n2=effAccFile->Get(h2name)->Clone("bb");

  n1->Add(n2);
  n1->Sumw2();
  n1->Divide(base);
  
  TH1* sig1=f->Get(h1name)->Clone(name);
  TH1* sig2=f->Get(h2name);

  sig1->Add(sig2);
  sig1->Divide(n1);

  delete n1;
  delete n2;

  return sig1;
}

#include "tdrstyle.C"
#include "zrapidityStandard.C"

void energyScale(TFile* base, TFile* vpe, TFile* vme,  TFile* vpf, TFile* vmf, const char* errfilename=0 ) {

  FILE* errfile=0;
  if (errfilename!=0) {
    errfile=fopen(errfilename,"wt");
    TDatime now;
    fprintf(errfile,"# Energy Scale Error # %s\n",now.AsString());
    fprintf(errfile,"#Bin Y_Min Y_Max Value Uncertainty\n");
  }

  gROOT->SetStyle("Plain");
  setTDRStyle();

  TH1* hb=extractDistribution("base",base,base);
  TH1* valep=extractDistribution("vpe",vpe,base);
  TH1* valem=extractDistribution("vme",vme,base);
  TH1* valfp=extractDistribution("vpf",vpf,base);
  TH1* valfm=extractDistribution("vmf",vmf,base);
  
  valep->Scale(hb->Integral(15,86)/valep->Integral(15,86));
  valem->Scale(hb->Integral(15,86)/valem->Integral(15,86));

  valep->Add(hb,-1);
  valem->Add(hb,-1);

  valep->Divide(hb);
  valem->Divide(hb);

  valfp->Scale(hb->Integral(15,86)/valfp->Integral(15,86));
  valfm->Scale(hb->Integral(15,86)/valfm->Integral(15,86));

  valfp->Add(hb,-1);
  valfm->Add(hb,-1);

  valfp->Divide(hb);
  valfm->Divide(hb);

  TH1* valm=valem->Clone("valm");
  TH1* valp=valep->Clone("valp");
  TH1* valave=valep->Clone("valave");

  for (int i=1; i<=valep->GetXaxis()->GetNbins(); i++) {
    valm->SetBinContent(i,sqrt(pow(valem.GetBinContent(i),2)+
			      pow(valfm.GetBinContent(i),2)));
    valp->SetBinContent(i,sqrt(pow(valep.GetBinContent(i),2)+
			      pow(valfp.GetBinContent(i),2)));
    valave->SetBinContent(i,0.5*(
				 valm->GetBinContent(i)+
				 valp->GetBinContent(i)
				 ));
    valave->SetBinError(i,0.5*sqrt(
			       pow(valm->GetBinError(i),2)+
			       pow(valp->GetBinError(i),2)
				 ));

    if (i>=15 && i<=86) {
      if (errfile!=0)
	fprintf(errfile,"%4d %5.2f %5.2f %5.3f %6.4f\n",i,
		valave->GetXaxis()->GetBinLowEdge(i),
		valave->GetXaxis()->GetBinUpEdge(i),
		valave->GetBinContent(i),
		valave->GetBinError(i)
		);
    }

  }

  TCanvas* c1=new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.15);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBottomMargin(0.10);

  //  valm->SetLineColor(kBlue);
  //  valp->SetLineColor(kRed);

  //  valm->SetLineWidth(2);
  //  valp->SetLineWidth(2);


  valave->SetLineWidth(2);
  valave->SetMaximum(0.05);
  valave->SetMinimum(0.0);
  
  valave->SetTitle(0);
  valave->SetStats(0);
  valave->GetXaxis()->SetRangeUser(-3.6,3.6);
  valave->GetXaxis()->CenterTitle(true);
  valave->GetYaxis()->SetTitle("Fractional Error from Energy Scale");
  valave->GetYaxis()->SetTitleOffset(1.8);
  valave->GetYaxis()->CenterTitle(true);
  valave->Draw("HIST");
  //  valm->Draw("HIST SAME");

  zrap_Prelim(0.6,0.9,0.6,0.87);
  c1->Print("energyScale.eps");

  if (errfile!=0) fclose(errfile);

}
