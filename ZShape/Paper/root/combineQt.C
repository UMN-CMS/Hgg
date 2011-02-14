#include "readStandardFile.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"

void combineQt(TFile* theory, int which) {
  TH1* truth=(TH1*)(theory->Get("mcEff/MUON-MUON/Acceptance/Z0_PtTL_masscut")->Clone("truth"));
  truth->SetDirectory(0);  



  setTDRStyle();
  TH1* elec_stat,*elec_syst;

  if (which==1) {
    elec_stat=readStandardFilePt("elec_stat","dsdqt_zee_final_avemig.csv",0);
    elec_syst=readStandardFilePt("elec_syst","dsdqt_zee_final_avemig.csv",1);
  } else {
    elec_stat=readStandardFilePt("elec_stat","dsdqt_zee_final_matrix.csv",0);
    elec_syst=readStandardFilePt("elec_syst","dsdqt_zee_final_matrix.csv",1);
  }
  TH1* mu_all=readStandardFilePt("mu_syst","dsdqt_muons.txt",0);

  TH1* elec_all=elec_syst->Clone("mu_all");
  //  TH1* mu_all=mu_syst->Clone("mu_all");
  TH1* comb_all=elec_syst->Clone("comb_all");

  TH1* pull_elec=new TH1F("pull_elec","pull_elec",20,-3,3);
  TH1* pull_mu=new TH1F("pull_mu","pull_mu",20,-3,3);

  TH1* pull_elec_qt=elec_stat->Clone("pull_elec_qt");
  TH1* pull_mu_qt=mu_all->Clone("pull_mu_qt");


  double td=0,tt=0;
  for (int i=1; i<=18; i++) {
    double total_err_elec=sqrt(pow(elec_stat->GetBinError(i),2)+pow(elec_syst->GetBinError(i),2));
    elec_all->SetBinError(i,total_err_elec);
    double bestEstimate,sumerr;
    double pelec=0,pmu=0;

    double total_err_mu=mu_all->GetBinError(i);
    sumerr=1.0/(total_err_elec*total_err_elec)+1.0/(total_err_mu*total_err_mu);
    
    bestEstimate=(elec_all->GetBinContent(i)/pow(total_err_elec,2)+mu_all->GetBinContent(i)/pow(total_err_mu,2))/sumerr;

    sumerr=sqrt(1.0/sumerr);
    
    pelec=(elec_all->GetBinContent(i)-bestEstimate)/total_err_elec;
    pull_elec->Fill(pelec);

    pull_elec_qt->SetBinContent(i,(elec_all->GetBinContent(i)/bestEstimate)-1);
    pull_elec_qt->SetBinError(i,total_err_elec/bestEstimate);
    
    pmu=((mu_all->GetBinContent(i)-bestEstimate)/total_err_mu);   
    pull_mu->Fill(pmu);

    pull_mu_qt->SetBinContent(i,(mu_all->GetBinContent(i)/bestEstimate)-1);
    pull_mu_qt->SetBinError(i,(total_err_mu)/bestEstimate);


    printf("%2d %.5e %.5e %.5e %.5e %.1f %.1f\n",i,elec_stat->GetBinContent(i),mu_all->GetBinContent(i),bestEstimate,sumerr,pelec,pmu);
    comb_all->SetBinContent(i,bestEstimate);
    comb_all->SetBinError(i,sumerr);
    td+=bestEstimate;
    tt+=truth->GetBinContent(i);
  }

  truth->Scale(1.0/(tt));

  for (int i=1; i<=18; i++) {
    truth->SetBinContent(i,truth->GetBinContent(i)/truth->GetBinWidth(i));
  }

  TH1* truth2=truth->Clone("truth2");
  
  for (int i=1; i<=18; i++) {
    truth2->SetBinContent(i,truth->GetBinContent(i)/comb_all->GetBinContent(i)-1);
    truth2->SetBinError(i,comb_all->GetBinError(i));
  }

  TCanvas* c1=new TCanvas("c1","c1",800,800);

  c1->SetBottomMargin(0.16);
  comb_all=zpt_rebinForPlot(comb_all);
  comb_all->SetMaximum(0.08);
  elec_all->GetXaxis()->SetRangeUser(0.7,500);
  comb_all->GetXaxis()->SetRangeUser(0.7,500);

  elec_all->SetMarkerStyle(25);
  mu_all->SetMarkerStyle(26);
  comb_all->SetMarkerStyle(21);


  comb_all->GetXaxis()->SetTitle("q_{T} [GeV/c]");
  comb_all->GetXaxis()->CenterTitle();
  comb_all->GetXaxis()->SetTitleOffset(1.2);

  comb_all->GetYaxis()->SetTitle("1/#sigma d#sigma/dq_{T}");
  comb_all->GetYaxis()->CenterTitle();

  comb_all->Draw("E");
  truth=zpt_rebinForPlot(truth);


  truth->Draw("SAME HIST");

  zrap_Prelim(0.75,0.9,0.75,0.82);
  zrap_Lumi(0.75,0.86,36);

  TLegend* tl=new TLegend(0.20,0.25,0.75,0.40);
  tl->AddEntry(comb_all,"Z#rightarrow e^{+}e^{-} and Z#rightarrow #mu^{+}#mu^{-} combined","P");
  tl->AddEntry(truth,"POWHEG + CT10 + Z2","L");

  tl->Draw();

  c1->SetLogx();

  c1->Print("dsdqt_combined_linear.png");
  c1->Print("dsdqt_combined_linear.eps");

  
  c1->SetLogy();

  c1->Print("dsdqt_combined_log.png");
  c1->Print("dsdqt_combined_log.eps");
  

  /*
  elec_all->Draw("E");
  mu_all->Draw("E SAME");
  comb_all->Draw("E SAME");
  */

  FILE* ftable=fopen("dsdqt_combined.tex","wt");
  FILE* fcsv=fopen("dsdqt_combined.csv","wt");
  fprintf(ftable,"\\begin{tabular}{|c|c||l|l|}\n\\hline\n");
  fprintf(ftable,"$q_{T,min}$ & $q_{T,max}$ & Measurement & Uncertainty  \\\\ \\hline\n");
  TDatime now;
  fprintf(fcsv,"# combined result electron/muon as of %s\n",now.AsString());
  fprintf(fcsv,"#bin qt_min qt_max value error\n");
  for (int i=1; i<=18; i++) {
    double bl=comb_all->GetXaxis()->GetBinLowEdge(i);
    double bh=comb_all->GetXaxis()->GetBinUpEdge(i);
    if (i==1) bl=0;
    int prec=int(-log10(comb_all->GetBinError(i)))+2;
    fprintf(ftable,"%7.1f & %7.1f & %7.*f & %7.*f \\\\ \n",
	    bl,bh,
	    prec,comb_all->GetBinContent(i),
	    prec,comb_all->GetBinError(i));
    fprintf(fcsv, "%4d %5.1f %5.1f %.9f %.9f\n",
	    i,bl,bh,
	    comb_all->GetBinContent(i),
	    comb_all->GetBinError(i));
  }
  fprintf(ftable,"\\hline\n\\end{tabular}\n");
  fclose(ftable);
  fclose(fcsv);

  TCanvas* c2=new TCanvas("c2","c2",1000,600);
  truth2=zpt_rebinForPlot(truth2);

  c2->Divide(2,1);
  c2->cd(1);
  pull_elec->Draw("HIST");
  c2->cd(2);
  pull_mu->Draw("HIST");

  TCanvas* c3=new TCanvas("c3","c3",800,400);
  c3->SetLeftMargin(0.14);
  c3->SetBottomMargin(0.16);

  truth2->SetMaximum(0.9);
  truth2->SetMinimum(-0.9);
  truth2->SetMarkerStyle(20);
  truth2->GetXaxis()->SetTitle("q_{T} [GeV/c]");
  truth2->GetXaxis()->SetRangeUser(0.7,500.0);
  truth2->GetXaxis()->CenterTitle();
  c3->SetLogx();

  truth2->GetYaxis()->SetTitle("q_{T,ll}/q_{T,combined} - 1");
  truth2->GetYaxis()->CenterTitle();  
  truth2->GetXaxis()->SetTitleOffset(1.2);
  truth2->GetYaxis()->SetTitleOffset(1.0);
  
  pull_mu_qt->SetMarkerStyle(24);
  
  truth2->GetXaxis()->SetRangeUser(0.7,500);

  truth2->Draw("HIST ");
  pull_elec_qt->Draw("E1 SAME");
  pull_mu_qt->Draw("E1 SAME");

  truth2->Draw("HIST SAME");

  TLegend* tlp=new TLegend(0.25,0.32,0.55,0.18);
  tlp->AddEntry(pull_elec_qt,"Z#rightarrowe^{+}e^{-}","P");
  tlp->AddEntry(pull_mu_qt,"Z#rightarrow#mu^{+}#mu^{-}","P");
  tlp->Draw();

  zrap_Prelim(0.35,0.9,0.35,0.82);
  zrap_Lumi(0.35,0.86,36);
  c3->Print("dsdqt_combined_deltas.eps");

}
