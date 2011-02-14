void extractTTbar(TFile* f0) {

  double normi=36.0*167/1164732; 

  TH1* hqt=f0->Get("ZFromData/ECAL80-ECAL95-MUO/C08-m(60,120)/Z0_Pt_masscut")->Clone("hqt");

  hqt->Sumw2();
  hqt->Scale(normi);

  FILE* record;
  record=fopen("background_ttbar_qt.csv","wt");
  fprintf(record,"#    ptBin     ptMin     ptMax   Value     Error\n");

  for (int i=1; i<=18; i++) {
    fprintf(record,"%8d %9.1f %9.1f",i,hqt->GetXaxis()->GetBinLowEdge(i), hqt->GetXaxis()->GetBinUpEdge(i));
    // total error includes 18% xsec and lumi uncertainty
    double terr=sqrt(pow(hqt->GetBinContent(i)*0.18,2)+pow(hqt->GetBinError(i),2));
    fprintf(record,"%9.2f %9.2f\n",hqt->GetBinContent(i),terr);
  }

  fclose(record);
  hqt->Draw();

  TH1* hye=f0->Get("ZFromData/ECAL80-ECAL95/C08-m(60,120)/Z0_Y_masscut")->Clone("hye");

  hye->Sumw2();
  hye->Scale(normi);

  FILE* record;
  record=fopen("background_ttbar_y_ee.csv","wt");
  fprintf(record,"#    ptBin     ptMin     ptMax   Value     Error\n");

  for (int i=1; i<=100; i++) {
    fprintf(record,"%8d %9.1f %9.1f",i,hye->GetXaxis()->GetBinLowEdge(i), hye->GetXaxis()->GetBinUpEdge(i));
    // total error includes 11% luminosity uncertainty
    double terr=sqrt(pow(hye->GetBinContent(i)*0.11,2)+pow(hye->GetBinError(i),2));
    fprintf(record,"%9.2f %9.2f\n",hye->GetBinContent(i),terr);
  }

  fclose(record);
  hye->Draw();

  TH1* hyf=f0->Get("ZFromData/ECAL80-HF/C08-m(60,120)/Z0_Y_masscut")->Clone("hyf");

  hyf->Sumw2();
  hyf->Scale(normi);

  FILE* record;
  record=fopen("background_ttbar_y_ef.csv","wt");
  fprintf(record,"#    ptBin     ptMin     ptMax   Value     Error\n");

  for (int i=1; i<=100; i++) {
    fprintf(record,"%8d %9.1f %9.1f",i,hyf->GetXaxis()->GetBinLowEdge(i), hyf->GetXaxis()->GetBinUpEdge(i));
    // total error includes 11% luminosity uncertainty
    double terr=sqrt(pow(hyf->GetBinContent(i)*0.11,2)+pow(hyf->GetBinError(i),2));
    fprintf(record,"%9.2f %9.2f\n",hyf->GetBinContent(i),terr);
  }

  fclose(record);
  hyf->Draw();

  



}
