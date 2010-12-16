#include <sstream>
#include <algorithm>
#include <list> 
#include "ZShape/EffAcc/interface/WgtProducer.h"
#include "Math/VectorUtil.h"
#include "TFile.h"
#include "TH2.h"

void WgtProducer::load(const std::string& fname){
  std::string noname="none";
  if (strcmp(fname.c_str(),noname.c_str()) != 0){
    
    TFile f(fname.c_str());
    TH2* c0= (TH2*)f.Get("gridRatio")->Clone("wgt_weights");
    
    weights_= c0;
    weights_->SetDirectory(0);
    nowgts_=false;
  }else{
    nowgts_=true;
  }
}

double WgtProducer::wgt(const double& pt, const double& y){

  if (!nowgts_){
    int bin=weights_->FindBin(fabs(y),pt);
    double rv;
    rv=weights_->GetBinContent(bin);
    return rv;
  }else{
    return 1.0;
  }
}
