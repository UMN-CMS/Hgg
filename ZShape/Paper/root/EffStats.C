// Laser Template
#include <iostream>
#include <map>
#include <iomanip>
#include <math.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <string>
#include <stdio.h>
#include <TText.h>
#include <fstream>
#include <TProfile.h>
#include <TObject.h>
#include "/home/haupt/MACROS/tdrstyle.C"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include <initializer_list>

using namespace std;


typedef struct{
  int color;
  int marker;
  double thickness;
  int linethick;
} myCust_t;

typedef std::map<std::string, myCust_t> RunType;
typedef RunType::iterator iRunType;

typedef std::map<std::string, RunType> RunTypeMap;
typedef RunTypeMap::iterator iRunTypeMap;


template <class T>
class myTClass {
    std::vector <T> Tvec;
    //typename std::vector<T>::iterator iTvec;
    std::string st,stlab,legtype;
    bool stroot,doIndi_;
    T baseT,summedT;
    double xmin, xmax, ymin, ymax;
    bool ratios, bayesratios, scale;
  public:
    myTClass (std::string mys)
      {st=mys; if (mys.find(".root") != string::npos) {stroot=true; std::cout << " root file " << mys << std::endl;} else{stroot=false; std::cout << "  not a root file " << mys << std::endl;}
      xmin=0; xmax=0; ymin=0; ymax=0;}
    void plot (std::string myfs, RunType &myrt, std::string ops);
    T getFiles(std::string myfs);
    void custimizePlots(void );   
    void doLegend(TLegend* myLeg, std::string legfops, std::string legops );
    void setXYLimits(double x1, double x2, double y1, double y2) {xmin=x1; xmax=x2; ymin=y1; ymax=y2;}
    void setPlotType(bool rt, bool br, bool sc) {ratios=rt; bayesratios=br; scale=sc; if (bayesratios) ratios=false;}
    void setSTlab(std::string stLab){stlab=stLab;}
    void setLegType(std::string leg) {legtype=leg;}
    void setIndiPlots(bool doindi) {doIndi_ = doindi;}
    void doIndiPlots(T myp, iRunType &ir);
    void doSummedQuad(void);
};

template <class T>
T myTClass<T>::getFiles (std::string myfs)
{
  TFile *myf = new TFile(myfs.c_str());
  T myp = (T) myf->Get(st.c_str());
  Tvec.push_back(myp);
  return myp;
}

template <class T>
void myTClass<T>::doLegend (TLegend* myLeg, std::string legfops, std::string legops)
{
   int i=0;
  for (auto iTvec = Tvec.begin(); iTvec != Tvec.end(); ++iTvec) 
  {
      if (i < 1)  myLeg->AddEntry((*iTvec),(*iTvec)->GetTitle(),legfops.c_str());
      else myLeg->AddEntry((*iTvec),(*iTvec)->GetTitle(),legops.c_str());
      i = 1;
  }
  //Simple little code to get out the stuff before the first space
  
}

template <class T>
void myTClass<T>::doSummedQuad (void)
{
   int i=0;
   summedT = (T) baseT->Clone("newSummedQuad"); 
   summedT->Reset();
   summedT->SetName("newSummedQuad");
   Int_t numBins = summedT->GetNbinsX();
   for (Int_t bin = 0 ; bin < numBins; ++bin)
   {
      double x2 = 0.0;
      for (auto iTvec = Tvec.begin(); iTvec != Tvec.end(); ++iTvec) 
      {
	double val = (*iTvec)->GetBinContent(bin);
	x2 += val*val;
      }
      double x = pow(x2,0.5);
      summedT->SetBinContent(bin,x);
   }
   
}

template <class T>
void myTClass<T>::doIndiPlots(T myp, iRunType &ir)
{
  if (!doIndi_) return;
  RooWorkspace* w = new RooWorkspace();
  RooRealVar* mass = new RooRealVar("mass","mass", 50, 140);
  RooDataHist *data = new RooDataHist("data", "data", *mass, myp); 
  w->import(*data);
  std::cout << " I was here " << std::endl;
  w->factory("Voigtian::signal(mass, mean[91.2,85,96], width[2.495], sigma[2.3,1,9])");
  w->factory("EXPR::backgroundPass('(RooMath::erfc((alpha - mass) * beta))*(( (mass-peak)*gamma < -70  )* (1e20) + ( (mass-peak)*gamma > 70 )* (0.0) + ( (mass-peak)*gamma >= -70 && (mass-peak)*gamma <= 70) * exp(-(mass-peak)*gamma))',mass,peak[91.1876],gamma[0.05, 0.0, 0.1],beta[0.051,0.0,0.2],alpha[62.0,35.0,105.0])");
  //w->factory("Voigtian::signal(mass, mean[91.2,85,96], width[2.495], sigma[7.0,1,9])");
  //w->factory("EXPR::backgroundPass('(RooMath::erfc((alpha - mass) * beta))*(( (mass-peak)*gamma < -70  )* (1e20) + ( (mass-peak)*gamma > 70 )* (0.0) + ( (mass-peak)*gamma >= -70 && (mass-peak)*gamma <= 70) * exp(-(mass-peak)*gamma))',mass,peak[91.1876],gamma[0.05, 0.0, 0.1],beta[0.051,0.0,0.2],alpha[62.0,35.0,75.0])");
  w->factory("SUM::pdfPass(numSignalAll[0.,1e10]*signal, numBackgroundPass[0.,1e10]*backgroundPass)");
  double totPassing = w->data("data")->sumEntries();
  std::cout << " total is " << totPassing << std::endl;
  w->var("numSignalAll")->setVal(0.95*totPassing);
  w->var("numBackgroundPass")->setVal(0.05*totPassing);
  bool quiet = true;
  w->pdf("pdfPass")->fitTo(*data, RooFit::Save(true), RooFit::Extended(true), RooFit::NumCPU(1), RooFit::PrintLevel(quiet?-1:1), RooFit::PrintEvalErrors(quiet?-1:1), RooFit::Warnings(!quiet));
  w->saveSnapshot("finalState",w->components());

  std::string fullabel = legtype + ((*ir).first);
  std::cout << " fulllabel " << fullabel << std::endl;

  
  TCanvas *mycc = new TCanvas("mycc","mycc",1200,700);
  RooPlot* frame = w->var("mass")->frame(RooFit::Title(fullabel.c_str())) ;
  w->data("data")->plotOn(frame);
  w->pdf("pdfPass")->plotOn(frame,RooFit::LineColor((*ir).second.color)) ;
  w->pdf("pdfPass")->plotOn(frame,RooFit::LineColor((*ir).second.color),RooFit::Components("backgroundPass"), RooFit::LineStyle(kDashed)) ;
  w->pdf("pdfPass")->paramOn(frame,RooFit::Label(fullabel.c_str()),RooFit::Format("NELU",RooFit::AutoPrecision(2)),RooFit::Layout(0.61, 0.98, 0.95));
  //w->pdf("pdfPass")->paramOn(frame, w->data("data"), (*ir).first.c_str(), 0, "NELU", 0.65, 0.98, 0.94);
  frame->GetYaxis()->SetTitleOffset(0.9);
  frame->Draw();

  //Print these out
  std::vector<std::string> mytype = {"png","C","root"};
  for (auto imytype=mytype.begin(); imytype!=mytype.end(); ++imytype)
  mycc->Print(Form("%s_%d.%s",stlab.c_str(),Tvec.size(),(*imytype).c_str()));
  delete mycc;
  delete w;
  delete mass;
  delete data;
  delete frame;
}

template <class T>
void myTClass<T>::plot (std::string myfs, RunType &myrt, std::string ops)
{
  T myp;
  if (myfs.find("SummedQuad") != string::npos ) {
    doSummedQuad();
    myp =  summedT;
  }
  else {
    TFile *myf = (stroot) ? new TFile(st.c_str()) : new TFile(myfs.c_str()) ;
    myf->cd();
    myp =  (stroot) ? ((T) myf->Get(myfs.c_str())) : ((T) myf->Get(st.c_str()));
    //T myp;
    //if (ops.find("SAME") == string::npos) myp =  (stroot) ? ((T) myf->Get(myfs.c_str())) : ((T) myf->Get("ZFromData/ECAL95-ECAL95/C07-HLT-GSF/Z0_Pt"));
    //else {myp =  (stroot) ? ((T) myf->Get(myfs.c_str())) : ((T) myf->Get(st.c_str()));}
    if (ops.find("SAME") == string::npos) {baseT = (T) myp->Clone("newbase"); baseT->Sumw2();} 
    if (ops.find("SAME") == string::npos) std::cout << " did I get here " << std::endl;
  }
  iRunType ir = myrt.begin();
  std::cout << " name is " << myp->GetName() << std::endl;
  myp->Sumw2();
  std::cout << " the title is " << (*myrt.begin()).first.c_str() << std::endl;

  //myp->GetXaxis()->SetRangeUser(4,140);
  if (myp){
          std::cout << " entries " << myp->GetEntries() << std::endl;
    if (myp->GetEntries() > 0) {
      doIndiPlots(myp,ir);
      if (scale) {
	std::cout << " Scale Factor is " << (double (baseT->Integral()))/( double ( myp->Integral())) << std::endl;
	std::cout << " base " << baseT->Integral() << " other " << myp->Integral() << std::endl;
	myp->Scale( (double (baseT->Integral()))/( double ( myp->Integral())) );
      }
      if (bayesratios)
      { 
	TGraphAsymmErrors *mypw = new TGraphAsymmErrors();
	mypw->BayesDivide(myp,baseT);
	mypw->GetYaxis()->SetRangeUser(0.0,3.00);
	//myp->GetXaxis()->SetRangeUser(110,140);
	mypw->GetXaxis()->SetTitle("Y_{Z0}");
	mypw->GetYaxis()->SetTitleOffset(0.9);
	mypw->SetTitle((*ir).first.c_str());
	mypw->SetMarkerColor((*ir).second.color);
	mypw->SetLineColor((*ir).second.color);
	mypw->SetMarkerStyle((*ir).second.marker);
	mypw->SetMarkerSize((*ir).second.thickness);
	mypw->SetLineStyle((*ir).second.thickness);
	mypw->SetLineWidth((*ir).second.linethick);
	mypw->Draw(ops.c_str());
      }
      if (ratios) myp->Divide(baseT);

      //myp->GetYaxis()->SetTitle("Fractional Error From Bin Correlations");
      //myp->GetYaxis()->SetTitle("Fractional Error From Efficiency Statistics");
      myp->GetYaxis()->SetTitle("Fractional Error ");
      myp->GetYaxis()->SetTitleSize(0.05);
      myp->GetXaxis()->CenterTitle();
      myp->GetYaxis()->CenterTitle();
      if (ymin > 0 || ymax > 0) myp->GetYaxis()->SetRangeUser(ymin,ymax);
      if (xmin > 0 || xmax > 0) myp->GetXaxis()->SetRangeUser(xmin,xmax);
      myp->GetYaxis()->SetTitleOffset(1.1);
      myp->SetTitle((*ir).first.c_str());
      myp->SetMarkerColor((*ir).second.color);
      myp->SetLineColor((*ir).second.color);
      myp->SetMarkerStyle((*ir).second.marker);
      myp->SetMarkerSize((*ir).second.thickness);
      myp->SetLineStyle((*ir).second.thickness);
      myp->SetLineWidth((*ir).second.linethick);
      if (!bayesratios) myp->Draw(ops.c_str());
      Tvec.push_back(myp);

    }
  }
}



int main () {
    
    //Map of run numbers to 'myCust_', like color, marker style
    setTDRStyle();
    
    MYSUPERCRUST
    
    RunTypeMap myRunMap;
   // myPair_t JRUN; JRUN.run=run; JRUN.myCust= RedCirc; myRunMap['myrun.root']=[]
   // myRunMap['myrun.root']['run']=RedCirc
   
    DARUNMAPCRUST
    
    std::string mylegstring = myLeg->GetHeader();
    size_t pos;
    pos = mylegstring.find(" ");
    mc.setLegType(mylegstring.substr(0,pos+1));
    ///LOOP over the map
    ///for each open the file, and plot the TObject
    //for the first one do not do 'same'
    //don't forget to get a TCanvas first.
    ///gStyle->SetPadLeftMargin(0.17);
    TCanvas *myc = new TCanvas("myc","myc",1400,800);
    
    myc->cd();
    myc->SetLogy(1);
    for (auto it = myRunMap.begin(); it != myRunMap.end(); it++)
    {
      std::string type = (*it).first;
      // rtype = (*it).second;
      std::cout <<" Working on " << type <<  std::endl;
      
      if ( it == myRunMap.begin()) mc.plot((*it).first,(*it).second,myops.c_str());
      else mc.plot((*it).first,(*it).second,my2ops.c_str());
      std::cout << " done with this one " << std::endl;
      
    }
    
    
    mc.doLegend(myLeg, legfops,legops );
    myLeg->SetFillColor(kWhite);
   // myLeg->SetTitle("EB-EE Normalized Mass Plots");
    myLeg->Draw();
    
    std::vector<std::string> mytype = {"png","C","root","eps"};
    //for (auto imytype=mytype.begin(); imytype!=mytype.end(); ++imytype) 
    for (auto imytype=mytype.begin(); imytype!=mytype.end(); ++imytype)
     myc->Print(Form("%s.%s",mystr.c_str(),(*imytype).c_str()));
    
    std::cout << " what is up here " << std::endl;
    return 0;
} //end of this silly macro



