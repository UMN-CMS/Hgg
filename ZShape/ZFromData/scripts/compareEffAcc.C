#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <THStack.h>
#include <TCanvas.h>
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <TLegend.h>
#include <TDirectory.h>
#include <TKey.h>
#include <time.h>
#include <stdio.h>
#include <TText.h>
#include <TLatex.h>
#include <fstream>
#include <iomanip>
#include <TArrayD.h>
#include <TDatime.h>
#include <TMath.h>

// This is my personal class I have used to write the efficiency times acceptance histograms
// This is essentionally a class that acts like a ROOT Macro
//------UPDATE----Jan-23rd-2009------
// Added the ability to compare the internal histograms.
//-----------------------------------

typedef std::map<int, std::string> effstr;
typedef effstr::iterator iteffstr;
typedef std::vector<std::string> vecstring;
typedef vecstring::iterator itvecstring;
typedef std::map<std::string, vecstring> mapstrings;
typedef mapstrings::iterator itmapstrings;
typedef std::map<int, int> colors;
typedef colors::iterator itcolors;
typedef std::vector<TH1F*> vechists;
typedef vechists::iterator itvechists;
typedef std::map<std::string, vechists> namedvechists;
typedef namedvechists::iterator itnamedvechists;
typedef std::vector<namedvechists> defhists;
typedef defhists::iterator itdefhists;


#define NUMEFFS 5
#define PHYSVARS 2

typedef struct{
   TH1F *eff;
   TH1F *den;
} PhysVar_t;

typedef struct{
  PhysVar_t var[PHYSVARS];
  TH2F *eff2d;
  TH2F *den2d;
} Histos_t;

typedef struct{
  vechists EffAcc;
  vechists ZRap;
  TH1F *ZRapMC;
  TH1F *ZRapMC_FD;
} McHistos_t;

typedef struct{
  vechists EffAcc;
  vechists ZRapAfter;
  TH1F *ZRapTotal;
  TH1F *ZMCRapTotal;
  TH1F *ZEffAccTotal;
  THStack *ZRapTotSk;
  THStack *ZEffAccTotSk;
  TH1F *ZResult;
} FinalHistos_t;

class EffAccHistos {

 public:
  EffAccHistos(std::string fromDataFile, std::string mcFile, vecstring ZDefVec, double scale);
  ~EffAccHistos();
  void getFiles(void);
  void getHistos(void);
  void printIndividualHistos(const char *ftype, bool withcolor = false);
  void printSumHistos(const char *ftype);
  void MakeResultFiles(const char *directory);

 private:

  std::string fromDataFile_;
  std::string mcFile_;
  vecstring ZDefVec_;
  effstr Effs_;
  colors Cols_;
  colors Mark_;
  Histos_t theHistos_[NUMEFFS];
  vechists FromDataHists_;
  McHistos_t MCHistos_;
  FinalHistos_t FinalHistos_;
  TFile *ZFDFile_;
  TFile *ZMCFile_;
  double scale_;
  const char *time_;

  defhists zfdDefHists_;
  defhists zmcDefHists_;
  
  double xlab_, ylabp_, ylabt_;
  double xlabc_, ylabpc_, ylabtc_;
  double xlabb_, ylabpb_, ylabtb_;
  double xlabr_, ylabpr_, ylabtr_;
  double xlabl_, ylabpl_, ylabtl_;

};


EffAccHistos::EffAccHistos(std::string fromDataFile, std::string mcFile, vecstring ZDefVec, double scale) :
  fromDataFile_(fromDataFile),
  mcFile_(mcFile),
  ZDefVec_(ZDefVec),
  scale_(1./scale)
{
   //File we initialize the strings for the efficiencies
   Effs_[0] = "GsfTrack";
   Effs_[1] = "Isolation";
   Effs_[2] = "ElectronId";
   Effs_[3] = "HLT";
   Effs_[4] = "HFElectronId";

  //Next we initialize the colors
   Cols_[5] = kBlack;
   Cols_[0] = kBlue;
   Cols_[1] = kMagenta;
   Cols_[2] = kRed;
   Cols_[3] = kCyan+1;
   Cols_[4] = 3;
   Cols_[6] = 40;
   Cols_[7] = 46;
   Cols_[8] = kMagenta;
   Cols_[9] = kBlue;
   Cols_[10] = kCyan+1;
   Cols_[11] = 3;
   Cols_[12] = 40;
   Cols_[13] = 46;
   Cols_[14] = kMagenta;
   Cols_[15] = kBlue;
   Cols_[16] = kCyan+1;
   Cols_[17] = 3;
   Cols_[18] = 40;
   Cols_[19] = 46;
 
   Mark_[0] = kOpenCircle;
   Mark_[1] = kOpenTriangleUp;
   Mark_[2] = kOpenCross;
   Mark_[3] = kOpenSquare;
   Mark_[4] = kOpenDiamond;
   Mark_[5] = kOpenStar;
   Mark_[6] = kFullCircle;
   Mark_[7] = kFullSquare;
   Mark_[8] = kFullTriangleUp;
   Mark_[9] = kFullStar;
   Mark_[10] = kFullTriangleDown;
   Mark_[11] = 33;
   Mark_[12] = 34;
   Mark_[13] = 31;
   Mark_[14] = 32;
   Mark_[15] = 2;
   Mark_[16] = 3;
   Mark_[17] = 5;
   Mark_[18] = 20;
   Mark_[19] = 21;


   //Cout some of the information for the record
   for (itvecstring mine = ZDefVec_.begin(); mine != ZDefVec_.end(); ++mine)
   {
      std::cout << " The Z Defs under consideration are " << (*mine) << std::endl; 
   }

   for (iteffstr it = Effs_.begin(); it != Effs_.end(); it++)
   {
      int type = (*it).first;
      // rtype = (*it).second;
      std::cout <<" Working on " << type << " " << (*it).second <<  std::endl;
   }

   //Now it is time to get the TIME
   time_t rawtime;
   struct tm * timeinfo;
   time (&rawtime);
   timeinfo = localtime (&rawtime);
   //time_ = asctime (timeinfo);
   //std::cout << " The time is " << time_ << std::endl;
   
   TDatime *mytime = new TDatime();
   time_ = mytime->AsString();
   
    std::cout << " The time is " << time_ << std::endl;
   
   xlab_=0.67;   ylabp_=0.90; ylabt_=0.87;
   xlabc_=0.56; ylabpc_=0.90; ylabtc_=0.87;
   xlabb_=0.88; ylabpb_=0.74; ylabtb_=0.71;
   xlabr_=0.88; ylabpr_=0.90; ylabtr_=0.87;
   xlabl_=0.32; ylabpl_=0.90; ylabtl_=0.87;

}

void
EffAccHistos::getFiles(void)
{

   ZFDFile_ =  new TFile(fromDataFile_.c_str());
   ZMCFile_ =  new TFile(mcFile_.c_str());

}


void
EffAccHistos::getHistos(void)
{
   std::cout << " get mc tree level tree " << std::endl;
   //MCHistos_.ZRapMC = (TH1F*)ZMCFile_->Get("mcEff/All/Z0_YTL"); //Changed to Z0_YTL rather than Z0_Y... to be the "truth"
   std::cout << " why am I here " << std::endl;
   //MCHistos_.ZRapMC = (TH1F*)ZMCFile_->Get("mcEff/AllInRange/Acceptance/Z0_YTL"); //Changed For ALL mass range  
   MCHistos_.ZRapMC = (TH1F*)ZMCFile_->Get("mcEff/AllInRange/C01-m(60,120)/Z0_YTL"); //Changed InRange Z0_YTL
   ////MCHistos_.ZRapMC = (TH1F*)ZMCFile_->Get("ZFromData/AllFirst/Z0_Y"); //Changed InRange Z0_YTL
   std::cout << " oh boy 0 " << std::endl;
   //MCHistos_.ZRapMC->Sumw2();
      std::cout << " oh boy 10 " << std::endl;
      if ( !MCHistos_.ZRapMC) std::cout << " why the ehllo not " << std::endl;
   MCHistos_.ZRapMC->GetXaxis()->CenterTitle(kTRUE);
      std::cout << " oh boy 20 " << std::endl;

   MCHistos_.ZRapMC->GetYaxis()->CenterTitle(kTRUE);
   MCHistos_.ZRapMC->GetYaxis()->SetTitleOffset(1.2);
   MCHistos_.ZRapMC->GetYaxis()->SetTitle(Form("Events/%0.2f Units of Rapidity",MCHistos_.ZRapMC->GetXaxis()->GetBinWidth(1)));

  std::cout << " got MC tree level tree " << std::endl;

   MCHistos_.ZRapMC_FD = (TH1F*)ZFDFile_->Get("ZFromData/AllFirst/Z0_Y");
   std::cout << " get data tree level tree " << std::endl;
   MCHistos_.ZRapMC_FD->Sumw2();
   MCHistos_.ZRapMC_FD->GetXaxis()->CenterTitle(kTRUE);
   MCHistos_.ZRapMC_FD->GetYaxis()->CenterTitle(kTRUE);
   MCHistos_.ZRapMC_FD->GetYaxis()->SetTitleOffset(1.2);
   MCHistos_.ZRapMC_FD->GetYaxis()->SetTitle(Form("Events/%0.2f Units of Rapidity",MCHistos_.ZRapMC_FD->GetXaxis()->GetBinWidth(1)));

  std::cout << " got data tree level tree " << std::endl;
   for (itvecstring mine = ZDefVec_.begin(); mine != ZDefVec_.end(); ++mine)
   {
      std::cout << " The Z Defs under consideration are " << (*mine) << std::endl;
      //FromDataHists_.push_back((TH1F*)ZFDFile_->Get(Form("ZFromData/%s/C07-HLT-EtaDet/Z0_Y",(*mine).c_str())));
      FromDataHists_.push_back((TH1F*)ZFDFile_->Get(Form("ZFromData/%s/C08-m(60,120)/Z0_Y",(*mine).c_str())));
      //MCHistos_.ZRap.push_back((TH1F*)ZMCFile_->Get(Form("mcEff/%s/C07-HLT-EtaDet/Z0_Y",(*mine).c_str())));
      ///MCHistos_.ZRap.push_back((TH1F*)ZMCFile_->Get(Form("ZFromData/%s/C08-m(60,120)/Z0_Y",(*mine).c_str())));
      MCHistos_.ZRap.push_back((TH1F*)ZMCFile_->Get(Form("mcEff/%s/C08-m(60,120)/Z0_Y",(*mine).c_str())));


      namedvechists ZFDhist;
      namedvechists ZMChist;
      //HERE I DO THE FILLING
      ///TDirectory *ZMCdir = ZMCFile_->GetDirectory(Form("ZFromData/%s",(*mine).c_str()));
      TDirectory *ZMCdir = ZMCFile_->GetDirectory(Form("mcEff/%s",(*mine).c_str()));
      TDirectory *ZFDdir = ZFDFile_->GetDirectory(Form("ZFromData/%s",(*mine).c_str()));
      TIter next(ZMCdir->GetListOfKeys());
      TKey * key;
      while(key=(TKey*) next() ){
          if ( strcmp(key->GetClassName(),"TDirectoryFile")==0 ) std::cout << "This is a " <<key->GetClassName() << " named " << key->GetName() << std::endl;
          else continue;
          TDirectory *ZMCdirtemp = ZMCdir->GetDirectory(key->GetName());
          //if (strstr(key->GetName(),"C08")) continue;
	  //if (strstr(key->GetName(),"C03")) continue;
          TDirectory *ZFDdirtemp = ZFDdir->GetDirectory(key->GetName());
          TIter nextd(ZMCdirtemp->GetListOfKeys());
          TKey * keyd;
          vechists singleHistsZMC;
          vechists singleHistsZFD;
          while(keyd=(TKey*) nextd() ){
             if ( strcmp(keyd->GetClassName(),"TH1F")==0 && !(strstr(keyd->GetName(),"TL")) ) std::cout << "This is a " <<keyd->GetClassName() << " named " << keyd->GetName() << std::endl;
             else continue;
             TH1F* zmctemphist = (TH1F*)ZMCdirtemp->Get(keyd->GetName());
             zmctemphist->GetXaxis()->CenterTitle(kTRUE);
             zmctemphist->GetYaxis()->CenterTitle(kTRUE);
			 const char *xtitlezmc = zmctemphist->GetXaxis()->GetTitle();
			 char *units = "Units of Rapidity";
			 if (strstr(xtitlezmc,"p_T") || strstr(xtitlezmc,"pt")) {units = "GeV";}
			 else if (strstr(xtitlezmc,"m")) {units = "GeV/c";}
			 else if (strstr(xtitlezmc,"phi")) {units = "radians";}
			 else if (strstr(xtitlezmc,"eta")) {units = "units of #eta";}
			 zmctemphist->GetYaxis()->SetTitleOffset(1.2);
			 zmctemphist->GetYaxis()->SetTitle(Form("Events/%0.2f %s",zmctemphist->GetXaxis()->GetBinWidth(1),units));
             //zmctemphist->Scale(scale_);
             singleHistsZMC.push_back(zmctemphist);
             TH1F* zfdtemphist = (TH1F*)ZFDdirtemp->Get(keyd->GetName());
             zfdtemphist->GetXaxis()->CenterTitle(kTRUE);
             zfdtemphist->GetYaxis()->CenterTitle(kTRUE);
			 const char *xtitlezfd = zmctemphist->GetXaxis()->GetTitle();
			 char *unitsfd = "Units of Rapidity";
			 if (strstr(xtitlezfd,"p_{T") ) {unitsfd = "GeV/c";}
			 else if (strstr(xtitlezfd,"m")) {unitsfd = "GeV/c^2";}
			 else if (strstr(xtitlezfd,"phi")) {unitsfd = "radians";}
			 else if (strstr(xtitlezfd,"eta")) {unitsfd = "units of #eta";}
			 zfdtemphist->GetYaxis()->SetTitleOffset(1.2);
             zfdtemphist->GetYaxis()->SetTitle(Form("Events/%0.2f %s",zfdtemphist->GetXaxis()->GetBinWidth(1),unitsfd));

             singleHistsZFD.push_back(zfdtemphist);
          }
          ZMChist[key->GetName()]=singleHistsZMC;
          ZFDhist[key->GetName()]=singleHistsZFD;
      }
      zfdDefHists_.push_back(ZFDhist);
      zmcDefHists_.push_back(ZMChist);

   }

  for (uint i =0; i < ZDefVec_.size(); ++i)
  {
     std::cout << " oh b 15 " << std::endl;
     TH1F *tempHist = (TH1F*) (MCHistos_.ZRap[i]->Clone(Form("Mytemp%d",i)));
     std::cout << " divide the first " << std::endl;
     MCHistos_.ZRap[i]->Sumw2();
     MCHistos_.ZRap[i]->GetXaxis()->CenterTitle(kTRUE);
     MCHistos_.ZRap[i]->GetYaxis()->CenterTitle(kTRUE);
     std::cout << " divide the first 0 " << std::endl;

     tempHist->Divide(MCHistos_.ZRap[i],MCHistos_.ZRapMC,1.0,1.0 );
std::cout << " divide the first 10 " << std::endl;
     tempHist->SetTitle(";Z0_Y;Efficiency X Acceptance");
     std::cout << " oh b 0 " << std::endl;
     tempHist->SetFillColor(Cols_[i]);
     tempHist->SetLineColor(Cols_[i]);
     tempHist->GetXaxis()->CenterTitle(kTRUE);
     tempHist->GetYaxis()->CenterTitle(kTRUE);
     //tempHist->Sumw2();
     MCHistos_.EffAcc.push_back(tempHist);
     std::cout << " oh b 10 " << std::endl;

  }
     std::cout << " oh b 20 " << std::endl;

   for (uint i =0; i < ZDefVec_.size(); ++i)
   {

     FromDataHists_[i]->SetLineColor(Cols_[i]);
     FromDataHists_[i]->GetXaxis()->CenterTitle(kTRUE);
     FromDataHists_[i]->GetYaxis()->CenterTitle(kTRUE);
     FromDataHists_[i]->SetFillColor(Cols_[i]);

     MCHistos_.ZRap[i]->SetLineColor(Cols_[i]);
     MCHistos_.ZRap[i]->GetXaxis()->CenterTitle(kTRUE);
     MCHistos_.ZRap[i]->GetYaxis()->CenterTitle(kTRUE);
     MCHistos_.ZRap[i]->SetFillColor(Cols_[i]);

     MCHistos_.EffAcc[i]->SetLineColor(Cols_[i]);
     MCHistos_.EffAcc[i]->GetXaxis()->CenterTitle(kTRUE);
     MCHistos_.EffAcc[i]->GetYaxis()->CenterTitle(kTRUE);
     MCHistos_.EffAcc[i]->SetFillColor(Cols_[i]);

   }

std::cout << " divide the first 50 " << std::endl;
  for (uint i =0; i < ZDefVec_.size(); ++i)
  {
     std::cout << " divide the first 60 " << std::endl;
     TH1F *tempHist = (TH1F*) (FromDataHists_[i]->Clone());
     std::cout << " divide the second " << std::endl;
     tempHist->Divide(FromDataHists_[i],MCHistos_.EffAcc[i],1.0,1.0 );
     tempHist->SetTitle(";Z0_Y;N_{selected}/(Efficiency X Acceptance)");
     tempHist->SetFillColor(Cols_[i]);
     tempHist->SetLineColor(Cols_[i]);
     tempHist->Sumw2();
     tempHist->GetXaxis()->CenterTitle(kTRUE);
     tempHist->GetYaxis()->CenterTitle(kTRUE);
     FinalHistos_.ZRapAfter.push_back(tempHist);
  }
std::cout << " divide the first 60 " << std::endl;
  FinalHistos_.ZRapTotal    = new TH1F("ZRapTotal","Z Rapidity;Y_{Z0};Events/0.1 units of rapidity",MCHistos_.ZRapMC->GetNbinsX(),-5.5,5.5);
  FinalHistos_.ZMCRapTotal  = new TH1F("ZMCRapTotal","Z Rapidity;Y_{Z0};Events/0.1 units of rapidity",MCHistos_.ZRapMC->GetNbinsX(),-5.5,5.5);
  FinalHistos_.ZEffAccTotal = new TH1F("EffAccTotal","EffxAcc Total;Y_{Z0};EffxAcc",MCHistos_.ZRapMC->GetNbinsX(),-5.5,5.5);
  FinalHistos_.ZRapTotSk = new THStack();
  FinalHistos_.ZEffAccTotSk= new THStack();
  FinalHistos_.ZRapTotal->GetXaxis()->CenterTitle(kTRUE);
  FinalHistos_.ZRapTotal->GetYaxis()->CenterTitle(kTRUE);
  FinalHistos_.ZMCRapTotal->GetXaxis()->CenterTitle(kTRUE);
  FinalHistos_.ZMCRapTotal->GetYaxis()->CenterTitle(kTRUE);
  FinalHistos_.ZEffAccTotal->GetXaxis()->CenterTitle(kTRUE);
  FinalHistos_.ZEffAccTotal->GetYaxis()->CenterTitle(kTRUE);

  for (uint i =0; i < ZDefVec_.size(); ++i)
  {
     FinalHistos_.ZRapTotSk->Add(FromDataHists_[i]);
     FinalHistos_.ZEffAccTotSk->Add(MCHistos_.EffAcc[i]);
     FinalHistos_.ZRapTotal->Add(FromDataHists_[i],1.0);
     FinalHistos_.ZMCRapTotal->Add(MCHistos_.ZRap[i],scale_); //do I need the scale factor here??
     FinalHistos_.ZEffAccTotal->Add(MCHistos_.EffAcc[i]);
  }

  FinalHistos_.ZResult = (TH1F*) MCHistos_.ZRapMC->Clone();
  FinalHistos_.ZResult->Divide(FinalHistos_.ZRapTotal,FinalHistos_.ZEffAccTotal,1.0,1.0 );
  FinalHistos_.ZResult->Sumw2();
}

void
EffAccHistos::printIndividualHistos(const char *ftype, bool withcolor)
{
  TText *plabel = new TText();
  plabel-> SetNDC();
  plabel -> SetTextFont(1);
  plabel -> SetTextColor(1);
  plabel -> SetTextSize(0.03);
  plabel -> SetTextAlign(22);
  plabel -> SetTextAngle(0);

  TText *tlabel = new TText();
  tlabel-> SetNDC();
  tlabel -> SetTextFont(1);
  tlabel -> SetTextColor(1);
  tlabel -> SetTextSize(0.02);
  tlabel -> SetTextAlign(22);
  tlabel -> SetTextAngle(0);

  TLatex *tllabel = new TLatex();
  tllabel-> SetNDC();
  tllabel -> SetTextFont(1);
  tllabel -> SetTextColor(1);
  tllabel -> SetTextSize(0.02);
  tllabel -> SetTextAlign(22);
  tllabel -> SetTextAngle(0);

  for (uint i =0; i < ZDefVec_.size(); ++i)
  {
     const char *rtype = ZDefVec_[i].c_str();
     TCanvas *tempCan = new TCanvas("tempcan","tempcan",800,600);
     FromDataHists_[i]->Draw("E1");
          plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");  //   tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
     tempCan->Print(Form("ZFromData_%s_Z0_Y.%s",rtype,ftype));

     MCHistos_.ZRap[i]->Draw();
          plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");     //tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
     tempCan->Print(Form("ZMC_%s_Z0_Y.%s",rtype,ftype));

     MCHistos_.EffAcc[i]->Draw("E1");
          plabel -> DrawText(xlab_, ylabp_, "CMS preliminary 2010");   //  tlabel -> DrawText(xlab_, ylabt_, Form("%s",time_));
     tempCan->Print(Form("ZMC_EffAcc_%s_Z0_Y.%s",rtype,ftype));

  }
  TCanvas *tempCan1 = new TCanvas("tempcan1","tempcan1",800,600);
  FinalHistos_.ZResult->SetMarkerStyle(kOpenCircle);
  FinalHistos_.ZResult->SetMarkerSize(0.9);
  FinalHistos_.ZResult->SetLineWidth(2);
  FinalHistos_.ZResult->Draw();
       plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");    // tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
  tempCan1->Print(Form("ZResult_Z0_Y.%s",ftype));

  FinalHistos_.ZEffAccTotal->Draw("hist");
       plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");    // tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
  tempCan1->Print(Form("ZEffAccTotal_Z0_Y.%s",ftype));

  FinalHistos_.ZRapTotal->Draw("pe1");
       plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");    // tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
  tempCan1->Print(Form("ZEffRapTotal_Z0_Y.%s",ftype));

  
  FinalHistos_.ZRapTotal->SetMarkerStyle(kOpenCircle);
  FinalHistos_.ZRapTotal->SetMarkerSize(0.9);
  FinalHistos_.ZRapTotal->SetLineWidth(2);
  if (FinalHistos_.ZRapTotal->GetMaximum() > FinalHistos_.ZMCRapTotal->GetMaximum())
    {
     FinalHistos_.ZRapTotal->Draw("pe");
     FinalHistos_.ZMCRapTotal->Draw("histsame");
    }
  else 
    {
     FinalHistos_.ZMCRapTotal->Draw("hist");
     FinalHistos_.ZRapTotal->Draw("samepe");
    }
  FinalHistos_.ZMCRapTotal->SetLineColor(kRed);
  
       plabel -> DrawText(xlabb_, ylabpb_, "CMS preliminary 2010");   
       TLegend *myLegj = new TLegend(.78,.80,.98,.95);
       myLegj->AddEntry(FinalHistos_.ZMCRapTotal,"Zee MC","l");
       myLegj->AddEntry(FinalHistos_.ZRapTotal,"7 TeV Data ","lp");
       myLegj->SetFillColor(kWhite);
       myLegj->Draw();
  tempCan1->Print(Form("ZEffMCRapTotal_Z0_Y.%s",ftype));
       myLegj->Delete();

  MCHistos_.ZRapMC->Draw("hist");
       plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");    // tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
  tempCan1->Print(Form("ZMCRAPFull_Z0_Y.%s",ftype));

  //MCHistos_.ZRapMC->Scale(1./8.117);
  FinalHistos_.ZResult->Draw("");
  MCHistos_.ZRapMC->Scale(scale_);
  MCHistos_.ZRapMC->SetLineColor(kRed);
  MCHistos_.ZRapMC->SetLineWidth(2);
  MCHistos_.ZRapMC->Draw("histsame");
  FinalHistos_.ZResult->Draw("same");
  
       plabel -> DrawText(xlabr_, ylabpr_, "preliminary 2010");     tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
  tempCan1->Print(Form("ZResult_With_MC_Z0_Y.%s",ftype));

  MCHistos_.ZRapMC_FD->SetLineColor(kBlue);
  MCHistos_.ZRapMC_FD->SetLineWidth(2);
  MCHistos_.ZRapMC_FD->Draw("histsame");
 
  tempCan1->Print(Form("ZResult_With_MC_With_FD_Z0_Y.%s",ftype));

  //Now I do the comparison histograms.
  for (uint i =0; i < ZDefVec_.size(); ++i)
  {
     const char *rtype = ZDefVec_[i].c_str();
     for (itnamedvechists zmcit = zmcDefHists_[i].begin(), zfdit = zfdDefHists_[i].begin(); zmcit != zmcDefHists_[i].end() || zfdit != zfdDefHists_[i].end(); ++zmcit, ++zfdit)
     {
        bool lastplot = false;
        const char *cuttype = (*zmcit).first.c_str();
        for(itvechists zmcvhit= (*zmcit).second.begin(), zfdvhit= (*zfdit).second.begin(); zmcvhit != (*zmcit).second.end() || zfdvhit != (*zfdit).second.end(); ++zmcvhit,++zfdvhit)
        {
           std::cout << " woohoo " << (*zmcvhit)->GetName() << std::endl;
	   (*zmcvhit)->Sumw2();
	   (*zmcvhit)->Scale(scale_); 
           //if (strstr((*zmcvhit)->GetName(),"Pt"))(*zmcvhit)->Rebin();
	   if ((*zfdvhit)->GetMaximum() > (*zmcvhit)->GetMaximum())(*zfdvhit)->Draw("P");
           else (*zmcvhit)->Draw("hist");
           if (strstr(cuttype,"C08")) lastplot = true;
           (*zmcvhit)->Draw("histsame");
           (*zmcvhit)->SetLineColor(kRed);
           (*zmcvhit)->SetLineWidth(2);
           (*zmcvhit)->SetFillColor(kWhite);
           //(*zmcvhit)->Sumw2();
           if (strstr((*zmcvhit)->GetName(),"P_t") || strstr((*zmcvhit)->GetName(),"Pt")) {(*zfdvhit)->GetXaxis()->SetRangeUser(0.,100.);(*zmcvhit)->GetXaxis()->SetRangeUser(0.,100.);}
           int placement = 0;
		   if (strstr((*zmcvhit)->GetName(),"Y")) {
		      if (strstr(rtype,"HF"))placement = 1;
			  else placement = 2;
		   }
		   else if (strstr((*zmcvhit)->GetName(),"m")) placement = 2;
		   else if (strstr((*zmcvhit)->GetName(),"eta")) {
		      if (strstr(rtype,"HF"))placement = 1;
			  else placement = 2;
           }
		   //(*zmcvhit)->Scale(scale_);
           (*zfdvhit)->Draw("sameeP");
           (*zfdvhit)->SetMarkerStyle(kOpenCircle);
           (*zfdvhit)->SetMarkerSize(0.9);
           (*zfdvhit)->SetMarkerColor(kBlack);
           (*zfdvhit)->SetLineColor(kBlack);
           (*zfdvhit)->SetLineWidth(2);
           if (strstr((*zmcvhit)->GetName(),"eta")) 
           {
              int electron = 0;
              if (strstr((*zmcvhit)->GetName(),"e1")) electron = 1;
              if (strstr((*zmcvhit)->GetName(),"e2")) electron = 2;
              if (electron > 0 ) (*zmcvhit)->GetXaxis()->SetTitle(Form("#eta_{d,e%d}",electron));
	      else (*zmcvhit)->GetXaxis()->SetTitle("#eta_{d,e}");
           }
           if (strstr((*zfdvhit)->GetName(),"eta")) 
           {
              int electron = 0;
              if (strstr((*zfdvhit)->GetName(),"e1")) electron = 1;
	      if (strstr((*zfdvhit)->GetName(),"e2")) electron = 2;
              if ( electron > 0 ) (*zfdvhit)->GetXaxis()->SetTitle(Form("#eta_{d,e%d}",electron));
	      else (*zfdvhit)->GetXaxis()->SetTitle("#eta_{d,e}");
           }
           if (placement ==0 )  { plabel -> DrawText(xlab_, ylabp_, "CMS  preliminary 2010");   }//  tlabel -> DrawText(xlab_, ylabt_, Form("%s",time_));}
           else  if (placement == 1 )  { plabel -> DrawText(xlabc_, ylabpc_, "CMS preliminary 2010");    }// tlabel -> DrawText(xlabc_, ylabtc_, Form("%s",time_));}
           else  if (placement == 2 )   {plabel -> DrawText(xlabb_, ylabpb_, "CMS preliminary 2010"); }//    tlabel -> DrawText(xlabb_, ylabtb_, Form("%s",time_));}
           else { plabel -> DrawText(xlab_, ylabp_, "CMS preliminary 2010");   }//  tlabel -> DrawText(xlab_, ylabt_, Form("%s",time_));}
 
            if (lastplot) {
              std::cout << " Did I get here " << std::endl;
              TLegend *myLeg1 = new TLegend(.78,.80,.98,.95);
			  if (placement == 1) {myLeg1->SetX1(0.46); myLeg1->SetX2(0.66);myLeg1->SetY1(0.70); myLeg1->SetY2(0.85);}
	       ///myLeg1->AddEntry((*zmcvhit),"Zee MC","l");
	       myLeg1->AddEntry((*zmcvhit),"Fast Sim","l");		  
               myLeg1->AddEntry((*zfdvhit),"7 TeV Data ","lp");
               myLeg1->SetFillColor(kWhite);
               myLeg1->Draw();
               tempCan1->Print(Form("Z_CompareFINAL_%s_%s.%s",rtype,(*zmcvhit)->GetName(),ftype));
               myLeg1->Delete();
            }
            TLegend *myLeg = new TLegend(.78,.78,.99,.91);
            myLeg->AddEntry((*zmcvhit),Form("ZMC %s %s",(*zmcvhit)->GetName(),cuttype),"l");
            myLeg->AddEntry((*zfdvhit),Form("ZFD %s %s",(*zfdvhit)->GetName(),cuttype),"lp");
            myLeg->SetFillColor(kWhite);
            myLeg->Draw();
            tempCan1->Print(Form("Z_Compare_%s_%s_%s.%s",rtype,cuttype,(*zmcvhit)->GetName(),ftype));
           
           
        }
     }
  }
}


void
EffAccHistos::printSumHistos(const char *ftype)
{
MCHistos_.ZRapMC->SetTitle("Z0_Y");

TText *plabel = new TText();
plabel-> SetNDC();
plabel -> SetTextFont(1);
plabel -> SetTextColor(1);
plabel -> SetTextSize(0.03);
plabel -> SetTextAlign(22);
plabel -> SetTextAngle(0);
  
TText *tlabel = new TText();
tlabel-> SetNDC();
tlabel -> SetTextFont(1);
tlabel -> SetTextColor(1);
tlabel -> SetTextSize(0.02);
tlabel -> SetTextAlign(22);
tlabel -> SetTextAngle(0);

TLatex *tllabel = new TLatex();
tllabel-> SetNDC();
//tllabel -> SetTextFont(1);
tllabel -> SetTextColor(1);
tllabel -> SetTextSize(0.03);
tllabel -> SetTextAlign(22);
tllabel -> SetTextAngle(0);

TCanvas *tempCan = new TCanvas("tempcan","tempcan",800,600);
TLegend *myLegl = new TLegend(.80,.45,.97,.75);
TLegend *myLegpl = new TLegend(.80,.65,.98,.95);
TLegend *myLegf = new TLegend(.80,.45,.97,.75);
myLegl->SetFillColor(kWhite);
myLegf->SetFillColor(kWhite);
for (uint i =0; i < ZDefVec_.size(); ++i)
  {
      const char *rtype = ZDefVec_[i].c_str();
      const char *myrtype = (strstr(rtype,"ECAL95-ECAL95-20")) ? ("WP95-WP95") :(rtype) ; 
      const char *myNrtype = (strstr(myrtype,"ECAL80-HF")) ? ("WP80-HF") :(myrtype) ; 
      myLegl->AddEntry(FromDataHists_[i],myNrtype,"l");
      myLegpl->AddEntry(FromDataHists_[i],myNrtype,"pl");
      myLegf->AddEntry(FromDataHists_[i],myNrtype,"f");
  }

FinalHistos_.ZEffAccTotSk->Draw("hist");
myLegf->Draw();
FinalHistos_.ZEffAccTotSk->GetXaxis()->SetTitle("Y_{Z0}");
     plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");   //  tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
tempCan->Print(Form("ZEffAccTotalStacked_Z0_Y.%s",ftype));


FinalHistos_.ZRapTotSk->Draw();
FinalHistos_.ZRapTotSk->GetXaxis()->SetTitle("Y_{Z0}");
myLegf->Draw();
     plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");    // tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
tempCan->Print(Form("ZRapTotalStacked_Z0_Y.%s",ftype));

if (FinalHistos_.ZRapTotSk->GetMaximum() > FinalHistos_.ZMCRapTotal->GetMaximum()) {
FinalHistos_.ZRapTotSk->Draw();
FinalHistos_.ZRapTotSk->GetXaxis()->SetTitle("Y_{Z0}");
FinalHistos_.ZMCRapTotal->Draw("samehist");
FinalHistos_.ZMCRapTotal->SetLineColor(kRed);
}
else 
{
FinalHistos_.ZMCRapTotal->Draw("hist");
FinalHistos_.ZMCRapTotal->SetLineColor(kRed);
FinalHistos_.ZRapTotSk->Draw("same");
FinalHistos_.ZRapTotSk->GetXaxis()->SetTitle("Y_{Z0}");
}
myLegf->SetTextSize(0.02);
myLegf->AddEntry(FinalHistos_.ZMCRapTotal,"Zee MC","l");
myLegf->Draw();
tempCan->Print(Form("ZRapMCTotalStacked_Z0_Y.%s",ftype));

for (uint i =0; i < ZDefVec_.size(); ++i)
  {
      FromDataHists_[i]->SetFillColor(kWhite);
      MCHistos_.EffAcc[i]->SetFillColor(kWhite);
      FromDataHists_[i]->SetMarkerStyle(Mark_[i]);
      FromDataHists_[i]->SetMarkerSize(0.9);
      FromDataHists_[i]->SetMarkerColor(Cols_[i]);
      FromDataHists_[i]->SetLineWidth(2);
  }



FinalHistos_.ZEffAccTotSk->Draw("histnostack");
myLegl->Draw();
     plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010"); //    tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
tempCan->Print(Form("ZEffAccTotalEach_Z0_Y.%s",ftype));
FinalHistos_.ZRapTotSk->Draw("e1p,nostack");
myLegl->Draw();
     plabel -> DrawText(xlabr_, ylabpr_, "CMS preliminary 2010");   //  tlabel -> DrawText(xlabr_, ylabtr_, Form("%s",time_));
tempCan->Print(Form("ZRapTotalEach_Z0_Y.%s",ftype));

FinalHistos_.ZResult->SetMarkerStyle(kOpenCircle);
FinalHistos_.ZResult->SetMarkerSize(0.9);
FinalHistos_.ZResult->SetLineWidth(2);
FinalHistos_.ZResult->Draw();
FinalHistos_.ZRapTotSk->Draw("same,nostack,e1p");
MCHistos_.ZRapMC->Scale(scale_);
MCHistos_.ZRapMC->SetLineColor(kRed);
MCHistos_.ZRapMC->SetLineWidth(2);
MCHistos_.ZRapMC->Draw("histsame");
FinalHistos_.ZResult->Draw("same");
myLegpl->AddEntry(FinalHistos_.ZResult,"Corrected Data","pl");
myLegpl->AddEntry(MCHistos_.ZRapMC,"POWHEG","pl");


myLegpl->Draw();
     plabel -> DrawText(xlabl_, ylabpl_, "CMS 7TEV preliminary 2010");     tllabel -> DrawLatex(xlabl_, ylabtl_, "35 pb^{-1}");
tempCan->Print(Form("ZRapCOOL_Z0_Y.%s",ftype));


//NOW I need to make the folded and the stat error distributions
Double_t mybins[200];
Double_t myerrs[200];
Double_t total = 0;
 
TH1F *myhc  = (TH1F* )FinalHistos_.ZResult->Clone("daresult");
TH1F *myhce = (TH1F*) FinalHistos_.ZResult->Clone("daresulte");
Int_t numBins = myhc->GetNbinsX();

TH1F *myhmc =  (TH1F*) MCHistos_.ZRapMC->Clone("daMC");
myhmc->SetLineColor(kRed);
myhmc->SetLineWidth(2);
myhmc->Scale(2./myhmc->Integral());

 for (Int_t bin = 0 ; bin < numBins/2+1 ; ++bin)
  {
    mybins[bin]=0; //simple emptying of the arrays.
    myerrs[bin]=0;
  }
  

 for (Int_t bin = 1 ; bin < numBins +1; ++bin)
   {
     //std::cout << " oh boy " << std::endl;
     Int_t mybin = bin - (numBins/2 +1) ;
     if (mybin >= 0 ) {mybin++;}
     if (mybin < 0 ) {mybin = TMath::Abs(mybin);}
     //std::cout << " oh boy 20  " << std::endl;
     
     Double_t myval = (Double_t ) myhc->GetBinContent(bin);
     Double_t mybe = myhc->GetBinError(bin);
     //std::cout << " oh boy 30 " << std::endl;
     
      if (mybe > 0 ) 
	{
	  mybins[mybin] += myval/(mybe*mybe);
	  myerrs[mybin] += 1./(mybe*mybe);
	  total += myval/2.;
	  //std::cout << " oh boy 35 " << std::endl;
	  std::cout << " mybin " << mybin << " myval " << myval << std::endl;
	  
	}  
      //std::cout << " oh boy 40 " << std::endl;
      
   }

 //here i need to scale TOTAL
  total = 12816;
  total = 14000;

  for (Int_t bin = 1 ; bin < numBins/2+1 ; ++bin)
    {
      if ( myerrs[bin] > 0 ) {
	Double_t daerr = TMath::Power(myerrs[bin],-0.5);
	
	Double_t dabin = mybins[bin]/myerrs[bin];

	std::cout << " da errorr " << daerr << " dabin " << dabin << " total " << total << std::endl;

	myhc->SetBinContent(bin + numBins/2 , dabin/total); //+ numBins/2 -1; because I don't make new hists

	myhc->SetBinError(bin + numBins/2 , daerr/total);

	myhce->SetBinContent(bin + numBins/2 , daerr/dabin);
	myhce->SetBinError(bin + numBins/2 , 0 );

      }
      else
	{
	  myhc->SetBinContent(bin, 0);
	  myhc->SetBinError(bin, 0);
	  myhce->SetBinContent(bin, 0);
	  //std::cout << " oh boy 59.2 " << std::endl;

	}
      //std::cout << " oh boy 59.5 " << std::endl;

    }


  myhc->SetMarkerStyle(kFullCircle);
  myhc->SetMarkerSize(1.1);

  //myhc->GetXaxis()->SetRange(10,35);
  myhc->GetYaxis()->SetTitle("1/#sigma (d#sigma/dY) / 0.1 Units of Rapdity");
  myhc->GetXaxis()->SetTitle("|Y_{Z0}|");
  myhc->GetXaxis()->SetRange(numBins/2+1,numBins);
  myhc->GetXaxis()->SetRangeUser(0,3.7);
  myhc->GetYaxis()->SetTitleOffset(1.3);

  myhc->Draw("pe");
  myhmc->Draw("same hist");
  myhc->Draw("same,pe");

  TLegend *legh = new TLegend(.72,.82,.98,.95);
  legh->AddEntry(myhc,"35 pb^{-1} Data ","p");
  legh->AddEntry(myhmc,"POWHEG","l");
  legh->SetFillColor(kWhite);
  legh->Draw();

  plabel -> DrawText(xlabl_, ylabpl_, "CMS preliminary 2010");     tllabel -> DrawLatex(xlabl_, ylabtl_, "36 pb^{-1}");
  tempCan->Print(Form("ZRapFOLD_Z0_Y.%s",ftype));
  
  //new plot
  myhce->SetMarkerStyle(kFullCircle);
  myhce->SetMarkerSize(1.1);
  myhce->Draw("p");
  tempCan->SetLogy(1);
  myhce->GetXaxis()->SetRangeUser(0,3.7);
  myhce->GetXaxis()->SetTitle("|Y_{Z0}|");
  myhce->GetYaxis()->SetRangeUser(0.01,1.);
  myhce->GetYaxis()->SetTitle("Stat Error Fraction / 0.1 Units of Y");
  myhce->GetYaxis()->SetTitleOffset(1.2);

  //pt = new TPaveText(0,0.7,2,1.0,"br");
  //pt->SetFillColor(19);
  //text = pt->AddText("35 pb^{-1} Data Statistical Error");
  //pt->SetFillColor(kWhite);
  //pt->Draw();

 plabel -> DrawText(xlabl_, ylabpl_, "CMS preliminary 2010");     tllabel -> DrawLatex(xlabl_, ylabtl_, "36 pb^{-1}");
  tempCan->Print(Form("ZRapFOLDSTAT_Z0_Y.%s",ftype));


//settting things back the way they were a few plots ago
MCHistos_.ZRapMC->Scale(1.0/scale_);
  tempCan->SetLogy(0);


 /*  
   for (int itype = 0; itype < PHYSVARS; ++itype)
   {
      TCanvas *tempCan = new TCanvas("

tempcan","tempcan");
      TLegend *myLeg = new TLegend(.75,.25,.97,.55);
      for (iteffstr it = Effs_.begin(); it != Effs_.end(); it++)
      {
         int type = (*it).first;
         const char *rtype = ((*it).second).c_str();
         theHistos_[type].var[itype].eff->Draw("same");
         theHistos_[type].var[itype].eff->SetLineColor(Cols_[type]);
         myLeg->AddEntry(theHistos_[type].var[itype].eff,rtype,"l");
         if (it == Effs_.begin()) theHistos_[type].var[itype].eff->SetMaximum(1.1);
         
      }
      myLeg->Draw();
      myLeg->SetFillColor(kWhite);
      tempCan->Print(Form("sbs_eff_%s.%s",PhysVec_[itype].c_str(),ftype));
   }
   for (int itype = 0; itype < PHYSVARS; ++itype)
   {
      TCanvas *tempCan = new TCanvas("tempcan","tempcan");
      TLegend *myLeg = new TLegend(.72,.72,.97,.93);
      for (iteffstr it = Effs_.begin(); it != Effs_.end(); it++)
      {
         int type = (*it).first;
         const char *rtype = ((*it).second).c_str();
         theHistos_[type].var[itype].den->Draw("same");
         theHistos_[type].var[itype].den->SetLineColor(Cols_[type]);
         myLeg->AddEntry(theHistos_[type].var[itype].den,rtype,"l");
         //if (it == Effs_.begin()) theHistos_[type].var[itype].den->SetMaximum(1.1);
      }
      myLeg->Draw();
      myLeg->SetFillColor(kWhite);
      tempCan->Print(Form("sbs_den_%s.%s",PhysVec_[itype].c_str(),ftype));
   }
*/
}

void EffAccHistos::MakeResultFiles(const char* directory)
{
  std::cout << " ok 0 " << std::endl;
  for (uint i =0; i < ZDefVec_.size(); ++i)
  {
      std::cout << " ok 10 " << std::endl;
     const char *rtype = ZDefVec_[i].c_str();
     ofstream outputfile;
	 std::cout << " ok 20 " << std::endl;
     outputfile.open (Form("%s/%s.txt",directory,rtype),std::ios::out);
     outputfile << "#put your comments after a '#'.\n";
     outputfile << "#Z def name: " << rtype <<"\n";
     outputfile << "#Bin"
	       << " " << std::setw(10)  <<  "YMin"
	       << " " << std::setw(10)  <<  "YMax"
	       << " " << std::setw(10)  << "Entries"
	       << "\n";
     
     if (outputfile.is_open())
     {
	    std::cout << " ok 30 " << std::endl;
		const TArrayD *bins = FromDataHists_[i]->GetXaxis()->GetXbins();
		//const Double_t mybins[50] = bins;
		std::cout << " ok 35 " << std::endl;
		//std::cout << "bins " << (FromDataHists_[i]->GetXaxis()->GetXbins()->GetArray())[0] << std::endl;
		//bins++;
		//std::cout << "bins+1 " << (FromDataHists_[i]->GetXaxis()->GetXbins()->GetArray())[1] << std::endl;
        for (int r=1; r<FromDataHists_[i]->GetNbinsX()+1; ++r) 
	    {  
		   std::cout << " ok 36 " << std::endl;
		   std::cout << " " << std::setw(10) << std::setprecision(4) << r << std::endl;
		   std::cout         << " " << std::setw(10) << std::setprecision(4) << FromDataHists_[i]->GetXaxis()->GetBinLowEdge(r)  << std::endl;
			std::cout  	     << " " << std::setw(10) << std::setprecision(4) <<  FromDataHists_[i]->GetXaxis()->GetBinUpEdge(r)  << std::endl;
		    //std::cout         << " " << std::setw(10) << std::setprecision(4) <<  bins->At(1) << std::endl;
			//std::cout  	     << " " << std::setw(10) << std::setprecision(4) <<   bins->At(r+1)<< std::endl;
		  	std::cout	     << " " << std::setw(10) << std::setprecision(4) <<  FromDataHists_[i]->GetBinContent(r)<< std::endl;
			std::cout	     << std::endl;
	       outputfile << r 
		             << " " << std::setw(10) << std::setprecision(4) <<  FromDataHists_[i]->GetXaxis()->GetBinLowEdge(r)
			  	     << " " << std::setw(10) << std::setprecision(4) <<  FromDataHists_[i]->GetXaxis()->GetBinUpEdge(r)
		  		     << " " << std::setw(10) << std::setprecision(4) <<  FromDataHists_[i]->GetBinContent(r)
				     << "\n";
	    }  
     }
  std::cout << " ok 40 " << std::endl;
  outputfile.close();

  }
}

EffAccHistos::~EffAccHistos()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

