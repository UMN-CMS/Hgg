#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <TLegend.h>
// This is my personal class I have used to write the efficiency histograms
// This is essentionally a class that acts like a ROOT Macro

typedef std::map<int, std::string> effstr;
typedef effstr::iterator iteffstr;
typedef std::vector<std::string> vecstring;
typedef vecstring::iterator itvecstring;
typedef std::map<int, int> colors;
typedef colors::iterator itcolors;

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



class EffHistos {

 public:
  EffHistos(vecstring FileVec, vecstring PhysVec);
  ~EffHistos();
  void getFiles(void);
  void getHistos(void);
  void printIndividualHistos(const char *ftype, bool withcolor = false);
  void printSumHistos(const char *ftype);

 private:

  vecstring FileVec_;
  vecstring PhysVec_;
  effstr Effs_;
  colors Cols_;
  TFile *theFiles_[NUMEFFS];
  Histos_t theHistos_[NUMEFFS];
};


EffHistos::EffHistos(vecstring FileVec, vecstring PhysVec) :
  FileVec_(FileVec),
  PhysVec_(PhysVec)
{
   //File we initialize the strings for the efficiencies
   Effs_[0] = "GsfTrack";
   Effs_[1] = "Isolation";
   Effs_[2] = "ElectronId";
   Effs_[3] = "HLT";
   Effs_[4] = "HFElectronId";

  //Next we initialize the colors
   Cols_[0] = kBlack;
   Cols_[1] = kRed;
   Cols_[2] = kMagenta;
   Cols_[3] = kBlue;
   Cols_[4] = kCyan+1;

   //Cout some of the information for the record
   for (itvecstring mine = FileVec_.begin(); mine != FileVec_.end(); ++mine)
   {
      std::cout << " The files are " << (*mine) << std::endl;
   }

   for (itvecstring mine = PhysVec_.begin(); mine != PhysVec_.end(); ++mine)
   {
      std::cout << " The Physics Variables are " << (*mine) << std::endl;
   }

   for (iteffstr it = Effs_.begin(); it != Effs_.end(); it++)
   {
      int type = (*it).first;
      // rtype = (*it).second;
      std::cout <<" Working on " << type << " " << (*it).second <<  std::endl;
   }
}

void
EffHistos::getFiles(void)
{
   for (iteffstr it = Effs_.begin(); it != Effs_.end(); it++)
   {
      int type = (*it).first;
      theFiles_[type]= new TFile((FileVec_[type]).c_str());
   }

}

void
EffHistos::getHistos(void)
{
   for (iteffstr it = Effs_.begin(); it != Effs_.end(); it++)
   {
      int type = (*it).first;
      for (int itype = 0; itype < PHYSVARS; ++itype)
      {
         theHistos_[type].var[itype].eff = 
                (TH1F*) theFiles_[type]->Get(Form("sbs_eff_%s",PhysVec_[itype].c_str()));
         theHistos_[type].var[itype].eff->GetXaxis()->SetTitle(PhysVec_[itype].c_str());
         theHistos_[type].var[itype].den = 
                (TH1F*) theFiles_[type]->Get(Form("sbs_den_%s",PhysVec_[itype].c_str()));
         theHistos_[type].var[itype].den->GetXaxis()->SetTitle(PhysVec_[itype].c_str());
      }
      theHistos_[type].eff2d = 
                (TH2F*) theFiles_[type]->Get(Form("sbs_eff_%s_%s",PhysVec_[1].c_str(),PhysVec_[0].c_str()));
      theHistos_[type].eff2d->GetXaxis()->SetTitle(PhysVec_[1].c_str());
      theHistos_[type].eff2d->GetYaxis()->SetTitle(PhysVec_[0].c_str());
      theHistos_[type].den2d = 
                (TH2F*) theFiles_[type]->Get(Form("sbs_den_%s_%s",PhysVec_[1].c_str(),PhysVec_[0].c_str()));
      theHistos_[type].den2d->GetXaxis()->SetTitle(PhysVec_[1].c_str());
      theHistos_[type].den2d->GetYaxis()->SetTitle(PhysVec_[0].c_str());
   }
}

void
EffHistos::printIndividualHistos(const char *ftype, bool withcolor)
{
   for (iteffstr it = Effs_.begin(); it != Effs_.end(); it++)
   {
      int type = (*it).first;
      const char *rtype = ((*it).second).c_str();
      TCanvas *tempCan = new TCanvas("tempcan","tempcan");

      for (int itype = 0; itype < PHYSVARS; ++itype)
      {
         theHistos_[type].var[itype].eff->Draw();
	 theHistos_[type].var[itype].eff->SetTitle(Form("%s %s", rtype,theHistos_[type].var[itype].eff->GetTitle() ));
         if (withcolor) theHistos_[type].var[itype].eff->SetLineColor(Cols_[type]);
         tempCan->Print(Form("sbs_eff_%s_%s.%s",PhysVec_[itype].c_str(),rtype,ftype));
         theHistos_[type].var[itype].den->Draw();
	 theHistos_[type].var[itype].den->SetTitle(Form("%s %s",rtype,theHistos_[type].var[itype].den->GetTitle() ));
         if (withcolor) theHistos_[type].var[itype].den->SetLineColor(Cols_[type]);
         tempCan->Print(Form("sbs_den_%s_%s.%s",PhysVec_[itype].c_str(),rtype,ftype));
      }
      theHistos_[type].eff2d->Draw("colz");
      theHistos_[type].eff2d->SetTitle(Form("%s %s",rtype,theHistos_[type].eff2d->GetTitle() ));
      if (withcolor) theHistos_[type].eff2d->SetLineColor(Cols_[type]);
      tempCan->Print(Form("sbs_eff_%s_%s_%s.%s",PhysVec_[1].c_str(),PhysVec_[0].c_str(),rtype,ftype));
      theHistos_[type].den2d->Draw("colz");
      theHistos_[type].den2d->SetTitle(Form("%s %s",rtype,theHistos_[type].den2d->GetTitle() ));
      if (withcolor) theHistos_[type].den2d->SetLineColor(Cols_[type]);
      tempCan->Print(Form("sbs_den_%s_%s_%s.%s",PhysVec_[1].c_str(),PhysVec_[0].c_str(),rtype,ftype));
   }

}


void
EffHistos::printSumHistos(const char *ftype)
{
   
   for (int itype = 0; itype < PHYSVARS; ++itype)
   {
      TCanvas *tempCan = new TCanvas("tempcan","tempcan");
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
}

EffHistos::~EffHistos()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

