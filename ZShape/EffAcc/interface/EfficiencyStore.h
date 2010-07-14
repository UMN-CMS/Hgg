#ifndef EFFSTORENEW_INC
#define EFFSTORENEW_INC

#include <string.h>
#include "ZShape/Base/interface/EffTableLoader.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TKey.h>

class EfficiencyStore {

 public:


  EfficiencyStore();
  //EfficiencyStore(TFile * file);
  EfficiencyStore(TFile *file, std::string EffName ="", std::string EffBinsFile ="" );
  EfficiencyStore(const std::string & textFileName);

  EffTableLoader* indexer() const { return efftable_; }

  //-----------------------------------------------------------------------------------//
  void setRootFile(TFile * file);  

  //-----------------------------------------------------------------------------------//
  void produceTxtFile(std::string &textFileName);
  
  //-----------------------------------------------------------------------------------//
  void produceHistograms(std::string rootFileName);
  
  //-----------------------------------------------------------------------------------//
  void produceHistograms(TFile * rootFile);


  
  TH1F* getValuesHisto1D()     {
    return values1DHisto_;
  };
  TH2F* getValuesHisto2D()     {return values2DHisto_;};
 
  TH1F* getSystPlusHisto1D()   {return systPlus1DHisto_;};
  TH2F* getSystPlusHisto2D()   {return systPlus2DHisto_;};

  TH1F* getSystMinusHisto1D()  {return systMinus1DHisto_;};
  TH2F* getSystMinusHisto2D()  {return systMinus2DHisto_;};
 
  TH1F* getDenominatorHisto1D(){return denominator1DHisto_;};
  TH2F* getDenominatorHisto2D(){return denominator2DHisto_;};

  int   getDimension(){return dimension_;};

  EffTableLoader* getEfftable(){return efftable_;};


 private:
  TFile * theFile_;
  std::string textFileName_;
  
  TH1F * values1DHisto_;     TH1F * systPlus1DHisto_;
  TH1F * systMinus1DHisto_;  TH1F * denominator1DHisto_;
  
  TH2F * values2DHisto_;     TH2F * systPlus2DHisto_;
  TH2F * systMinus2DHisto_;  TH2F * denominator2DHisto_;

  int dimension_;
  std::string efficiencyName_;
  std::string efficiencyBinsFileName_;

  std::string effName_;
  void produceTxtFile1D();
  void produce1DHistograms(TFile * rootFile);
  //  void produce2DHistograms(TFile * rootFile);
  void produceTxtFile2D();

  EffTableLoader* efftable_;
 
  void doInit(void);
    

};


#endif
