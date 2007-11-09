#ifndef EFFSTORE_INC
#define EFFSTORE_INC

#include <string.h>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TKey.h>

class EfficiencyStore {

 public:


  EfficiencyStore();
  //EfficiencyStore(TFile * file);
  EfficiencyStore(TFile *file, std::string EffName ="" , std::string PhysVar ="" );
  EfficiencyStore(std::string & textFileName);  
  

  //-----------------------------------------------------------------------------------//
  void setRootFile(TFile * file);  

  //-----------------------------------------------------------------------------------//
  void produceTxtFile(std::string &textFileName);
  
  //-----------------------------------------------------------------------------------//
  void produceHistograms(std::string rootFileName);
  
  //-----------------------------------------------------------------------------------//
  void produceHistograms(TFile * rootFile);


  
  TH1F* getValuesHisto1D()     {
    std::cout << "inside the getter method values1DHisto_: " << values1DHisto_ << std::endl;
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

 private:
  TFile * theFile_;
  std::string textFileName_;
  
  TH1F * values1DHisto_;     TH1F * systPlus1DHisto_;
  TH1F * systMinus1DHisto_;  TH1F * denominator1DHisto_;
  
  TH2F * values2DHisto_;     TH2F * systPlus2DHisto_;
  TH2F * systMinus2DHisto_;  TH2F * denominator2DHisto_;

  int dimension_;
  std::string efficiencyName_;

  std::string effName_;
  std::string physVar_;
  void produceTxtFile1D();
  void produce1DHistograms(TFile * rootFile);
  //  void produce2DHistograms(TFile * rootFile);
  void produceTxtFile2D();
    

};


#endif
