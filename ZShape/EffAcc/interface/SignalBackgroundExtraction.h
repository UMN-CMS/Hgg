#ifndef SIGBKGEXTNEW_INC
#define SIGBKGEXTNEW_INC

#include <string.h>
#include "ZShape/Base/interface/EffTableLoader.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TKey.h>
#include <map>

typedef struct{
  Double_t val;
  Double_t max;
  Double_t min; 
  bool con;
} Vars_t;

typedef std::map<std::string, Vars_t > varstore;
typedef varstore::iterator varit;
typedef std::map<int, varstore > varmap;
typedef varmap::iterator varmapit;



class SignalBackgroundExtraction {

 public:
 

  SignalBackgroundExtraction();
  //SignalBackgroundExtraction(TFile * file);
  SignalBackgroundExtraction(TFile *file, std::string EffName ="", std::string EffType = "", bool doSig  = true, bool doBkg = true, bool doSimpleOut = true, double masslo = 70.0, double masshi = 110.0);
  //SignalBackgroundExtraction(TFile *file, std::string EffName ="");
  SignalBackgroundExtraction(const std::string & sigtextFileName = "0", const std::string & bkgtextFileName = "0" );

  EffTableLoader* indexersig() const { return sigefftable_; }
  EffTableLoader* indexerbkg() const { return bkgefftable_; }


  //-----------------------------------------------------------------------------------//
  void setRootFile(TFile * file);  

  //-----------------------------------------------------------------------------------//
  void vecToHist(void);

  //-----------------------------------------------------------------------------------//
  void histToVec(void);

  //-----------------------------------------------------------------------------------//
  void produceTxtFile(std::string &textFileName, std::string sbtype);
  
  //-----------------------------------------------------------------------------------//
  void produceHistograms(std::string rootFileName);
  
  //-----------------------------------------------------------------------------------//
  void produceHistograms(TFile * rootFile, bool issig);

  //-----------------------------------------------------------------------------------//
  void setVarMinMax(double vmin, double vmax);

  //-----------------------------------------------------------------------------------//
  void produceConfig(std::string myconfig, std::string mytemplate,std::string sbtype);

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

  EffTableLoader* getSigEfftable(){return sigefftable_;};
  EffTableLoader* getBkgEfftable(){return bkgefftable_;};


 private:
  TFile * theFile_;
  std::string sigtextFileName_, bkgtextFileName_, textFileName_ ;
  
  TH1F * values1DHisto_;     TH1F * systPlus1DHisto_;
  TH1F * systMinus1DHisto_;  TH1F * denominator1DHisto_;
  
  TH2F * values2DHisto_;     TH2F * systPlus2DHisto_;
  TH2F * systMinus2DHisto_;  TH2F * denominator2DHisto_;

  std::vector<float> bin1Min_;
  std::vector<float> bin1Max_;
  std::vector<float> bin2Min_;
  std::vector<float> bin2Max_;
  std::vector<float> values1D_;
  std::vector<float> systPlus1D_;
  std::vector<float> systMinus1D_;
  std::vector<float> denominator1D_;
  std::vector<std::string> signameVec_;
  std::vector<std::string> bkgnameVec_;
  std::vector<int> ptbin_;
  std::vector<int> etabin_;
  std::vector<int> rapbin_;

  double masslo_, masshi_;

  int dimension_;
  std::string efficiencyName_;
  std::string efficiencyBinsFileName_;

  std::string effName_;
  std::string effType_;
  double vmax_,vmin_;
  Double_t eff_, effSystLo_, effSystHi_;
  bool isPt_,isBoth_,doSig_,doBkg_,isRap_,doSimpleOut_;
  void produceTxtFile1D(std::string sbtype);
  void produce1DHistograms(TFile * rootFile, bool issig);
  //  void produce2DHistograms(TFile * rootFile);
  void produceTxtFile2D();

  EffTableLoader *sigefftable_, *bkgefftable_;
 
  void doInit(void);



  varmap SigVars_;
  varmap BkgVars_;
    

};


#endif
