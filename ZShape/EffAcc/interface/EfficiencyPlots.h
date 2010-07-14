#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
// Used for 2D efficiency reading/writing 
#include "ZShape/Base/interface/EffTableLoader.h"
// This is my personal class I have used to write the efficiency histograms
// This is essentionally a class that acts like a ROOT Macro

typedef std::map<int, std::string> effstr;
typedef effstr::iterator iteffstr;
typedef std::vector<std::string> vecstring;
typedef vecstring::iterator itvecstring;
typedef std::map<int, int> colors;
typedef colors::iterator itcolors;





class EfficiencyPlots {

 public:
  EfficiencyPlots(char *txtFile, char *cutType, int binStart, char *binsFile, int effDom);
  ~EfficiencyPlots();
  void makeHistos(void);
  void printIndividualHistos(const char *ftype, bool withcolor = false);
  void printSumHistos(const char *ftype);

 private:
 
  char *txtFile_; 
  char *cutType_;

  effstr Effs_;
  colors Cols_;
  colors Marker_;

  EffTableLoader *efftable_, *binsTable_; // This holds the efficiency bins information
  int numbins_, binStart_;
  char *binsFile_;
  int effDom_;
  int vareff_;

  char *time_;

  std::string physTypeDom_,physTypeVar_,physUnitsVar_, physUnitsDom_;
  
  std::vector<TGraphAsymmErrors* > effgraphs_;

  std::string trueeff;

  TStyle *mystyle_;
  
};

