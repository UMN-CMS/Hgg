#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

#include "ZShape/EffAcc/interface/EfficiencyPlots.h"

int main( int argc, char **argv ) {

  std::cout << "\nThis program converts efficiencies txt files to plots." << std::endl;
  char *eff;
  char *bins;
  char *binsfile;
  int domeff=0;

  if(argc<7 ){
    std::cout << "\nlacking arguments. Usage:    efficiencyTxt2Plots <path/efficiencyFileName.root> <path/effFileName.txt> <eff name> <effbins>\n" << std::endl;
    return -1;
  }
  else if (argc == 7){
    eff = argv[3];
    bins= argv[4];
    binsfile =  argv[5];
    domeff = atoi(argv[6]);
  }
  else if (argc >7){
    std::cout << "too many arguments.   Usage: efficiencyTxt2Plots <path/efficiencyFileName.root> <path/effFile.txt> Efficiency  EfficiencyBinsFile\n" << std::endl;
    return -1;
  }

  char * effRootFileName = argv[1];
  char * effTxtFileName = argv[2];
  std::cout << "root file containing efficiency is: " << effRootFileName << std::endl;
  std::cout << "txt file for histo efficiency is: " << effTxtFileName << std::endl;
  std::cout << "txt file holding bin info is: " << bins << std::endl;

  //TFile * effFile = new TFile(effRootFileName);

  //if (effFile){     
  //  std::cout << "This is the content of the root file:\n\n" << std::endl;
  //  effFile->ls();
  //}
  //else
   // {
    //  std::cout << "in efficienncyFileHisto program, file: " << effRootFileName << " not found" << std::endl;      
     // return -1;
   // }
  
  
  // created with root file in constructor
  int numbins = atoi(bins);
  std::cout << " about to initiate the plots" << std::endl;  
  EfficiencyPlots Txt2Plots(effTxtFileName,eff,numbins,binsfile,domeff);
  std::cout << " Just initiated the plots" << std::endl;
  std::string textFileName(effTxtFileName);
 std::cout << " about to make the histos" << std::endl;
  Txt2Plots.makeHistos();
  std::cout << " Just made the histos" << std::endl;
           

 // Txt2Plots.printIndividualHistos("png");
std::cout << " Just printed indiviuals" << std::endl;

  Txt2Plots.printSumHistos("png");
std::cout << " Just printed sumed" << std::endl;
 Txt2Plots.printSumHistos("eps");

}

