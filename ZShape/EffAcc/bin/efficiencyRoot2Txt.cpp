#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

#include "ZShape/EffAcc/interface/EfficiencyStore.h"

int main( int argc, char **argv ) {

  std::cout << "\nThis program converts efficiencies from Root to txt format." << std::endl;
  char *eff ="";
  char *bins="";

  if(argc<5 ){
    std::cout << "\nlacking arguments. Usage:    efficiencyRoot2Txt <path/efficiencyFileName.root> <path/effFileName.txt> <eff name> <effbins>\n" << std::endl;
    return -1;
  }
  else if (argc == 5){
    eff = argv[3];
    bins= argv[4];
  }
  else if (argc >5){
    std::cout << "too many arguments.   Usage: efficiencyFileHisto <path/efficiencyFileName.root> <path/effFile.txt> Efficiency  EfficiencyBinsFile\n" << std::endl;
    return -1;
  }

  char * effRootFileName = argv[1];
  char * effTxtFileName = argv[2];
  std::cout << "root file containing efficiency is: " << effRootFileName << std::endl;
  std::cout << "txt file for histo efficiency is: " << effTxtFileName << std::endl;
  std::cout << "txt file holding bin info is: " << bins << std::endl;

  TFile * effFile = new TFile(effRootFileName);

  if (effFile){     
    std::cout << "This is the content of the root file:\n\n" << std::endl;
    effFile->ls();
  }
  else
    {
      std::cout << "in efficienncyFileHisto program, file: " << effRootFileName << " not found" << std::endl;      
      return -1;
    }
  
  
  // created with root file in constructor
  EfficiencyStore Root2TxtObject(effFile, eff, bins);

  std::string textFileName(effTxtFileName);
  
  Root2TxtObject.produceTxtFile(textFileName);

}

