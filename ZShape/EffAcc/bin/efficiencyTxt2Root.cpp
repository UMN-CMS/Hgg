#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

#include "ZShape/EffAcc/interface/EfficiencyStore.h"

int main( int argc, char **argv ) {

  std::cout << "\nThis program converts efficiencies from txt to Root format." << std::endl;

  if(argc<3 ){
    std::cout << "\nlacking arguments. Usage:    efficiencyTxt2Root <path/effFileName.txt> <path/efficiencyFileName.root>\n"
	      << std::endl;
    return -1;
  }

  else if (argc >3){
    std::cout << "too many arguments.   Usage: efficiencyFileHisto <path/effFileName.txt> <path/efficiencyFileName.root>\n"
	      << std::endl;
    return -1;
  }

  char * effTxtFileName = argv[1];
  char * effRootFileName = argv[2];
  std::cout << "txt file containing efficiency is: " << effTxtFileName << std::endl;
  std::cout << "root file for histo efficiency is: " << effRootFileName << std::endl;

  // created with text file in constructor
  std::string textInputName(effTxtFileName);

  EfficiencyStore Txt2RootObject( textInputName );
  // the method to produce txt file is called with name of product as argument
  Txt2RootObject.produceHistograms( std::string(effRootFileName) );

}

