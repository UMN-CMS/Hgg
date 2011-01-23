#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

#include "ZShape/EffAcc/interface/SignalBackgroundExtraction.h"

int main( int argc, char **argv ) {

  std::cout << "\nThis program converts efficiencies from Root to txt format." << std::endl;
  char* eff ;
  char* bins;
  float min=-100;
  float max=100; 
  bool sig = true;
  bool bkg = true;

  if(argc<5 ){
    std::cout << "\nlacking arguments. Usage:    efficiencyRoot2Txt <path/efficiencyFileName.root> <path/effFileName.txt> <eff name> <efftype> Min Max doSig doBkg\n" << std::endl;
    return -1;
  }
  else if (argc == 9){
    eff = argv[3];
    bins= argv[4];
    min = atof(argv[5]);
    max = atof(argv[6]);
    sig = atof(argv[7]);
    bkg = atof(argv[8]);
  }
  else if (argc == 7 || argc == 5){
    eff = argv[3];
    bins= argv[4];
    sig = 1;
    bkg = 1;
    if (argc == 7)
      {
	min = atof(argv[5]);
	max = atof(argv[6]);	
      }
  }
  else if (argc >9){
    std::cout << "too many arguments.   Usage: efficiencyFileHisto <path/efficiencyFileName.root> <path/effFile.txt> Efficiency  EffType Min Max doSig doBkg \n" << std::endl;
    return -1;
  }
  else if (argc == 8 || argc == 6){
    std::cout << "one more argument.   Usage: efficiencyFileHisto <path/efficiencyFileName.root> <path/effFile.txt> Efficiency   EffType Min Max DoSig DoBkg\n" << std::endl;
    return -1;
  }


  std::string effRootFileName = argv[1];
  std::string effTxtFileName = argv[2];

  std::string effTxtFileNameSig = "Sig"+effTxtFileName;
  std::string effTxtFileNameBkg = "Bkg"+effTxtFileName;
  
  std::cout << "root file containing efficiency is: " << effRootFileName << std::endl;
  std::cout << "txt file for signal histo efficiency is: " << effTxtFileNameSig << std::endl;
  std::cout << "txt file for backround histo efficiency is: " << effTxtFileNameBkg << std::endl;
  std::cout << "Type of Efficicnecy (eta,pt,rap) is: " << bins << std::endl;
  std::cout << "Min and Max for other variable is: " << min << " " << max << std::endl;

  TFile * effFile = new TFile(effRootFileName.c_str());

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
  SignalBackgroundExtraction Root2TxtObject(effFile, eff,bins,sig,bkg,true,60,120);

  //std::string textFileName(effTxtFileName.c_str());
  Root2TxtObject.setVarMinMax(min,max); //I should only do this is there are actually 7 args....
  Root2TxtObject.produceTxtFile(effTxtFileNameSig,"signal");
  Root2TxtObject.produceTxtFile(effTxtFileNameBkg,"background");


}

