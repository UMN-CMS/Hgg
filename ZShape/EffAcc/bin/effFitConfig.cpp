#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

#include "ZShape/EffAcc/interface/EfficiencyFit.h"

int main( int argc, char **argv ) {

  std::cout << "\nThis program converts templates to proper values." << std::endl;


  if(argc<7 ){
    std::cout << "\nlacking arguments. Usage: SigFile BkgFile SigTemplate BkgTemplate SigOut Bkg Out  \n" << std::endl;
    return -1;
  }
  else if (argc >7){
    std::cout << "too many arguments.   Usage: efficiencyFileHisto SigFile BkgFile SigTemplate BkgTemplate SigOut Bkg Out \n" << std::endl;
    return -1;
  }



  std::string SigFileName = argv[1];
  std::string BkgFileName = argv[2];

  std::string SigTemplateFileName = argv[3];
  std::string BkgTemplateFileName = argv[4];

  std::string SigOutputFileName = argv[5];
  std::string BkgOutputFileName = argv[6];


  
  EfficiencyFit Root2TxtObject(SigFileName.c_str(),BkgFileName.c_str());

  Root2TxtObject.produceConfig(SigOutputFileName.c_str(),SigTemplateFileName.c_str(),"signal");
  Root2TxtObject.produceConfig(BkgOutputFileName.c_str(),BkgTemplateFileName.c_str(),"background");


}

