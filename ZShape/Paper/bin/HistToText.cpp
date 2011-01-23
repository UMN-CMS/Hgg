#include <iostream>
#include <map>
#include <iomanip>
#include <math.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <string>
#include <stdio.h>
#include <TText.h>
#include <fstream>
#include <TProfile.h>
#include <TObject.h>
#include <TKey.h>
#include <TPad.h>
#include <TLegend.h>
#include <TClass.h>
#include <TBox.h>
#include <TPave.h>
#include <TFrame.h>
#include <TWbox.h>
#include <TBox.h>
#include <TLegendEntry.h>
#include <TAttCanvas.h>


int main (int argc, char **argv) {
  if(argc<4  || argc > 4){
    std::cout << "\nlacking arguments. Usage:   HistToText <path/file.root> <histname> <textoutput.txt>\n" << std::endl;
    std::cout << "\n\n\n";
    std::cout << "1. If a TKey found in the root file is a TH1F then it just get the TH1F from the file <histname>. Then the <textoutput.txt> will have the proper format directly from the hist. " << std::endl;

    std::cout << "\n2. If a TKey found in the file is a TCanvas the it gets the TCanvas and takes the histogram <histname> FROM the TCanvas and output the <textoutput.txt> in the proper format." << std::endl;

    std::cout << "\n3. If a TKey found in the file is a TDirectoryFile (a TDirectory) then you can just take directly the histogram <histname> from the file and then gives the proper output <textoutput.txt>." << std::endl;

    std::cout << "1. Allows simple histogram in the main part of the rootfile." << std::endl;

    std::cout <<"\n  HistToText MyRoot.root MyHist MyOutput.txt " << std::endl;

    std::cout << "\n2. Allows taking the histogram from a TCanvas (this time MyHist is a member of a TCanvas)" << std::endl;

    std::cout <<"\n HistToText MyRoot.root MyHist MyOutput.txt " << std::endl;
 
    std::cout << "\n3. Is a generalization of 1 where you allow other class types and get the histogram from a directory. " << std::endl;

    std::cout << "\nHistToText MyRoot.root Mydir1/mydir2/MyHist MyOutput.txt " << std::endl;

    return -1;
  }

   
  char *filename = argv[1];
  char *histname = argv[2];
  char *outtext = argv[3];
  
  TFile *miofile = new TFile(filename,"read");
  
  TIter next(miofile->GetListOfKeys());
  TKey * key;
  TH1 *myh=0;
  std::cout << " did I get here ? " << std::endl;
    while((key=((TKey*) next()) )){
      std::cout << " is is a " << key->GetClassName() << std::endl;
    if (strcmp(key->GetClassName(),"TCanvas")==0)
      {
	  std::cout << " it is a canvas " << std::endl;
	  TCanvas *c1 = (TCanvas*)miofile->Get(key->GetName());
	  myh = (TH1*) c1->GetPrimitive(histname);
	  break;
      }
    if (strcmp(key->GetClassName(),"TH1F")==0)
      {
	  std::cout << " it is a TH1F " << std::endl;
	  myh = (TH1*) miofile->Get(histname);
	  break;
      }
    if (strcmp(key->GetClassName(),"TDirectoryFile")==0)
      {
	  std::cout << " it is a TDirectory " << std::endl;
	  myh = (TH1*) miofile->Get(histname);
	  break;
      }
    }
  std::cout << " histname " << histname << std::endl;
  //std::cout << " hist size " << myh->GetEntries() << std::endl;
    //Now the fun stuff... open the text file and output the histogram stuff.
  ofstream theOutFile;
  theOutFile.open (outtext,std::ios::out);
  theOutFile << "# \n";
  theOutFile << "# " <<  std::setw(10) <<"Bin " << std::setw(10) << " Y_Min " <<  std::setw(10) <<" Y_Max " <<  std::setw(10) <<" Value " <<  std::setw(10) << " Uncertainty" << "\n";
  for (int k=1; k<=myh->GetNbinsX(); k++) {
    double xlow = myh->GetXaxis()->GetBinLowEdge(k);
    double xhi = myh->GetXaxis()->GetBinUpEdge(k);
    double val = myh->GetBinContent(k);
    double err = myh->GetBinError(k);

    theOutFile << std::setw(10) << k << " " << std::setw(10) << xlow << " " << std::setw(10)<< xhi << " " << std::setw(10)<< val << " " << std::setw(10)<< err << "\n";
    //h->SetBinEntries(k,1);
  }

  
  theOutFile.close();




}
