#include <fstream>
#include <iomanip>

#include <string.h>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TKey.h>

#include "ZShape/EffAcc/interface/EfficiencyStore.h"
#include "PhysicsTools/TagAndProbe/interface/EffTableLoader.h"

EfficiencyStore::EfficiencyStore()
{
  theFile_ =0;
  efftable_ =0;
  doInit();

  textFileName_="";
  
}

EfficiencyStore::EfficiencyStore(TFile * file, std::string EffName, std::string EffBinsFile)
{ 
  doInit();
  effName_ = EffName;
  efficiencyBinsFileName_ = EffBinsFile;
  efftable_ = new EffTableLoader(EffBinsFile);
  std::cout << "class EfficiencyStore created with root file: " << file->GetName() << std::endl;
  setRootFile(file);
  textFileName_="";
}


EfficiencyStore::EfficiencyStore(const std::string & textFileName)
{ 
  theFile_ =0;
  doInit();
  textFileName_ = textFileName;
  efftable_ = new EffTableLoader(textFileName_);
  std::cout << "class EfficiencyStore created with text file: " << textFileName_ << std::endl;
  produceHistograms(0);
}


void EfficiencyStore::doInit(void)
{
  values1DHisto_=0;   systPlus1DHisto_=0;
  systMinus1DHisto_=0;  denominator1DHisto_=0;
  values2DHisto_=0;  systPlus2DHisto_=0;
  systMinus2DHisto_=0; denominator2DHisto_=0;

}
  

//-----------------------------------------------------------------------------------//
  
void EfficiencyStore::setRootFile(TFile * file){
  
  theFile_ = file;
  if ( (!theFile_) &&
       theFile_->IsOpen() == kFALSE
       )
    {
      std::cout << "invalid root file. Returning." << std::endl;
    }

  std::cout<<"object EfficiencyStore opened file: " << theFile_->GetName() << std::endl;
    
    
  // getting list of keys of all objects in the file
  TIter next(theFile_->GetListOfKeys());
  TKey * key;
    
  values1DHisto_=0;     systPlus1DHisto_=0;
  systMinus1DHisto_=0;  denominator1DHisto_=0;
    
  values2DHisto_=0;     systPlus2DHisto_=0;
  systMinus2DHisto_=0;  denominator2DHisto_=0;
    
    
  // looping on objects in the file aiming at the 4 needed histograms 
  
    
  while(key=((TKey*) next()) ){
      
    //    std::cout << "object in file has name: " << key->GetName() << std::endl;
    TH1F * tmp1;    TH2F * tmp2;
    if(  strstr(key->GetName(),"eff") &&   strcmp(key->GetClassName(),"TH1F")==0 )
         { 
   	   tmp1 = (TH1F*) theFile_->Get(key->GetName());
	   values1DHisto_ = (TH1F *) tmp1->Clone();
	   values1DHisto_->SetDirectory(0);
	   //	values1DHisto_ = (TH1F*) theFile_->Get(key->GetName());
	   //	   std::cout << " This Object is values_ " <<  key->GetName() << std::endl;
         }
       else if(  strstr(key->GetName(),"systematicPlus_") &&  strcmp(key->GetClassName(),"TH1F")==0) 
         {
	   tmp1 = (TH1F*) theFile_->Get(key->GetName());
	   systPlus1DHisto_ = (TH1F *) tmp1->Clone();
	   systPlus1DHisto_->SetDirectory(0);
	   //	   std::cout << " This Object is sytemtaicPlus_ " <<  key->GetName() << std::endl;        
         }
       else if(  strstr(key->GetName(),"systematicMinus_") &&   strcmp(key->GetClassName(),"TH1F")==0) 
         {
	   tmp1 = (TH1F*) theFile_->Get(key->GetName());
	   systMinus1DHisto_ = (TH1F *) tmp1->Clone();
	   systMinus1DHisto_->SetDirectory(0);
	   //	   std::cout << " This Object is systematicMinus_ " <<  key->GetName() << std::endl;
         }
       else if(  strstr(key->GetName(),"den") &&   strcmp(key->GetClassName(),"TH1F")==0)
         {
	   tmp1 = (TH1F*) theFile_->Get(key->GetName());
	   denominator1DHisto_ = 	(TH1F *) tmp1->Clone();
	   denominator1DHisto_->SetDirectory(0);
	   //	   std::cout << " This Object is denominator_ " <<  key->GetName() << std::endl;
         }
      
    if(  strstr(key->GetName(),"eff") && strcmp(key->GetClassName(),"TH2F")==0)
         {
	   tmp2 = (TH2F*) theFile_->Get(key->GetName());
	   values2DHisto_  = (TH2F *) tmp2->Clone();
	   values2DHisto_->SetDirectory(0);
         }
    else if(  strstr(key->GetName(),"systematicPlus_") && strcmp(key->GetClassName(),"TH2F")==0)
         {
	   tmp2 = (TH2F*) theFile_->Get(key->GetName());
	   systPlus2DHisto_ =  (TH2F *) tmp2->Clone();
	   systPlus2DHisto_->SetDirectory(0);
         }
    else if(  strstr(key->GetName(),"systematicMinus_") && strcmp(key->GetClassName(),"TH2F")==0)
         {
	   tmp2 = (TH2F*) theFile_->Get(key->GetName());
	   systMinus2DHisto_ =  (TH2F *) tmp2->Clone();
	   systMinus2DHisto_->SetDirectory(0);
         }
    else if(  strstr(key->GetName(),"den") && strcmp(key->GetClassName(),"TH2F")==0)
         {
	   tmp2 = (TH2F*) theFile_->Get(key->GetName());
	   denominator2DHisto_ =  (TH2F *) tmp2->Clone();
	   denominator2DHisto_->SetDirectory(0);
	 }
  }
  std::cout << "leaving setRootFile  values1DHisto_ :" << values1DHisto_ << std::endl;
  
}
  
  

//-----------------------------------------------------------------------------------//
void EfficiencyStore::produceTxtFile(std::string &textFileName){
      
  // check if there are 4 histograms of the same type and with the same number of bins;
  bool histosAre1D =false;
  if(  ( values1DHisto_ && systPlus1DHisto_ &&  systMinus1DHisto_  &&  denominator1DHisto_ ) ) histosAre1D = true;
  if(  ( values1DHisto_ && denominator1DHisto_ ) ) histosAre1D = true; //To allow for now systematic histos
  bool histosAre2D =false;
  if(  ( values2DHisto_ && systPlus2DHisto_ &&  systMinus2DHisto_  &&  denominator2DHisto_ ) ) histosAre2D = true;
  if(  ( values2DHisto_ &&  denominator2DHisto_ ) ) histosAre2D = true; //To allow for now systematic histos

  if ( (!histosAre1D)   && (!histosAre2D))
    {
      std::cout << "there need be 4 1-d or 2-d histograms which were not found in file: " 
		<< theFile_->GetName() << ". Text file cannot be produced" << std::endl;
      return;
    }



  textFileName_ = textFileName;

  if(histosAre1D ) produceTxtFile1D();
  if(histosAre2D ) produceTxtFile2D();

}


  
//-----------------------------------------------------------------------------------//
void EfficiencyStore::produceHistograms(std::string rootFileName){
  TFile * rootFile = new TFile(rootFileName.c_str(),"RECREATE");
  produceHistograms(rootFile);
}    


  
//-----------------------------------------------------------------------------------//
void EfficiencyStore::produceHistograms(TFile * rootFile){
    
  // check if the text file is available as needed for input
  if (textFileName_ == ""){
    std::cout << "no inpunt txt file give, cannot produce histograms." <<std::endl;
    return;
  }
    
  dimension_=0;
  char str [1000];    char effName[1000];
  ifstream is(textFileName_.c_str());

  while ( 	is.getline(str,1000)){
    // std::cout << str <<"\n";
    if ( strstr(str,"#dimension:") )           sscanf(str, "#dimension:%d",  &dimension_);
    if ( strstr(str,"#efficiency name:") )  sscanf(str, "#efficiency name: %s",  effName);
  }
  is.close();
  efficiencyName_ = std::string(effName);

  // call dedicated private method depending on the dimension of the histogram
  //  std::cout<<"histogram dimension is: " << dimension_ << std::endl;
  if(dimension_==1) produce1DHistograms(rootFile);
  //  else if(dimension_==2) produce2DHistograms(rootFile);
  else  {
    std::cout<<"*** histogram dimension is: " << dimension_ << " but only 1 and 2 are valid.\n" <<std::endl;
  }
    
  return;

}

  


//-----------------------------------------------------------------------------------//
void EfficiencyStore::produceTxtFile1D(){

  int numBinvalues1D = values1DHisto_ ->GetXaxis()->GetNbins();
  ///int numBinsystPlus1D = systPlus1DHisto_ ->GetXaxis()->GetNbins();
  ///int numBinsystMinus1D = systMinus1DHisto_->GetXaxis()->GetNbins();
  int numBindenominator1D =denominator1DHisto_->GetXaxis()->GetNbins();

  if (!
      numBinvalues1D==
	  //numBinsystPlus1D==
      // numBinsystMinus1D==
      numBindenominator1D
      )
    {
      std::cout<<"The four 1d histograms do not have the same number of bins. Cannot produce 1d text file." << std::endl;
      return;
    }    


  // assuming implicitely that if number of bins is the same
  // also the binning is the same

  std::vector<float> binMin;
  std::vector<float> binMax;
  std::vector<float> values1D;
  std::vector<float> systPlus1D;
  std::vector<float> systMinus1D;
  std::vector<float> denominator1D;

  for (int u =1; u<(1+numBinvalues1D); u++) 
    {
      binMin.push_back( values1DHisto_->GetBinLowEdge(u) );
      binMax.push_back( values1DHisto_->GetBinLowEdge(u) + values1DHisto_->GetBinWidth(u) );
      values1D.push_back( values1DHisto_->GetBinContent(u) );
      //systPlus1D.push_back( systPlus1DHisto_->GetBinContent(u) );
      //systMinus1D.push_back( systMinus1DHisto_->GetBinContent(u) );
      systPlus1D.push_back( 1.0 ); //To allow for no systematic histos
      systMinus1D.push_back( -1.0 ); //To allow for no systematic histos
      denominator1D.push_back( denominator1DHisto_->GetBinContent(u) );
    }


  std::string theHistoTitle(values1DHisto_ ->GetTitle());
  int v = theHistoTitle.find( std::string("_") );
  std::string theEffTitle = theHistoTitle.substr( (v+1) );
  //std::string theEffTitle = theHistoTitle.substr(0 ,v );
  

  ofstream the1DEffFile;
  the1DEffFile.open (textFileName_.c_str(),std::ios::out);

  the1DEffFile << "#put your comments after a '#'.\n";
  the1DEffFile << "#efficiency name: " << effName_  <<"\n";
  the1DEffFile << "#dimension: 1 \n\n\n";
  the1DEffFile << " " << std::setw(10) << "#PtMin"
	       << " " << std::setw(10) << "PtMax"
	       << " " << std::setw(10)  << "EtaMin"
	       << " " << std::setw(10)  <<  "EtaMax"
	       << " " << std::setw(10)  << "NumParms"
	       << " " << std::setw(10)  << "eff"
	       << " " << std::setw(10) << "statp"
	       << " " << std::setw(10) << "statm"
	       << " " << std::setw(10) << "den"
	       << "\n";
     
  if (the1DEffFile.is_open())
    {
	
      for (int r=0; r<numBinvalues1D; r++) 
	{
          std::vector<std::pair<float, float> > EffInfo = efftable_->GetCellInfo(r);
	  the1DEffFile << " " << std::setw(10) << std::setprecision(4) <<  EffInfo[0].first 
		       << " " << std::setw(10) << std::setprecision(4) <<  EffInfo[0].second
 		       << " " << std::setw(10) << std::setprecision(4) <<  EffInfo[1].first 
		       << " " << std::setw(10) << std::setprecision(4) <<  EffInfo[1].second
		       << " " << std::setw(10) << std::setprecision(4) <<  4
		       << " " << std::setw(10) << std::setprecision(4) <<  values1D[r]
		       << " " << std::setw(10) << std::setprecision(4) << systPlus1D[r]
		       << " " << std::setw(10) << std::setprecision(4) << systMinus1D[r]
		       << " " << std::setw(10) << std::setprecision(5) << denominator1D[r]
		       << "\n";
	}
    }

  the1DEffFile.close();
    
}

  
  

void EfficiencyStore::produce1DHistograms(TFile * rootFile){

  int binCounter = efftable_->size();

  if (rootFile!=0) rootFile->cd();
  std::string title = std::string("values_")+efficiencyName_ + "EtaDet_Pt"; /// JUST A HACK for the MOMENT 
  TH1F * valuesHisto = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("systematicMinus_")+efficiencyName_;
  TH1F * systematicMHisto = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("systematicPlus_")+efficiencyName_;
  TH1F * systematicPHisto = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("denominator_")+efficiencyName_;
  TH1F * denominatorHisto = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);

  if (rootFile==0) {
    valuesHisto->SetDirectory(0);
    systematicMHisto->SetDirectory(0);
    systematicPHisto->SetDirectory(0);
    denominatorHisto->SetDirectory(0);
    values1DHisto_=valuesHisto;
    systPlus1DHisto_=systematicMHisto;
    systMinus1DHisto_=systematicPHisto;
    denominator1DHisto_=denominatorHisto;
  }


  for (int u=1; u<=binCounter; u++)
    {
      std::vector<float> EffBinParms = efftable_->correctionEff(u-1);
      valuesHisto->SetBinContent(u, EffBinParms[1]);
      systematicPHisto->SetBinContent(u, EffBinParms[2]);
      systematicMHisto->SetBinContent(u, EffBinParms[3]);
      denominatorHisto->SetBinContent(u, EffBinParms[4]);
    }

  if (rootFile!=0) {
    valuesHisto->Write();systematicMHisto ->Write();
    systematicPHisto->Write(); denominatorHisto->Write();
    rootFile->Close();
  }
}
  

//void produce2DHistograms(TFile * rootFile){
//  std::cout<<"produce2DHistograms not yet implemented" << std::endl;
//}
  
void EfficiencyStore::produceTxtFile2D(){
  std::cout<<"produceTxtFile2D not yet implemented" << std::endl;
}
