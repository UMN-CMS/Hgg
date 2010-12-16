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
#include "ZShape/Base/interface/EffTableLoader.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"


EfficiencyStore::EfficiencyStore()
{
  theFile_ =0;
  efftable_ =0;
  doInit();

  textFileName_="";
  
}

/*
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
*/

EfficiencyStore::EfficiencyStore(TFile * file, std::string EffName, std::string EffType)
{ 
  doInit();
  effName_ = EffName;
  effType_ = EffType;
  std::cout << "class EfficiencyStore created with root file: " << file->GetName() << std::endl;
  isEta_ = strstr(EffType.c_str(),"eta");
  isBoth_ = (bool (strstr(EffType.c_str(),"eta")) && bool (strstr(EffType.c_str(),"pt")));
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
    
  while((key=((TKey*) next()) )){
    if (strcmp(key->GetClassName(),"TDirectoryFile")==0)
      {
	TDirectory *mdir = (TDirectory*) theFile_->Get(key->GetName());
	TIter next1(mdir->GetListOfKeys());
	TKey * key1;
	while((key1=((TKey*) next1())) ){
	  //at this step I need to check if I want eta or pt or whatever...
	  //std::cout << " key1 " << key1->GetName() << " type " << key1->GetClassName() << std::endl;
	  if (strcmp(key1->GetClassName(),"TDirectoryFile")==0 && strcmp(key1->GetName(),effType_.c_str())==0 )
	    {	    
	      TDirectory *mdir1 = (TDirectory*) mdir->Get(key1->GetName());
	      TIter next2(mdir1->GetListOfKeys());
	      TKey * key2;
	      
	      while((key2=((TKey*) next2())) ){
		//std::cout << " key2 " << key2->GetName() << " type " << key2->GetClassName() << std::endl;
		if (strcmp(key2->GetClassName(),"TDirectoryFile")==0 && strstr(key2->GetName(),"bin") )
		  {
		    TDirectory *mdir2 = (TDirectory*) mdir1->Get(key2->GetName());
		    TIter next3(mdir2->GetListOfKeys());
		    TKey * key3;
		    
		    while((key3=((TKey*) next3())) ){
		      //std::cout << " key3 " << key3->GetName() << " type " << key3->GetClassName() << std::endl;
		      //here is where I get the w's!!!
		      if ( strcmp(key3->GetName(), "w")==0 )
			{
			  std::cout << " Found " << key3->GetName() << " in " << key2->GetName() << " in "<< key1->GetName() << " in "<< key->GetName() << std::endl;
 			  RooWorkspace* w = (RooWorkspace*) mdir2->Get(key3->GetName());
			  //w->Print();
			  RooRealVar *et  = (isEta_ && !isBoth_) ? w->var("probe_sc_eta"): w->var("probe_sc_et") ; //dear anyone reading this, it is EITHER et or eta!!
			  RooRealVar *eta = isBoth_ ? w->var("probe_sc_eta"): 0 ;
			  //Here I will one day need to do 'two-D'... that is why I have the bin1 and bin2
			  RooRealVar *eff = w->var("efficiency");
			  RooRealVar *nS  = w->var("numSignalAll");
			  if (! et ) std::cout << " no et? " << std::endl;
			  if (! eff ) std::cout << " no eff? " << std::endl;
			  if (! nS ) std::cout << " no nS? " << std::endl;
			  if (! nS || ! et ) continue; 
			  int hbinnum =et->getBinning().binNumber(et->getVal());
			  Double_t hval = et->getBinning().binHigh(hbinnum);
			  Double_t lval = et->getBinning().binLow(hbinnum);

			  //std::cout << "Efficiency: " << eff->getVal() << " Den: " << nS->getVal() <<  std::endl;
			  std::cout << "Efficiency: " << eff->getVal() << " Den: " << nS->getVal() << " et min: " << lval << " et max: " <<   hval << " Eff Err Hi: " <<eff->getErrorHi() << " Eff Err Lo: " <<eff->getErrorLo() << std::endl;  
			  //crazy hack below... look into why I needed to do something like this...
			  //MAYBE for more events using the default errors will be just 'fine'.
			  Double_t errorHi = eff->getErrorHi() < 0.0001 ? (1.0 - eff->getVal()) : eff->getErrorHi() ;
			  std::cout << " new error high is " << errorHi << std::endl;
			  //do the bins different based on eta or pt
			  if (isBoth_)
			  {
			    int hbinnumeta =eta->getBinning().binNumber(eta->getVal());
			    Double_t hvaleta = eta->getBinning().binHigh(hbinnumeta);
			    Double_t lvaleta = eta->getBinning().binLow(hbinnumeta); 
			    bin1Min_.push_back( lval );     bin1Max_.push_back( hval);
			    bin2Min_.push_back( lvaleta );  bin2Max_.push_back( hvaleta );
			  }
			  else 
			  {
			    if (isEta_) {bin2Min_.push_back( lval );  bin2Max_.push_back( hval );}
			    else { bin1Min_.push_back( lval );  bin1Max_.push_back( hval);}			  
			  }
			  values1D_.push_back( eff->getVal());
			  denominator1D_.push_back( nS->getVal());
			  //systPlus1D_.push_back(errorHi);
			  //systMinus1D_.push_back(eff->getErrorLo());

			  //HERE I actually need to do a much better systematic treatment
			  // Use eff = (Npass)/(Npass+Nfail)
			  // Which we can write as....
			  // eff = (Npasstotal- BPass)/(Ntotal - BPass - BFail)

			  //MY PROBLEM HERE IS THAT I CALCULARE BPass and BFail for the entire 
			  // MASS region. Thus, the numbers can vary but not affect the efficiency. 

			  // The problem is that the subtraction is implicit but can vary based on numerics.
			  // SO for the time being the systematic due to the efficiency background
			  //    subtraction will be over estimated...
			  //FIRST I NEED To DEFINE A RESTRICTED MASS RANGE
			  
			  w->var("mass")->setRange("tight",80,100); 
			  //Here we assume 100% of the +- fails in this range can affect the efficiency systematically
			  
			  
			  RooRealVar *numBackgroundPass = w->var("numBackgroundPass");
			  Double_t NPScale =  w->pdf("backgroundPass")->createIntegral(*w->var("mass"),RooFit::NormSet(*w->var("mass")),RooFit::Range("tight"))->getVal();
			  Double_t NBP   = numBackgroundPass->getVal();
			  Double_t NBPLo = numBackgroundPass->getErrorLo()*NPScale;
			  Double_t NBPHi = numBackgroundPass->getErrorHi()*NPScale;

			  RooRealVar *numBackgroundFail = w->var("numBackgroundFail");
			  Double_t NFScale =  w->pdf("backgroundFail")->createIntegral(*w->var("mass"),RooFit::NormSet(*w->var("mass")),RooFit::Range("tight"))->getVal();
			  Double_t NBF   = numBackgroundFail->getVal();
			  Double_t NBFLo = numBackgroundFail->getErrorLo()*NFScale;
			  Double_t NBFHi = numBackgroundFail->getErrorHi()*NFScale;

			  Double_t NSig = (nS->getVal());
			  Double_t EFF =  eff->getVal();

			  Double_t NT =  NSig;
			  Double_t NP =  NSig*EFF;

			  setEffSysts(EFF,NT,NP,NBPLo,NBPHi,NBFLo,NBFHi);
			  systPlus1D_.push_back(effSystHi_);
			  systMinus1D_.push_back(effSystLo_);



			  // std::cout << " has a range " << et->hasRange() << std::endl; 
			  //std::cout << " ave et is " << et->getVal() << std::endl;
			  //std::cout << " numbins " << et->numBins() << std::endl;
			  //et->printExtras(std::cout);
			  //int highbin =et->getBinning().binNumber(et->getVal());
			  //std::cout << std::endl << " bin info " << et->getBinning().binHigh(highbin) << " low " << et->getBinning().binLow(highbin) << std::endl;
			}
		    }
		  }
	      }
	    }
	}
      }
  }
    
  vecToHist();
  std::cout << "leaving setRootFile  values1DHisto_ :" << values1DHisto_ << std::endl;
  
}


//-----------------------------------------------------------------------------------//
void EfficiencyStore::vecToHist(void){
  //Simple internal program converts the vector quantities into the histogram ones.
  //Maybe I should keep some 'bools' around to test if the vec's or hist's are around.....

  int binCounter = values1D_.size();

  efficiencyName_ = "empty";

  std::string title = std::string("values_")+efficiencyName_ + "EtaDet_Pt"; /// JUST A HACK for the MOMENT 
  values1DHisto_ = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("systematicMinus_")+efficiencyName_;
  systMinus1DHisto_ = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("systematicPlus_")+efficiencyName_;
  systPlus1DHisto_ = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("denominator_")+efficiencyName_;
  denominator1DHisto_ = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);

  for (int i=0; i<binCounter; ++i)
    {
      values1DHisto_->SetBinContent(i, values1D_[i] );
      systPlus1DHisto_->SetBinContent(i, systPlus1D_[i]);
      systMinus1DHisto_->SetBinContent(i, systMinus1D_[i]);
      denominator1DHisto_->SetBinContent(i, denominator1D_[i]);
    }
  
}

//-----------------------------------------------------------------------------------//
void EfficiencyStore::histToVec(void){
  
}


//-----------------------------------------------------------------------------------//
void EfficiencyStore::setVarMinMax(double vmin, double vmax){
  vmax_=vmax;
  vmin_=vmin;
  if (isBoth_) return;
  for (uint j = 0; j < values1D_.size(); ++j)
    {
      if (isEta_) 
	{
	  bin1Min_.push_back(vmin);
	  bin1Max_.push_back(vmax);
	}
      else
	{
	  bin2Min_.push_back(vmin);
	  bin2Max_.push_back(vmax);
	}
    }
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
 void EfficiencyStore::setEffSysts(Double_t EFF, Double_t NT, Double_t NP, Double_t NBPLo, Double_t NBPHi, Double_t NBFLo, Double_t NBFHi){
  //This function is designed to find the lowest and highest systematic changes to the efficiency
  
  Double_t EFFBPLo = (NT+NBPLo) > 0 ? (NP+NBPLo)/(NT+NBPLo) : EFF;
  Double_t EFFBPHi = (NT+NBPHi) > 0 ? (NP+NBPHi)/(NT+NBPHi) : EFF;
  Double_t EFFBFLo = (NT+NBFLo) > 0 ? (NP)/(NT+NBFLo) : EFF;
  Double_t EFFBFHi = (NT+NBFHi) > 0 ? (NP)/(NT+NBFHi) : EFF;

  Double_t EFFLO = EFF;
  Double_t EFFHI = EFF;
  
  if ( EFFBPLo < EFFLO ) EFFLO = EFFBPLo ; 
  else if ( EFFBPLo > EFFHI ) EFFHI = EFFBPLo ; 
  
  if ( EFFBPHi < EFFLO ) EFFLO = EFFBPHi ; 
  else if ( EFFBPHi > EFFHI ) EFFHI = EFFBPHi ; 
  
  if ( EFFBFLo < EFFLO ) EFFLO = EFFBFLo ; 
  else if ( EFFBFLo > EFFHI ) EFFHI = EFFBFLo ; 
  
  if ( EFFBFHi < EFFLO ) EFFLO = EFFBFHi ; 
  else if ( EFFBFHi > EFFHI ) EFFHI = EFFBFHi ;

  if ( EFFHI > 1.0 ) EFFHI = 1.0;
  if ( EFFLO < 0.0 ) EFFLO = 0.0;
  
  eff_ = EFF;
  effSystLo_ = EFFLO;
  effSystHi_ = EFFHI;
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

  //int numBinvalues1D = values1DHisto_ ->GetXaxis()->GetNbins();
  //int numBinsystPlus1D = systPlus1DHisto_ ->GetXaxis()->GetNbins();
  //int numBinsystMinus1D = systMinus1DHisto_->GetXaxis()->GetNbins();
  //int numBindenominator1D =denominator1DHisto_->GetXaxis()->GetNbins();

  //if (!
  //    (numBinvalues1D==
  //    numBinsystPlus1D==
  //    numBinsystMinus1D==
  //     numBindenominator1D)
  //    )
  //  {
  //    std::cout<<"The four 1d histograms do not have the same number of bins. Cannot produce 1d text file." << std::endl;
  //    return;
  //  }    




  //std::string theHistoTitle(values1DHisto_ ->GetTitle());
  //int v = theHistoTitle.find( std::string("_") );
  //std::string theEffTitle = theHistoTitle.substr( (v+1) );
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
	       << " " << std::setw(10) << "systp"
	       << " " << std::setw(10) << "systm"
	       << " " << std::setw(10) << "den"
	       << "\n";
     
  if (the1DEffFile.is_open())
    {
	
      for (uint r=0; r< bin1Min_.size(); r++) 
	{
          //std::vector<std::pair<float, float> > EffInfo = efftable_->GetCellInfo(r);
	  the1DEffFile << " " << std::setw(10) << std::setprecision(4) <<  bin1Min_[r]  
		       << " " << std::setw(10) << std::setprecision(4) <<  bin1Max_[r]
 		       << " " << std::setw(10) << std::setprecision(4) <<  bin2Min_[r]
		       << " " << std::setw(10) << std::setprecision(4) <<  bin2Max_[r]
		       << " " << std::setw(10) << std::setprecision(4) <<  4
		       << " " << std::setw(10) << std::setprecision(4) <<  values1D_[r]
		       << " " << std::setw(10) << std::setprecision(4) << systPlus1D_[r]
		       << " " << std::setw(10) << std::setprecision(4) << systMinus1D_[r]
		       << " " << std::setw(10) << std::setprecision(5) << denominator1D_[r]
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
