#include "ZShape/EffAcc/interface/EfficiencyCut.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

TRandom* EfficiencyCut::randomGenerator_=0; 

EfficiencyCut::EfficiencyCut ( TH1F * histo, EffTableLoader* ind) : theIndexer(ind) {

  if (randomGenerator_==0) randomGenerator_=new TRandom3(0);
  
  if(! histo)
    {
      std::cout << "null TF1F* pointer in EfficiencyCut constructor, returning."<< std::endl;
      return;
    }

  std::cout << "building EfficiencyCut - constructor - num entries is: " << histo->GetEntries() << std::endl;
  //  std::cerr << "theTmpHisto_ is: " << theTmpHisto_<< std::endl;

  char name[1024];
  snprintf(name,1023,"%s ",histo->GetName());
  theClonedEffHisto_ = (TH1F *) histo->Clone(name);
  theClonedEffHisto_ ->SetDirectory(0);
}
 
bool EfficiencyCut::passesCut(int index, float level) const
{
   
  int theBin =  index + 1;
  if (theBin <= 0 ) 
    {
      return false; // underflow
    }
  if (theBin ==  ( theClonedEffHisto_->GetNbinsX()  +1)  ) 
    {

      return false; // underflow
    }


  float theEfficiency = theClonedEffHisto_->GetBinContent(theBin);

 
  float randNum;
  if (level<0) randNum=randomGenerator_->Uniform(0.,1.);
  else randNum=level;
  lastRandomLevel_=randNum;

    if ( randNum < theEfficiency ) 
       {
	   //std::cout << "EfficiencyCut passescut in cut: " <<  theClonedEffHisto_->GetTitle() <<  " which falls in bin: " << theBin << " and the cut was passed" 
	  //            " (randNum: " << randNum << " eff: " << theEfficiency << ")" <<std::endl;
	 // std::cout << " Returning TRUE " << std::endl;
	 return true;
       }
   else{
           //  std::cout << "iEfficiencyCut passescut n cut: " <<  theClonedEffHisto_->GetTitle() << " which falls in bin: " << theBin << " and the cut was not passed" 
         //	" (randNum: " << randNum << " eff: " << theEfficiency << ")" << std::endl;
           //   std::cout << " Returning FALSE " << std::endl;
         return false;
    }

  std::cout << " Returning FALSE " << std::endl;
  return false;

} 


bool EfficiencyCut::passesCut( const ZShapeElectron& elec) const {
  return passesCut(indexOf(elec));
}

bool EfficiencyCut::passesCut( const ZShapeElectron& elec, float level) const {
  return passesCut(indexOf(elec),level);
}


int EfficiencyCut::indexOf( const ZShapeElectron& elec) const {
  return theIndexer->GetBandIndex(elec.p4_.Pt(),elec.detEta_);
}

