#include "ZShape/EffAcc/interface/EfficiencyCut.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"


EfficiencyCut::EfficiencyCut ( void ) {
  std::cerr << "in EfficiencyCut void constructor" << std::endl ;
}

EfficiencyCut::EfficiencyCut ( TH1F * histo ):theTmpHisto_(histo) {
  
  if(! histo)
    {
      std::cout << "null TF1F* pointer in EfficiencyCut constructor, returning."<< std::endl;
      return;
    }

  std::cout << "building EfficiencyCut - constructor - num entries is: " << theTmpHisto_->GetEntries() << std::endl;
  std::cerr << "theTmpHisto_ is: " << theTmpHisto_<< std::endl;
  
  theClonedEffHisto_ = (TH1F *) theTmpHisto_->Clone("clonedEfficiencyHisto");
  //   std::cout << "cloned EfficiencyCut obj created with histo named: " 
  // 	    <<   theClonedEffHisto_ ->GetTitle() 
  // 	    << " and pointer: " << theClonedEffHisto_ << std::endl;
  theClonedEffHisto_ ->SetDirectory(0);
}


bool EfficiencyCut::passesCut(float variable)
{
  //  std::cerr << "in EfficiencyCut passesCut, original histo is: " <<   theTmpHisto_ ->GetTitle() << " and pointer: " <<  theTmpHisto_ << std::endl;
  std::cerr << "in EfficiencyCut passesCut, cloned histo is: " <<   theClonedEffHisto_ ->GetTitle() << " and pointer: " << theClonedEffHisto_ << std::endl;
   
  int theBin =  theClonedEffHisto_ -> Fill(variable,0);
  if (theBin == -1 ) 
    {
      std::cout << "in EfficiencyCut passesCut piff underflow" << std::endl;
      return false; // underflow
    }
  if (theBin ==  ( theClonedEffHisto_->GetNbinsX()  +1)  ) 
    {
      std::cout << "in EfficiencyCut passesCut overflow" << std::endl;
      return false; // underflow
    }


  float theEfficiency = theClonedEffHisto_->GetBinContent(theBin);

  TRandom3 randomNum(0); // put this back
  float randNum = randomNum.Uniform(0., 1.); // put this back

  if ( randomNum.Uniform(0., 1.) < theEfficiency ) 
     if ( randNum < theEfficiency ) 
       {
	 //  	std::cout << "EfficiencyCut passescut in cut: " <<  theClonedEffHisto_->GetTitle() << " variable is: " <<
	 //  	  variable << " which falls in bin: " << theBin << " and the cut was passed" 
	 //  	  " (randNum: " << randNum << " eff: " << theEfficiency << ")" <<std::endl;
	 return true;
       }
     else{
       //        std::cout << "iEfficiencyCut passescut n cut: " <<  theClonedEffHisto_->GetTitle() << " variable is: " <<
       //  	variable << " which falls in bin: " << theBin << " and the cut was not passed" 
       //  	" (randNum: " << randNum << " eff: " << theEfficiency << ")" << std::endl;
       return false;
     }


  return true;

} 
 
