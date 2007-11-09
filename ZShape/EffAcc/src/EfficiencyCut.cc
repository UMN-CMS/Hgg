#include "ZShape/EffAcc/interface/EfficiencyCut.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


EfficiencyCut::EfficiencyCut ( void ) {
  std::cerr << "in EfficiencyCut void constructor" << std::endl ;
}

EfficiencyCut::EfficiencyCut ( TH1F * histo ) {
  
  if(! histo)
    {
      std::cout << "null TF1F* pointer in EfficiencyCut constructor, returning."<< std::endl;
      return;
    }

  std::cout << "building EfficiencyCut - constructor - num entries is: " << histo->GetEntries() << std::endl;
  //  std::cerr << "theTmpHisto_ is: " << theTmpHisto_<< std::endl;
  
  theClonedEffHisto_ = (TH1F *) histo->Clone("clonedEfficiencyHisto");
  //   std::cout << "cloned EfficiencyCut obj created with histo named: " 
  // 	    <<   theClonedEffHisto_ ->GetTitle() 
  // 	    << " and pointer: " << theClonedEffHisto_ << std::endl;
  theClonedEffHisto_ ->SetDirectory(0);
  std::string variable = histo->GetName();
  if (strstr(variable.c_str(),"eta")      || strstr(variable.c_str(),"Eta") ) theCutVariable_=cv_Eta;
  else if (strstr(variable.c_str(),"phi") || strstr(variable.c_str(),"Phi") ) theCutVariable_=cv_Phi;
  else if (strstr(variable.c_str(),"pt")  || strstr(variable.c_str(),"et")  ||
	   strstr(variable.c_str(),"p_t") || strstr(variable.c_str(),"e_t")||
           strstr(variable.c_str(),"Pt") )
    theCutVariable_=cv_Pt;
  else {
    edm::LogWarning("ZShape") << "Unknown variable in efficiency cut: '" << variable << "'";
  }
  //std::cout << "For CutVariable " <<  variable << " the type is " << theCutVariable_ << std::endl;
}

bool EfficiencyCut::passesCut(float variable) const
{
  //  std::cerr << "in EfficiencyCut passesCut, original histo is: " <<   theTmpHisto_ ->GetTitle() << " and pointer: " <<  theTmpHisto_ << std::endl;
  //  std::cerr << "in EfficiencyCut passesCut, cloned histo is: " <<   theClonedEffHisto_ ->GetTitle() << " and pointer: " << theClonedEffHisto_ << std::endl;
   
  int theBin =  theClonedEffHisto_ -> FindBin(variable);
  if (theBin <= 0 ) 
    {
      std::cout << "in EfficiencyCut passesCut piff underflow : " << variable << std::endl;
      return false; // underflow
    }
  if (theBin ==  ( theClonedEffHisto_->GetNbinsX()  +1)  ) 
    {
      std::cout << "in EfficiencyCut passesCut overflow : " << variable << std::endl;
      return false; // underflow
    }


  float theEfficiency = theClonedEffHisto_->GetBinContent(theBin)
;
  TRandom3 randomNum(0); // put this back
  float randNum = randomNum.Uniform(0., 1.); // put this back
  //std::cout << "EffCut is " << theEfficiency << " in " << theClonedEffHisto_->GetTitle() << " variable is: " << variable << " which falls in bin: " << theBin << std::endl;
  if ( randomNum.Uniform(0., 1.) < theEfficiency ) 
     if ( randNum < theEfficiency ) 
       {
	 //  	std::cout << "EfficiencyCut passescut in cut: " <<  theClonedEffHisto_->GetTitle() << " variable is: " <<
	 //  	  variable << " which falls in bin: " << theBin << " and the cut was passed" 
	 // 	  " (randNum: " << randNum << " eff: " << theEfficiency << ")" <<std::endl;
	 return true;
       }
     else{
       //     std::cout << "iEfficiencyCut passescut n cut: " <<  theClonedEffHisto_->GetTitle() << " variable is: " <<
       //	variable << " which falls in bin: " << theBin << " and the cut was not passed" 
       //	" (randNum: " << randNum << " eff: " << theEfficiency << ")" << std::endl;
       return false;
     }


  return true;

} 
 

bool EfficiencyCut::passesCut( const ::math::PtEtaPhiMLorentzVector& elec) const {
  switch (theCutVariable_) {
  case (cv_Eta) : return passesCut(elec.eta());
  case (cv_Phi) : return passesCut(elec.phi());
  case (cv_Pt) : return passesCut(elec.Pt());
  case (cv_Energy) : return passesCut(elec.E());
  default: return false;
  }
}
