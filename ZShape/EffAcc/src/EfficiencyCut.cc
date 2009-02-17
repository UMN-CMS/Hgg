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

  char name[1024];
  snprintf(name,1023,"%s ",histo->GetName());
  theClonedEffHisto_ = (TH1F *) histo->Clone(name);
  //   std::cout << "cloned EfficiencyCut obj created with histo named: " 
  // 	    <<   theClonedEffHisto_ ->GetTitle() 
  // 	    << " and pointer: " << theClonedEffHisto_ << std::endl;
  theClonedEffHisto_ ->SetDirectory(0);
  std::string variable = histo->GetName();
  if (strstr(variable.c_str(),"DetEta")      || strstr(variable.c_str(),"EtaDet") ) theCutVariable_=cv_DetEta;
  else if (strstr(variable.c_str(),"eta")      || strstr(variable.c_str(),"Eta") ) theCutVariable_=cv_Eta;
  else if (strstr(variable.c_str(),"phi") || strstr(variable.c_str(),"Phi") ) theCutVariable_=cv_Phi;
  else if (strstr(variable.c_str(),"pt")  || strstr(variable.c_str(),"et")  ||
	   strstr(variable.c_str(),"p_t") || strstr(variable.c_str(),"e_t")||
           strstr(variable.c_str(),"Pt") )
    theCutVariable_=cv_Pt;
  else {
    edm::LogWarning("ZShape") << "Unknown variable in efficiency cut: '" << variable << "'";
  }
  if ( (theCutVariable_== cv_DetEta) && (strstr(variable.c_str(),"Pt"))) theCutVariable_=cv_DetEta_Pt;
  std::cout << "For CutVariable " <<  variable << " the type is " << theCutVariable_ << std::endl;
}

bool EfficiencyCut::passesCut(float variable) const
{
  //  std::cerr << "in EfficiencyCut passesCut, original histo is: " <<   theTmpHisto_ ->GetTitle() << " and pointer: " <<  theTmpHisto_ << std::endl;
  //  std::cerr << "in EfficiencyCut passesCut, cloned histo is: " <<   theClonedEffHisto_ ->GetTitle() << " and pointer: " << theClonedEffHisto_ << std::endl;
   
  int theBin =  theClonedEffHisto_ -> FindBin(variable);
  if (theBin <= 0 ) 
    {
      /*
      std::cout << "in " << theClonedEffHisto_->GetName() <<" passesCut piff underflow : " << variable;
      switch (theCutVariable_) {
      case(cv_Eta) : std::cout << " eta"; break;
      case(cv_Phi) : std::cout << " phi"; break;
      case(cv_Pt) : std::cout << " pt"; break;
      case(cv_Energy) : std::cout << " energy"; break;
      }
      std::cout << std::endl;
      */
      return false; // underflow
    }
  if (theBin ==  ( theClonedEffHisto_->GetNbinsX()  +1)  ) 
    {
      /*
      std::cout << "in " << theClonedEffHisto_->GetName() << " passesCut overflow : " << variable;
      switch (theCutVariable_) {
      case(cv_Eta) : std::cout << " eta"; break;
      case(cv_Phi) : std::cout << " phi"; break;
      case(cv_Pt) : std::cout << " pt"; break;
      case(cv_Energy) : std::cout << " energy"; break;
      }
      std::cout << std::endl;
      */
      return false; // underflow
    }


  float theEfficiency = theClonedEffHisto_->GetBinContent(theBin);

  TRandom3 randomNum(0); // put this back
  float randNum = randomNum.Uniform(0., 1.); // put this back
  //std::cout << "EffCut is " << theEfficiency << " in " << theClonedEffHisto_->GetTitle() << " variable is: " << variable << " which falls in bin: " << theBin << std::endl;
  // GF: I do not understand this: why comparing rand with theEfficiency twice?
  //if ( randomNum.Uniform(0., 1.) < theEfficiency ) 
    if ( randNum < theEfficiency ) 
       {
	 //  std::cout << "EfficiencyCut passescut in cut: " <<  theClonedEffHisto_->GetTitle() << " variable is: " <<
	 // variable << " which falls in bin: " << theBin << " and the cut was passed" 
	 //             " (randNum: " << randNum << " eff: " << theEfficiency << ")" <<std::endl;
	 // std::cout << " Returning TRUE " << std::endl;
	 return true;
       }
   else{
         //     std::cout << "iEfficiencyCut passescut n cut: " <<  theClonedEffHisto_->GetTitle() << " variable is: " <<
         //	variable << " which falls in bin: " << theBin << " and the cut was not passed" 
         //	" (randNum: " << randNum << " eff: " << theEfficiency << ")" << std::endl;
         //     std::cout << " Returning FALSE " << std::endl;
         return false;
    }

  std::cout << " Returning FALSE " << std::endl;
  return false;

} 
 
bool EfficiencyCut::passesCut(int index) const
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

  TRandom3 randomNum(0); // put this back
  float randNum = randomNum.Uniform(0., 1.); // put this back

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
  switch (theCutVariable_) {
  case (cv_DetEta) : return passesCut(float(elec.detEta_));
  case (cv_Eta) : return passesCut(float(elec.p4_.eta()));
  case (cv_Phi) : return passesCut(float(elec.p4_.phi()));
  case (cv_Pt) : return passesCut(float(elec.p4_.Pt()));
  case (cv_Energy) : return passesCut(float(elec.p4_.E()));
  default: return false;
  }
}


bool EfficiencyCut::passesCut( const ZShapeElectron& elec, EffTableLoader* effTable) const {
  switch (theCutVariable_) {
  case (cv_DetEta_Pt) : {
       int index = effTable->GetBandIndex(elec.p4_.Pt(),elec.detEta_);
       //std::cout << " Index " << index << " Pt " << elec.p4_.Pt() << " DetEta " << elec.detEta_ << std::endl;
       return passesCut(index);}
  default: return passesCut(elec);
  }
}

