#ifndef EFFCUT_INC
#define EFFCUT_INC

#include "TH1.h"

class EfficiencyCut
{
  public :

  EfficiencyCut( void );
  EfficiencyCut( TH1F * histo );

  bool passesCut( float variable );

 private:
  TH1F * theClonedEffHisto_;
  TH1F * theTmpHisto_;

};
#endif
