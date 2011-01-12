#ifndef ZShape_Base_Systematics_h_included
#define ZShape_Base_Systematics_h_included

#include "ZShapeElectron.h"

namespace zshape {
  
  class EnergyScale {
  public:
    EnergyScale(int ECAL, int HF); 
    void rescale(ZShapeElectron& electron);
  private:
    int shiftECAL_,shiftHF_;
  };

}

#endif // ZShape_Base_Systematics_h_included
