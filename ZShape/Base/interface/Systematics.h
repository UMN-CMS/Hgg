#ifndef ZShape_Base_Systematics_h_included
#define ZShape_Base_Systematics_h_included

#include "ZShapeElectron.h"

namespace zshape {
  
  class EnergyScale {
  public:
    EnergyScale(bool positive); 
    void rescale(ZShapeElectron& electron);
  private:
    bool upshift_;
  };

}

#endif // ZShape_Base_Systematics_h_included
