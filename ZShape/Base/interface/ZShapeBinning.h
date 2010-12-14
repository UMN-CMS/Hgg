#ifndef ZShapeBinning_h_included
#define ZShapeBinning_h_included

namespace zshape {

  static const int    y_bins=100;
  static const double y_min=-5.0;
  static const double y_max= 5.0;

  static const int    pt_bins=19;
  static const double pt_binning[]={ 0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };

}


#endif // ZShapeBinning_h_included
