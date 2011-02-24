#ifndef ZShapeBinning_h_included
#define ZShapeBinning_h_included

namespace zshape {

  static const int    y_bins=100;
  static const double y_min=-5.0;
  static const double y_max= 5.0;

  static const int    pt_bins=19;
  static const double pt_binning[]={ 0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };
  //static const double pt_binning[]={ 0, 4.0, 8.0, 12.0, 16.0, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };

  static const int e_eta_bins=14;
  static const double e_eta_binning[]={-5.,-3., -2.5, -2., -1.5, -1., -0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0};
  
  static const int e_pt_bins=7;
  static const double e_pt_binning[]={20.,25., 30., 35., 40., 50., 120., 500.};
  
}


#endif // ZShapeBinning_h_included
