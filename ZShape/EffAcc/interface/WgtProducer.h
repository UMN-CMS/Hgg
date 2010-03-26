#ifndef _WGTPRODUCER_H
#define _WGTPRODUCER_H

// system include files
#include <string>

// user include files

class TH2;

//
// class decleration
//

class WgtProducer {
public:
  
  WgtProducer() {}
  void load(const std::string& filename);
  double wgt(const double& pt, const double& y);

private:

  TH2* weights_;
  bool nowgts_;
 
};

#endif

