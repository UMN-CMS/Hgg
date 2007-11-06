#ifndef ZSHAPE_BASE_ZSHAPEZDEF_H
#define ZSHAPE_BASE_ZSHAPEZDEF_H 1

#include <vector>
#include <string>
class ZShapeEvent;

/** \class ZShapeZDef
  *
  * Basically, a Z-definition is a list of cuts (by name) which must pass.
  * These cuts apply to the first electron, the second electron, and then
  * the di-electron. 
  * 
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class ZShapeZDef {
public:
  enum CriterionType { crit_Z=0, crit_E1=1, crit_E2=2 };
  /// add a criterion
  void addCriterion(CriterionType ctype, const std::string& cit);
  /// get the number of criteria
  int criteriaCount(CriterionType ctype) const;
  /// get the n-th criterion
  const std::string& criteria(CriterionType ctype, int n) const;
  /// does this event pass all criteria of this Z Definition?
  bool pass(const ZShapeEvent& evt, bool* swapPairing=0) const;
  /// does this event pass all criteria of this Z Definition up-to-and-including the n-th for each item?
  bool pass(const ZShapeEvent& evt, int ne1, int ne2=1000, int nz=1000, bool* swapPairing=0) const;
private:
  bool pass(const ZShapeEvent::ElectronStruct& e1, const ZShapeEvent::ElectronStruct& e2, int ne1, int ne2, int nz) const;
  std::vector<std::string> m_requirements[3];
};

#endif
