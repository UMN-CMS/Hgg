#include "ZShape/Base/interface/ZShapeEvent.h"
#include "ZShape/Base/interface/ZShapeZDef.h"

void ZShapeZDef::addCriterion(ZShapeZDef::CriterionType ct, const std::string& crit) {
  m_requirements[ct].push_back(crit);
}

int ZShapeZDef::criteriaCount(ZShapeZDef::CriterionType ct) const {
  return (int)(m_requirements[ct].size());
}

static const std::string BAD("Overflow Error");
const std::string& ZShapeZDef::criteria(ZShapeZDef::CriterionType ct, int n) const {
  return (n<0 || n>=criteriaCount(ct))?(BAD):(m_requirements[ct][n]);
}

bool ZShapeZDef::pass(const ZShapeEvent& evt, bool* swapPairing) const {
  return pass(evt,m_requirements[0].size(),m_requirements[1].size(),m_requirements[2].size(),swapPairing);
}

bool ZShapeZDef::pass(const ZShapeEvent& evt, int ne1, int ne2, int nz,bool* swapPairing) const {
  // try straight-up assignment first
  bool p=pass(evt.elec(0),evt.elec(1),ne1,ne2,nz);
  if (!p) { // try swapping electrons!
    p=pass(evt.elec(1),evt.elec(0),ne1,ne2,nz);
    if (swapPairing!=0) *swapPairing=true;
  } else if (swapPairing!=0) *swapPairing=false;
  return p;
}

bool ZShapeZDef::pass(const ZShapeElectron& e1, const ZShapeElectron& e2, int ne1, int ne2, int nz) const {
  bool ok=true;
  // concrete assignment
  for (unsigned int i=0; ok && int(i)<ne1 && i<m_requirements[crit_E1].size(); ++i) {
    if (!e1.passed(m_requirements[crit_E1][i])) ok=false;
  }
  for (unsigned int i=0; ok && int(i)<ne2 && i<m_requirements[crit_E2].size(); ++i) {
    if (!e2.passed(m_requirements[crit_E2][i])) ok=false;
  }
  /*
  for (unsigned int i=0; ok && i<ne1 && i<m_requirements[0].size(); ++i) {
    if (!e1.passed(m_requirements[0][i])) ok=false;
  }
  */
  return ok;
}
