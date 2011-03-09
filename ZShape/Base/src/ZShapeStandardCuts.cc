#include "ZShape/Base/interface/ZShapeStandardCuts.h"

void ZShapeStandardCuts::acceptanceCuts(ZShapeElectron& elec) const {
  elec.cutResult("ACC(EB)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EB));
  elec.cutResult("ACC(EE)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE));
  elec.cutResult("ACC(MUON)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_MUON));
  elec.cutResult("ACC(EE+MUON)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_muon));
  elec.cutResult("ACC(EE+TRK)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_tracker));
  elec.cutResult("ACC(ECAL)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ECAL));
  elec.cutResult("ACC(ECAL+TRK)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ECAL) && (!theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_notracker)));
  elec.cutResult("ACC(ECAL+MUON)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ECAL_muon));
  elec.cutResult("ACC(ECAL-TRK)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ECAL) && (theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_notracker)));
  elec.cutResult("ACC(HF)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_HF));
  elec.cutResult("ACC(ANY)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ANY));
  elec.cutResult("ACC(ANY-EEnotrk)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ANY) && (!theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_notracker)));
}

void ZShapeStandardCuts::ptCuts(ZShapeElectron& elec) const {
  double pt=elec.p4_.Pt();
  elec.cutResult("PT0",pt>=-50000);
  elec.cutResult("PT5",pt>=5);
  elec.cutResult("PT10",pt>=10);
  elec.cutResult("PT15",pt>=15);
  elec.cutResult("PT17",pt>=17);
  elec.cutResult("PT20",pt>=20);
  elec.cutResult("PT25",pt>=25);
  elec.cutResult("PT30",pt>=30);
}

void ZShapeStandardCuts::gen_ptCuts(reco::GenParticle& gelec,ZShapeElectron& elec) const {
  double pt=gelec.polarP4().Pt();
  elec.cutResult("GPT0",pt>=-50000);
  elec.cutResult("GPT5",pt>=5);
  elec.cutResult("GPT10",pt>=10);
  elec.cutResult("GPT15",pt>=15);
  elec.cutResult("GPT17",pt>=17);
  elec.cutResult("GPT20",pt>=20);
  elec.cutResult("GPT25",pt>=25);
  elec.cutResult("GPT30",pt>=30);
}


void ZShapeStandardCuts::dummyCuts(ZShapeElectron& elec) const {
  elec.cutResult("DummyPass",true);
  elec.cutResult("DummyFail",false);
}
