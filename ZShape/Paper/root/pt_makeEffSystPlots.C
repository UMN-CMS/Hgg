{
#include "plotEffStat.C"
//#include "makefinal.C"

TFile *my95    = new TFile("data/TFEcalEffSystECAL80-ECAL95Def_WP95-EFFSYST_001.root");
TFile *my95hf  = new TFile("data/TFEcalEffSystECAL80-ECAL95Def_WP95-EFFSYST_001.root");
TFile *myo95p  = new TFile("WP95P.root","RECREATE");
TFile *myo95m  = new TFile("WP95M.root","RECREATE");

TFile *my80    = new TFile("data/TFEcalEffSystECAL80-ECAL95Def_WP80-EFFSYST_001.root");
TFile *my80hf  = new TFile("data/TFEcalEffSystECAL80-HFDef_WP80-EFFSYST_001.root");
TFile *myo80p  = new TFile("WP80P.root","RECREATE");
TFile *myo80m  = new TFile("WP80M.root","RECREATE");


TFile *myGSF   = new TFile("data/TFEcalEffSystECAL80-ECAL95Def_GsfTrack-EtaDet-EFFSYST_001.root");
TFile *myGSFhf = new TFile("data/TFEcalEffSystECAL80-HFDef_GsfTrack-EtaDet-EFFSYST_001.root");
TFile *myoGSFp = new TFile("GSFP.root","RECREATE");
TFile *myoGSFm = new TFile("GSFM.root","RECREATE");


TFile *myHFEID   = new TFile("data/TFEcalEffSystECAL80-HFDef_HFElectronId-EtaDet-EFFSYST_001.root");
TFile *myHFEIDhf = new TFile("data/TFEcalEffSystECAL80-HFDef_HFElectronId-EtaDet-EFFSYST_001.root");
TFile *myoHFEIDp = new TFile("HFEIDP.root","RECREATE");
TFile *myoHFEIDm = new TFile("HFEIDM.root","RECREATE");


TFile *myHFSC   = new TFile("data/TFEcalEffSystECAL80-HFDef_HFSuperCluster-Et-EFFSYST_001.root");
TFile *myHFSChf = new TFile("data/TFEcalEffSystECAL80-HFDef_HFSuperCluster-Et-EFFSYST_001.root");
TFile *myoHFSCp = new TFile("HFSCP.root","RECREATE");
TFile *myoHFSCm = new TFile("HFSCM.root","RECREATE");


TFile *mySC    = new TFile("data/TFEcalEffSystECAL80-ECAL95Def_Supercluster-Eta-EFFSYST_001.root");
TFile *mySChf  = new TFile("data/TFEcalEffSystECAL80-HFDef_Supercluster-Eta-EFFSYST_001.root");
TFile *myoSCp   = new TFile("SCP.root","RECREATE");
TFile *myoSCm   = new TFile("SCM.root","RECREATE");


TFile *myHLT   = new TFile("data/TFEcalEffSystECAL80-ECAL95Def_HLT-GSF-EFFSYST_001.root");
TFile *myHLThf = new TFile("data/TFEcalEffSystECAL80-HFDef_HLT-GSF-EFFSYST_001.root");
TFile *myoHLTp = new TFile("HLTP.root","RECREATE");
TFile *myoHLTm = new TFile("HLTM.root","RECREATE");


TFile *myNT    = new TFile("data/TFEcalEffSystECAL80-NTLooseDef_NTTightElectronId-EtaDet-EFFSYST_001.root");
TFile *myNThf  = new TFile("data/TFEcalEffSystECAL80-NTLooseDef_NTTightElectronId-EtaDet-EFFSYST_001.root");
TFile *myoNTp   = new TFile("NTP.root","RECREATE");
TFile *myoNTm   = new TFile("NTM.root","RECREATE");


plotEffSyst(my95,my95hf,myo95p, "Z0_Pt_masscut", "Plus");
myo95p->Write();
plotEffSyst(my95,my95hf,myo95m, "Z0_Pt_masscut", "Minus");
myo95m->Write();
plotEffSyst(my80,my80hf, myo80p, "Z0_Pt_masscut", "Plus");
myo80p->Write();
plotEffSyst(my80,my80hf, myo80m, "Z0_Pt_masscut", "Minus");
myo80m->Write();
plotEffSyst(myGSF,myGSFhf, myoGSFp, "Z0_Pt_masscut", "Plus");
myoGSFp->Write();
plotEffSyst(myGSF,myGSFhf, myoGSFm, "Z0_Pt_masscut", "Minus");
myoGSFm->Write();
plotEffSyst(myHFEID,myHFEIDhf,myoHFEIDp, "Z0_Pt_masscut", "Plus");
myoHFEIDp->Write();
plotEffSyst(myHFEID,myHFEIDhf,myoHFEIDm, "Z0_Pt_masscut", "Minus");
myoHFEIDm->Write();
plotEffSyst(myHFSC,myHFSChf,myoHFSCp, "Z0_Pt_masscut", "Plus");
myoHFSCp->Write();
plotEffSyst(myHFSC,myHFSChf,myoHFSCm, "Z0_Pt_masscut", "Minus");
myoHFSCm->Write();
plotEffSyst(mySC,mySChf,myoSCp, "Z0_Pt_masscut", "Plus");
myoSCp->Write();
plotEffSyst(mySC,mySChf,myoSCm, "Z0_Pt_masscut", "Minus");
myoSCm->Write();
plotEffSyst(myHLT,myHLThf,myoHLTp, "Z0_Pt_masscut", "Plus");
myoHLTp->Write();
plotEffSyst(myHLT,myHLThf,myoHLTm, "Z0_Pt_masscut", "Minus");
myoHLTm->Write();
plotEffSyst(myNT,myNThf,myoNTp, "Z0_Pt_masscut", "Plus");
myoNTp->Write();
plotEffSyst(myNT,myNThf,myoNTm, "Z0_Pt_masscut", "Minus");
myoNTm->Write();
myoNTm->Close();





}

