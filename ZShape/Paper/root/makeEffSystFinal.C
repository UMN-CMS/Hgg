
{
#include "plotEffStat.C"
//#include "makefinal.C"


TFile *myo95p  = new TFile("WP95P.root");
TFile *myo95m  = new TFile("WP95M.root");
TFile *myo80p  = new TFile("WP80P.root");
TFile *myo80m  = new TFile("WP80M.root");
TFile *myoGSFp = new TFile("GSFP.root");
TFile *myoGSFm = new TFile("GSFM.root");
TFile *myoHFEIDp = new TFile("HFEIDP.root");
TFile *myoHFEIDm = new TFile("HFEIDM.root");
TFile *myoHFSCp = new TFile("HFSCP.root");
TFile *myoHFSCm = new TFile("HFSCM.root");
TFile *myoSCp   = new TFile("SCP.root");
TFile *myoSCm   = new TFile("SCM.root");
TFile *myoHLTp = new TFile("HLTP.root");
TFile *myoHLTm = new TFile("HLTM.root");
TFile *myoNTp   = new TFile("NTP.root");
TFile *myoNTm   = new TFile("NTM.root");

TFile *myo95  = new TFile("WP95.root","RECREATE");
TFile *myo80  = new TFile("WP80.root","RECREATE");
TFile *myoGSF = new TFile("GSF.root","RECREATE");
TFile *myoHFEID = new TFile("HFEID.root","RECREATE");
TFile *myoHFSC = new TFile("HFSC.root","RECREATE");
TFile *myoSC   = new TFile("SC.root","RECREATE");
TFile *myoHLT = new TFile("HLT.root","RECREATE");
TFile *myoNT   = new TFile("NT.root","RECREATE");

makefinal(myo95p, myo95m, myo95);
myo95->Write();
makefinal(myo80p, myo80m, myo80);
myo80->Write();
makefinal(myoGSFp, myoGSFm, myoGSF);
myoGSF->Write();
makefinal(myoHFEIDp, myoHFEIDm, myoHFEID);
myoHFEID->Write();
makefinal(myoHFSCp, myoHFSCm, myoHFSC);
myoHFSC->Write();
makefinal(myoSCp, myoSCm, myoSC);
myoSC->Write();
makefinal(myoHLTp, myoHLTm, myoHLT);
myoHLT->Write();
makefinal(myoNTp, myoNTm, myoNT);
myoNT->Write();
myoNT->Close();

}
