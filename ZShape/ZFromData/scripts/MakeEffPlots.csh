#!/bin/tcsh -f
####################
# Input and outpur directory arguments
#
#
#####################

set outdir=$1

set MTP="/data/cms004/oradata/TP/NEW"
set isoFile=$MTP/test_electroneff_GsfToIso_EtaPt.root
set gsfFile=$MTP/test_electroneff_ScToGsf_EtaPt.root
set eidFile=$MTP/test_electroneff_IsoToId_EtaPt.root
set hltFile=$MTP/test_electroneff_IdToHLT_EtaPt.root
set hfiFile=$MTP/test_electroneff_HFSCToHFId_EtaPt.root

root -l -b <<!
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(10);
std::vector<std::string> filevec;
filevec.push_back("$gsfFile");
filevec.push_back("$isoFile");
filevec.push_back("$eidFile");
filevec.push_back("$hltFile");
filevec.push_back("$hfiFile");

std::vector<std::string> physvec;
physvec.push_back("EtaDet");
physvec.push_back("Pt");

//.L  makeEffPlots.C++
//makeeffplots(filevec)

.L makeEffClass.C++

EffHistos *myEffHistos = new EffHistos(filevec,physvec);
myEffHistos->getFiles();
myEffHistos->getHistos();
myEffHistos->printSumHistos("png");
myEffHistos->printIndividualHistos("png");

!

echo "now we are making the webpages"
if (! -d $outdir) then
mkdir ${outdir}
endif

cp *png ${outdir}/.
./makeWebPage.csh ${outdir}

echo "That's all folks"

#End of Efficiency Script
