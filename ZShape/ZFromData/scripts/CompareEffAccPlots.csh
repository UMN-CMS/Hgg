#!/bin/tcsh -f
########################################
# Input and outpur directory arguments
#
#
########################################

set outdir=$1
set filetype=$2
set deffile="/home/grad/haupt/CMSSW/ZShape/CMSSW_3_6_5/src/ZShape/EffAcc/data"
set mytype=$3

root -l -b <<!
#include <map>
typedef std::vector<std::string> zdef;
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(10);
.L tdrstyle.C
setTDRStyle();

zdef zdefvec;

std::cout << " Im Trying here " << std::endl;
//std::map<std::string, zdef> zdefvecvec;  
std::cout << " Im Trying here -1 " << std::endl;
std::cout << " My type " << $mytype << std::endl;

int mytype = $mytype;
std::cout << " My type " << mytype << std::endl;
//zdefvec.push_back("Tight-EB-Loose-ECAL");
//zdefvec.push_back("Tight-ECAL-Loose-EE");
//zdefvec.push_back("Tight-ECAL-Loose-ECAL"); 
//zdefvec.push_back("Tight-ECAL-HF");
//zdefvec.push_back("ECAL95-ECAL95");
//zdefvec.push_back("ECAL90-ECAL90");
//zdefvec.push_back("ECAL85-ECAL85");
//zdefvec.push_back("ECAL80-ECAL80");
//zdefvec.push_back("ECAL70-ECAL70");
//zdefvec.push_back("ECAL60-ECAL60");
//zdefvec.push_back("ECAL95-HF");
//zdefvec.push_back("ECAL90-HF");
//zdefvec.push_back("ECAL85-HF");
//zdefvec.push_back("ECAL80-HF");
//zdefvec.push_back("ECAL70-HF");
//zdefvec.push_back("ECAL60-HF");
//zdefvec.push_back("ECAL95-HFTight");
//zdefvec.push_back("ECAL90-HFTight");
//zdefvec.push_back("ECAL85-HFTight");
//zdefvec.push_back("ECAL80-HFTight");
//zdefvec.push_back("ECAL70-HFTight");
//zdefvec.push_back("ECAL60-HFTight");
//zdefvec.push_back("Golden-EB-EB");
//zdefvec.push_back("ECAL80ISO-HFTight");
//zdefvec.push_back("ECAL80EID-HFTight");
//zdefvec.push_back("ECAL80ISOCONV-HFTight");
//zdefvec.push_back("ECAL80EIDCONV-HFTight");
//zdefvec.push_back("ECAL95EB-ECAL95EB");
//zdefvec.push_back("ECAL95EE-ECAL95EE");
//zdefvec.push_back("ECAL95EB-ECAL95EE");
//zdefvec.push_back("ECAL80-HFTight-20");
//zdefvec.push_back("ECAL95-ECAL95-20");
std::cout << " My type again " << mytype << std::endl;

if ( mytype == 0 ) {std::cout << " Im Trying here 0 " << std::endl;zdefvec.push_back("ECAL80-HF");zdefvec.push_back("ECAL95-ECAL95");}

if ( mytype == 1  ) {std::cout << " Im Trying here 10 " << std::endl;std::cout << " in 1 " << std::endl;zdefvec.push_back("ECAL95-ECAL95");}  

if ( mytype == 2 ) {zdefvec.push_back("ECAL80-ECAL80");}
if ( mytype == 3 ){zdefvec.push_back("ECAL95EB-ECAL95EB");zdefvec.push_back("ECAL95EE-ECAL95EE");zdefvec.push_back("ECAL95EB-ECAL95EE");} 

if ( $mytype == 4 ){zdefvec.push_back("ECAL80EB-ECAL80EB");zdefvec.push_back("ECAL80EE-ECAL80EE");zdefvec.push_back("ECAL80EB-ECAL80EE");}

if ( $mytype == 14 ){zdefvec.push_back("ECAL80EB-ECAL95EB");zdefvec.push_back("ECAL80EE-ECAL95EE");zdefvec.push_back("ECAL80EB-ECAL80EE");zdefvec.push_back("ECAL80-HF");}

if ( $mytype == 114 ){zdefvec.push_back("ECAL80EB-ECAL95EB");zdefvec.push_back("ECAL80EE-ECAL95EE");zdefvec.push_back("ECAL80EB-ECAL80EE");}

if ( $mytype == 5 ){zdefvec.push_back("ECAL80-ECAL95");zdefvec.push_back("ECAL80-HFTight-20");}

if ( $mytype == 6 ){zdefvec.push_back("ECAL80-ECAL95");zdefvec.push_back("ECAL80-HF");}

if ( $mytype == 16 ){zdefvec.push_back("ECAL80-ECAL95");zdefvec.push_back("ECAL80-HF");zdefvec.push_back("ECAL80-NTTight");}

if ( $mytype == 161 ){zdefvec.push_back("ECAL80EB-ECAL95EB");zdefvec.push_back("ECAL80EE-ECAL95EE");zdefvec.push_back("ECAL80EB-ECAL80EE");zdefvec.push_back("ECAL80-HF");zdefvec.push_back("ECAL80-NTTight");}

if ( $mytype == 160 ){zdefvec.push_back("ECAL80-NTTight");}



if ( $mytype == 7 ){defvec.push_back("ECAL95-ECAL95");zdefvec.push_back("ECAL85-HFTight");}

if ( $mytype == 8 ){zdefvec.push_back("ECAL95-ECAL95");zdefvec.push_back("ECAL85-HF");}

if ( $mytype == 9 ){zdefvec.push_back("ECAL95-ECAL95");zdefvec.push_back("ECAL90-HFTight");}

if ( $mytype == 10 ){zdefvec.push_back("ECAL95-ECAL95");zdefvec.push_back("ECAL90-HF");}

if ( $mytype == 11 ){zdefvec.push_back("ECAL95-ECAL95");zdefvec.push_back("ECAL95-HFTight");}

if ( $mytype == 12 ){zdefvec.push_back("ECAL95-ECAL95");zdefvec.push_back("ECAL95-HF");}

if ( $mytype == 13 ){zdefvec.push_back("ECAL95-HF");zdefvec.push_back("ECAL90-HF");zdefvec.push_back("ECAL85-HF");zdefvec.push_back("ECAL80-HF");zdefvec.push_back("ECAL95-HFTight");zdefvec.push_back("ECAL90-HFTight");zdefvec.push_back("ECAL85-HFTight");zdefvec.push_back("ECAL80-HFTight");}


std::cout << " DID I GET HERE " << std::endl;

gStyle->SetTitleAlign(22);

.L compareEffAcc.C+
//81.17
//8.0593665
//16.118733
//6.88835 for M20.... (cross fingers)...
//13.7767
//4.6675 For original Summer09
//4.88 or 5.0057
//2.4153439 for the combined
//3.9 For the full full but some must not have gone, so 4.031
//////
// 49852.5 for 1
// 498.525 for 100
// 253.1 for 197
//2932.5 for 17
//2077.2 for 24
//
//EffAccHistos *myEffAccHistos = new EffAccHistos("/data/whybee0c/user/haupt/Electrons/TNPTREE10/ZeeFromDataJeremyTEMPLATE/histoJason10_only1hlt_TRIAL_ZFULL.root","/data/whybee0c/user/haupt/Electrons/TNPTREE10/ZeeMCFromDataJeremyTEMPLATE/ZFD.Dec7thNEW.root",zdefvec, 108.1259);
EffAccHistos *myEffAccHistos = new EffAccHistos("/local/cms/user/haupt/Electrons/TNPTREE10/MYMCZFD.root","/local/cms/user/haupt/Electrons/TNPTREE10/testProdZEfficiencySmeared_10MTEMPLATE/MCHist.root",zdefvec, 51.48); //1935 //1650 //141.0 for OLD Z2TUNE //380 for 10 mil powheg //1785 for 100 mil
//2960 for the 4 inv pb of Nov 2010
//WAS 3724
//2208 for 4.5
//ZMC 1930 for 5.3
//Z2Tune 2086 for 15 
//2851 for 3 in-pb
myEffAccHistos->getFiles();
myEffAccHistos->getHistos();
myEffAccHistos->printSumHistos("${filetype}");
myEffAccHistos->printIndividualHistos("${filetype}");
myEffAccHistos->MakeResultFiles("${deffile}");

!

##echo "now we are making the webpages"
if (! -d $outdir) then
mkdir ${outdir}
endif

if (! -d $outdir/Compare ) then
mkdir ${outdir}/Compare
endif

if (! -d $outdir/FinalCompare ) then
mkdir ${outdir}/FinalCompare
endif

cp Z* ${outdir}/.

mv ${outdir}/Z_Compare* ${outdir}/Compare/.
./makeEffAccWebPage.csh ${outdir} ${filetype}


echo "That's all folks"
rm Z*${filetype}

##End of Efficiency*Acceptance Script
