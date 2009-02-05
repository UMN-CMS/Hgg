#!/bin/tcsh -f
########################################
# Input and outpur directory arguments
#
#
########################################

set outdir=$1

root -l -b <<!
#include <map.h>
typedef std::vector<std::string> zdef;
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(10);
zdef zdefvec;

std::cout << " I'm Trying here " << std::endl;
std::map<std::string, zdef> zdefvecvec;
//zdefvec.push_back("Tight-EB-Loose-ECAL");
zdefvec.push_back("Tight-ECAL-Loose-ECAL"); 
zdefvec.push_back("Tight-ECAL-HF");
//zdefvec.push_back("Golden-EB-EB");

.L compareEffAcc.C++
//81.17
EffAccHistos *myEffAccHistos = new EffAccHistos("/data/cms004/oradata/ZEE/histoJeremy09_only1hlt_TRIAL.root","/data/cms004/oradata/ZMC/histo_10M_part.root",zdefvec, 8.117);
myEffAccHistos->getFiles();
myEffAccHistos->getHistos();
myEffAccHistos->printSumHistos("png");
myEffAccHistos->printIndividualHistos("png");

!

##echo "now we are making the webpages"
if (! -d $outdir) then
mkdir ${outdir}
endif

if (! -d $outdir/Compare ) then
mkdir ${outdir}/Compare
endif

cp Z_Compare*png ${outdir}/Compare/.

cp Z*png ${outdir}/.
./makeEffAccWebPage.csh ${outdir}





echo "That's all folks"
rm Z*png

##End of Efficiency*Acceptance Script
