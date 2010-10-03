#!/bin/tcsh -f
set numb=""
if ($#argv > 0) then
    echo " ha ha, we have two args "
    set numb="_${1}"
endif
echo $numb

set TP="/home/grad/haupt/CMSSW/ZShape/CMSSW_3_6_3_patch2/src/PhysicsTools/TagAndProbe/test/CRAB"


efficiencyRoot2Txt $TP/testEff_GSFToISO.root Iso-Pt$numb.txt Iso pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToISO.root Iso-Eta$numb.txt Iso-Eta eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToISO80.root Iso80-Pt$numb.txt Iso80 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToISO80.root Iso80-Eta$numb.txt Iso80 eta 20 120
efficiencyRoot2Txt $TP/testEff_SCToGSF.root GsfTrack-EtaDet$numb.txt GsfTrack eta 20 120
efficiencyRoot2Txt $TP/testEff_SCToGSF.root GsfTrack-Pt$numb.txt GsfTrack pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_ISOToID.root ElectronId-EtaDet$numb.txt ElectronId eta 20 120
efficiencyRoot2Txt $TP/testEff_ISOToID.root ElectronId-Pt$numb.txt ElectronId pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_ISOToID80.root ElectronId80-EtaDet$numb.txt ElectronId80 eta 20 120
efficiencyRoot2Txt $TP/testEff_ISOToID80.root ElectronId80-Pt$numb.txt ElectronId80 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_IDToHLT.root HLT-EtaDet$numb.txt HLT eta 20 120
efficiencyRoot2Txt $TP/testEff_IDToHLT.root HLT-Pt$numb.txt HLT pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_HFSCToEID.root HFEId-EtaDet$numb.txt HFEId eta 20 120
efficiencyRoot2Txt $TP/testEff_HFSCToEID.root HFEId-Pt$numb.txt HFEId pt -5 5

#efficiencyRoot2Txt $TP/test_electroneff_MCtoSC_EtaPt$numb.root Supercluster-Eta$numb.txt SuperCluster $typegs

#end



#end of script








