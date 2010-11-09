#!/bin/tcsh -f
set numb=""
if ($#argv > 0) then
    echo " ha ha, we have two args "
    set numb="_${1}"
endif
echo $numb

#set TP="/home/grad/haupt/CMSSW/ZShape/CMSSW_3_8_4_patch3/src/PhysicsTools/TagAndProbe/test/CRAB"
set TP="/data/whybee1a/user/haupt/Electrons/TNPTREE10"

efficiencyRoot2Txt $TP/testEff_GSFToISO.root Iso-Pt$numb.txt Iso pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToISO.root Iso-Eta$numb.txt Iso-Eta eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToISO.root Iso-PtEta$numb.txt Iso-PtEta pt_eta 
efficiencyRoot2Txt $TP/testEff_GSFToISO80.root Iso80-Pt$numb.txt Iso80 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToISO80.root Iso80-Eta$numb.txt Iso80 eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToISO80.root Iso80-PtEta$numb.txt Iso80 pt_eta 
efficiencyRoot2Txt $TP/testEff_SCToGSF.root GsfTrack-EtaDet$numb.txt GsfTrack eta 20 120
efficiencyRoot2Txt $TP/testEff_SCToGSF.root GsfTrack-Pt$numb.txt GsfTrack pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_SCToGSF.root GsfTrack-PtEta$numb.txt GsfTrack pt_eta 
efficiencyRoot2Txt $TP/testEff_ISOToID.root ElectronId-EtaDet$numb.txt ElectronId eta 20 120
efficiencyRoot2Txt $TP/testEff_ISOToID.root ElectronId-Pt$numb.txt ElectronId pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_ISOToID.root ElectronId-PtEta$numb.txt ElectronId pt_eta
efficiencyRoot2Txt $TP/testEff_ISOToID80.root ElectronId80-EtaDet$numb.txt ElectronId80 eta 20 120
efficiencyRoot2Txt $TP/testEff_ISOToID80.root ElectronId80-Pt$numb.txt ElectronId80 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_ISOToID80.root ElectronId80-PtEta$numb.txt ElectronId80 pt_eta
efficiencyRoot2Txt $TP/testEff_IDToHLT.root HLT-EtaDet$numb.txt HLT eta 20 120
efficiencyRoot2Txt $TP/testEff_IDToHLT.root HLT-Pt$numb.txt HLT pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_IDToHLT.root HLT-PtEta$numb.txt HLT pt_eta
efficiencyRoot2Txt $TP/testEff_HFSCToEID.root HFEId-EtaDet$numb.txt HFEId eta 20 120
efficiencyRoot2Txt $TP/testEff_HFSCToEID.root HFEId-Pt$numb.txt HFEId pt -5 5
efficiencyRoot2Txt $TP/testEff_HFSCToEID.root HFEId-PtEta$numb.txt HFEId pt_eta
efficiencyRoot2Txt $TP/testEff_HFSCToEIDT.root HFEIdT-EtaDet$numb.txt HFEIdT eta 20 120
efficiencyRoot2Txt $TP/testEff_HFSCToEIDT.root HFEIdT-Pt$numb.txt HFEIdT pt -5 5
efficiencyRoot2Txt $TP/testEff_HFSCToEIDT.root HFEIdT-PtEta$numb.txt HFEIdT pt_eta
efficiencyRoot2Txt $TP/testEff_GSFToWP80.root WP80-Pt$numb.txt WP80 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToWP80.root WP80-Eta$numb.txt WP80 eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToWP80.root WP80-PtEta$numb.txt WP80 pt_eta
efficiencyRoot2Txt $TP/testEff_GSFToWP95.root WP95-Pt$numb.txt WP95 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToWP95.root WP95-Eta$numb.txt WP95 eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToWP95.root WP95-PtEta$numb.txt WP95 pt_eta
efficiencyRoot2Txt $TP/testEff_GSFToWP90.root WP90-Pt$numb.txt WP90 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToWP90.root WP90-Eta$numb.txt WP90 eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToWP90.root WP90-PtEta$numb.txt WP90 pt_eta
efficiencyRoot2Txt $TP/testEff_GSFToWP85.root WP85-Pt$numb.txt WP85 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToWP85.root WP85-Eta$numb.txt WP85 eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToWP85.root WP85-PtEta$numb.txt WP85 pt_eta
efficiencyRoot2Txt $TP/testEff_GSFToWP70.root WP70-Pt$numb.txt WP70 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToWP70.root WP70-Eta$numb.txt WP70 eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToWP70.root WP70-PtEta$numb.txt WP70 pt_eta
efficiencyRoot2Txt $TP/testEff_GSFToWP60.root WP60-Pt$numb.txt WP60 pt -2.5 2.5
efficiencyRoot2Txt $TP/testEff_GSFToWP60.root WP60-Eta$numb.txt WP60 eta 20 120
efficiencyRoot2Txt $TP/testEff_GSFToWP60.root WP60-PtEta$numb.txt WP60 pt_eta
efficiencyRoot2Txt $TP/testEff_SCToNTEID.root NTEId-EtaDet$numb.txt NTEId eta 20 120
efficiencyRoot2Txt $TP/testEff_SCToNTEID.root NTEId-Pt$numb.txt NTEId pt -3 3
efficiencyRoot2Txt $TP/testEff_SCToNTEID.root NTEId-PtEta$numb.txt NTEId pt_eta
efficiencyRoot2Txt $TP/testEff_SCToNTEIDT.root NTEIdT-EtaDet$numb.txt NTEIdT eta 20 120
efficiencyRoot2Txt $TP/testEff_SCToNTEIDT.root NTEIdT-Pt$numb.txt NTEIdT pt -3 3
efficiencyRoot2Txt $TP/testEff_SCToNTEIDT.root NTEIdT-PtEta$numb.txt NTEIdT pt_eta




#efficiencyRoot2Txt $TP/test_electroneff_MCtoSC_EtaPt$numb.root Supercluster-Eta$numb.txt SuperCluster $typegs

#end



#end of script








