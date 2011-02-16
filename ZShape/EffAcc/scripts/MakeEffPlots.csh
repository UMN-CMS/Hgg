#!/bin/tcsh -f
set numb=""
set mnumb=""
if ($#argv > 0) then
    echo " ha ha, we have two args "
    set numb="_${1}"
    set mnumb="${1}"
endif
echo $numb

set TP="/local/cms/user/haupt/Electrons/TNPTREE10"
set effs="/home/grad/haupt/CMSSW/ZShape/CMSSW_3_9_5_patch2/src/ZShape/EffAcc/scripts"
#foreach type ( EtaDet Pt )
set type="$TP/zzEffBinsFile.txt"
set typehf="$TP/zzEffBinsFileHF.txt"
set typehf2="$TP/zzEffBinsFileHF2.txt"
set typehf4="$TP/zzEffBinsFileHF4.txt"

set typept="$TP/zzEffBinsFilePt.txt"
set typepteta="$TP/zzEffBinsFilePtEta.txt"
set typeetapt="$TP/zzEffBinsFileEtaPt.txt"
set typent="$TP/zzEffBinsFileNT.txt"
set typesc="$TP/zzEffBinsFileSC.txt"
set typesceta="$TP/zzEffBinsFileSCEta.txt"




efficiencyTxt2Plots NOTHING $effs/Iso-Pt$numb.txt Isolation$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/Iso-Eta$numb.txt IsolationEta$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/Iso-PtEta$numb.txt IsolationPtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/Iso-PtEta$numb.txt IsolationEtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/GsfTrack-EtaDet$numb.txt GsfTrack$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/GsfTrack-Pt$numb.txt GsfTrackPt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/GsfTrack-PtEta$numb.txt GsfTrackPtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/GsfTrack-PtEta$numb.txt GsfTrackEtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/ElectronId-EtaDet$numb.txt ElectronId$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/ElectronId-Pt$numb.txt ElectronIdPt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/ElectronId-PtEta$numb.txt ElectronIdPtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/ElectronId-PtEta$numb.txt ElectronIdEtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/HLT-EtaDet$numb.txt HLT$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/HLT-Pt$numb.txt HLTPt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/HLT-PtEta$numb.txt HLTPtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/HLT-PtEta$numb.txt HLTEtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/HLT80-EtaDet$numb.txt HLT80$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/HLT80-Pt$numb.txt HLT80Pt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/HLT80-PtEta$numb.txt HLT80PtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/HLT80-PtEta$numb.txt HLT80EtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/HFEId-Pt$numb.txt HFEId$numb 3 $typehf2 0
efficiencyTxt2Plots NOTHING $effs/HFEId-PtEta$numb.txt HFEIdEtaPt$numb 3 $typehf4 0
efficiencyTxt2Plots NOTHING $effs/HFEIdT-Pt$numb.txt HFEIdT$numb 3 $typehf2 0
efficiencyTxt2Plots NOTHING $effs/Iso80-Pt$numb.txt Isolation80$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/Iso80-Eta$numb.txt IsolationEta80$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/Iso80-PtEta$numb.txt IsolationPtEta80$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/Iso80-PtEta$numb.txt IsolationEtaPt80$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/ElectronId80-EtaDet$numb.txt ElectronId80$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/ElectronId80-Pt$numb.txt ElectronId80Pt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/ElectronId80-PtEta$numb.txt ElectronId80PtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/ElectronId80-PtEta$numb.txt ElectronId80EtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/WP80-Pt$numb.txt WP80Pt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/WP80-Eta$numb.txt WP80Eta$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/WP80-PtEta$numb.txt WP80PtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/WP80-PtEta$numb.txt WP80EtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/WP95-Pt$numb.txt WP95Pt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/WP95-Eta$numb.txt WP95Eta$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/WP95-PtEta$numb.txt WP95PtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/WP95-PtEta$numb.txt WP95EtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/WP90-Pt$numb.txt WP90Pt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/WP90-Eta$numb.txt WP90Eta$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/WP90-PtEta$numb.txt WP90PtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/WP90-PtEta$numb.txt WP90EtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/WP85-Pt$numb.txt WP85Pt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/WP85-Eta$numb.txt WP85Eta$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/WP85-PtEta$numb.txt WP85PtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/WP85-PtEta$numb.txt WP85EtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/WP70-Pt$numb.txt WP70Pt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/WP70-Eta$numb.txt WP70Eta$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/WP70-PtEta$numb.txt WP70PtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/WP70-PtEta$numb.txt WP70EtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/WP60-Pt$numb.txt WP60Pt$numb 3 $typept 0
efficiencyTxt2Plots NOTHING $effs/WP60-Eta$numb.txt WP60Eta$numb 3 $type 1
efficiencyTxt2Plots NOTHING $effs/WP60-PtEta$numb.txt WP60PtEta$numb 3 $typepteta 1
efficiencyTxt2Plots NOTHING $effs/WP60-PtEta$numb.txt WP60EtaPt$numb 3 $typeetapt 0
efficiencyTxt2Plots NOTHING $effs/NTEId-EtaDet$numb.txt NTEId$numb 3 $typent 1
efficiencyTxt2Plots NOTHING $effs/NTEIdT-EtaDet$numb.txt NTEIdT$numb 3 $typent 1

##These are only needed if doing the MC stuff
efficiencyTxt2Plots NOTHING $effs/SC-EtaDet$numb.txt SCEta$numb 3 $typesceta 1
efficiencyTxt2Plots NOTHING $effs/SC-PtEta$numb.txt SCPtEta$numb 3 $typesc 1



#efficiencyTxt2Plots NOTHING $effs/Supercluster-Eta$numb.txt SuperCluster$numb 3 $type 0


#efficiencyRoot2Txt $TP/test_electroneff_GsfToIso_EtaPt$numb.root Iso-Pt$numb.txt Iso $typept
#efficiencyRoot2Txt $TP/test_electroneff_ScToGsf_EtaPt$numb.root GsfTrack-EtaDet$numb.txt GsfTrack $typegs
#efficiencyRoot2Txt $TP/test_electroneff_IsoToId_EtaPt$numb.root ElectronId-EtaDet$numb.txt ElectronId $type
#efficiencyRoot2Txt $TP/test_electroneff_IdToHLT_EtaPt$numb.root HLT-EtaDet$numb.txt HLT $type
#efficiencyRoot2Txt $TP/test_electroneff_HFSCToHFId_EtaPt$numb.root HFEId-EtaDet$numb.txt HFEId $typehf
#efficiencyRoot2Txt $TP/test_electroneff_MCtoSC_EtaPt$numb.root Supercluster-Eta$numb.txt SuperCluster $typegs

#end



#end of script








