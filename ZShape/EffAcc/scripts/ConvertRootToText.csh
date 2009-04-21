#!/bin/tcsh -f
set numb=""
if ($#argv > 0) then
    echo " ha ha, we have two args "
    set numb="_${1}"
endif
echo $numb

set TP="/data/whybee0c/user/haupt/Electrons/TNPTREE"
#foreach type ( EtaDet Pt )
set type="$TP/EffBinsFile.txt"
#set typehf="$TP/EffBinsFileHFEta.txt"
set typehf="$TP/EffBinsFileHF.txt"
set typept="$TP/EffBinsFilePt.txt"
set typegs="$TP/EffBinsFilePtEta.txt"
set typegsgs="$TP/EffBinsFilePtEtaGS.txt"

efficiencyRoot2Txt $TP/test_electroneff_GsfToIso_EtaPt$numb.root Iso-Pt$numb.txt Iso $typept
efficiencyRoot2Txt $TP/test_electroneff_ScToGsf_EtaPt$numb.root GsfTrack-EtaDet$numb.txt GsfTrack $typegs
efficiencyRoot2Txt $TP/test_electroneff_IsoToId_EtaPt$numb.root ElectronId-EtaDet$numb.txt ElectronId $type
efficiencyRoot2Txt $TP/test_electroneff_IdToHLT_EtaPt$numb.root HLT-EtaDet$numb.txt HLT $type
efficiencyRoot2Txt $TP/test_electroneff_HFSCToHFId_EtaPt$numb.root HFEId-EtaDet$numb.txt HFEId $typehf
efficiencyRoot2Txt $TP/test_electroneff_MCtoSC_EtaPt$numb.root Supercluster-Eta$numb.txt SuperCluster $typegs

#end



#end of script








