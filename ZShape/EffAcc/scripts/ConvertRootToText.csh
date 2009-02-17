#!/bin/tcsh -f

set TP="/data/cms004/oradata/TP/NEW"
#foreach type ( EtaDet Pt )
set type="EffBinsFile.txt"
set typehf="EffBinsFileHF.txt"

efficiencyRoot2Txt $TP/test_electroneff_GsfToIso_EtaPt.root Iso-Pt.txt Iso $type
efficiencyRoot2Txt $TP/test_electroneff_ScToGsf_EtaPt.root GsfTrack-EtaDet.txt GsfTrack $type
efficiencyRoot2Txt $TP/test_electroneff_IsoToId_EtaPt.root ElectronId-EtaDet.txt ElectronId $type
efficiencyRoot2Txt $TP/test_electroneff_IdToHLT_EtaPt.root HLT-EtaDet.txt HLT $type
efficiencyRoot2Txt $TP/test_electroneff_HFSCToHFId_EtaPt.root HFEId-EtaDet.txt HFEId $typehf

#end



#end of script








