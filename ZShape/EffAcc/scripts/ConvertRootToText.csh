#!/bin/tcsh -f

set TP="/data/cms004/oradata/TP/NEW"
foreach type ( EtaDet Pt )

efficiencyRoot2Txt $TP/test_electroneff_GsfToIso_EtaPt.root Iso-${type}.txt Iso $type
efficiencyRoot2Txt $TP/test_electroneff_ScToGsf_EtaPt.root GsfTrack-${type}.txt GsfTrack $type
efficiencyRoot2Txt $TP/test_electroneff_IsoToId_EtaPt.root ElectronId-${type}.txt ElectronId $type
efficiencyRoot2Txt $TP/test_electroneff_IdToHLT_EtaPt.root HLT-${type}.txt HLT $type
efficiencyRoot2Txt $TP/test_electroneff_HFSCToHFId_EtaPt.root HFEId-${type}.txt HFEId $type

end



#end of script








