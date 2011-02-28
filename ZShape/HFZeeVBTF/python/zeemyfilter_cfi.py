import FWCore.ParameterSet.Config as cms

# to access values of EldId cuts
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi

zeemyfilter = cms.EDFilter('ZeeMyFilter',
                           acceptedElectronIDs = cms.vint32( 7 ),
                           ECALid1 = cms.string('simpleEleId80relIso'),
                           ECALid2 = cms.string('simpleEleId95relIso'),
                           minEtECAL = cms.double(20),
                           PtZ_min = cms.double(90),
                           PtZ_max = cms.double(600),
                           myName = cms.string('highPt'),
                           DoLog = cms.bool(True)
                    )
