import FWCore.ParameterSet.Config as cms

process = cms.Process("Calib1")

process.load("FWCore.MessageService.MessageLogger_cfi")

from PhysicsTools.PatAlgos.tools.coreTools import *
## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## global tag for data
process.GlobalTag.globaltag = 'START38_V7::All' ## needed for CMSSW_3_8_0 due to changes in the DB access for JEC ## process.GlobalTag.globaltag = cms.string('GR_R_35X_V8B::All')


from PhysicsTools.PatAlgos.tools.metTools import *
removeMCMatching(process, ['All'])

process.selectedPatElectrons.cut = 'pt > 10. & abs(eta) < 12. & trackIso < 5.1 & ecalIso < 5.0 & hcalIso < 3.4 & scSigmaIEtaIEta < 0.028 '

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(110000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(         )
)

process.calib = cms.EDFilter('HFZCalib',
                             hfClusterShapes = cms.untracked.InputTag("hfEMClusters"),
                             hfRecoEcalCandidate = cms.untracked.InputTag("hfRecoEcalCandidate"),
                             hfHits = cms.untracked.InputTag("hfreco"),
                             selectedPatElectrons = cms.untracked.string('patElectrons')
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string('HFZCalib_from_data.root')
)


process.p = cms.Path(process.patDefaultSequence*process.calib)
#process.p = cms.Path(process.patDefaultSequence*process.hfEMClusteringSequence*process.demo)