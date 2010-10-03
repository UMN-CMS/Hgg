import FWCore.ParameterSet.Config as cms

ZIntoElectronsEventProducer = cms.EDProducer("ZEventProducer",
    deltaR = cms.double(0.2),
    minET = cms.double(0.05),
    src = cms.untracked.InputTag('generator'),
    needMatch3 = cms.untracked.bool(True)
)
