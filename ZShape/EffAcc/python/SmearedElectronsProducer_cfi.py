# produces cones around electrons type 1
# reconstructing electrons and provides smearing

import FWCore.ParameterSet.Config as cms

SmearedElectronsProducer = cms.EDProducer("ZSmearingProducer",
                                          zsrc = cms.untracked.InputTag("ZIntoElectronsEventProducer","ZEventParticles"),
                                          doSmearing = cms.untracked.bool(True),
                                          quiet = cms.untracked.bool(True),
                                          zElectronsCone = cms.double(0.08),
                                          randomSeed = cms.int32(134),
                                          HF = cms.PSet(
                                             stocastic = cms.double(2.765),
                                             constant = cms.double(0.056)
                                          )
)
