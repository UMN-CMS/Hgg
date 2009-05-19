# produces cones around electrons type 1
# reconstructing electrons and provides smearing

import FWCore.ParameterSet.Config as cms

SmearedElectronsProducer = cms.EDProducer("ZSmearingProducer",
    zsrc = cms.untracked.InputTag("ZIntoElectronsEventProducer","ZEventParticles"),
    doSmearing = cms.untracked.bool(True),
    quiet = cms.untracked.bool(True),
    zElectronsCone = cms.double(0.08),
    randomSeed = cms.int32(134),
    SmearTable = cms.FileInPath('ZShape/EffAcc/data/CBBinsFile.txt'),
    HF = cms.PSet(
       stocastic = cms.double(2.8554),
       constant = cms.double(0.054798)
        #tuned by kevin on 04-16-09
    ),
    EB = cms.PSet(
       p0 = cms.double(0.166),
       p1 = cms.double(-0.176),
       p2 = cms.double(0.793),
       c = cms.double(0.011),
       alpha = cms.double(1.1),
       mean = cms.double(1.00),
       n = cms.double(3.0) 	   
    ),
     EE = cms.PSet(
       p0 = cms.double(6.29),
       p1 = cms.double(-4.53),
       p2 = cms.double(0.95),
       alpha = cms.double(1.2),
       mean = cms.double(1.013),
       n = cms.double(3.0) 	
    )
                                          

)
