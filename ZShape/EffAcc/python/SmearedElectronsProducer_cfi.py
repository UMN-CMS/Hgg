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
       stocastic = cms.double(2.61),
       constant = cms.double(0.050),
       mean = cms.double(0.997)
        #tuned by kevin on 04-16-09, mean added by J. Mans 05-22-09
    ),
    EB = cms.PSet(
       #tuned by j. Mans on 06-08-09
       p0 = cms.double(0.040),
       p1 = cms.double(1.06),
       p2 = cms.double(4.78),
       c = cms.double(0.004),
       alpha = cms.double(1.93),
       mean = cms.double(0.995),
       n = cms.double(3.0) 	   
    ),
     EE = cms.PSet(
       p0 = cms.double(1.1),
       p1 = cms.double(0.72),
       p2 = cms.double(0.15),
       alpha = cms.double(1.7),
       mean = cms.double(0.994),
       n = cms.double(3.0) 	
    )
                                          

)
