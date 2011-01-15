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
       stocastic = cms.double(3.5),
       constant = cms.double(0.085),
       mean = cms.double(1.04346)
        #tuned by kevin on 04-16-09, mean added by J. Mans 05-22-09
    ),
    EB = cms.PSet(
       #tuned by j. Mans on 06-08-09
       p0 = cms.double(0.0416507) , 
       p1 = cms.double(1.06),
       p2 = cms.double(4.78),
       c = cms.double(0.009), #0.004
       alpha = cms.double(1.379),
       mean = cms.double(0.999956), #back to 1.0 from 0.983
       n = cms.double(2.0) 	   
    ),
     EE = cms.PSet(
       p0 = cms.double(1.62647),
       p1 = cms.double(0.72),
       p2 = cms.double(0.15), #0.15
       alpha = cms.double(1.61528),
       mean = cms.double(0.977413),
       n = cms.double(2.0) 	
    )
                                          

)

