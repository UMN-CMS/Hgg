import FWCore.ParameterSet.Config as cms

process = cms.Process("Zefficiency")
process.TimerService = cms.Service("TimerService")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("")
)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_10M_partBUILDINGTTEST.root')
)

process.load("ZShape.EffAcc.SmearedElectronsProducer_cfi")
process.load("ZShape.EffAcc.ZEfficiencyStd_cfi")

# errors used are associated to the parameter vaues of SmearedElectronsProducer_cfi.py V00-05-08 :
#a-EB-down
#process.SmearedElectronsProducer.EB.p0 = 0.0391
#a-EB-up
#process.SmearedElectronsProducer.EB.p0 = 0.0429
#a-EE-down
#process.SmearedElectronsProducer.EE.p0 = 1.407
#a-EE-up
#process.SmearedElectronsProducer.EE.p0 = 1.541
#alpha-EB-down
#process.SmearedElectronsProducer.EB.alpha = 1.359
#alpha-EB-up
#process.SmearedElectronsProducer.EB.alpha = 1.413
#alpha-EE-down
#process.SmearedElectronsProducer.EE.alpha = 1.573
#alpha-EE-up
#process.SmearedElectronsProducer.EE.alpha = 1.697
#c-EB-down
#process.SmearedElectronsProducer.EB.c = 0.0081
#c-EB-up
#process.SmearedElectronsProducer.EB.c = 0.0099
#c-HFm-down
#process.SmearedElectronsProducer.HF.constantm = 0.127
#c-HFm-up
#process.SmearedElectronsProducer.HF.constantm = 0.151
#c-HFp-down
#process.SmearedElectronsProducer.HF.constantp = 0.168
#c-HFp-up
#process.SmearedElectronsProducer.HF.constantp = 0.19


process.p = cms.Path(process.SmearedElectronsProducer
                     + process.mcEff
                     )

