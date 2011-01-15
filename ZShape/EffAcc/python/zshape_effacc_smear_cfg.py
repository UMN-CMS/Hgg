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
    fileNames = cms.untracked.vstring("" )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_10M_partBUILDINGTTEST.root')
)

process.load("ZShape.EffAcc.SmearedElectronsProducer_cfi")
process.load("ZShape.EffAcc.ZEfficiencyStd_cfi")

# errors used are associated to the parameter vaues of SmearedElectronsProducer_cfi.py V00-05-03 :
#a-EB-down
#process.SmearedElectronsProducer.EB.p0 = 0.04143
#a-EB-up
#process.SmearedElectronsProducer.EB.p0 = 0.04197
#a-EE-down
#process.SmearedElectronsProducer.EE.p0 = 1.623
#a-EE-up
#process.SmearedElectronsProducer.EE.p0 = 1.629
#a-HF-down
#process.SmearedElectronsProducer.HF.stocastic = 3.432
#a-HF-up
#process.SmearedElectronsProducer.HF.stocastic = 3.48
#alpha-EB-down
#process.SmearedElectronsProducer.EB.alpha = 1.37853
#alpha-EB-up
#process.SmearedElectronsProducer.EB.alpha = 1.37947
#alpha-EE-down
#process.SmearedElectronsProducer.EE.alpha = 1.607
#alpha-EE-up
#process.SmearedElectronsProducer.EE.alpha = 1.623
#c-EB-down
#process.SmearedElectronsProducer.EB.c = 0.0088
#c-EB-up
#process.SmearedElectronsProducer.EB.c = 0.0092
#c-HF-down
#process.SmearedElectronsProducer.HF.constant = 0.086
#c-HF-up
#process.SmearedElectronsProducer.HF.constant = 0.088


# ad-hoc assumption of a +-5% variation on each parameter 
#a-EB-Vdown
#process.SmearedElectronsProducer.EB.p0 = 0.0437332
#a-EB-Vup
#process.SmearedElectronsProducer.EB.p0 = 0.0395682
#a-EE-Vdown
#process.SmearedElectronsProducer.EE.p0 = 1.6224
#a-EE-Vup
#process.SmearedElectronsProducer.EE.p0 = 1.70779
#a-HF-Vdown
#process.SmearedElectronsProducer.HF.stocastic = 3.44762
#a-HF-Vup
#process.SmearedElectronsProducer.HF.stocastic = 3.62907
#alpha-EB-Vdown
#process.SmearedElectronsProducer.EB.alpha = 1.31005
#alpha-EB-Vup
#process.SmearedElectronsProducer.EB.alpha = 1.44795
#alpha-EE-Vdown
#process.SmearedElectronsProducer.EE.alpha = 1.53452
#alpha-EE-Vup
#process.SmearedElectronsProducer.EE.alpha = 1.69604
#c-EB-Vdown
#process.SmearedElectronsProducer.EB.c = 0.00855
#c-EB-Vup
#process.SmearedElectronsProducer.EB.c = 0.00945
#c-HF-Vdown
#process.SmearedElectronsProducer.HF.constant = 0.0821845
#c-HF-Vup
#process.SmearedElectronsProducer.HF.constant = 0.0908355



process.p = cms.Path(process.SmearedElectronsProducer
                     + process.mcEff
                     )
