import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")


# initialize  MessageLogger

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.limit = 25
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.TimerService = cms.Service("TimerService",
    useCPUtime = cms.untracked.bool(True)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
     moduleSeeds = cms.PSet(
         VtxSmeared = cms.untracked.uint32(123456787),
     ),
     # This is to initialize the random engine of the source
     sourceSeed = cms.untracked.uint32(123456787)
 )



process.load("ZShape.EffAcc.PythiaZee10TeV_cfi")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(3000))
process.PythiaSource.pythiaPylistVerbosity = 1
process.PythiaSource.maxEventsToPrint = 1


process.load("Configuration.StandardSequences.VtxSmearedGauss_cff")
process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")
process.load("ZShape.EffAcc.ZEventProducer_cfi")

process.dumpEv = cms.EDAnalyzer("EventContentAnalyzer")

process.p1 = cms.Path(
          process.VtxSmeared
        * process.ZIntoElectronsEventProducer
        * process.dumpEv
)


process.SlimEvent = cms.OutputModule("PoolOutputModule",
                                     fileName = cms.untracked.string('myfile_pythia.root'),
                                     outputCommands = cms.untracked.vstring(
                                      "drop *",
                                      "keep *_ZIntoElectronsEventProducer_*_*"
   )
)

process.outpath = cms.EndPath(process.SlimEvent)

