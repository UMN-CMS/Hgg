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
         VtxSmeared = cms.untracked.uint32(123456786),
	 generator = cms.untracked.uint32(123456786)		 
     ),
     # This is to initialize the random engine of the source
     sourceSeed = cms.untracked.uint32(123456786)
 )

#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#    theSource = cms.PSet(
#        initialSeed = cms.untracked.uiHAHAnt32(123456789),
#        engineName = cms.untracked.string('HepJamesRandom')
#    ),
#    generator = cms.PSet(
#        initialSeed = cms.untracked.uiHAHAnt32(123456789),
#        engineName = cms.untracked.string('HepJamesRandom')
#    )
#    ,saveFileName = cms.untracked.string('')
# )

process.load("ZShape.EffAcc.PythiaZee7TeV_cfi")
#process.load("ZShape.EffAcc.PythiaZee10TeV_cfi")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(300000))
process.generator.pythiaPylistVerbosity = 1
process.generator.maxEventsToPrint = 1


process.load("Configuration.StandardSequences.VtxSmearedGauss_cff")
process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")
process.load("ZShape.EffAcc.ZEventProducer_cfi")

process.dumpEv = cms.EDAnalyzer("EventContentAnalyzer")

process.p1 = cms.Path(
	  process.generator
        * process.VtxSmeared
        * process.ZIntoElectronsEventProducer
     #  * process.dumpEv
)


process.SlimEvent = cms.OutputModule("PoolOutputModule",
                                     fileName = cms.untracked.string('myfile_pythia.root'),
                                     outputCommands = cms.untracked.vstring(
                                      "drop *",
                                      "keep *_ZIntoElectronsEventProducer_*_*"
   )
)

process.outpath = cms.EndPath(process.SlimEvent)

