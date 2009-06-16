import FWCore.ParameterSet.Config as cms

process = cms.Process("RAPSKIMMER")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0005/32EB7039-DE1A-DE11-97EE-003048724749.root')
)


process.maxEvents = cms.untracked.PSet(
   #input = cms.untracked.int32(-1)
   #input = cms.untracked.int32(123200)
   input = cms.untracked.int32(22760)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TimerService = cms.Service("TimerService",
    useCPUtime = cms.untracked.bool(True)
)

process.load("ZShape.ZFromData.BackgroundSkimmer_cfi") 
	
		  
process.BackgroundSkimmer = cms.EDFilter("BackgroundSkimmer",
	allProbeCandTags = cms.untracked.VInputTag('tpMapEnd' ),
)


#process.p1 = cms.Path(process.hfEMClusteringSequence * process.lepton_cands+process.TPEdm)
process.p1 = cms.Path(process.lepton_cands+process.electron_sequence+process.BackgroundSkimmer)
		
process.outpath = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(   #'drop *', 
        'keep *'),
    fileName = cms.untracked.string('Ohkj.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    )

)

process.the_end = cms.EndPath(process.outpath)

