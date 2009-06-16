import FWCore.ParameterSet.Config as cms

process = cms.Process("ZFROMDATA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0005/1E2DDC37-EC1A-DE11-BB33-0030487F933D.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0005/32BEC84F-DE1A-DE11-9602-0030487EB003.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0005/32EB7039-DE1A-DE11-97EE-003048724749.root',
#'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0005/46ECF224-E81A-DE11-99A0-0030487D7B79.root',
#'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0005/8AC19B36-EC1A-DE11-99EA-0030487E4B8D.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0010/7EC42E56-841F-DE11-89EA-0030487D5E67.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0010/886441B8-831F-DE11-9D72-0030487D5E67.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0010/98754C5E-841F-DE11-979B-0030487F1665.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0010/DADEAD9A-861F-DE11-94B6-0030487F929F.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0011/4A154264-911F-DE11-A563-0030487277D1.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0012/389C8F78-971F-DE11-AF31-0030487F1A49.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0012/7C8E7674-971F-DE11-9D4E-0030487F171D.root',
'file:/data/whybee0b/phedex/store/mc/Summer08/Zee_M20/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0012/B06EC480-971F-DE11-AD3D-0030487F16BF.root'			
	)
)


process.maxEvents = cms.untracked.PSet(
   #input = cms.untracked.int32(-1)
   #input = cms.untracked.int32(123200)
   input = cms.untracked.int32(22760)
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('InclusiveCheckerHisto.root')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TimerService = cms.Service("TimerService",
    useCPUtime = cms.untracked.bool(True)
)

process.load("ZShape.ZFromData.InclusiveCheck_cfi") 
	
		  
process.ZFromData = cms.EDAnalyzer("InclusiveCheck",
    quiet = cms.untracked.bool(False),
    outHistogramsFile = cms.untracked.string('base_Jan23rd_1.root'),
    ExtraFromDataHistos = cms.untracked.bool(True),

	dRMatchCut = cms.untracked.double(0.5),
	dPtMatchCut = cms.untracked.double(1.5),
	
	#MCTag = cms.untracked.InputTag('thepreSuperClusters'),
	#SCTag = cms.untracked.InputTag('theSuperClusters'),
        MCTag = cms.untracked.InputTag('theSuperClusters'),
	SCTag = cms.untracked.InputTag('theGsfElectrons'),
	#SCTag = cms.untracked.InputTag('electrons'),
        #SCTag = cms.untracked.InputTag('pixelMatchGsfElectrons'),
		 
	writeHistoBeforeEndJob = cms.untracked.bool(False),
		 
	EffBinsFile =  cms.untracked.string('EffBinsFilePtEta.txt'),
	scEtCut = cms.untracked.double(5.0),
		 
	CBBinsFile = cms.untracked.string('CBBinsFile.txt'),
	FromFile   = cms.untracked.bool(False),
	InputFile  = cms.untracked.string('histoMCOLD.root'),
		 
	#allProbeCandTags = cms.untracked.VInputTag('theHLTInc',      'theGsfRap','theIsolation','theId','theHLT','tpMapEnd' ),
	allProbeCandTags = cms.untracked.VInputTag('theHLTInc','theGsfElectrons','theIsolation','theId','theHLT', 'tpMapEnd'),
	NumMatch         = cms.untracked.vint32(       0      ,         1       ,       1      ,    1  ,   0    ,          0)

	#Here is where the bulk of everything goes, into the Z definitions
)

process.outpath = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(   #'drop *', 
        'keep *_*_*_Zefficiency'),
    fileName = cms.untracked.string('Ohkj.root')
)

#process.p1 = cms.Path(process.hfEMClusteringSequence * process.lepton_cands+process.TPEdm)
process.p1 = cms.Path(process.lepton_cands+process.electron_sequence+process.ZFromData)
process.the_end = cms.EndPath(process.outpath)

