import FWCore.ParameterSet.Config as cms

process = cms.Process("ZFROMDATA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_1.root',
       'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_2.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_3.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_4.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_5.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_6.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_7.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_8.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_9.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_10.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_11.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_12.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_13.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_14.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_15.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_16.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_17.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_18.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_19.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_20.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_21.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_22.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_23.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_24.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_25.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_26.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_27.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_28.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_29.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_30.root'
	)
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('histo.root')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TimerService = cms.Service("TimerService",
    useCPUtime = cms.untracked.bool(True)
)

process.load("ZShape.ZFromData.ZFromDataElectrons_cfi") 

process.ZFromData = cms.EDAnalyzer("ZFromData",
    quiet = cms.untracked.bool(False),
    outHistogramsFile = cms.untracked.string('base_DAug4th.root'),
    ExtraFromDataHistos = cms.untracked.bool(True),

	dRMatchCut = cms.untracked.double(0.3),
	dPtMatchCut = cms.untracked.double(0.3),

	Effs = cms.untracked.VPSet(
	  cms.PSet(
	         name = cms.untracked.string('Supercluster-Eta'),
	         effFile = cms.untracked.string('Supercluster-Eta.txt')
	  ),
	  cms.PSet(
	         name = cms.untracked.string('GsfTrack-Eta'),
	         effFile = cms.untracked.string('GsfTrack-Eta.txt')
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Iso-Pt'),
	         effFile = cms.untracked.string('Iso-Pt.txt')
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ElectronId-Eta'),
	         effFile = cms.untracked.string('ElectronId-Eta.txt')
	  ),
	  cms.PSet(
	         name = cms.untracked.string('HLT-Eta'),
	         effFile = cms.untracked.string('HLT-Eta.txt')
	  )),
	  
	TagProbeProducer = cms.untracked.InputTag('tpMapGsfElectrons'),
	
	CutNames          = cms.untracked.vstring("Supercluster-Eta", "GsfTrack-Eta",   "Iso-Pt",        "ElectronId-Eta", "HLT-Eta"),
	allProbeCandTags  = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId")),
	passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHLT")),
    ExactMatch        = cms.untracked.vint32(0,0,0,0,0),
	
	writeHistoBeforeEndJob = cms.untracked.bool(False),

	#Here is where the bulk of everything goes, into the Z definitions
	ZDefs = cms.untracked.VPSet(
	  cms.PSet(
	         name = cms.untracked.string('Golden-EB-EB'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-Eta","Iso-Pt","ElectronId-Eta", "HLT-Eta"),
			 e2 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-Eta","Iso-Pt","ElectronId-Eta", "HLT-Eta"),						 
			 z = cms.untracked.vstring()
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-EB-Loose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
										"GsfTrack-Eta","Iso-Pt","ElectronId-Eta","HLT-Eta"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
										"GsfTrack-Eta","Iso-Pt","ElectronId-Eta"),						 
			 z = cms.untracked.vstring()
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Maximal'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-Eta","Iso-Pt","ElectronId-Eta", "HLT-Eta"),
			 e2 = cms.untracked.vstring("ACC(ANY)","PT10"),						 
			 z = cms.untracked.vstring()
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-Eta","Iso-Pt","ElectronId-Eta", "HLT-Eta"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-Eta","Iso-Pt","ElectronId-Eta"),						 
			 z = cms.untracked.vstring()
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-VeryLoose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-Eta","Iso-Pt","ElectronId-Eta", "HLT-Eta"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta"),						 
			 z = cms.untracked.vstring()
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-AlmostVeryLoose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-Eta","Iso-Pt","ElectronId-Eta", "HLT-Eta"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20"),						 
			 z = cms.untracked.vstring()
	  )
	)
)

process.outpath = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(   #'drop *', 
        'keep *_*_*_Zefficiency'),
    fileName = cms.untracked.string('OutDoisdMar12.root')
)

#process.p1 = cms.Path(process.hfEMClusteringSequence * process.lepton_cands+process.TPEdm)
process.p1 = cms.Path(process.lepton_cands+process.ZFromData)
process.the_end = cms.EndPath(process.outpath)

