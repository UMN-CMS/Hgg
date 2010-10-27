import FWCore.ParameterSet.Config as cms

process = cms.Process("ZFROMDATA")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
#process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')
process.GlobalTag.globaltag = 'GR10_P_V10::All'





process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( )                      
)

#process.load("ZShape.ZFromData.GoodLumis_cfi")

from ZShape.ZFromData.GoodLumis_cfi import *
#from ZShape.ZFromData.GoodLumis_July16ReReco_cfi import *
#from ZShape.ZFromData.GoodLumis_June14thReReco_cfi import *

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
    GoodLumis1+GoodLumis2+
    GoodLumis3+
    GoodLumis4+
    GoodLumis5
)

#print(process.source.lumisToProcess)

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('histoJeremy09_only1hlt_TRIAL.root')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    fileMode = cms.untracked.string('NOMERGE')
)

process.TimerService = cms.Service("TimerService",
    useCPUtime = cms.untracked.bool(True)
)

process.load("ZShape.ZFromData.ZFromDataElectrons_cfi") 

#process.theHLT.hltTag = cms.untracked.InputTag("HLT_Photon15_L1R","","HLT")
#process.theHLTGsf.hltTag = cms.untracked.InputTag("HLT_Photon15_L1R","","HLT")
#process.theHLT.hltTag = cms.untracked.InputTag("HLT_Photon15_Cleaned_L1R","","HLT")
#process.theHLTGsf.hltTag = cms.untracked.InputTag("HLT_Photon15_Cleaned_L1R","","HLT")
#process.theHLT.hltTag = cms.untracked.InputTag("HLT_Ele15_SW_L1R","","HLT")
#process.theHLTGsf.hltTag = cms.untracked.InputTag("HLT_Ele15_SW_L1R","","HLT")
#process.theHLT.hltTag = cms.untracked.InputTag("HLT_Ele15_SW_CaloEleId_L1R","","HLT")
#process.theHLTGsf.hltTag = cms.untracked.InputTag("HLT_Ele15_SW_CaloEleId_L1R","","HLT")


process.ZFromData = cms.EDAnalyzer("ZFromData",
    quiet = cms.untracked.bool(False),
    src = cms.untracked.InputTag('generator'),
    outHistogramsFile = cms.untracked.string('base_Jan23rd_1.root'),
    ExtraFromDataHistos = cms.untracked.bool(False),
    doMC = cms.untracked.bool(False),                               
    WeightsFile = cms.untracked.string('none'),                               

	dRMatchCut = cms.untracked.double(0.3),
	dPtMatchCut = cms.untracked.double(0.6),

	Effs = cms.untracked.VPSet(
	  cms.PSet(
	         name = cms.untracked.FileInPath('ZShape/EffAcc/data/Supercluster-Eta.txt'),
	         effFile = cms.untracked.string('Supercluster-Eta.txt')
	#  ),
	#  cms.PSet(
	#         name = cms.untracked.string('GsfTrack-EtaDet'),
	#         effFile = cms.untracked.string('GsfTrack-EtaDet.txt')
	#  ),
	#  cms.PSet(
	#         name = cms.untracked.string('Iso-Pt'),
	#         effFile = cms.untracked.string('Iso-Pt.txt')
	#  ),
	#  cms.PSet(
	#         name = cms.untracked.string('ElectronId-EtaDet'),
	#         effFile = cms.untracked.string('ElectronId-EtaDet.txt')
	#  ),
	#  cms.PSet(
	#         name = cms.untracked.string('HLT-EtaDet'),
	#         effFile = cms.untracked.string('HLT-EtaDet.txt')
	  )),
	  
	#TagProbeProducer = cms.untracked.InputTag('tpMapGsfElectrons'),
	#TagProbeProducer = cms.untracked.InputTag('tpMapSuperClusters'),
	#TagProbeProducer = cms.untracked.InputTag('tpMapGsfAndHF'),
        TagProbeProducer = cms.untracked.InputTag('tpMapWP95AndHF'),
	GsfProducer = cms.untracked.InputTag('theGsfElectrons'),
                                   
        CutNames          = cms.untracked.vstring("Supercluster-Eta", "GsfTrack-EtaDet",   "Iso-Pt",        "ElectronId-EtaDet", "HLT-EtaDet", "HFElectronId-EtaDet", "HFSuperCluster-Et","HFTightElectronId-EtaDet","EID95","ISO95","EID90","ISO90","EID85","ISO85","EID80","ISO80","EID70","ISO70","EID60","ISO60","HLT-GSF","ISO80Only","ISO80Conv","EID80Only","EID80Conv","WP95","WP90","WP85","WP80","NTLooseElectronId-EtaDet","NTTightElectronId-EtaDet" ),
	allProbeCandTags  = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHFSuperClusters"), cms.InputTag("theHFSuperClusters"),cms.InputTag("theHFSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("NTElecLoose"),cms.InputTag("NTElecLoose")),
	passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHLT"), cms.InputTag("HFElectronID"), cms.InputTag("theHFSuperClusters"), cms.InputTag("HFElectronIDTight"), cms.InputTag("ElectronID95"), cms.InputTag("Iso95"), cms.InputTag("ElectronID90"), cms.InputTag("Iso90"), cms.InputTag("ElectronID85"), cms.InputTag("Iso85"), cms.InputTag("ElectronID80"), cms.InputTag("Iso80"), cms.InputTag("ElectronID70"), cms.InputTag("Iso70"), cms.InputTag("ElectronID60"), cms.InputTag("Iso60"), cms.InputTag("theHLTGsf"), cms.InputTag("Iso80Only"), cms.InputTag("Iso80WConv"), cms.InputTag("ElectronID80Only"), cms.InputTag("ElectronID80WConv"), cms.InputTag("ElectronID95"),cms.InputTag("ElectronID90"),cms.InputTag("ElectronID85"),cms.InputTag("ElectronID80"),cms.InputTag("NTElecLoose"),cms.InputTag("NTElecTight")),
        ExactMatch        = cms.untracked.vint32(1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),

	writeHistoBeforeEndJob = cms.untracked.bool(False),

	#Here is where the bulk of everything goes, into the Z definitions
	ZDefs = cms.untracked.VPSet(
	  cms.PSet(
	         name = cms.untracked.string('Golden-EB-EB'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT10",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT10",
                                         #"GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),	
					 "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),						 
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-EB-Loose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
										"GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet","HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
										"GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),						 
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Maximal'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ANY)","PT10"),						 
			 Z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),		 
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-VeryLoose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta"),		 
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-AlmostVeryLoose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20"),			 
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Tight-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-EEHighEta'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20","EEEdge-Eta"),			 
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
            name = cms.untracked.string('AllInRange'),
            Z = cms.untracked.vstring('m(70,110)'),
            e1 = cms.untracked.vstring('PT0'),
            e2 = cms.untracked.vstring('PT0')
        ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-EE'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
				         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),			 
			 Z = cms.untracked.vstring('m(70,110)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL95-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL95-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL95-ECAL95'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","PT20","WP95"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (True)
          ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL90-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO90","EID90", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL90-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO90","EID90", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL90-ECAL90'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO90","EID90", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO90","EID90"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (True)
          ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL85-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO85","EID85", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL85-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO85","EID85", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL85-ECAL85'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO85","EID85", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO85","EID85"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (True)
          ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL80-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL80-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL80-ECAL80'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","PT20","WP80"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (True)
          ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL70-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO70","EID70", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL70-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO70","EID70", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL70-ECAL70'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO70","EID70", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO70","EID70"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (True)
          ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL60-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO60","EID60", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL60-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO60","EID60", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('ECAL60-ECAL60'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO60","EID60", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO60","EID60"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (True)
          ),
          cms.PSet(
	         name = cms.untracked.string('ECAL80ISO-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO80Only","DummyPass", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
	         name = cms.untracked.string('ECAL80EID-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","EID80Only","DummyPass", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
	         name = cms.untracked.string('ECAL80ISOCONV-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO80Conv","DummyPass", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
	         name = cms.untracked.string('ECAL80EIDCONV-HFTight'),
	         e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","EID80Conv","DummyPass", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
                 name = cms.untracked.string('ECAL80-HFTight-15'),
                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT15",
                                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT15","HFTightElectronId-EtaDet"),
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
                 name = cms.untracked.string('ECAL80-HFTight-17'),
                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT17",
                                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
                 name = cms.untracked.string('ECAL80-HFTight-20'),
                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
                 name = cms.untracked.string('ECAL80-HFTight-25'),
                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT25",
                                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
                 name = cms.untracked.string('ECAL80-HFTight-30'),
                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT30",
                                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
                 name = cms.untracked.string('ECAL95-ECAL95-15'),
                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT15",
                                            "GsfTrack-EtaDet","ISO95","EID95", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT15",
                                            "GsfTrack-EtaDet","ISO95","EID95"),
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (True)
          ),
          cms.PSet(
                 name = cms.untracked.string('ECAL95-ECAL95-17'),
                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT17",
                                            "GsfTrack-EtaDet","ISO95","EID95", "HLT-GSF"),
                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT17",
                                            "GsfTrack-EtaDet","ISO95","EID95"),
                 Z = cms.untracked.vstring('m(70,110)'),
                 ptorder = cms.untracked.bool (True)
          ),
          cms.PSet(
                name = cms.untracked.string('ECAL95-ECAL95-20'),
                e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                           "GsfTrack-EtaDet","ISO95","EID95", "HLT-GSF"),
                e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                           "GsfTrack-EtaDet","ISO95","EID95"),
                Z = cms.untracked.vstring('m(70,110)'),
                ptorder = cms.untracked.bool (True)
          ),
          cms.PSet(
                name = cms.untracked.string('ECAL95-ECAL95-25'),
                e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT25",
                                           "GsfTrack-EtaDet","ISO95","EID95", "HLT-GSF"),
                e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT25",
                                           "GsfTrack-EtaDet","ISO95","EID95"),
                Z = cms.untracked.vstring('m(70,110)'),
                ptorder = cms.untracked.bool (True)
          ),
          cms.PSet(
                name = cms.untracked.string('ECAL95-ECAL95-30'),
                e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT30",
                                           "GsfTrack-EtaDet","ISO95","EID95", "HLT-GSF"),
                e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT30",
                                           "GsfTrack-EtaDet","ISO95","EID95"),
                Z = cms.untracked.vstring('m(70,110)'),
                ptorder = cms.untracked.bool (True)
          ),
          cms.PSet(
                name = cms.untracked.string('ECAL95EE-ECAL95EE'),
                e1 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
                                           "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
                e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
                                           "GsfTrack-EtaDet","PT20","WP95"),
                Z = cms.untracked.vstring('m(70,110)'),
                ptorder = cms.untracked.bool (True)
          ),
          cms.PSet(
                name = cms.untracked.string('ECAL95EB-ECAL95EB'),
                e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                           "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
                e2 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                           "GsfTrack-EtaDet","PT20","WP95"),
                Z = cms.untracked.vstring('m(70,110)'),
                ptorder = cms.untracked.bool (True)
          ),
          cms.PSet(
                name = cms.untracked.string('ECAL95EB-ECAL95EE'),
                e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                           "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
                e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
                                           "GsfTrack-EtaDet","PT20","WP95"),
                Z = cms.untracked.vstring('m(70,110)'),
                ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
            name = cms.untracked.string('ECAL80-NTLoose'),
            e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                       "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
            e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT10","Supercluster-Eta","PT20","NTLooseElectronId-EtaDet"),		
            Z = cms.untracked.vstring('m(70,110)'),
            ptorder = cms.untracked.bool (False)
          ),
          cms.PSet(
            name = cms.untracked.string('ECAL80-NTTight'),
            e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                       "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
            e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT10","Supercluster-Eta","PT20","NTTightElectronId-EtaDet"),		
            Z = cms.untracked.vstring('m(70,110)'),
            ptorder = cms.untracked.bool (False)
          ),
	  cms.PSet(
            name = cms.untracked.string('ECAL95-NTLoose'),
            e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                       "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
            e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT10","Supercluster-Eta","PT20","NTLooseElectronId-EtaDet"),		
            Z = cms.untracked.vstring('m(70,110)'),
            ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
            name = cms.untracked.string('ECAL95-NTTight'),
            e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                       "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
            e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT10","Supercluster-Eta","PT20","NTTightElectronId-EtaDet"),		
            Z = cms.untracked.vstring('m(70,110)'),
            ptorder = cms.untracked.bool (False))          
	)
)


#process.load("ISpy/Analyzers/ISpy_Producer_cff")
#process.p2 = cms.Path(process.iSpy_sequence)

#process.p1 = cms.Path(process.hfEMClusteringSequence * process.lepton_cands+process.TPEdm)
process.p1 = cms.Path(process.lepton_cands+process.ZFromData)

