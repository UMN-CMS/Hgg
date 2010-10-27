import FWCore.ParameterSet.Config as cms

process = cms.Process("Zefficiency")
process.TimerService = cms.Service("TimerService")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
    )

process.dumpEv = cms.EDAnalyzer("EventContentAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("" )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_10M_partBUILDINGTTEST.root')
)

process.load("ZShape.EffAcc.SmearedElectronsProducer_cfi")

process.mcEff = cms.EDAnalyzer("ZEfficiencyCalculator",
    Effs = cms.untracked.VPSet(
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/HLT-Pt.txt'),
            name = cms.untracked.string('HLT-EtaDet')
        ), 
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/HLT-Pt.txt'),
            name = cms.untracked.string('HLT-GSF')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/ElectronId-Pt.txt'),
            name = cms.untracked.string('ElectronId-EtaDet')
        ), 
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/Iso-PtEta.txt'),
            name = cms.untracked.string('Iso-Pt')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/ElectronId-Pt.txt'),
            name = cms.untracked.string('EID95')
        ), 
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/Iso-PtEta.txt'),
            name = cms.untracked.string('ISO95')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/ElectronId80-Pt.txt'),
            name = cms.untracked.string('EID80')
        ), 
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/Iso80-PtEta.txt'),
            name = cms.untracked.string('ISO80')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/GsfTrack-Pt.txt'),
            name = cms.untracked.string('GsfTrack-EtaDet')
        ), 
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/Supercluster-Eta.txt'),
            name = cms.untracked.string('Supercluster-Eta')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/HFEId-Pt.txt'), 
            name = cms.untracked.string('HFElectronId-EtaDet')
        ),
         cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/HFEIdT-Pt.txt'), 
            name = cms.untracked.string('HFTightElectronId-EtaDet')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/HFSC-Et_910.txt'), 
            name = cms.untracked.string('HFSuperCluster-Et')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/NTEId-Pt.txt'),
            name = cms.untracked.string('NTLooseElectronId')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/NTEIdT-Pt.txt'),
            name = cms.untracked.string('NTTightElectronId')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/WP95-PtEta.txt'),
            name = cms.untracked.string('WP95')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/WP90-PtEta.txt'),
            name = cms.untracked.string('WP90')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/WP85-PtEta.txt'),
            name = cms.untracked.string('WP85')
        ),
        cms.PSet(
            effFile = cms.FileInPath('ZShape/EffAcc/data/WP80-PtEta.txt'),
            name = cms.untracked.string('WP80')
        )                                    
    ),
    writeHistoBeforeEndJob = cms.untracked.bool(False),
    zsrc = cms.untracked.InputTag("SmearedElectronsProducer","ZEventParticles"),
    zElectronsCollection = cms.untracked.InputTag("SmearedElectronsProducer","ZEventParticles"),
    zTreeLevelElectronsCollection = cms.untracked.InputTag("ZIntoElectronsEventProducer","ZEventEle3"), quiet = cms.untracked.bool(True),
    zElectronsCone = cms.double(0.00),
    ZDefs = cms.untracked.VPSet(cms.PSet(
        name = cms.untracked.string('Golden-EB-EB'),
        Z = cms.untracked.vstring('m(70,110)'),
        e1 = cms.untracked.vstring(
            'ACC(EB)', 
            'PT10', 
            'Supercluster-Eta', 
            'PT10', 
            'GsfTrack-EtaDet', 
            'Iso-Pt', 
            'ElectronId-EtaDet', 
            'HLT-EtaDet'),
        e2 = cms.untracked.vstring('ACC(EB)', 
            'PT10', 
            'Supercluster-Eta', 
            'PT10', 
            'GsfTrack-EtaDet', 
            'Iso-Pt', 
            'ElectronId-EtaDet') #, 
            #'HLT-EtaDet')
    ), 
        cms.PSet(
            name = cms.untracked.string('Tight-EB-Loose-ECAL'),
            Z = cms.untracked.vstring('m(70,110)'),
            e1 = cms.untracked.vstring('ACC(EB)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-EtaDet', 
                'Iso-Pt', 
                'ElectronId-EtaDet', 
                'HLT-EtaDet'),
            e2 = cms.untracked.vstring('ACC(ECAL+TRK)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-EtaDet', 
                'Iso-Pt', 
                'ElectronId-EtaDet')
        ), 
        cms.PSet(
            name = cms.untracked.string('Tight-ECAL-Loose-ECAL'),
            Z = cms.untracked.vstring('m(70,110)'),
            e1 = cms.untracked.vstring('ACC(ECAL+TRK)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-EtaDet', 
                'Iso-Pt', 
                'ElectronId-EtaDet', 
                'HLT-EtaDet'),
            e2 = cms.untracked.vstring('ACC(ECAL+TRK)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-EtaDet', 
                'Iso-Pt', 
                'ElectronId-EtaDet')
        ), 
        cms.PSet(
            name = cms.untracked.string('Tight-ECAL-HF'),
            Z = cms.untracked.vstring('m(70,110)'),
            e1 = cms.untracked.vstring('ACC(ECAL+TRK)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-EtaDet', 
                'Iso-Pt', 
                'ElectronId-EtaDet', 
                'HLT-EtaDet'),
            e2 = cms.untracked.vstring('ACC(HF)', 
                'PT10', 
	        'HFSuperCluster-Et',
                'PT20',
                'HFElectronId-EtaDet')
        ),
	cms.PSet(
            name = cms.untracked.string('AllInRange'),
            Z = cms.untracked.vstring('m(70,110)'),
            e1 = cms.untracked.vstring('PT0'),
            e2 = cms.untracked.vstring('PT0')
        ),
	cms.PSet(
            name = cms.untracked.string('AllIn'),
            Z = cms.untracked.vstring('m(35,10000)'),
            e1 = cms.untracked.vstring('PT0'),
            e2 = cms.untracked.vstring('PT0')
        ),					
	cms.PSet(
            name = cms.untracked.string('Tight-ECAL-Loose-EE'),
            e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                                       "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
            e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
                                       "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),			 
            Z = cms.untracked.vstring('m(70,110)')
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
            name = cms.untracked.string('Maximal'),
            Z = cms.untracked.vstring('m(70,110)'),
            e1 = cms.untracked.vstring('ACC(ECAL+TRK)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-Eta', 
                'Iso-Pt', 
                'ElectronId-Eta', 
                'HLT-Eta'),
            e2 = cms.untracked.vstring('ACC(ANY)', 
                'PT10')
        ))
)

process.p = cms.Path(process.SmearedElectronsProducer
                     + process.mcEff
                     #+ process.dumpEv
                     )
 
