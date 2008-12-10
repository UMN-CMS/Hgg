import FWCore.ParameterSet.Config as cms

process = cms.Process("Zefficiency")
process.TimerService = cms.Service("TimerService")

process.load('FWCore/MessageService/MessageLogger_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(10000)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000000)
    )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_1.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_2.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_3.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_4.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_5.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_6.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_7.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_8.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_9.root',
    'file:/data/whybee0c/user/jmmans/Zee_shape_10TeV_GEN/Zee_shape_10TeV_GEN_10.root'
                                      )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_1M.root')
)

process.mcEff = cms.EDFilter("ZEfficiencyCalculator",
    Effs = cms.untracked.VPSet(cms.PSet(
        effFile = cms.untracked.string('HLT-EtaDet.txt'),
        name = cms.untracked.string('HLT-Eta')
    ), 
        cms.PSet(
            effFile = cms.untracked.string('ElectronId-EtaDet.txt'),
            name = cms.untracked.string('ElectronId-Eta')
        ), 
        cms.PSet(
            effFile = cms.untracked.string('Iso-Pt.txt'),
            name = cms.untracked.string('Iso-Pt')
        ), 
        cms.PSet(
            effFile = cms.untracked.string('GsfTrack-EtaDet.txt'),
            name = cms.untracked.string('GsfTrack-Eta')
        ), 
        cms.PSet(
            effFile = cms.untracked.string('Supercluster-Eta.txt'),
            name = cms.untracked.string('Supercluster-Eta')
        ),
        cms.PSet(
            effFile = cms.untracked.string('HFEId-EtaDet.txt'),
            name = cms.untracked.string('HFElectronId-Eta')
        )

                               ),
    writeHistoBeforeEndJob = cms.untracked.bool(False),
    zsrc = cms.untracked.InputTag("ZIntoElectronsEventProducer","ZEventParticles"),
    quiet = cms.untracked.bool(True),
    zElectronsCone = cms.double(0.08),
    ZDefs = cms.untracked.VPSet(cms.PSet(
        name = cms.untracked.string('Golden-EB-EB'),
        Z = cms.untracked.vstring(),
        e1 = cms.untracked.vstring(
            'ACC(EB)', 
            'PT10', 
            'Supercluster-Eta', 
            'PT20', 
            'GsfTrack-Eta', 
            'Iso-Pt', 
            'ElectronId-Eta', 
            'HLT-Eta'),
        e2 = cms.untracked.vstring('ACC(EB)', 
            'PT10', 
            'Supercluster-Eta', 
            'PT20', 
            'GsfTrack-Eta', 
            'Iso-Pt', 
            'ElectronId-Eta', 
            'HLT-Eta')
    ), 
        cms.PSet(
            name = cms.untracked.string('Tight-EB-Loose-ECAL'),
            Z = cms.untracked.vstring(),
            e1 = cms.untracked.vstring('ACC(EB)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-Eta', 
                'Iso-Pt', 
                'ElectronId-Eta', 
                'HLT-Eta'),
            e2 = cms.untracked.vstring('ACC(ECAL)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-Eta', 
                'Iso-Pt', 
                'ElectronId-Eta')
        ), 
        cms.PSet(
            name = cms.untracked.string('Tight-ECAL-Loose-ECAL'),
            Z = cms.untracked.vstring(),
            e1 = cms.untracked.vstring('ACC(ECAL+TRK)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-Eta', 
                'Iso-Pt', 
                'ElectronId-Eta', 
                'HLT-Eta'),
            e2 = cms.untracked.vstring('ACC(ECAL)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-Eta', 
                'Iso-Pt', 
                'ElectronId-Eta')
        ), 
        cms.PSet(
            name = cms.untracked.string('Tight-ECAL-HF'),
            Z = cms.untracked.vstring(),
            e1 = cms.untracked.vstring('ACC(ECAL+TRK)', 
                'PT10', 
                'Supercluster-Eta', 
                'PT20', 
                'GsfTrack-Eta', 
                'Iso-Pt', 
                'ElectronId-Eta', 
                'HLT-Eta'),
            e2 = cms.untracked.vstring('ACC(HF)', 
                'PT10', 
                '', 
                'PT20', 
                '',
                '',
                'HF-ElectronId-Eta')
        ), 
        cms.PSet(
            name = cms.untracked.string('Acceptance-Tracker'),
            Z = cms.untracked.vstring(),
            e1 = cms.untracked.vstring('ACC(ECAL+TRK)', 
                'PT10', 
                'PT20','PT20'),
            e2 = cms.untracked.vstring('ACC(ANY-EEnotrk)', 
                'PT10', 
                'PT10','PT20')
        ), 
        cms.PSet(
            name = cms.untracked.string('Maximal'),
            Z = cms.untracked.vstring(),
            e1 = cms.untracked.vstring('ACC(ECAL)', 
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

process.p = cms.Path(process.mcEff)


