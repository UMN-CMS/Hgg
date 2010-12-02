import FWCore.ParameterSet.Config as cms

process = cms.Process("Res2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
                            fileMode = cms.untracked.string('NOMERGE')
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use

    fileNames=cms.untracked.vstring(    )
                            ,
                            inputCommands = cms.untracked.vstring('keep *'
                                                                  ,'drop *_hfRecoEcalCandidate_*_RECO' # drop hfRecoEcalCandidate as remade in this process
                                                                  )
                            )



# # in this file I collect the lumi sections intervals I am interested in
# from ZShape.HFZeeVBTF.myGoodlumi_cfi import *
# # this is the equivalent of the JSON file 
# process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#     #                                                                  propmtRecoJsonAug4
#     # + propmtRecoJsonAug4
#     fullPromptRecoUpTo144114
#     )


# see explanations in : https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideGoodLumiSectionsJSONFile#cmsRun
# in 36x need tags: V01-00-02 FWCore/PythonUtilities     and      V00-04-00 PhysicsTools/PythonAnalysis

import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
#myLumis = LumiList.LumiList(filename = '/home/hep/franzoni/cmssw/CMSSW_3_6_1_patch3_Z_HF_update7/src/ZShape/HFZeeVBTF/fileLists/').getCMSSWString().split(',')
myLumis = LumiList.LumiList(filename = '/home/hep/franzoni/cmssw/CMSSW_3_6_1_patch3_Z_HF_update7/src/ZShape/HFZeeVBTF/fileLists/Cert_132440-148058_7TeV_StreamExpress_Collisions10_JSON.txt').getCMSSWString().split(',')
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)


## global tag for data
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = 'START38_V7::All' ## needed for CMSSW_3_8_0 due to changes in the DB access for JEC ## process.GlobalTag.globaltag = cms.string('GR_R_35X_V8B::All')

process.es_ascii = cms.ESSource("HcalTextCalibrations",
                                                       input = cms.VPSet(
           cms.PSet(
               object = cms.string('RespCorrs'),
               #file = cms.FileInPath('ZShape/HFZeeVBTF/test/prelim_correctionsHBHEHF.txt')
               file = cms.FileInPath('Work/HFRescaler/data/prelim_correctionsHBHEHF.txt')
                          )
                  )
                                                       )
process.prefer("es_ascii")

process.hfrecalib=cms.EDProducer('HFRescaler',
                                                                 input = cms.InputTag('hfreco'),
                                                                 invert = cms.bool(False)
                                 )
process.load("RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff")
## feed the re-calibrated HF rechits into a new clustering
process.hfEMClusters.hits=cms.InputTag("hfrecalib")
# 0 = no corrections applied; 1=corrections applie
process.hfEMClusters.correctionType=0




process.load("RecoEgamma.EgammaHFProducers.hfRecoEcalCandidate_cfi")
process.hfRecoEcalCandidate.intercept2DCut=0.3
process.hfRecoEcalCandidate.e9e25Cut      =0.94
# 0.94 is the same as default in the HF cluster producer

process.TFileService = cms.Service("TFileService",
       fileName = cms.string( )
)

# to access values of EldId cuts
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi

from ZShape.HFZeeVBTF.hfzeevbtf_cfi import *

process.IdIsoRej                     = demo.clone()
process.IdIsoRej.myName              = cms.string('HFZeeVBTF-IdIsoRej')
process.IdIsoRej.acceptedElectronIDs = cms.vint32( 7 )
process.IdIsoRej.DoLog               = cms.bool(True)
#if I don't RErerun the clustering in this job, there will be a problem
#process.IdIsoRej.hFElectrons         = cms.InputTag("hfRecoEcalCandidate","","hfzeevbtfAnalysFromSkim")

# EWK analyzer: different Z definitions with one ECAL electron and one HF electron
# ---> this is the instance to run AFTER the Wenu EWK filter

process.IdIso                     = demo.clone()
process.IdIso.myName              = cms.string('HFZeeVBTF-IdIso')
process.IdIso.acceptedElectronIDs = cms.vint32( 3, 7 )
#if I don't RErerun the clustering in this job, there will be a problem
#process.IdIso.hFElectrons         = cms.InputTag("hfRecoEcalCandidate","","hfzeevbtfAnalysFromSkim")

process.IdRej                     = demo.clone()
process.IdRej.myName              = cms.string('HFZeeVBTF-IdRej')
process.IdRej.acceptedElectronIDs = cms.vint32( 5, 7 )
#if I don't RErerun the clustering in this job, there will be a problem
#process.IdRej.hFElectrons         = cms.InputTag("hfRecoEcalCandidate","","hfzeevbtfAnalysFromSkim")

process.IsoRej                     = demo.clone()
process.IsoRej.myName              = cms.string('HFZeeVBTF-IsoRej')
process.IsoRej.acceptedElectronIDs = cms.vint32( 6, 7 )
#if I don't RErerun the clustering in this job, there will be a problem
#process.IsoRej.hFElectrons         = cms.InputTag("hfRecoEcalCandidate","","hfzeevbtfAnalysFromSkim")

process.Id                     = demo.clone()
process.Id.myName              = cms.string('HFZeeVBTF-Id')
process.Id.acceptedElectronIDs = cms.vint32( 1, 3, 5, 7 )
#if I don't RErerun the clustering in this job, there will be a problem
#process.Id.hFElectrons         = cms.InputTag("hfRecoEcalCandidate","","hfzeevbtfAnalysFromSkim")

process.Iso                     = demo.clone()
process.Iso.myName              = cms.string('HFZeeVBTF-Iso')
process.Iso.acceptedElectronIDs = cms.vint32( 2, 3, 6, 7 )
#if I don't RErerun the clustering in this job, there will be a problem
#process.Iso.hFElectrons         = cms.InputTag("hfRecoEcalCandidate","","hfzeevbtfAnalysFromSkim")

process.Rej                     = demo.clone()
process.Rej.myName              = cms.string('HFZeeVBTF-pincoPonco-Rej')
process.Rej.acceptedElectronIDs = cms.vint32( 4, 5, 6, 7 )
process.Rej.hFElectrons         = cms.InputTag("hfRecoEcalCandidate","","hfzeevbtfAnalysFromSkim")

from RecoEgamma.EgammaHFProducers.hfRecoEcalCandidate_cfi import *
process.hfRecoEcalCandidateLoose            = hfRecoEcalCandidate.clone()
process.hfRecoEcalCandidateLoose.hfclusters = cms.InputTag("hfEMClusters")
process.hfRecoEcalCandidateLoose.intercept2DCut=-100# this is to avoid JUST completely the usage of the esel
process.hfRecoEcalCandidateLoose.e9e25Cut      =0.94

process.IdIsoRejHFIsoOnly                     = demo.clone()
process.IdIsoRejHFIsoOnly.myName              = cms.string('HFZeeVBTF-IdIsoRejHFIsonly')
process.IdIsoRejHFIsoOnly.acceptedElectronIDs = cms.vint32( 7 )
# the following params are, respectively: e9e25_loose, e9e25_tight,  var2d_loose, var2d_tight,  eCOREe9_loose, eCOREe9_tight,  eSeL_loose, eSeL_tight;
process.IdIsoRejHFIsoOnly.hFselParams         = cms.vdouble(0.90, 0.94,      -9999, -9999,    0.7, 0.85,     9999, 9999)
process.IdIsoRejHFIsoOnly.DoLog               = cms.bool(True)
process.IdIsoRejHFIsoOnly.hFElectrons         = cms.InputTag("hfRecoEcalCandidateLoose","","hfzeevbtfAnalysFromSkim")

#  this module is run at the beginning of the job sequence just to count the number of events
#  the analysis has been running on
process.demoBefCuts                     = demo.clone()
process.demoBefCuts.myName              = cms.string('demoBefCuts')
process.demoBefCuts.acceptedElectronIDs = cms.vint32( 7 )
process.demoBefCuts.DoLog               = cms.bool(False)
#process.demoBefCuts.hFElectrons         = cms.InputTag("hfRecoEcalCandidateLoose","","hfzeevbtfAnalysFromSkim")
#if I don't RErerun the clustering in this job, there will be a problem
#process.demoBefCuts.hFElectrons   = cms.InputTag("hfRecoEcalCandidate","","")

from ZShape.HFZeeVBTF.elewithmet_cfi import *
#myDesiredMetCollection = "caloMET"
myDesiredMetCollection = "pfMET"
#myDesiredMetCollection = "tcMET"

process.metEleIdIsoRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIdIsoRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"",""),
#                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIdIso = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 3, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIdIso'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
#                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                  metCollectionTag = cms.InputTag(myDesiredMetCollection,"",""),
                                 )


process.metEleIdRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 5, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIdRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
#                                metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"",""),
                                 )


process.metEleIsoRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIsoRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
#                                metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"",""),
                                 )


process.metEleId = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 1, 3, 5 ,7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleId'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
#                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                  metCollectionTag = cms.InputTag(myDesiredMetCollection,"",""),
                                 )


process.metEleIso = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 2, 3, 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIso'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
#                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"",""),
                                 )


process.metEleRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 4, 5, 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
#                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"",""),
                                 )


import FWCore.Modules.printContent_cfi
process.dumpEv = FWCore.Modules.printContent_cfi.printContent.clone()

process.p = cms.Path(
    # process.EG_1e28 *
    # process.hfrecalib *
    # process.hfEMClusteringSequence * 
    process.hfRecoEcalCandidate *     # redo the candidates, with the locally defined HF id cuts
    process.hfRecoEcalCandidateLoose  # never remove this because a module will need its product (which are not originally present in RECO)
    # * process.dumpEv
    process.demoBefCuts # this one is just to count events (needed in MC!)
    * process.IdIso
    * process.IdRej
    * process.IsoRej
    * process.Id
    * process.Iso
    * process.Rej
    * process.IdIsoRej
    * process.IdIsoRejHFIsoOnly
    * process.metEleIdIsoRej
    * process.metEleIdIso
    * process.metEleIdRej
    * process.metEleIsoRej
    * process.metEleId
    * process.metEleIso
    * process.metEleRej
    )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
