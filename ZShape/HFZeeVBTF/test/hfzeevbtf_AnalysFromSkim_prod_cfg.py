import FWCore.ParameterSet.Config as cms

process = cms.Process("hfzeevbtfAnalysFromSkim")

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
    fileNames=cms.untracked.vstring( )
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


process.load("RecoEgamma.EgammaHFProducers.hfRecoEcalCandidate_cfi")
process.hfRecoEcalCandidate.intercept2DCut=0.3
process.hfRecoEcalCandidate.e9e25Cut      =0.94
# 0.94 is the same as default in the HF cluster producer
# 0: correction is off; 1: correction is activated
process.hfRecoEcalCandidate.correctionType=0

process.TFileService = cms.Service("TFileService",
       fileName = cms.string('')
)

# to access values of EldId cuts
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi

from ZShape.HFZeeVBTF.hfzeevbtf_cfi import *

process.IdIsoRej                     = demo.clone()
process.IdIsoRej.myName              = cms.string('HFZeeVBTF-IdIsoRej')
process.IdIsoRej.acceptedElectronIDs = cms.vint32( 7 )
process.IdIsoRej.DoLog               = cms.bool(True)


# EWK analyzer: different Z definitions with one ECAL electron and one HF electron
# ---> this is the instance to run AFTER the Wenu EWK filter
process.IdIso                     = demo.clone()
process.IdIso.myName              = cms.string('HFZeeVBTF-IdIso')
process.IdIso.acceptedElectronIDs = cms.vint32( 3, 7 )


process.IdRej                     = demo.clone()
process.IdRej.myName              = cms.string('HFZeeVBTF-IdRej')
process.IdRej.acceptedElectronIDs = cms.vint32( 5, 7 )


process.IsoRej                     = demo.clone()
process.IsoRej.myName              = cms.string('HFZeeVBTF-IsoRej')
process.IsoRej.acceptedElectronIDs = cms.vint32( 6, 7 )


process.Id                     = demo.clone()
process.Id.myName              = cms.string('HFZeeVBTF-Id')
process.Id.acceptedElectronIDs = cms.vint32( 1, 3, 5, 7 )


process.Iso                     = demo.clone()
process.Iso.myName              = cms.string('HFZeeVBTF-Iso')
process.Iso.acceptedElectronIDs = cms.vint32( 2, 3, 6, 7 )


process.Rej                     = demo.clone()
process.Rej.myName              = cms.string('HFZeeVBTF-pincoPonco-Rej')
process.Rej.acceptedElectronIDs = cms.vint32( 4, 5, 6, 7 )


from RecoEgamma.EgammaHFProducers.hfRecoEcalCandidate_cfi import *
process.hfRecoEcalCandidateLoose            = hfRecoEcalCandidate.clone()
process.hfRecoEcalCandidateLoose.hfclusters = cms.InputTag("hfEMClusters")
process.hfRecoEcalCandidateLoose.intercept2DCut=-100# this is to avoid JUST completely the usage of the esel
process.hfRecoEcalCandidateLoose.e9e25Cut      =0.94
# 0: correction is off; 1: correction is activated
process.hfRecoEcalCandidateLoose.correctionType=0

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
process.demoBefCuts.hFElectrons   = cms.InputTag("hfRecoEcalCandidate","","")

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
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIdIso = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 3, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIdIso'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIdRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 5, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIdRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIsoRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIsoRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleId = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 1, 3, 5 ,7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleId'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIso = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 2, 3, 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIso'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 4, 5, 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


import FWCore.Modules.printContent_cfi
process.dumpEv = FWCore.Modules.printContent_cfi.printContent.clone()

process.p = cms.Path(
      # process.EG_1e28 *
    process.hfRecoEcalCandidate
    * process.hfRecoEcalCandidateLoose
    * process.demoBefCuts # this one is just to count events (needed in MC!)
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
    #* process.dumpEv
    )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
