import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.genParticles_cfi import *

#from PhysicsTools.TagAndProbe.tag_probe_electron_cfi import *

from RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff import *

from Geometry.CMSCommonData.cmsIdealGeometryXML_cff import *

from Geometry.CaloEventSetup.CaloGeometry_cff import *

from Configuration.EventContent.EventContent_cff import *

#  electron isolation  ################
#
from RecoEgamma.EgammaIsolationAlgos.eleIsoFromDepsModules_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *

hfRecoEcalCandidate.Correct = True
##hfRecoEcalCandidate.e9e25Cut = 0
hfRecoEcalCandidate.intercept2DCut = 0.32

hfSuperClusterCandidate = hfRecoEcalCandidate.clone()
hfSuperClusterCandidate.e9e25Cut = 0
hfSuperClusterCandidate.intercept2DCut = -99


#  Calculate efficiency for *SuperClusters passing as GsfElectron* 
#
#  Tag           =  isolated GsfElectron with Robust ID, passing HLT, and  
#                    within the fiducial volume of ECAL
#  Probe --> Passing Probe   =  
#  SC --> GsfElectron --> isolation --> id --> Trigger
#############################################################



##   ____                         ____ _           _            
##  / ___| _   _ _ __   ___ _ __ / ___| |_   _ ___| |_ ___ _ __ 
##  \___ \| | | | '_ \ / _ \ '__| |   | | | | / __| __/ _ \ '__|
##   ___) | |_| | |_) |  __/ |  | |___| | |_| \__ \ ||  __/ |   
##  |____/ \__,_| .__/ \___|_|   \____|_|\__,_|___/\__\___|_|   
##  
    ############################################# 
    ############  SuperClusters  ################
    ############################################# 

HybridSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
    src = cms.InputTag("correctedHybridSuperClusters","", "RECO"),
    particleType = cms.string('gamma')
)
EBSuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("HybridSuperClusters"),
    cut = cms.string('abs( eta ) < 1.4442')
)

EndcapSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
    src = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","", "RECO"),
    particleType = cms.string('gamma')
)
EESuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("EndcapSuperClusters"),
    cut = cms.string('abs( eta ) > 1.560 & abs( eta ) < 2.5')
)

allSuperClusters = cms.EDFilter("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("EBSuperClusters"), cms.InputTag("EESuperClusters"),cms.InputTag("theHFSuperClusters"))
)




theHFSuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("hfSuperClusterCandidate"),
    cut = cms.string('et > 5.0')
)


#superClusters = cms.EDFilter("SuperClusterMerger",
#   src = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters","", "RECO"),
#                       cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","", "RECO")
 #                      #cms.InputTag("theHFSuperClusters","", "RECO")
#                       ) 
#)

#superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
#   src = cms.InputTag("superClusters"),
#   particleType = cms.string('gamma'),
#)


# My final selection of superCluster candidates 
theSuperClusters = cms.EDFilter("CandViewSelector",
    #src = cms.InputTag("superClusterCands"),
    src = cms.InputTag("allSuperClusters"),
    cut = cms.string('et  > 5.0 && abs(eta)<2.5 && !(1.4442< abs(eta) <1.560)'),
    filter = cms.bool(True) #LOOK UP WHAT THIS DOES
)

sc_sequence = cms.Sequence( ( hfSuperClusterCandidate * theHFSuperClusters) * HybridSuperClusters * EBSuperClusters + EndcapSuperClusters * EESuperClusters * allSuperClusters * theSuperClusters)



##    ____      __ _____ _           _                   
##   / ___|___ / _| ____| | ___  ___| |_ _ __ ___  _ __  
##  | |  _/ __| |_|  _| | |/ _ \/ __| __| '__/ _ \| '_ \ 
##  | |_| \__ \  _| |___| |  __/ (__| |_| | | (_) | | | |
##   \____|___/_| |_____|_|\___|\___|\__|_|  \___/|_| |_|
##  
    ############################################# 
    ############  GsfElectrons  ################
    #############################################
 
theGsfElectrons = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string('(abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.560) && (( caloEnergy * sin( caloPosition.theta ) )  > 5.0)')
)


HFElectronID = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("hfRecoEcalCandidate"),
    cut = cms.string('et > 5.0')
)

theGsfHf = cms.EDFilter("CandViewMerger",
    #src = cms.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theHFSuperClusters"), cms.InputTag("theEEHFGapSuperClusters"))
    src = cms.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theHFSuperClusters"))
)


##     ___           _       _   _             
##    |_ _|___  ___ | | __ _| |_(_) ___  _ __  
##     | |/ __|/ _ \| |/ _` | __| |/ _ \| '_ \ 
##     | |\__ \ (_) | | (_| | |_| | (_) | | | |
##    |___|___/\___/|_|\__,_|\__|_|\___/|_| |_|

    
## #  Isolation  ################

theIsolation = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string(theGsfElectrons.cut.value() +
         " && (( isEB && (dr04TkSumPt<7.2) && (dr04EcalRecHitSumEt<5.7) && (dr04HcalTowerSumEt<8.1)"
         " && (deltaEtaSuperClusterTrackAtVtx<0.0071) && (deltaPhiSuperClusterTrackAtVtx<1.0)"
         " && (sigmaEtaEta<0.01))"
         " || (isEE && (dr04TkSumPt<5.1) && (dr04EcalRecHitSumEt<5.0) && (dr04HcalTowerSumEt<3.4)"
         " && (deltaEtaSuperClusterTrackAtVtx<0.0066) && (deltaPhiSuperClusterTrackAtVtx<1.0)"
         " && (sigmaEtaEta<0.028)))")
)

##    _____ _           _                     ___    _ 
##   | ____| | ___  ___| |_ _ __ ___  _ __   |_ _|__| |
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \   | |/ _` |
##   | |___| |  __/ (__| |_| | | (_) | | | |  | | (_| |
##   |_____|_|\___|\___|\__|_|  \___/|_| |_| |___\__,_|
##   

# Electron ID  ######
theId = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string(theIsolation.cut.value() +
          " && ( (isEB && sigmaIetaIeta<0.01 && deltaEtaSuperClusterTrackAtVtx<0.0071)"
          "|| (isEE && sigmaIetaIeta<0.028 && deltaEtaSuperClusterTrackAtVtx<0.0066) )")   
)


##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
##   

# Trigger  ##################
theHLT = cms.EDProducer("trgMatchedGsfElectronProducer",                     
    InputProducer = cms.InputTag("theId"),                          
    #hltTag = cms.untracked.InputTag("HLT_Ele15_SW_LooseTrackIso_L1R","","HLT"),
    hltTag = cms.untracked.InputTag("HLT_L1SingleEG8","","HLT"),    
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
)


electron_sequence = cms.Sequence(theGsfElectrons * theGsfHf * theIsolation * theId * theHLT * HFElectronID )



 


##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
#  Tag & probe selection ######
    ################################################################ 
    ############  All Tag / Probe Association Maps  ###############
    ################################################################ 
    #// Remember that tag will always be "theHLT" collection.
    #//
    #// Probe can be one of the following collections: 
    #// "theSuperClusters", "theGsfElectrons", "theIsolation", "theId".
    #// 
    #// Passing Probe can be one of the following collections: 
    #// "theGsfElectrons", "theIsolation", "theId", "theHLT".
    #//




tpMapSuperClusters = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("theSuperClusters theSuperClusters"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("30 < mass < 120"),
)


tpMapGsfElectrons = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("theGsfElectrons theGsfElectrons"), # charge coniugate states are implied
    cut   = cms.string("60 < mass < 120"),
    checkCharge = cms.bool(False),
)


tpMapIsolation = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("theHLT theIsolation"), # charge coniugate states are implied
    cut   = cms.string("60 < mass < 120"),
    checkCharge = cms.bool(False),
)

tpMapId =  cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("theHLT theId"), # charge coniugate states are implied
    cut   = cms.string("60 < mass < 120"),
    checkCharge = cms.bool(False),
)

tpMapHFSuperClusters = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("theHLT theHFSuperClusters"), # charge coniugate states are implied
    cut   = cms.string("30 < mass < 120"),
    checkCharge = cms.bool(False),
)



tpMapGsfAndHF =  cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("theHLT theGsfHf"), # charge coniugate states are implied
    cut   = cms.string("30 < mass < 120"),
    checkCharge = cms.bool(False),
)



tpMap_sequence = cms.Sequence( tpMapSuperClusters + tpMapGsfElectrons + tpMapIsolation + tpMapId + tpMapHFSuperClusters + tpMapGsfAndHF)
#tpMap_sequence = cms.Sequence( tpMapGsfElectrons + tpMapIsolation + tpMapId + tpMapHFSuperClusters + tpMapGsfAndHF)

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        


    ################################################################ 
    ############  All Truth-matched collections  ###################
    ################################################################ 
    #// find generator particles matching by DeltaR


SuperClustersMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theSuperClusters"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

GsfElectronsMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theGsfElectrons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

IsolationMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theIsolation"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

IdMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theId"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

HLTMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theHLT"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

HFSCMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theHFSuperClusters"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

HFIDMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("HFElectronID"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)


truthMatch_sequence = cms.Sequence( SuperClustersMatch + GsfElectronsMatch + IsolationMatch + IdMatch + HLTMatch + HFSCMatch + HFIDMatch)
#truthMatch_sequence = cms.Sequence( GsfElectronsMatch + IsolationMatch + IdMatch + HLTMatch + HFSCMatch + HFIDMatch)
#lepton_cands = cms.Sequence(genParticles * hfEMClusteringSequence * sc_sequence * electron_sequence * tpMap_sequence * truthMatch_sequence)

lepton_cands = cms.Sequence(hfEMClusteringSequence * sc_sequence * electron_sequence * tpMap_sequence )


   

