import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.genParticles_cfi import *

#from PhysicsTools.TagAndProbe.tag_probe_electron_cfi import *

from RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff import *

from Geometry.CMSCommonData.cmsIdealGeometryXML_cff import *

from Geometry.CaloEventSetup.CaloGeometry_cff import *

from Configuration.EventContent.EventContent_cff import *

hfRecoEcalCandidate.Correct = True
hfRecoEcalCandidate.e9e25Cut = 0
hfRecoEcalCandidate.intercept2DCut = 0.32

hfSuperClusterCandidate = hfRecoEcalCandidate.clone()
hfSuperClusterCandidate.e9e25Cut = 0
hfSuperClusterCandidate.intercept2DCut = -99


    ############################################# 
    ############  GsfElectrons  ################
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


theHFSuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("hfSuperClusterCandidate"),
    cut = cms.string('et > 15.0')
)


allSuperClusters = cms.EDFilter("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("EBSuperClusters"), cms.InputTag("EESuperClusters"),cms.InputTag("theHFSuperClusters"))
)


theSuperClusters = cms.EDFilter("CandViewSelector",
    #src = cms.InputTag("superClusterCands"),
    src = cms.InputTag("allSuperClusters"),
    cut = cms.string('et  > 15.0 && abs(eta)<2.8 && !(1.4442< abs(eta) <1.560)'),
    filter = cms.bool(True) #LOOK UP WHAT THIS DOES
)

####NOW I DO THE MAP
tpMapSuperClusters = cms.EDFilter("CandViewShallowCloneCombiner",
    decay = cms.string("theSuperClusters allSuperClusters"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("0 < mass < 140"),
    filter = cms.bool(True) #LOOK UP WHAT THIS DOES
)




sc_sequence = cms.Sequence( ( hfSuperClusterCandidate * theHFSuperClusters) * HybridSuperClusters * EBSuperClusters + EndcapSuperClusters * EESuperClusters * allSuperClusters * theSuperClusters * tpMapSuperClusters  )

theGsfElectrons = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string('(abs(superCluster.eta)<2.8) && !(1.4442<abs(superCluster.eta)<1.560) && (( caloEnergy * sin( caloPosition.theta ) )  > 15.0)')
)

HFElectronID = cms.EDFilter("CandViewSelector",
                                src = cms.InputTag("hfRecoEcalCandidate"),
                                cut = cms.string('et > 15.0')
                            )


theGsfHf = cms.EDFilter("CandViewMerger",
    #src = cms.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theHFSuperClusters"), cms.InputTag("theEEHFGapSuperClusters"))
    src = cms.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theHFSuperClusters"))
)

tpMapGsfElectrons = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("theGsfElectrons allSuperClusters"), # charge coniugate states are implied
    cut   = cms.string("45 < mass < 130"),
    checkCharge = cms.bool(False),
)

theHLTGsfHf = cms.EDFilter("CandViewMerger",
    #src = cms.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theHFSuperClusters"), cms.InputTag("theEEHFGapSuperClusters"))
    src = cms.VInputTag(cms.InputTag("ElectronID95"), cms.InputTag("HFElectronID"))
)

tpMapHLTElectrons = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("theHLT theHLTGsfHf"), # charge coniugate states are implied
    cut   = cms.string("60 < mass < 120"),
    checkCharge = cms.bool(False),
)


#tpMapGsfElectrons = cms.EDProducer("CandViewShallowCloneCombiner",
#                                       decay = cms.string("theGsfElectrons allSuperClusters"), # charge coniugate states are implied
#                                       cut   = cms.string("30 < mass < 130"),
#                                       checkCharge = cms.bool(False),
#                                   )

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







from ElectroWeakAnalysis.ZEE.simpleEleIdSequence_cff import *
from RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff import *
from PhysicsTools.PatAlgos.patSequences_cff import *

patElectronIDs = cms.Sequence(simpleEleIdSequence)

#electronEcalRecHitIsolationLcone.ecalBarrelRecHitProducer = cms.InputTag("reducedEcalRecHitsEB")
#electronEcalRecHitIsolationScone.ecalBarrelRecHitProducer = cms.InputTag("reducedEcalRecHitsEB")
#electronEcalRecHitIsolationLcone.ecalEndcapRecHitProducer = cms.InputTag("reducedEcalRecHitsEE")
#electronEcalRecHitIsolationScone.ecalEndcapRecHitProducer = cms.InputTag("reducedEcalRecHitsEE")
#electronEcalRecHitIsolationLcone.ecalBarrelRecHitCollection = cms.InputTag("")
#electronEcalRecHitIsolationScone.ecalBarrelRecHitCollection = cms.InputTag("")
#electronEcalRecHitIsolationLcone.ecalEndcapRecHitCollection = cms.InputTag("")
#electronEcalRecHitIsolationScone.ecalEndcapRecHitCollection = cms.InputTag("")
patElectronIsolation = cms.Sequence(egammaIsolationSequence)


makePatElectrons = cms.Sequence(patElectronIDs*patElectronIsolation*patElectrons)

patElectrons.addElectronID = cms.bool(True)
patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),
)
patElectrons.addGenMatch = cms.bool(False)
patElectrons.embedGenMatch = cms.bool(False)
patElectrons.isoDeposits = cms.PSet()
patElectrons.userIsolation = cms.PSet()

Iso95 =  cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectrons"),
    cut = cms.string('pt > 15.0 &&  (electronID("simpleEleId95relIso") > 6-0.1 && electronID("simpleEleId95relIso") < 7+.1)')                
)

ElectronID95 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectrons"),
    cut = cms.string('pt > 15.0 && ( electronID("simpleEleId95relIso") > 7-0.1 && electronID("simpleEleId95relIso") < 7+.1)')                
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
                  hltTag = cms.untracked.InputTag("HLT_Photon15_L1R","","HLT"),
                  triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
                        )


tpMapHFZs = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("theHLT HFElectronID"), # charge coniugate states are implied
                                       cut   = cms.string("70 < mass < 110"),
                                       checkCharge = cms.bool(False),
                                   )


#electron_sequence = cms.Sequence(electrons * theGsfElectrons * theGsfHf * theIsolation * eidRobust * theId * theHLT * HFElectronID )

#electron_sequence = cms.Sequence(electrons * theGsfElectrons *  HFElectronID * theGsfHf * tpMapEnd )

electron_sequence = cms.Sequence(theGsfElectrons * theGsfHf * tpMapGsfElectrons )

hfelectron_sequence = cms.Sequence(theGsfElectrons * theGsfHf * theIsolation * theId * theHLT * HFElectronID * tpMapHFZs )

felectron_sequence = cms.Sequence(theGsfElectrons * theGsfHf * makePatElectrons * Iso95 * ElectronID95 * theIsolation * theId * theHLT * HFElectronID *  theHLTGsfHf * tpMapHLTElectrons  )

#lepton_cands = cms.Sequence(genParticles * hfEMClusteringSequence * sc_sequence * electron_sequence * tpMap_sequence * truthMatch_sequence)
#lepton_cands = cms.Sequence(genParticles *  ZIntoElectronsEventProducer * SmearedElectronsProducer * hfEMClusteringSequence * sc_sequence * electrons * theGsfElectrons)
lepton_cands = cms.Sequence( hfEMClusteringSequence * sc_sequence )


   

