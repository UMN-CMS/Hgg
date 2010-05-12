import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.genParticles_cfi import *

#from PhysicsTools.TagAndProbe.tag_probe_electron_cfi import *

from RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff import *

from Geometry.CMSCommonData.cmsIdealGeometryXML_cff import *

from Geometry.CaloEventSetup.CaloGeometry_cff import *

from Configuration.EventContent.EventContent_cff import *

hfRecoEcalCandidate.Correct = True
##hfRecoEcalCandidate.e9e25Cut = 0
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
    cut = cms.string('et > 5.0')
)


allSuperClusters = cms.EDFilter("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("EBSuperClusters"), cms.InputTag("EESuperClusters"),cms.InputTag("theHFSuperClusters"))
)


theSuperClusters = cms.EDFilter("CandViewSelector",
    #src = cms.InputTag("superClusterCands"),
    src = cms.InputTag("allSuperClusters"),
    cut = cms.string('et  > 5.0 && abs(eta)<2.8 && !(1.4442< abs(eta) <1.560)'),
    filter = cms.bool(True) #LOOK UP WHAT THIS DOES
)

####NOW I DO THE MAP
tpMapSuperClusters = cms.EDFilter("CandViewShallowCloneCombiner",
    decay = cms.string("theSuperClusters allSuperClusters"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("30 < mass < 140"),
    filter = cms.bool(True) #LOOK UP WHAT THIS DOES
)




sc_sequence = cms.Sequence( ( hfSuperClusterCandidate * theHFSuperClusters) * HybridSuperClusters * EBSuperClusters + EndcapSuperClusters * EESuperClusters * allSuperClusters * theSuperClusters * tpMapSuperClusters  )



#electron_sequence = cms.Sequence(electrons * theGsfElectrons * theGsfHf * theIsolation * eidRobust * theId * theHLT * HFElectronID )

#electron_sequence = cms.Sequence(electrons * theGsfElectrons *  HFElectronID * theGsfHf * tpMapEnd )




#lepton_cands = cms.Sequence(genParticles * hfEMClusteringSequence * sc_sequence * electron_sequence * tpMap_sequence * truthMatch_sequence)
#lepton_cands = cms.Sequence(genParticles *  ZIntoElectronsEventProducer * SmearedElectronsProducer * hfEMClusteringSequence * sc_sequence * electrons * theGsfElectrons)
lepton_cands = cms.Sequence( hfEMClusteringSequence * sc_sequence )


   

