import FWCore.ParameterSet.Config as cms

from Validation.HGCalValidation.CaloParticleSelectionForEfficiency_cfi import *

from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociation
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociation

from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels, ticlIterLabelsMerge

labelTst = [cms.InputTag("ticlTracksters"+iteration) for iteration in ticlIterLabelsMerge]
labelTst.extend([cms.InputTag("ticlSimTracksters", "fromCPs"), cms.InputTag("ticlSimTracksters")])
#I will add the general simTracksters collection later, be careful with the tree below. 
#labelTst.extend([cms.InputTag("ticlSimTracksters", "fromCPs")])
lcInputMask = [cms.InputTag("ticlTracksters"+iteration) for iteration in ticlIterLabels]
lcInputMask.extend([cms.InputTag("ticlSimTracksters", "fromCPs"), cms.InputTag("ticlSimTracksters")])
#lcInputMask.extend([cms.InputTag("ticlSimTracksters", "fromCPs")])

HGCalAnalysis = cms.EDAnalyzer(
    "HGCalAnalysis",
    # selection of CP for evaluation of efficiency #
    CaloParticleSelectionForEfficiency,

    label_lcl = layerClusterCaloParticleAssociation.label_lc,
    label_tst = cms.VInputTag(labelTst),
    label_simTS = cms.InputTag("ticlSimTracksters"),
    label_simTSFromCP = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_ticlSeedGlobal = cms.InputTag("ticlSeedingGlobal"),
    label_ticlSeedTrk = cms.InputTag("ticlSeedingTrk"),    
    associator = cms.untracked.InputTag("layerClusterCaloParticleAssociationProducer"),
    associatorSim = cms.untracked.InputTag("layerClusterSimClusterAssociationProducer"),
    #General event info
    SaveGeneralInfo = cms.untracked.bool(True),
    HGCEEuncalibRecHitCollection = cms.InputTag('HGCalUncalibRecHit:HGCEEUncalibRecHits'),
    HGCHEFuncalibRecHitCollection = cms.InputTag('HGCalUncalibRecHit:HGCHEFUncalibRecHits'),
    HGCHEBuncalibRecHitCollection = cms.InputTag('HGCalUncalibRecHit:HGCHEBUncalibRecHits'),
    #RecHitsRaw tree
    doRecHitsRawTree = cms.untracked.bool(False),
    #RecHits tree
    doRecHitsTree = cms.untracked.bool(True),
    #CaloParticle tree
    doCaloParticleTree = cms.untracked.bool(False),
    #SimCluster tree
    doSimClustersTree = cms.untracked.bool(True),
    doSimClustersFromCPs = cms.untracked.bool(False),
    label_SimClusters = cms.InputTag("SimClusters"),
    label_SimClustersLevel = cms.InputTag("ClusterLevel"),
    #Layer Cluster tree
    doLayerClustersTree = cms.untracked.bool(False),
    label_layerClusterPlots = cms.InputTag("hgcalLayerClusters"),
    label_LCToCPLinking = cms.InputTag("LCToCP_association"),
    #Select caloParticles for efficiency or pass through
    doCaloParticleSelection = cms.untracked.bool(True),
    #Trackster and SimTrakcsters trees
    doTrackstersPlots = cms.untracked.bool(False),
    doOnlyTrackstersMerge = cms.untracked.bool(True),
    doEdges = cms.untracked.bool(False),
    label_TS = cms.InputTag("Morphology"),
    label_TSToCPLinking = cms.InputTag("TSToCP_linking"),
    label_TSToSTSPR = cms.InputTag("TSToSTS_patternRecognition"),
    #SimTrackster from CPs tree
    #doSimTrackstersFromCPsPlots = cms.untracked.bool(True),
    ### sim input configuration ###
    label_cp_effic = layerClusterCaloParticleAssociation.label_cp,
    label_cp_fake = cms.InputTag("mix","MergedCaloTruth"),
    #simClusters
    label_scl = layerClusterSimClusterAssociation.label_scl,
    simVertices = cms.InputTag("g4SimHits"),
    LayerClustersInputMask = cms.VInputTag(lcInputMask),
    #311: K0, 130: K0_short, 310: K0_long
    pdgIdCPs = cms.vint32(11, -11, 13, -13, 22, 111, 211, -211, 321, -321, 311, 130, 310),
    #Total number of layers of HGCal that we want to monitor
    #Could get this also from HGCalImagingAlgo::maxlayer but better to get it from here
    totallayers_to_monitor = cms.int32(52),
    #The cumulative material budget in front of each layer. To be more specific, it
    #is the material budget just in front of the active material (not including it).
    #This file is created using the official material budget code.
    cummatbudinxo = cms.FileInPath('HGCalValidator/HGCalAnalysis/data/D88.cumulative.xo'),    
    trees = cms.untracked.vstring(["RecHitsRawFromHitMap","SimClusters","CaloParticles","LayerClusters","ticlTrackstersTrkEM","ticlTrackstersEM","ticlTrackstersTrk","ticlTrackstersHAD","ticlTrackstersMerge","ticlSimTracksters_fromCPs","ticlSimTracksters"]),
    createTree = cms.untracked.bool(True)
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(HGCalAnalysis,
    label_cp_fake = "mixData:MergedCaloTruth"
)

from Configuration.Eras.Modifier_phase2_hgcalV10_cff import phase2_hgcalV10
phase2_hgcalV10.toModify(HGCalAnalysis, totallayers_to_monitor = cms.int32(50))

from Configuration.Eras.Modifier_phase2_hgcalV16_cff import phase2_hgcalV16
phase2_hgcalV16.toModify(HGCalAnalysis, totallayers_to_monitor = cms.int32(47))


