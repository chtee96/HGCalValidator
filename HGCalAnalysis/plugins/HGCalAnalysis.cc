//
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"
#include "DataFormats/HGCalReco/interface/TICLSeedingRegion.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Validation/HGCalValidation/interface/CaloParticleSelector.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"

#include <boost/algorithm/string.hpp> 

#include "HGCalValidator/HGCalAnalysis/interface/HGCalInfo.h"

#include <TFile.h>

#include <cmath>
#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
using namespace edm;

class HGCalAnalysis : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit HGCalAnalysis(const edm::ParameterSet &);
  ~HGCalAnalysis() override {}

private:

  void beginJob() override;
  void endJob() override;
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void analyze(edm::Event const &, edm::EventSetup const &) override;
  void endRun(edm::Run const &, edm::EventSetup const &) override;
  void cpParametersAndSelection(std::vector<CaloParticle> const& cPeff,
                                std::vector<SimVertex> const& simVertices,
                                std::vector<size_t>& selected_cPeff,
                                unsigned int layers,
                                std::unordered_map<DetId, const HGCRecHit*> const&) const;


private:

  //event
  hgcal_validation::eventInfo fEventInfo; 
  //rechits raw from hit map
  hgcal_validation::recHitRawInfo fRecHitRawInfo;
  //rechits 
  hgcal_validation::recHitInfo fRecHitInfo, fRecHitInfo_sc, fRecHitInfo_lc;
  //simClusters
  hgcal_validation::simClustersInfo fSimClustersInfo;
  //layerClusters
  hgcal_validation::layerClustersInfo fLayerClustersInfo;
  //caloParticles
  hgcal_validation::caloParticlesInfo fCaloParticlesInfo;
  //Tracksters
  hgcal_validation::trackstersInfo fTrackstersInfo;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::InputTag label_lcl;
  std::vector<edm::InputTag> label_tst;
  edm::InputTag label_simTS, label_simTSFromCP;
  edm::InputTag label_ticlSeedGlobal;
  edm::InputTag label_ticlSeedTrk;
  edm::InputTag associator_;
  edm::InputTag associatorSim_;
  const bool SaveGeneralInfo_;
  const bool doRecHitsRawTree_;
  const bool doRecHitsTree_;
  const bool doCaloParticleTree_;
  const bool doSimClustersTree_;
  const bool doSimClustersFromCPs_;
  edm::InputTag label_SimClustersPlots_, label_SimClustersLevel_;
  const bool doLayerClustersTree_;
  edm::InputTag label_layerClustersPlots_, label_LCToCPLinking_;
  const bool doCaloParticleSelection_;
  const bool doTrackstersPlots_;
  const bool doOnlyTrackstersMerge_;
  const bool doEdges_;
  edm::InputTag label_TS_, label_TSToCPLinking_, label_TSToSTSPR_;
  std::vector<edm::InputTag> label_clustersmask;
  //const bool doSimTrackstersFromCPsPlots_;

  const edm::FileInPath cummatbudinxo_;
  std::vector<std::string> trees_;
  bool createTree_;

  std::vector<edm::EDGetTokenT<reco::CaloClusterCollection>> labelToken;
  edm::EDGetTokenT<std::vector<SimCluster>> simClusters_;
  edm::EDGetTokenT<reco::CaloClusterCollection> layerclusters_;
  std::vector<edm::EDGetTokenT<ticl::TracksterCollection>> label_tstTokens;
  edm::EDGetTokenT<ticl::TracksterCollection> simTracksters_;
  edm::EDGetTokenT<ticl::TracksterCollection> simTrackstersFromCPs_;
  edm::EDGetTokenT<ticl::TracksterCollection> simTracksters_fromCPs_;
  edm::EDGetTokenT<std::map<uint, std::vector<uint>>> simTrackstersMap_;
  edm::EDGetTokenT<std::vector<TICLSeedingRegion>> ticlSeedingGlobalToken_;
  edm::EDGetTokenT<std::vector<TICLSeedingRegion>> ticlSeedingTrkToken_;
  edm::EDGetTokenT<std::vector<CaloParticle>> label_cp_effic;
  edm::EDGetTokenT<std::vector<CaloParticle>> label_cp_fake;
  edm::EDGetTokenT<std::vector<SimVertex>> simVertices_;
  std::vector<edm::EDGetTokenT<std::vector<float>>> clustersMaskTokens_;
  edm::EDGetTokenT<HGCeeUncalibratedRecHitCollection> eeUncalibRecHitCollection_;
  edm::EDGetTokenT<HGChefUncalibratedRecHitCollection> hefUncalibRecHitCollection_;
  edm::EDGetTokenT<HGChebUncalibratedRecHitCollection> hebUncalibRecHitCollection_;
  edm::EDGetTokenT<std::unordered_map<DetId, const HGCRecHit*>> hitMap_;
  edm::EDGetTokenT<hgcal::RecoToSimCollection> associatorMapRtS;
  edm::EDGetTokenT<hgcal::SimToRecoCollection> associatorMapStR;
  CaloParticleSelector cpSelector;
  edm::EDGetTokenT<hgcal::SimToRecoCollectionWithSimClusters> associatorMapSimtR;
  edm::EDGetTokenT<hgcal::RecoToSimCollectionWithSimClusters> associatorMapRtSim;

  std::shared_ptr<hgcal::RecHitTools> tools_;
  std::map<double, double> cummatbudg;
  std::vector<int> particles_to_monitor_;
  unsigned totallayers_to_monitor_;

  //One tree for each object we monitor
  std::unordered_map<std::string, TTree * > fTree;

};

HGCalAnalysis::HGCalAnalysis(const edm::ParameterSet& pset) : 
  caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
  label_lcl(pset.getParameter<edm::InputTag>("label_lcl")),
  label_tst(pset.getParameter<std::vector<edm::InputTag>>("label_tst")),
  label_simTS(pset.getParameter<edm::InputTag>("label_simTS")),
  label_simTSFromCP(pset.getParameter<edm::InputTag>("label_simTSFromCP")),
  label_ticlSeedGlobal(pset.getParameter<edm::InputTag>("label_ticlSeedGlobal")),
  label_ticlSeedTrk(pset.getParameter<edm::InputTag>("label_ticlSeedTrk")),
  associator_(pset.getUntrackedParameter<edm::InputTag>("associator")),
  associatorSim_(pset.getUntrackedParameter<edm::InputTag>("associatorSim")),
  SaveGeneralInfo_(pset.getUntrackedParameter<bool>("SaveGeneralInfo")),
  doRecHitsRawTree_(pset.getUntrackedParameter<bool>("doRecHitsRawTree")),
  doRecHitsTree_(pset.getUntrackedParameter<bool>("doRecHitsTree")),
  doCaloParticleTree_(pset.getUntrackedParameter<bool>("doCaloParticleTree")),
  doSimClustersTree_(pset.getUntrackedParameter<bool>("doSimClustersTree")),
  doSimClustersFromCPs_(pset.getUntrackedParameter<bool>("doSimClustersFromCPs")),
  label_SimClustersPlots_(pset.getParameter<edm::InputTag>("label_SimClusters")),
  label_SimClustersLevel_(pset.getParameter<edm::InputTag>("label_SimClustersLevel")),
  doLayerClustersTree_(pset.getUntrackedParameter<bool>("doLayerClustersTree")),
  label_layerClustersPlots_(pset.getParameter<edm::InputTag>("label_layerClusterPlots")),
  label_LCToCPLinking_(pset.getParameter<edm::InputTag>("label_LCToCPLinking")),
  doCaloParticleSelection_(pset.getUntrackedParameter<bool>("doCaloParticleSelection")),
  doTrackstersPlots_(pset.getUntrackedParameter<bool>("doTrackstersPlots")),
  doOnlyTrackstersMerge_(pset.getUntrackedParameter<bool>("doOnlyTrackstersMerge")),
  doEdges_(pset.getUntrackedParameter<bool>("doEdges")),
  label_TS_(pset.getParameter<edm::InputTag>("label_TS")),
  label_TSToCPLinking_(pset.getParameter<edm::InputTag>("label_TSToCPLinking")),
  label_TSToSTSPR_(pset.getParameter<edm::InputTag>("label_TSToSTSPR")),
  label_clustersmask(pset.getParameter<std::vector<edm::InputTag>>("LayerClustersInputMask")),
//  doSimTrackstersPlots_(pset.getUntrackedParameter<bool>("doSimTrackstersPlots")),
//  doSimTrackstersFromCPsPlots_(pset.getUntrackedParameter<bool>("doSimTrackstersFromCPsPlots")),
  cummatbudinxo_(pset.getParameter<edm::FileInPath>("cummatbudinxo")),
  trees_(pset.getUntrackedParameter<std::vector<std::string> >("trees")),
  createTree_(pset.getUntrackedParameter<bool>("createTree"))
{

  usesResource(TFileService::kSharedResource);

  //In this way we can easily generalize to associations between other objects also.
  const edm::InputTag& label_cp_effic_tag = pset.getParameter<edm::InputTag>("label_cp_effic");
  const edm::InputTag& label_cp_fake_tag = pset.getParameter<edm::InputTag>("label_cp_fake");

  label_cp_effic = consumes<std::vector<CaloParticle>>(label_cp_effic_tag);
  label_cp_fake = consumes<std::vector<CaloParticle>>(label_cp_fake_tag);

  simVertices_ = consumes<std::vector<SimVertex>>(pset.getParameter<edm::InputTag>("simVertices"));

  for (auto& itag : label_clustersmask) {
    clustersMaskTokens_.push_back(consumes<std::vector<float>>(itag));
  }

  associatorMapSimtR = consumes<hgcal::SimToRecoCollectionWithSimClusters>(associatorSim_);
  associatorMapRtSim = consumes<hgcal::RecoToSimCollectionWithSimClusters>(associatorSim_);

  cpSelector = CaloParticleSelector(pset.getParameter<double>("ptMinCP"),
                                    pset.getParameter<double>("ptMaxCP"),
                                    pset.getParameter<double>("minRapidityCP"),
                                    pset.getParameter<double>("maxRapidityCP"),
                                    pset.getParameter<double>("lipCP"),
                                    pset.getParameter<double>("tipCP"),
                                    pset.getParameter<int>("minHitCP"),
                                    pset.getParameter<int>("maxSimClustersCP"),
                                    pset.getParameter<bool>("signalOnlyCP"),
                                    pset.getParameter<bool>("intimeOnlyCP"),
                                    pset.getParameter<bool>("chargedOnlyCP"),
                                    pset.getParameter<bool>("stableOnlyCP"),
                                    pset.getParameter<bool>("notConvertedOnlyCP"),
                                    pset.getParameter<std::vector<int>>("pdgIdCP"));

  eeUncalibRecHitCollection_ =  consumes<HGCeeUncalibratedRecHitCollection>(pset.getParameter<edm::InputTag>("HGCEEuncalibRecHitCollection"));
  hefUncalibRecHitCollection_ = consumes<HGChefUncalibratedRecHitCollection>(pset.getParameter<edm::InputTag>("HGCHEFuncalibRecHitCollection"));
  hebUncalibRecHitCollection_ = consumes<HGChebUncalibratedRecHitCollection>(pset.getParameter<edm::InputTag>("HGCHEBuncalibRecHitCollection"));

  simTrackstersMap_ = consumes<std::map<uint, std::vector<uint>>>(edm::InputTag("ticlSimTracksters"));

  hitMap_ = consumes<std::unordered_map<DetId, const HGCRecHit*>>(edm::InputTag("hgcalRecHitMapProducer"));

  simClusters_ = consumes<std::vector<SimCluster>>(pset.getParameter<edm::InputTag>("label_scl"));

  layerclusters_ = consumes<reco::CaloClusterCollection>(label_lcl);

  for (auto& itag : label_tst) {
    label_tstTokens.push_back(consumes<ticl::TracksterCollection>(itag));
  }

  simTracksters_ = consumes<ticl::TracksterCollection>(label_simTS);
  simTracksters_fromCPs_ = consumes<ticl::TracksterCollection>(label_simTSFromCP);

  ticlSeedingGlobalToken_ = consumes<std::vector<TICLSeedingRegion>>(label_ticlSeedGlobal);
  ticlSeedingTrkToken_ = consumes<std::vector<TICLSeedingRegion>>(label_ticlSeedTrk);
  associatorMapRtS = consumes<hgcal::RecoToSimCollection>(associator_);
  associatorMapStR = consumes<hgcal::SimToRecoCollection>(associator_);

  tools_.reset(new hgcal::RecHitTools());

  particles_to_monitor_ = pset.getParameter<std::vector<int>>("pdgIdCPs");
  totallayers_to_monitor_ = pset.getParameter<int>("totallayers_to_monitor");

  //For the material budget file here
  std::ifstream fmb(cummatbudinxo_.fullPath().c_str());
  double thelay = 0.;
  double mbg = 0.;
  for (unsigned ilayer = 1; ilayer <= 2 * totallayers_to_monitor_; ++ilayer) {
    fmb >> thelay >> mbg;
    cummatbudg.insert(std::pair<double, double>(thelay, mbg));
  }

  fmb.close();

}

void HGCalAnalysis::beginJob() {

  edm::Service<TFileService> fs;
  if (createTree_) {

    for (const auto& mtr : trees_) {
      //Be careful here: The trees_ are defined in the config file and should be
      //in agreement with the branches we will create below.
      //"RecHitsRawFromHitMap","SimClusters","CaloParticles","LayerClusters","Tracksters"
      fTree[mtr] = fs->make<TTree>(mtr.data(), mtr.data()); 
      //We want the event and rechit branches in all trees
      hgcal_validation::createEventBranches(fTree[mtr], fEventInfo);
      //hgcal_validation::createRecHitBranches(fTree[mtr], fRecHitInfo);
    }
    //make the structure we want in each tree
    hgcal_validation::createRecHitRawBranches(fTree["RecHitsRawFromHitMap"], fRecHitRawInfo);
    hgcal_validation::createSimClustersBranches(fTree["SimClusters"], fSimClustersInfo);
    hgcal_validation::createCaloParticlesBranches(fTree["CaloParticles"], fCaloParticlesInfo, fSimClustersInfo);
    hgcal_validation::createLayerClustersBranches(fTree["LayerClusters"], fLayerClustersInfo);
    //Here the trees should be in agreement with the created in the tree_ config or a
    //segmentation fault will be produced
    for (const auto& itag : label_tst) {
      hgcal_validation::createTrackstersBranches(fTree[itag.label()], fTrackstersInfo);
      //Add to the trackster branches the relevant LayerCluster and SimCluster ones. 
      hgcal_validation::createSimClustersBranches(fTree[itag.label()], fSimClustersInfo);
      hgcal_validation::createLayerClustersBranches(fTree[itag.label()], fLayerClustersInfo);

    }
    
  }

}

void HGCalAnalysis::analyze(const edm::Event &event, const edm::EventSetup &setup) {
 
  edm::ESHandle<CaloGeometry> geom = setup.getHandle(caloGeomToken_);
  tools_->setGeometry(*geom);

  edm::Handle<ticl::TracksterCollection> simTracksterHandle;
  event.getByToken(simTracksters_, simTracksterHandle);
  ticl::TracksterCollection const& simTracksters = *simTracksterHandle;

  edm::Handle<ticl::TracksterCollection> simTracksterFromCPHandle;
  event.getByToken(simTracksters_fromCPs_, simTracksterFromCPHandle);
  ticl::TracksterCollection const& simTrackstersFromCPs = *simTracksterFromCPHandle;

  edm::Handle<std::map<uint, std::vector<uint>>> simTrackstersMapHandle;
  event.getByToken(simTrackstersMap_, simTrackstersMapHandle);
  const std::map<uint, std::vector<uint>> cpToSc_SimTrackstersMap = *simTrackstersMapHandle;

  edm::Handle<std::vector<CaloParticle>> caloParticleHandle;
  event.getByToken(label_cp_effic, caloParticleHandle);
  std::vector<CaloParticle> const& caloParticles = *caloParticleHandle;

  edm::Handle<std::vector<SimVertex>> simVerticesHandle;
  event.getByToken(simVertices_, simVerticesHandle);
  std::vector<SimVertex> const& simVertices = *simVerticesHandle;

  //Apart from the hitMap with rechits we need the uncalibrated rechits
  //which may have duplicates.
  const auto& ee_hits = event.get(eeUncalibRecHitCollection_);
  const auto& fh_hits = event.get(hefUncalibRecHitCollection_);
  const auto& bh_hits = event.get(hebUncalibRecHitCollection_);
  std::unordered_map<DetId, std::vector<const HGCUncalibratedRecHit *> > UnCalibHitMap;
  UnCalibHitMap.clear();
    
  for (const auto& hit : ee_hits) {
    // std::cout << hit.id().rawId() << std::endl;
    // std::cout << hit.amplitude() << std::endl;
    UnCalibHitMap[hit.id()].emplace_back(&hit);
  }
  for (const auto& hit : fh_hits) { UnCalibHitMap[hit.id()].emplace_back(&hit); }
  for (const auto& hit : bh_hits) { UnCalibHitMap[hit.id()].emplace_back(&hit); }

  // for (const auto& hit : ee_hits) {
  //   if (UnCalibHitMap[hit.id()].size() >=2) {std::cout << "SIZE " << UnCalibHitMap[hit.id()].size() << std::endl;;}
  // }
  edm::Handle<std::unordered_map<DetId, const HGCRecHit*>> hitMapHandle;
  event.getByToken(hitMap_, hitMapHandle);
  const std::unordered_map<DetId, const HGCRecHit*>* hitMap = &*hitMapHandle;

  std::vector<size_t> cPIndices;
  //Consider CaloParticles coming from the hard scatterer
  //excluding the PU contribution and save the indices.
  removeCPFromPU(caloParticles, cPIndices);
  
  std::vector<size_t> selected_cPeff;
  cpParametersAndSelection(caloParticles, simVertices, selected_cPeff, totallayers_to_monitor_, *hitMap);

  //---------------------------------------------------------------------------------------------------
  //Event 
  //---------------------------------------------------------------------------------------------------
  if(SaveGeneralInfo_){
    hgcal_validation::initEventInfo(fEventInfo);
    hgcal_validation::fillEventInfo(fEventInfo,event);
  }

  //---------------------------------------------------------------------------------------------------
  //RecHitsRaw 
  //---------------------------------------------------------------------------------------------------
  if (doRecHitsRawTree_){
    hgcal_validation::initRecHitRawInfo(fRecHitRawInfo);
    hgcal_validation::fillRecHitRawInfo(fRecHitRawInfo, *hitMap, tools_, totallayers_to_monitor_);
    fTree["RecHitsRawFromHitMap"]->Fill();
  }

  //---------------------------------------------------------------------------------------------------
  //SimClusters
  //---------------------------------------------------------------------------------------------------
  if (doSimClustersTree_) {
    edm::Handle<std::vector<SimCluster>> simClustersHandle;
    event.getByToken(simClusters_, simClustersHandle);
    std::vector<SimCluster> const& simClusters = *simClustersHandle;
    
    std::vector<SimCluster> simClustersFromCPs;
    simClustersFromCPs.clear();
    for (const auto& cpId : cPIndices) {
      const SimClusterRefVector& simClusterRefVector = caloParticles[cpId].simClusters();
      for (const auto& it_sc : simClusterRefVector) {
      	SimCluster sc = (*(it_sc));
	simClustersFromCPs.push_back(sc);
      }
    }//end of loop over caloParticles


    //To check also the unmatched rechits
    hgcal_validation::initRecHitRawInfo(fRecHitRawInfo);
    hgcal_validation::fillRecHitRawInfo(fRecHitRawInfo, *hitMap, tools_, totallayers_to_monitor_);

    //This is for the matched ones
    hgcal_validation::initRecHitInfo(fRecHitInfo);
    //Now to the SimClusters
    hgcal_validation::initSimClustersInfo(fSimClustersInfo);
    if (doSimClustersFromCPs_){
      hgcal_validation::fillSimClustersInfo(fSimClustersInfo, fRecHitInfo, simClustersFromCPs, UnCalibHitMap, *hitMap, tools_, totallayers_to_monitor_, cummatbudg);
    } else{
      hgcal_validation::fillSimClustersInfo(fSimClustersInfo, fRecHitInfo, simClusters, UnCalibHitMap, *hitMap, tools_, totallayers_to_monitor_, cummatbudg);
    }
    fTree["SimClusters"]->Fill();

  }//end of simClusters

  //---------------------------------------------------------------------------------------------------
  //LayerClusters
  //---------------------------------------------------------------------------------------------------
  if (doLayerClustersTree_) {
    edm::Handle<reco::CaloClusterCollection> clusterHandle;
    event.getByToken(layerclusters_, clusterHandle);
    const reco::CaloClusterCollection& clusters = *clusterHandle;

    if (doRecHitsTree_){hgcal_validation::initRecHitInfo(fRecHitInfo);}
    hgcal_validation::initLayerClustersInfo(fLayerClustersInfo);
    hgcal_validation::fillLayerClustersInfo(fLayerClustersInfo, fRecHitInfo, clusters, 9999, 9999, false, UnCalibHitMap, *hitMap, tools_, totallayers_to_monitor_, cummatbudg);
    fTree["LayerClusters"]->Fill();

  }

  //---------------------------------------------------------------------------------------------------
  //CaloParticles
  //---------------------------------------------------------------------------------------------------
  if (doCaloParticleTree_) {

    //We will also need the SimClusters, so we initialize them too. 
    hgcal_validation::initSimClustersInfo(fSimClustersInfo);
    hgcal_validation::initCaloParticlesInfo(fCaloParticlesInfo);
    hgcal_validation::fillCaloParticlesInfo(fCaloParticlesInfo, fSimClustersInfo, caloParticles, cPIndices, *hitMap, tools_, totallayers_to_monitor_);
    fTree["CaloParticles"]->Fill();
  }

  //---------------------------------------------------------------------------------------------------
  //Tracksters
  //---------------------------------------------------------------------------------------------------
  if (doTrackstersPlots_) {

    edm::Handle<std::vector<TICLSeedingRegion>> ticlSeedingGlobalH;
    event.getByToken(ticlSeedingGlobalToken_, ticlSeedingGlobalH);
    auto const& ticlSeedingGlobal = *ticlSeedingGlobalH.product();

    edm::Handle<std::vector<TICLSeedingRegion>> ticlSeedingTrkH;
    event.getByToken(ticlSeedingTrkToken_, ticlSeedingTrkH);
    auto const& ticlSeedingTrk = *ticlSeedingTrkH.product();

    edm::Handle<std::vector<SimCluster>> simClustersHandle;
    event.getByToken(simClusters_, simClustersHandle);
    std::vector<SimCluster> const& simClusters = *simClustersHandle;

    edm::Handle<reco::CaloClusterCollection> clusterHandle;
    event.getByToken(layerclusters_, clusterHandle);
    const reco::CaloClusterCollection& clusters = *clusterHandle;

    for (unsigned int wml = 0; wml < label_tstTokens.size(); wml++) {

      //Will work only with merge collection for now
      if ( doOnlyTrackstersMerge_ && label_tst[wml].label() != "ticlTrackstersMerge" ){ continue; }
      
      edm::Handle<ticl::TracksterCollection> tracksterHandle;
      event.getByToken(label_tstTokens[wml], tracksterHandle);
      const ticl::TracksterCollection& tracksters = *tracksterHandle;

      hgcal_validation::initRecHitInfo(fRecHitInfo_sc);
      hgcal_validation::initRecHitInfo(fRecHitInfo_lc);
      hgcal_validation::initSimClustersInfo(fSimClustersInfo);
      hgcal_validation::initLayerClustersInfo(fLayerClustersInfo);
      hgcal_validation::initTrackstersInfo(fTrackstersInfo);

      // Pattern recognition
      // hgcal_validation::prepare_tracksters_to_SimTracksters(fTrackstersInfo,
      // 							    tracksters,
      // 							    clusters,
      // 							    simTracksters,
      // 							    1,
      // 							    simTrackstersFromCPs,
      // 							    cpToSc_SimTrackstersMap,
      // 							    simClusters,
      // 							    caloParticleHandle.id(),
      // 							    caloParticles,
      // 							    cPIndices,
      // 							    selected_cPeff,
      // 							    *hitMap,
      // 							    totallayers_to_monitor_);


      hgcal_validation::fillTrackstersInfo(fTrackstersInfo, fSimClustersInfo, fLayerClustersInfo, fRecHitInfo_sc, fRecHitInfo_lc, tracksters, simClusters, clusters, ticlSeedingGlobal, ticlSeedingTrk, UnCalibHitMap, *hitMap, tools_, doEdges_, totallayers_to_monitor_, cummatbudg);

      fTree[label_tst[wml].label()]->Fill();
      //fTree["LayerClusters"]->Fill();
      // fTree["SimClusters"]->Fill();

    }//end of loop over Trackster input labels
  }  


}

void HGCalAnalysis::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {}

void HGCalAnalysis::endRun(edm::Run const &iEvent, edm::EventSetup const &) {}

void HGCalAnalysis::endJob() {}

void HGCalAnalysis::cpParametersAndSelection(std::vector<CaloParticle> const& cPeff,
					     std::vector<SimVertex> const& simVertices,
					     std::vector<size_t>& selected_cPeff,
					     unsigned int layers,
					     std::unordered_map<DetId, const HGCRecHit*> const& hitMap) const {
  selected_cPeff.reserve(cPeff.size());

  size_t j = 0;
  for (auto const& caloParticle : cPeff) {
    int id = caloParticle.pdgId();

    if (!doCaloParticleSelection_ || (doCaloParticleSelection_ && cpSelector(caloParticle, simVertices))) {
      selected_cPeff.push_back(j);
    }
    ++j;
  }  //end of loop over caloparticles
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalAnalysis);
