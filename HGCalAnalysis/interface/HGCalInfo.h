#ifndef Validation_HGCalValidation_HGCalInfo_h
#define Validation_HGCalValidation_HGCalInfo_h

#include <vector>

#include <TTree.h>

// Forward declarations
class EncodedEventId;

const double ScoreCutTStoSTSFakeMerge_[] = {0.6, 1.e-09};  //FLT_MIN
const double ScoreCutSTStoTSPurDup_[] = {0.2, 1.e-11};     //FLT_MIN

namespace hgcal_validation
{

  //---------------------------------------------------------------------------------------------------
  //Some useful structs for later when checking associators. 
  //---------------------------------------------------------------------------------------------------
  struct detIdInfoInCluster {
    bool operator==(const detIdInfoInCluster& o) const { return clusterId == o.clusterId; };
    long unsigned int clusterId;
    float fraction;
  };

  struct detIdInfoInTrackster {
    bool operator==(const detIdInfoInTrackster& o) const { return tracksterId == o.tracksterId; };
    unsigned int tracksterId;
    long unsigned int clusterId;
    float fraction;
  };

  struct caloParticleOnLayer {
    unsigned int caloParticleId;
    float energy = 0;
    std::vector<std::pair<DetId, float>> hits_and_fractions;
    std::unordered_map<unsigned int, std::pair<float, float>> layerClusterIdToEnergyAndScore;
  };

  //---------------------------------------------------------------------------------------------------
  //Event 
  //---------------------------------------------------------------------------------------------------
  struct eventInfo {
    Long64_t run;
    Long64_t LS;
    Long64_t evnum; 
    Long64_t processid;
    Long64_t bx;
    Long64_t orbit;
  };

  void initEventInfo(eventInfo &evInfo) {
    evInfo.run       = (Long64_t) -99999.99;
    evInfo.LS        = (Long64_t) -99999.99;
    evInfo.evnum     = (Long64_t) -99999.99;
    evInfo.processid = (Long64_t) -99999.99;
    evInfo.bx        = (Long64_t) -99999.99;
    evInfo.orbit     = (Long64_t) -99999.99;
  }

   void createEventBranches(TTree *fTree, eventInfo &evInfo) {
     fTree->Branch("run",&evInfo.run);
     fTree->Branch("lumi",&evInfo.LS);
     fTree->Branch("event",&evInfo.evnum);
     fTree->Branch("processid",&evInfo.processid);
     fTree->Branch("bx",&evInfo.bx);
     fTree->Branch("orbit",&evInfo.orbit);
   }

  void fillEventInfo(eventInfo &evInfo, const edm::Event& iEvent) {
    evInfo.run   = iEvent.id().run();
    evInfo.LS    = iEvent.id().luminosityBlock();
    evInfo.evnum = iEvent.id().event();
    evInfo.bx    = iEvent.bunchCrossing();
    evInfo.orbit = iEvent.orbitNumber();
  }

  //---------------------------------------------------------------------------------------------------
  //RecHitsRaw from hitMap
  //---------------------------------------------------------------------------------------------------
  struct recHitRawInfo {
    std::vector<float> rechit_raw_eta;
    std::vector<float> rechit_raw_phi;
    std::vector<float> rechit_raw_pt;
    std::vector<float> rechit_raw_energy;
    std::vector<float> rechit_raw_x;
    std::vector<float> rechit_raw_y;
    std::vector<float> rechit_raw_z;
    std::vector<float> rechit_raw_time;
    std::vector<float> rechit_raw_thickness;
    std::vector<int> rechit_raw_layer;
    std::vector<int> rechit_raw_wafer_u;
    std::vector<int> rechit_raw_wafer_v;
    std::vector<int> rechit_raw_cell_u;
    std::vector<int> rechit_raw_cell_v;
    std::vector<unsigned int> rechit_raw_detid;
    std::vector<bool> rechit_raw_isHalf;
    std::vector<int> rechit_raw_flags;
    std::vector<float> rechit_raw_radius;
  };

  void initRecHitRawInfo(recHitRawInfo &rhInfo) {
    rhInfo.rechit_raw_eta.clear();
    rhInfo.rechit_raw_phi.clear();
    rhInfo.rechit_raw_pt.clear();
    rhInfo.rechit_raw_energy.clear();
    rhInfo.rechit_raw_x.clear();
    rhInfo.rechit_raw_y.clear();
    rhInfo.rechit_raw_z.clear();
    rhInfo.rechit_raw_time.clear();
    rhInfo.rechit_raw_thickness.clear();
    rhInfo.rechit_raw_layer.clear();
    rhInfo.rechit_raw_wafer_u.clear();
    rhInfo.rechit_raw_wafer_v.clear();
    rhInfo.rechit_raw_cell_u.clear();
    rhInfo.rechit_raw_cell_v.clear();
    rhInfo.rechit_raw_detid.clear();
    rhInfo.rechit_raw_isHalf.clear();
    rhInfo.rechit_raw_flags.clear();
    rhInfo.rechit_raw_radius.clear();
  }

  void createRecHitRawBranches(TTree *fTree, recHitRawInfo &rhInfo) {
    fTree->Branch("rechit_raw_eta", &rhInfo.rechit_raw_eta);
    fTree->Branch("rechit_raw_phi", &rhInfo.rechit_raw_phi);
    fTree->Branch("rechit_raw_pt", &rhInfo.rechit_raw_pt);
    fTree->Branch("rechit_raw_energy", &rhInfo.rechit_raw_energy);
    fTree->Branch("rechit_raw_x", &rhInfo.rechit_raw_x);
    fTree->Branch("rechit_raw_y", &rhInfo.rechit_raw_y);
    fTree->Branch("rechit_raw_z", &rhInfo.rechit_raw_z);
    fTree->Branch("rechit_raw_time", &rhInfo.rechit_raw_time);
    fTree->Branch("rechit_raw_thickness", &rhInfo.rechit_raw_thickness);
    fTree->Branch("rechit_raw_layer", &rhInfo.rechit_raw_layer);
    fTree->Branch("rechit_raw_wafer_u", &rhInfo.rechit_raw_wafer_u);
    fTree->Branch("rechit_raw_wafer_v", &rhInfo.rechit_raw_wafer_v);
    fTree->Branch("rechit_raw_cell_u", &rhInfo.rechit_raw_cell_u);
    fTree->Branch("rechit_raw_cell_v", &rhInfo.rechit_raw_cell_v);
    fTree->Branch("rechit_raw_detid", &rhInfo.rechit_raw_detid);
    fTree->Branch("rechit_raw_isHalf", &rhInfo.rechit_raw_isHalf);
    fTree->Branch("rechit_raw_flags", &rhInfo.rechit_raw_flags);
    fTree->Branch("rechit_raw_radius", &rhInfo.rechit_raw_radius);

  }

  void fillRecHitRawInfo(recHitRawInfo &rhInfo,
			 std::unordered_map<DetId, const HGCRecHit*> const& hitMap,
			 std::shared_ptr<hgcal::RecHitTools> recHitTools,
			 unsigned int layers) {
    
    for(const auto &i : hitMap ){
      const DetId rh_detid = i.first;
      const HGCRecHit* hit = i.second;

      const GlobalPoint position = recHitTools->getPosition(rh_detid);
      const double energy = hit->energy();

      int layer =
          recHitTools->getLayerWithOffset(rh_detid) + layers * ((recHitTools->zside(rh_detid) + 1) >> 1) - 1;
      
      //Will continue with the same convention as before for raw hits flags. 
      int flags = 0x3;
    
      std::pair<int, int> wafer;
      std::pair<int, int> cell;
      double thickness;

      if (rh_detid.det() == DetId::Forward || rh_detid.det() == DetId::HGCalEE || rh_detid.det() == DetId::HGCalHSi) {
	thickness = recHitTools->getSiThickness(rh_detid);
	wafer = recHitTools->getWafer(rh_detid);
	cell = recHitTools->getCell(rh_detid);
      } else {
	thickness = std::numeric_limits<std::float_t>::max();
	wafer = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
	cell = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
      }

      const bool isHalfCell = recHitTools->isHalfCell(rh_detid);
      const double eta = recHitTools->getEta(position);
      const double phi = recHitTools->getPhi(position);
      const double pt = recHitTools->getPt(position,energy);
      const double radius = 
      ((rh_detid.det() == DetId::Forward || rh_detid.det() == DetId::HGCalEE || rh_detid.det() == DetId::HGCalHSi) ? recHitTools->getRadiusToSide(rh_detid) : -1.);

      rhInfo.rechit_raw_eta.push_back(eta);
      rhInfo.rechit_raw_phi.push_back(phi);
      rhInfo.rechit_raw_pt.push_back(pt);
      rhInfo.rechit_raw_energy.push_back(energy);
      rhInfo.rechit_raw_layer.push_back(layer);
      rhInfo.rechit_raw_wafer_u.push_back(wafer.first);
      rhInfo.rechit_raw_wafer_v.push_back(wafer.second);
      rhInfo.rechit_raw_cell_u.push_back(cell.first);
      rhInfo.rechit_raw_cell_v.push_back(cell.second);
      rhInfo.rechit_raw_detid.push_back(rh_detid);
      rhInfo.rechit_raw_x.push_back(position.x());
      rhInfo.rechit_raw_y.push_back(position.y());
      rhInfo.rechit_raw_z.push_back(position.z());
      rhInfo.rechit_raw_time.push_back(hit->time());
      rhInfo.rechit_raw_thickness.push_back(thickness);
      rhInfo.rechit_raw_isHalf.push_back(isHalfCell);
      rhInfo.rechit_raw_flags.push_back(flags);
      rhInfo.rechit_raw_radius.push_back(radius);


    } //end of loop over RecHits in hitMap
  } //end of fillRecHitInfo

  //---------------------------------------------------------------------------------------------------
  //RecHits from the relevant related object (Layercluster, CaloParticle, Trackster, SimTrackster).
  //Depending on what object we monitor we will mark the relevant id with the object index in the 
  //collection and all others marked with -1. 
  //---------------------------------------------------------------------------------------------------
  struct recHitInfo {
    std::vector<float> rechit_eta;
    std::vector<float> rechit_phi;
    std::vector<float> rechit_pt;
    std::vector<float> rechit_energy;
    std::vector<float> rechit_recostructable_energy;
    std::vector<float> rechit_SoN;
    std::vector<float> rechit_uncalib_energy;
    std::vector<float> rechit_matbudget;
    std::vector<float> rechit_x;
    std::vector<float> rechit_y;
    std::vector<float> rechit_z;
    std::vector<float> rechit_time;
    std::vector<float> rechit_thickness;
    std::vector<int> rechit_layer;
    std::vector<int> rechit_wafer_u;
    std::vector<int> rechit_wafer_v;
    std::vector<int> rechit_cell_u;
    std::vector<int> rechit_cell_v;
    std::vector<unsigned int> rechit_detid;
    std::vector<bool> rechit_isHalf;
    std::vector<int> rechit_flags;
    std::vector<int> rechit_layerclusterid;
    std::vector<int> rechit_simclusterid;
    /* std::vector<int> rechit_caloparticleid; */
    std::vector<int> rechit_tracksterid;
    /* std::vector<int> rechit_simtracksterid; */
    std::vector<float> rechit_radius;
  };

  void initRecHitInfo(recHitInfo &rhInfo) {
    rhInfo.rechit_eta.clear();
    rhInfo.rechit_phi.clear();
    rhInfo.rechit_pt.clear();
    rhInfo.rechit_energy.clear();
    rhInfo.rechit_recostructable_energy.clear();
    rhInfo.rechit_SoN.clear();
    rhInfo.rechit_uncalib_energy.clear();
    rhInfo.rechit_matbudget.clear();
    rhInfo.rechit_x.clear();
    rhInfo.rechit_y.clear();
    rhInfo.rechit_z.clear();
    rhInfo.rechit_time.clear();
    rhInfo.rechit_thickness.clear();
    rhInfo.rechit_layer.clear();
    rhInfo.rechit_wafer_u.clear();
    rhInfo.rechit_wafer_v.clear();
    rhInfo.rechit_cell_u.clear();
    rhInfo.rechit_cell_v.clear();
    rhInfo.rechit_detid.clear();
    rhInfo.rechit_isHalf.clear();
    rhInfo.rechit_flags.clear();
    rhInfo.rechit_layerclusterid.clear();
    rhInfo.rechit_simclusterid.clear();
    rhInfo.rechit_tracksterid.clear();
    rhInfo.rechit_radius.clear();
  }

  void createRecHitBranches(TTree *fTree, recHitInfo &rhInfo) {
    fTree->Branch("rechit_eta", &rhInfo.rechit_eta);
    fTree->Branch("rechit_phi", &rhInfo.rechit_phi);
    fTree->Branch("rechit_pt", &rhInfo.rechit_pt);
    fTree->Branch("rechit_energy", &rhInfo.rechit_energy);
    fTree->Branch("rechit_recostructable_energy", &rhInfo.rechit_recostructable_energy);
    fTree->Branch("rechit_SoN", &rhInfo.rechit_SoN);
    fTree->Branch("rechit_uncalib_energy", &rhInfo.rechit_uncalib_energy);
    fTree->Branch("rechit_matbudget", &rhInfo.rechit_matbudget);
    fTree->Branch("rechit_x", &rhInfo.rechit_x);
    fTree->Branch("rechit_y", &rhInfo.rechit_y);
    fTree->Branch("rechit_z", &rhInfo.rechit_z);
    fTree->Branch("rechit_time", &rhInfo.rechit_time);
    fTree->Branch("rechit_thickness", &rhInfo.rechit_thickness);
    fTree->Branch("rechit_layer", &rhInfo.rechit_layer);
    fTree->Branch("rechit_wafer_u", &rhInfo.rechit_wafer_u);
    fTree->Branch("rechit_wafer_v", &rhInfo.rechit_wafer_v);
    fTree->Branch("rechit_cell_u", &rhInfo.rechit_cell_u);
    fTree->Branch("rechit_cell_v", &rhInfo.rechit_cell_v);
    fTree->Branch("rechit_detid", &rhInfo.rechit_detid);
    fTree->Branch("rechit_isHalf", &rhInfo.rechit_isHalf);
    fTree->Branch("rechit_flags", &rhInfo.rechit_flags);
    fTree->Branch("rechit_layerclusterid", &rhInfo.rechit_layerclusterid);
    fTree->Branch("rechit_tracksterid", &rhInfo.rechit_tracksterid);
    fTree->Branch("rechit_simclusterid", &rhInfo.rechit_simclusterid);
    fTree->Branch("rechit_radius", &rhInfo.rechit_radius);
  }

  //---------------------------------------------------------------------------------------------------
  //SimClusters
  //---------------------------------------------------------------------------------------------------
  struct simClustersInfo {

    // The PDG ID of the first associated gen particle. If there are no
    // gen particles associated then it returns type() from the first SimTrack.
    std::vector<int> simcluster_pdgId;
    // EncodedEventId is taken from the first SimTrack only, but there shouldn't be any
    // SimTracks from different crossings in the SimCluster. 
    std::vector<EncodedEventId> simcluster_eventId;
    // The index of the simCluster in the event collection
    std::vector<unsigned int> simcluster_id;
    //particleId
    std::vector<uint64_t> simcluster_particleId;
    // Electric charge. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_charge;
    // Magnitude of momentum vector. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_p;
    // Energy. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_energy;
    // Transverse energy. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_et;
    // Mass. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_mass;
    // Transverse mass. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_mt;
    // Transverse momentum. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_pt;
    // Momentum azimuthal angle. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_phi;
    // Momentum polar angle. Note this is taken from the first SimTrack only.
    std::vector<float> simcluster_theta;
    // Momentum pseudorapidity. Note this is taken from the simtrack before the calorimeter
    std::vector<float> simcluster_eta;
    // Rapidity. Note this is taken from the simtrack before the calorimeter
    std::vector<float> simcluster_rapidity;
    // Returns status() from the first gen particle, or -99 if there are no gen particles attached.
    std::vector<int> simcluster_status;
    // is long lived?
    std::vector<bool> simcluster_longLived;
    // Gives the total number of SimHits, in the cluster 
    std::vector<int> simcluster_numberOfSimHits;
    // returns the accumulated sim energy in the cluster from PCaloHit energy
    std::vector<float> simcluster_simEnergy;
    // simHit detids related to simCluster 
    std::vector<std::vector<uint32_t>> simcluster_hits;
    // Detector the simHit of the simcluster is related 
    std::vector<std::vector<int>> simcluster_hits_dets;
    // simCluster indices with at least one matched simhit with rechit via detid
    std::vector<std::vector<uint32_t>> simcluster_matched_hits;
    // simHits/recHits cell thickness
    std::vector<std::vector<int>> simcluster_hits_thickness;
    // fractions of the related matched recHits
    std::vector<std::vector<float>> simcluster_fractions;
    // layers the simcluster expandes to, could go up to 100 depending of geometry
    // Will use a set because there will be many duplicates per simCluster coming 
    // from the corresponding simhits
    std::vector<std::vector<unsigned int>> simcluster_layers;
    // wafer_u the simhit detid belongs too
    std::vector<std::vector<int>> simcluster_wafers_u;
    // wafer_v the simhit detid belongs too
    std::vector<std::vector<int>> simcluster_wafers_v;
    // cell_u the simhit detid belongs too
    std::vector<std::vector<int>> simcluster_cells_u;
    // cell_v the simhit detid belongs too
    std::vector<std::vector<int>> simcluster_cells_v;
    // cell type the simhit detid belongs too
    std::vector<std::vector<int>> simcluster_cells_type;
    // z side of the simhit
    std::vector<std::vector<int>> simcluster_cells_zside;
    // rechit info
    std::vector<std::vector<float>> simcluster_rechit_eta;
    std::vector<std::vector<float>> simcluster_rechit_phi;
    std::vector<std::vector<float>> simcluster_rechit_pt;
    std::vector<std::vector<float>> simcluster_rechit_energy;
    std::vector<std::vector<float>> simcluster_rechit_uncalib_energy;
    std::vector<std::vector<float>> simcluster_rechit_matbudget;
    std::vector<std::vector<float>> simcluster_rechit_recostructable_energy;
    std::vector<std::vector<float>> simcluster_rechit_SoN;
    std::vector<std::vector<float>> simcluster_rechit_x;
    std::vector<std::vector<float>> simcluster_rechit_y;
    std::vector<std::vector<float>> simcluster_rechit_z;
    std::vector<std::vector<float>> simcluster_rechit_time;
    std::vector<std::vector<int>> simcluster_rechit_simclusterid;


  };

  void initSimClustersInfo(simClustersInfo &scInfo) {
    scInfo.simcluster_pdgId.clear();
    scInfo.simcluster_eventId.clear();
    scInfo.simcluster_id.clear();
    scInfo.simcluster_particleId.clear();
    scInfo.simcluster_charge.clear();
    scInfo.simcluster_p.clear();
    scInfo.simcluster_energy.clear();
    scInfo.simcluster_et.clear();
    scInfo.simcluster_mass.clear();
    scInfo.simcluster_mt.clear();
    scInfo.simcluster_pt.clear();
    scInfo.simcluster_phi.clear();
    scInfo.simcluster_theta.clear();
    scInfo.simcluster_eta.clear();
    scInfo.simcluster_rapidity.clear();
    scInfo.simcluster_status.clear();
    scInfo.simcluster_longLived.clear();
    scInfo.simcluster_numberOfSimHits.clear();
    scInfo.simcluster_simEnergy.clear();
    scInfo.simcluster_hits.clear();
    scInfo.simcluster_hits_dets.clear();
    scInfo.simcluster_matched_hits.clear();
    scInfo.simcluster_hits_thickness.clear();
    scInfo.simcluster_fractions.clear();
    scInfo.simcluster_layers.clear();
    scInfo.simcluster_wafers_u.clear();
    scInfo.simcluster_wafers_v.clear();
    scInfo.simcluster_cells_u.clear();
    scInfo.simcluster_cells_v.clear();
    scInfo.simcluster_cells_type.clear();
    scInfo.simcluster_cells_zside.clear();
    scInfo.simcluster_rechit_eta.clear();
    scInfo.simcluster_rechit_phi.clear();
    scInfo.simcluster_rechit_pt.clear();
    scInfo.simcluster_rechit_energy.clear();
    scInfo.simcluster_rechit_uncalib_energy.clear();
    scInfo.simcluster_rechit_matbudget.clear();
    scInfo.simcluster_rechit_recostructable_energy.clear();
    scInfo.simcluster_rechit_SoN.clear();
    scInfo.simcluster_rechit_x.clear();
    scInfo.simcluster_rechit_y.clear();
    scInfo.simcluster_rechit_z.clear();
    scInfo.simcluster_rechit_time.clear();
    scInfo.simcluster_rechit_simclusterid.clear();


  }

  void createSimClustersBranches(TTree *fTree, simClustersInfo &scInfo) {
    fTree->Branch("simcluster_pdgId",&scInfo.simcluster_pdgId);
    fTree->Branch("simcluster_eventId",&scInfo.simcluster_eventId);
    fTree->Branch("simcluster_id",&scInfo.simcluster_id);
    fTree->Branch("simcluster_particleId",&scInfo.simcluster_particleId);
    fTree->Branch("simcluster_charge",&scInfo.simcluster_charge);
    fTree->Branch("simcluster_simcluster_p",&scInfo.simcluster_p);
    fTree->Branch("simcluster_energy",&scInfo.simcluster_energy);
    fTree->Branch("simcluster_et",&scInfo.simcluster_et);
    fTree->Branch("simcluster_mass",&scInfo.simcluster_mass);
    fTree->Branch("simcluster_mt",&scInfo.simcluster_mt);
    fTree->Branch("simcluster_pt",&scInfo.simcluster_pt);
    fTree->Branch("simcluster_phi",&scInfo.simcluster_phi);
    fTree->Branch("simcluster_theta",&scInfo.simcluster_theta);
    fTree->Branch("simcluster_eta",&scInfo.simcluster_eta);
    fTree->Branch("simcluster_rapidity",&scInfo.simcluster_rapidity);
    fTree->Branch("simcluster_status",&scInfo.simcluster_status);
    fTree->Branch("simcluster_longLived",&scInfo.simcluster_longLived);
    fTree->Branch("simcluster_numberOfSimHits",&scInfo.simcluster_numberOfSimHits);
    fTree->Branch("simcluster_simEnergy",&scInfo.simcluster_simEnergy);
    fTree->Branch("simcluster_hits",&scInfo.simcluster_hits);
    fTree->Branch("simcluster_hits_dets",&scInfo.simcluster_hits_dets);
    fTree->Branch("simcluster_matched_hits",&scInfo.simcluster_matched_hits);
    fTree->Branch("simcluster_hits_thickness",&scInfo.simcluster_hits_thickness);
    fTree->Branch("simcluster_fractions",&scInfo.simcluster_fractions);
    fTree->Branch("simcluster_layers",&scInfo.simcluster_layers);
    fTree->Branch("simcluster_wafers_u",&scInfo.simcluster_wafers_u);
    fTree->Branch("simcluster_wafers_v",&scInfo.simcluster_wafers_v);
    fTree->Branch("simcluster_cells_u",&scInfo.simcluster_cells_u);
    fTree->Branch("simcluster_cells_v",&scInfo.simcluster_cells_v);
    fTree->Branch("simcluster_cells_type",&scInfo.simcluster_cells_type);
    fTree->Branch("simcluster_cells_zside",&scInfo.simcluster_cells_zside);
    fTree->Branch("simcluster_rechit_eta", &scInfo.simcluster_rechit_eta);
    fTree->Branch("simcluster_rechit_phi", &scInfo.simcluster_rechit_phi);
    fTree->Branch("simcluster_rechit_pt", &scInfo.simcluster_rechit_pt);
    fTree->Branch("simcluster_rechit_energy", &scInfo.simcluster_rechit_energy);
    fTree->Branch("simcluster_rechit_uncalib_energy", &scInfo.simcluster_rechit_uncalib_energy);
    fTree->Branch("simcluster_rechit_matbudget", &scInfo.simcluster_rechit_matbudget);
    fTree->Branch("simcluster_rechit_recostructable_energy", &scInfo.simcluster_rechit_recostructable_energy);
    fTree->Branch("simcluster_rechit_SoN", &scInfo.simcluster_rechit_SoN);
    fTree->Branch("simcluster_rechit_x", &scInfo.simcluster_rechit_x);
    fTree->Branch("simcluster_rechit_y", &scInfo.simcluster_rechit_y);
    fTree->Branch("simcluster_rechit_z", &scInfo.simcluster_rechit_z);
    fTree->Branch("simcluster_rechit_time", &scInfo.simcluster_rechit_time);
    fTree->Branch("simcluster_rechit_simclusterid", &scInfo.simcluster_rechit_simclusterid);

  }

  void fillSimClustersInfo(simClustersInfo &scInfo, 
			   recHitInfo &rhInfo,
			   std::vector<SimCluster> const& simClusters,
			   std::unordered_map<DetId, std::vector<const HGCUncalibratedRecHit *> > const& UnCalibHitMap,
			   std::unordered_map<DetId, const HGCRecHit*> const& hitMap,
			   std::shared_ptr<hgcal::RecHitTools> recHitTools, 
			   unsigned int layers,
			   std::map<double, double> cummatbudg){


    std::vector<uint32_t> hits;
    std::vector<int> hits_dets;
    std::vector<uint32_t> matched_hits;
    std::vector<int> hits_thickness;
    std::vector<float> fractions;
    std::vector<unsigned int> simcluster_layers;
    std::vector<int> wafers_u;
    std::vector<int> wafers_v;
    std::vector<int> cells_u;
    std::vector<int> cells_v;
    std::vector<int> cells_type;
    std::vector<int> cells_zside;

    //loop through simClusters
    for (unsigned int ic = 0; ic < simClusters.size(); ++ic) {
      const auto& sc = simClusters[ic];
      const auto& hitsAndFractions = sc.hits_and_fractions();

      //Clear here once for each simCluster before looping in simHits. 
      hits.clear();
      hits_dets.clear();
      matched_hits.clear();
      hits_thickness.clear();
      fractions.clear();
      simcluster_layers.clear();
      wafers_u.clear();
      wafers_v.clear();
      cells_u.clear();
      cells_v.clear();
      cells_type.clear();
      cells_zside.clear();
      initRecHitInfo(rhInfo);

      //For the hits thickness of the simcluster.
      double thickness = 0.;

      //loop through hits of the simCluster
      for (const auto& hAndF : hitsAndFractions) {
	const DetId sh_detid = hAndF.first;

	//The layer the cluster belongs to. As mentioned in the mapping above, it takes into account -z and +z.
	int layerid =
	  recHitTools->getLayerWithOffset(sh_detid) + layers * ((recHitTools->zside(sh_detid) + 1) >> 1) - 1;

	//We keep the detid of all relevant simhits of simCluster
	hits.push_back(sh_detid.rawId());
	//std::cout << sh_detid.rawId() << std::endl;
	//We keep the detid.det, the detector of all relevant simhits of simCluster
	hits_dets.push_back(sh_detid.det());
	//We keep the fractions of all relevant simhits of simCluster
	fractions.push_back(hAndF.second);
	// The layers the current simCluster expands. We will save it per hit to compare and know where 
	// each hit is later on the analysis of the tree. 
	simcluster_layers.push_back(layerid);
 
	if (sh_detid.det() == DetId::Forward || sh_detid.det() == DetId::HGCalEE || sh_detid.det() == DetId::HGCalHSi) {
	  HGCSiliconDetId thedetId = HGCSiliconDetId(sh_detid);
	  thickness = recHitTools->getSiThickness(sh_detid);
	  std::pair<int, int> this_wafer = recHitTools->getWafer(sh_detid);
	  std::pair<int, int> this_cell = recHitTools->getCell(sh_detid);
	  wafers_u.push_back(this_wafer.first);
	  wafers_v.push_back(this_wafer.second);
	  cells_u.push_back(this_cell.first);
	  cells_v.push_back(this_cell.second);
	  cells_type.push_back(thedetId.type());
	  cells_zside.push_back(thedetId.zside());
	  /* std::cout << "layerid " << layerid << " sh_detid.layer() " << thedetId.layer() << " recHitTools->getLayerWithOffset(sh_detid) " << recHitTools->getLayerWithOffset(sh_detid) << std::endl; */
	} else {
	  HGCScintillatorDetId thedetId = HGCScintillatorDetId(sh_detid);
	  thickness = std::numeric_limits<std::float_t>::max();
	  wafers_u.push_back(std::numeric_limits<unsigned int>::max());
	  wafers_v.push_back(std::numeric_limits<unsigned int>::max());
	  cells_u.push_back(std::numeric_limits<unsigned int>::max());
	  cells_v.push_back(std::numeric_limits<unsigned int>::max());
	  cells_type.push_back(thedetId.type());
	  cells_zside.push_back(thedetId.zside());
	  /* std::cout << "layerid " << layerid << " sh_detid.layer() " << thedetId.layer() << " recHitTools->getLayerWithOffset(sh_detid) " << recHitTools->getLayerWithOffset(sh_detid) << std::endl; */
	}

	hits_thickness.push_back(thickness);

	//At this point we want to check if the simhit has a corresponding rechit and 
	//save all the relevant information.

	//Check if the current simHit has a related recHit via detid
	std::unordered_map<DetId, const HGCRecHit*>::const_iterator itcheck = hitMap.find(sh_detid);
	//If there is a corresponding recHit save the detid 
	if (itcheck != hitMap.end()) {

	  const HGCRecHit* hit = itcheck->second;

	  const GlobalPoint position = recHitTools->getPosition(sh_detid);
	  const double energy = hit->energy();
	  const double recostructable_energy = hit->energy() * hAndF.second;
	  const double SoN = hit->signalOverSigmaNoise();

	  // For the uncalibrated energy, which is in MIPs.
	  std::unordered_map<DetId, std::vector<const HGCUncalibratedRecHit *>>::const_iterator uncalib_check = UnCalibHitMap.find(sh_detid);
	  const double uncalib_energy = uncalib_check->second.front()->amplitude();


	  matched_hits.push_back(sh_detid.rawId());
	  //std::cout << "1 " << sh_detid.rawId() << std::endl;
	  /* std::cout << recHitTools->getLayerWithOffset(itcheck->first) + layers * ((recHitTools->zside(itcheck->first) + 1) >> 1) - 1 << std::endl; */
	  /* std::cout << "lastLayerEE() " << recHitTools->lastLayerEE() << " lastLayerFH() " << recHitTools->lastLayerFH() << " firstLayerBH() " << recHitTools->firstLayerBH() << " lastLayerBH() " << recHitTools->lastLayerBH() << " recHitTools->zside " << recHitTools->zside(itcheck->first) << std::endl; */

	  rhInfo.rechit_eta.push_back(recHitTools->getEta(position));
	  rhInfo.rechit_phi.push_back(recHitTools->getPhi(position));
	  rhInfo.rechit_pt.push_back(recHitTools->getPt(position,energy));
	  rhInfo.rechit_energy.push_back(energy);
	  rhInfo.rechit_uncalib_energy.push_back(uncalib_energy);
	  rhInfo.rechit_matbudget.push_back(cummatbudg[(double)layerid + 1.]);
	  rhInfo.rechit_recostructable_energy.push_back(recostructable_energy);
	  rhInfo.rechit_SoN.push_back(SoN);
	  rhInfo.rechit_x.push_back(position.x());
	  rhInfo.rechit_y.push_back(position.y());
	  rhInfo.rechit_z.push_back(position.z());
	  rhInfo.rechit_time.push_back(hit->time());
	  rhInfo.rechit_simclusterid.push_back(ic);

	} else{
	  //This simhit wasn't reconstructed. Save the info with as simclusterid -1, so that we know. 
	  //std::cout << "2 " << sh_detid.rawId() << std::endl;
	  matched_hits.push_back(9999);

	  rhInfo.rechit_eta.push_back(-9999.);
	  rhInfo.rechit_phi.push_back(-9999.);
	  rhInfo.rechit_pt.push_back(-9999.);
	  rhInfo.rechit_energy.push_back(-9999.);
	  rhInfo.rechit_uncalib_energy.push_back(-9999.);
	  rhInfo.rechit_matbudget.push_back(-9999.);
	  rhInfo.rechit_recostructable_energy.push_back(-9999.);
	  rhInfo.rechit_SoN.push_back(-9999.);
	  rhInfo.rechit_x.push_back(-9999.);
	  rhInfo.rechit_y.push_back(-9999.);
	  rhInfo.rechit_z.push_back(-9999.);
	  rhInfo.rechit_time.push_back(-9999.);
	  rhInfo.rechit_simclusterid.push_back(-1);
	  
	}




      }//end of loop through simhits of the current simCluster


      scInfo.simcluster_pdgId.push_back(sc.pdgId());
      scInfo.simcluster_eventId.push_back(sc.eventId());
      scInfo.simcluster_id.push_back(ic);
      scInfo.simcluster_particleId.push_back(sc.particleId());
      scInfo.simcluster_charge.push_back(sc.charge());
      scInfo.simcluster_p.push_back(sc.p());
      scInfo.simcluster_energy.push_back(sc.energy());
      scInfo.simcluster_et.push_back(sc.et());
      scInfo.simcluster_mass.push_back(sc.mass());
      scInfo.simcluster_mt.push_back(sc.mt());
      scInfo.simcluster_pt.push_back(sc.pt());
      scInfo.simcluster_phi.push_back(sc.phi());
      scInfo.simcluster_theta.push_back(sc.theta());
      scInfo.simcluster_eta.push_back(sc.eta());
      scInfo.simcluster_rapidity.push_back(sc.rapidity());
      scInfo.simcluster_status.push_back(sc.status());
      scInfo.simcluster_longLived.push_back(sc.longLived());
      scInfo.simcluster_numberOfSimHits.push_back(sc.numberOfSimHits());
      scInfo.simcluster_simEnergy.push_back(sc.simEnergy());
      scInfo.simcluster_hits.push_back(hits);
      scInfo.simcluster_hits_dets.push_back(hits_dets);
      scInfo.simcluster_hits_thickness.push_back(hits_thickness);
      scInfo.simcluster_matched_hits.push_back(matched_hits);
      scInfo.simcluster_fractions.push_back(fractions);
      scInfo.simcluster_layers.push_back(simcluster_layers);
      scInfo.simcluster_wafers_u.push_back(wafers_u);
      scInfo.simcluster_wafers_v.push_back(wafers_v);
      scInfo.simcluster_cells_u.push_back(cells_u);
      scInfo.simcluster_cells_v.push_back(cells_v);
      scInfo.simcluster_cells_type.push_back(cells_type);
      scInfo.simcluster_cells_zside.push_back(cells_zside);
      scInfo.simcluster_rechit_eta.push_back(rhInfo.rechit_eta);
      scInfo.simcluster_rechit_phi.push_back(rhInfo.rechit_phi);
      scInfo.simcluster_rechit_pt.push_back(rhInfo.rechit_pt);
      scInfo.simcluster_rechit_energy.push_back(rhInfo.rechit_energy);
      scInfo.simcluster_rechit_uncalib_energy.push_back(rhInfo.rechit_uncalib_energy);
      scInfo.simcluster_rechit_matbudget.push_back(rhInfo.rechit_matbudget);
      scInfo.simcluster_rechit_recostructable_energy.push_back(rhInfo.rechit_recostructable_energy);
      scInfo.simcluster_rechit_SoN.push_back(rhInfo.rechit_SoN);
      scInfo.simcluster_rechit_x.push_back(rhInfo.rechit_x);
      scInfo.simcluster_rechit_y.push_back(rhInfo.rechit_y);
      scInfo.simcluster_rechit_z.push_back(rhInfo.rechit_z);
      scInfo.simcluster_rechit_time.push_back(rhInfo.rechit_time);
      scInfo.simcluster_rechit_simclusterid.push_back(rhInfo.rechit_simclusterid);

    }// end of loop through simClusters


 

  }// end of fillSimClustersInfo

  //---------------------------------------------------------------------------------------------------
  //LayerClusters
  //---------------------------------------------------------------------------------------------------
  struct layerClustersInfo {
    std::vector<unsigned int> layerCluster_id;
    std::vector<unsigned int> layerCluster_trackster_id;
    std::vector<float> layerCluster_eta;
    std::vector<float> layerCluster_phi;
    std::vector<float> layerCluster_pt;
    std::vector<float> layerCluster_energy;
    std::vector<float> layerCluster_x;
    std::vector<float> layerCluster_y;
    std::vector<float> layerCluster_z;
    std::vector<int> layerCluster_layer;
    std::vector<int> layerCluster_nhitCore;
    std::vector<int> layerCluster_nhitAll;
    std::vector<float> layerCluster_matbudget;
    std::vector<std::vector<unsigned int>> layerCluster_rechits;
    std::vector<int> layerCluster_rechitSeed;
    std::vector<std::vector<float>> layerCluster_rechit_eta;
    std::vector<std::vector<float>> layerCluster_rechit_phi;
    std::vector<std::vector<float>> layerCluster_rechit_pt;
    std::vector<std::vector<float>> layerCluster_rechit_energy;
    std::vector<std::vector<float>> layerCluster_rechit_uncalib_energy;
    std::vector<std::vector<float>> layerCluster_rechit_x;
    std::vector<std::vector<float>> layerCluster_rechit_y;
    std::vector<std::vector<float>> layerCluster_rechit_z;
    std::vector<std::vector<float>> layerCluster_rechit_time;
    std::vector<std::vector<float>> layerCluster_rechit_thickness;
    std::vector<std::vector<int>> layerCluster_rechit_layer;
    std::vector<std::vector<int>> layerCluster_rechit_wafer_u;
    std::vector<std::vector<int>> layerCluster_rechit_wafer_v;
    std::vector<std::vector<int>> layerCluster_rechit_cell_u;
    std::vector<std::vector<int>> layerCluster_rechit_cell_v;
    std::vector<std::vector<unsigned int>> layerCluster_rechit_detid;
    std::vector<std::vector<bool>> layerCluster_rechit_isHalf;
    std::vector<std::vector<int>> layerCluster_rechit_flags;
    std::vector<std::vector<int>> layerCluster_rechit_layerclusterid;
    /* std::vector<std::vector<int>> layerCluster_rechit_simclusterid; */
    /* std::vector<std::vector<int>> layerCluster_rechit_caloparticleid; */
    std::vector<std::vector<int>> layerCluster_rechit_tracksterid; 
    /* std::vector<std::vector<int>> layerCluster_rechit_simtracksterid; */
    std::vector<std::vector<float>> layerCluster_rechit_radius;




  };

  void initLayerClustersInfo(layerClustersInfo &lcInfo) {
    lcInfo.layerCluster_id.clear();
    lcInfo.layerCluster_trackster_id.clear();
    lcInfo.layerCluster_eta.clear();
    lcInfo.layerCluster_phi.clear();
    lcInfo.layerCluster_pt.clear();
    lcInfo.layerCluster_energy.clear();
    lcInfo.layerCluster_x.clear();
    lcInfo.layerCluster_y.clear();
    lcInfo.layerCluster_z.clear();
    lcInfo.layerCluster_layer.clear();
    lcInfo.layerCluster_nhitCore.clear();
    lcInfo.layerCluster_nhitAll.clear();
    lcInfo.layerCluster_matbudget.clear();
    lcInfo.layerCluster_rechits.clear();
    lcInfo.layerCluster_rechitSeed.clear();
    lcInfo.layerCluster_rechit_eta.clear();
    lcInfo.layerCluster_rechit_phi.clear();
    lcInfo.layerCluster_rechit_pt.clear();
    lcInfo.layerCluster_rechit_energy.clear();
    lcInfo.layerCluster_rechit_uncalib_energy.clear();
    lcInfo.layerCluster_rechit_x.clear();
    lcInfo.layerCluster_rechit_y.clear();
    lcInfo.layerCluster_rechit_z.clear();
    lcInfo.layerCluster_rechit_time.clear();
    lcInfo.layerCluster_rechit_thickness.clear();
    lcInfo.layerCluster_rechit_layer.clear();
    lcInfo.layerCluster_rechit_wafer_u.clear();
    lcInfo.layerCluster_rechit_wafer_v.clear();
    lcInfo.layerCluster_rechit_cell_u.clear();
    lcInfo.layerCluster_rechit_cell_v.clear();
    lcInfo.layerCluster_rechit_detid.clear();
    lcInfo.layerCluster_rechit_isHalf.clear();
    lcInfo.layerCluster_rechit_flags.clear();
    lcInfo.layerCluster_rechit_layerclusterid.clear();
    /* lcInfo.layerCluster_rechit_simclusterid.clear(); */
    /* lcInfo.layerCluster_rechit_caloparticleid.clear(); */
    lcInfo.layerCluster_rechit_tracksterid.clear();
    /* lcInfo.layerCluster_rechit_simtracksterid.clear(); */
    lcInfo.layerCluster_rechit_radius.clear();



  }

  void createLayerClustersBranches(TTree *fTree, layerClustersInfo &lcInfo) {
    fTree->Branch("layerCluster_id", &lcInfo.layerCluster_id);
    fTree->Branch("layerCluster_trackster_id", &lcInfo.layerCluster_trackster_id);
    fTree->Branch("layerCluster_eta", &lcInfo.layerCluster_eta);
    fTree->Branch("layerCluster_phi", &lcInfo.layerCluster_phi);
    fTree->Branch("layerCluster_pt", &lcInfo.layerCluster_pt);
    fTree->Branch("layerCluster_energy", &lcInfo.layerCluster_energy);
    fTree->Branch("layerCluster_x", &lcInfo.layerCluster_x);
    fTree->Branch("layerCluster_y", &lcInfo.layerCluster_y);
    fTree->Branch("layerCluster_z", &lcInfo.layerCluster_z);
    fTree->Branch("layerCluster_layer", &lcInfo.layerCluster_layer);
    fTree->Branch("layerCluster_nhitCore", &lcInfo.layerCluster_nhitCore);
    fTree->Branch("layerCluster_nhitAll", &lcInfo.layerCluster_nhitAll);
    fTree->Branch("layerCluster_matbudget", &lcInfo.layerCluster_matbudget);
    fTree->Branch("layerCluster_rechits", &lcInfo.layerCluster_rechits);
    fTree->Branch("layerCluster_rechitSeed", &lcInfo.layerCluster_rechitSeed);
    fTree->Branch("layerCluster_rechit_eta", &lcInfo.layerCluster_rechit_eta);
    fTree->Branch("layerCluster_rechit_phi", &lcInfo.layerCluster_rechit_phi);
    fTree->Branch("layerCluster_rechit_pt", &lcInfo.layerCluster_rechit_pt);
    fTree->Branch("layerCluster_rechit_energy", &lcInfo.layerCluster_rechit_energy);
    fTree->Branch("layerCluster_rechit_uncalib_energy", &lcInfo.layerCluster_rechit_uncalib_energy);
    fTree->Branch("layerCluster_rechit_x", &lcInfo.layerCluster_rechit_x);
    fTree->Branch("layerCluster_rechit_y", &lcInfo.layerCluster_rechit_y);
    fTree->Branch("layerCluster_rechit_z", &lcInfo.layerCluster_rechit_z);
    fTree->Branch("layerCluster_rechit_time", &lcInfo.layerCluster_rechit_time);
    fTree->Branch("layerCluster_rechit_thickness", &lcInfo.layerCluster_rechit_thickness);
    fTree->Branch("layerCluster_rechit_layer", &lcInfo.layerCluster_rechit_layer);
    fTree->Branch("layerCluster_rechit_wafer_u", &lcInfo.layerCluster_rechit_wafer_u);
    fTree->Branch("layerCluster_rechit_wafer_v", &lcInfo.layerCluster_rechit_wafer_v);
    fTree->Branch("layerCluster_rechit_cell_u", &lcInfo.layerCluster_rechit_cell_u);
    fTree->Branch("layerCluster_rechit_cell_v", &lcInfo.layerCluster_rechit_cell_v);
    fTree->Branch("layerCluster_rechit_detid", &lcInfo.layerCluster_rechit_detid);
    fTree->Branch("layerCluster_rechit_isHalf", &lcInfo.layerCluster_rechit_isHalf);
    fTree->Branch("layerCluster_rechit_flags", &lcInfo.layerCluster_rechit_flags);
    fTree->Branch("layerCluster_rechit_layerclusterid", &lcInfo.layerCluster_rechit_layerclusterid);
    /* fTree->Branch("layerCluster_", &lcInfo.layerCluster_rechit_simclusterid); */
    /* fTree->Branch("layerCluster_", &lcInfo.layerCluster_rechit_caloparticleid); */
    fTree->Branch("layerCluster_rechit_tracksterid", &lcInfo.layerCluster_rechit_tracksterid);
    /* fTree->Branch("layerCluster_", &lcInfo.layerCluster_rechit_simtracksterid); */
    fTree->Branch("layerCluster_rechit_radius", &lcInfo.layerCluster_rechit_radius);


  }

  void fillLayerClustersInfo(layerClustersInfo &lcInfo, 
			     recHitInfo &rhInfo,
			     const reco::CaloClusterCollection& clusters,
			     unsigned int lcIdfromTrackster,
			     unsigned int trId,
			     bool fromTrackster,
			     std::unordered_map<DetId, std::vector<const HGCUncalibratedRecHit *> > const& UnCalibHitMap,
			     std::unordered_map<DetId, const HGCRecHit*> const& hitMap,
			     std::shared_ptr<hgcal::RecHitTools> recHitTools, 
			     unsigned int layers,
			     std::map<double, double> cummatbudg){

    
    auto nLayerClusters = clusters.size();
    std::vector<unsigned int> rhIndices;

    for (unsigned int lcId = 0; lcId < nLayerClusters; ++lcId) {
      if (fromTrackster && lcId != lcIdfromTrackster){continue;}
      const std::vector<std::pair<DetId, float>>& hits_and_fractions = clusters[lcId].hitsAndFractions();
      unsigned int numberOfHitsInLC = hits_and_fractions.size();
 
      const auto firstHitDetId = hits_and_fractions[0].first;
      int lcLayerId =
        recHitTools->getLayerWithOffset(firstHitDetId) + layers * ((recHitTools->zside(firstHitDetId) + 1) >> 1) - 1;

      //Counter of the hits in a layerCluster with fraction 1. 
      int ncoreHit = 0;
      //looking for the rechit seed for this layerCluster
      int rhSeed = 0;
      float maxEnergy = -1.;

      unsigned int rechit_index = 0;
      rhIndices.clear();
      initRecHitInfo(rhInfo);

      for (unsigned int hitId = 0; hitId < numberOfHitsInLC; hitId++) {
	DetId rh_detid = hits_and_fractions[hitId].first;
	const auto rhFraction = hits_and_fractions[hitId].second;

	std::unordered_map<DetId, const HGCRecHit*>::const_iterator itcheck = hitMap.find(rh_detid);
	const HGCRecHit* hit = itcheck->second;

	const GlobalPoint position = recHitTools->getPosition(rh_detid);
	const double energy = hit->energy();

	// For the uncalibrated energy, which is in MIPs.
 	std::unordered_map<DetId, std::vector<const HGCUncalibratedRecHit *>>::const_iterator uncalib_check = UnCalibHitMap.find(rh_detid);
	const double uncalib_energy = uncalib_check->second.front()->amplitude();
	
	//Want to check if each layer gives the correct layerid against the first hit that is attributed to 
	//the layerCluster id. 
	int rhLayerId =
          recHitTools->getLayerWithOffset(rh_detid) + layers * ((recHitTools->zside(rh_detid) + 1) >> 1) - 1;
  
	//Will continue with the same convention as before for hits flags. 
	int flags = 0x0;
	if (rhFraction > 0. && rhFraction < 1.){
	  flags = 0x1;
	} else if (rhFraction < 0.){
	  //We shouldn't be here in LC rechits. 
	  flags = 0x3;
	} else if (rhFraction == 0.){
	  flags = 0x2;
	}

	ncoreHit += int(rhFraction);
	if (hit->energy() > maxEnergy) {
          rhSeed = rechit_index;
          maxEnergy = hit->energy();
        }
        /* rhIndices.push_back(rechit_index); */

	std::pair<int, int> wafer;
	std::pair<int, int> cell;
	double thickness;

	if (rh_detid.det() == DetId::Forward || rh_detid.det() == DetId::HGCalEE || rh_detid.det() == DetId::HGCalHSi) {
	  thickness = recHitTools->getSiThickness(rh_detid);
	  wafer = recHitTools->getWafer(rh_detid);
	  cell = recHitTools->getCell(rh_detid);
	} else {
	  thickness = std::numeric_limits<std::float_t>::max();
	  wafer = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
	  cell = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
	}
	const bool isHalfCell = recHitTools->isHalfCell(rh_detid);
	const double radius = 
	  ((rh_detid.det() == DetId::Forward || rh_detid.det() == DetId::HGCalEE || rh_detid.det() == DetId::HGCalHSi) ? recHitTools->getRadiusToSide(rh_detid) : -1.);

	rhInfo.rechit_eta.push_back(recHitTools->getEta(position));
	rhInfo.rechit_phi.push_back(recHitTools->getPhi(position));
	rhInfo.rechit_pt.push_back(recHitTools->getPt(position,energy));
	rhInfo.rechit_energy.push_back(energy);
	rhInfo.rechit_uncalib_energy.push_back(uncalib_energy);
	rhInfo.rechit_layer.push_back(lcLayerId);
	rhInfo.rechit_wafer_u.push_back(wafer.first);
	rhInfo.rechit_wafer_v.push_back(wafer.second);
	rhInfo.rechit_cell_u.push_back(cell.first);
	rhInfo.rechit_cell_v.push_back(cell.second);
	rhInfo.rechit_detid.push_back(rh_detid);
	rhInfo.rechit_x.push_back(position.x());
	rhInfo.rechit_y.push_back(position.y());
	rhInfo.rechit_z.push_back(position.z());
	rhInfo.rechit_time.push_back(hit->time());
	rhInfo.rechit_thickness.push_back(thickness);
	rhInfo.rechit_isHalf.push_back(isHalfCell);
	rhInfo.rechit_flags.push_back(flags);
	rhInfo.rechit_radius.push_back(radius);
	rhInfo.rechit_layerclusterid.push_back(rhLayerId);
	rhInfo.rechit_tracksterid.push_back(trId);

	++rechit_index;

      }  // end loop over hits on a LayerCluster

      double pt = clusters[lcId].energy() / cosh(clusters[lcId].eta());
      lcInfo.layerCluster_id.push_back(lcId);
      lcInfo.layerCluster_trackster_id.push_back(trId);
      lcInfo.layerCluster_eta.push_back(clusters[lcId].eta());
      lcInfo.layerCluster_phi.push_back(clusters[lcId].phi());
      lcInfo.layerCluster_pt.push_back(pt);
      lcInfo.layerCluster_energy.push_back(clusters[lcId].energy());
      lcInfo.layerCluster_x.push_back(clusters[lcId].x());
      lcInfo.layerCluster_y.push_back(clusters[lcId].y());
      lcInfo.layerCluster_z.push_back(clusters[lcId].z());
      lcInfo.layerCluster_layer.push_back(lcLayerId);
      lcInfo.layerCluster_nhitCore.push_back(ncoreHit);
      lcInfo.layerCluster_nhitAll.push_back(numberOfHitsInLC);
      //lcLayerId starts from 0, while in the input file we start from 1. 
      lcInfo.layerCluster_matbudget.push_back(cummatbudg[(double)lcLayerId + 1.]);
      lcInfo.layerCluster_rechitSeed.push_back(rhSeed);
      lcInfo.layerCluster_rechits.push_back(rhIndices);
      lcInfo.layerCluster_rechit_eta.push_back(rhInfo.rechit_eta);
      lcInfo.layerCluster_rechit_phi.push_back(rhInfo.rechit_phi);
      lcInfo.layerCluster_rechit_pt.push_back(rhInfo.rechit_pt);
      lcInfo.layerCluster_rechit_energy.push_back(rhInfo.rechit_energy);
      lcInfo.layerCluster_rechit_uncalib_energy.push_back(rhInfo.rechit_uncalib_energy);
      lcInfo.layerCluster_rechit_layer.push_back(rhInfo.rechit_layer);
      lcInfo.layerCluster_rechit_wafer_u.push_back(rhInfo.rechit_wafer_u);
      lcInfo.layerCluster_rechit_wafer_v.push_back(rhInfo.rechit_wafer_v);
      lcInfo.layerCluster_rechit_cell_u.push_back(rhInfo.rechit_cell_u);
      lcInfo.layerCluster_rechit_cell_v.push_back(rhInfo.rechit_cell_v);
      lcInfo.layerCluster_rechit_detid.push_back(rhInfo.rechit_detid);
      lcInfo.layerCluster_rechit_x.push_back(rhInfo.rechit_x);
      lcInfo.layerCluster_rechit_y.push_back(rhInfo.rechit_y);
      lcInfo.layerCluster_rechit_z.push_back(rhInfo.rechit_z);
      lcInfo.layerCluster_rechit_time.push_back(rhInfo.rechit_time);
      lcInfo.layerCluster_rechit_thickness.push_back(rhInfo.rechit_thickness);
      lcInfo.layerCluster_rechit_isHalf.push_back(rhInfo.rechit_isHalf);
      lcInfo.layerCluster_rechit_flags.push_back(rhInfo.rechit_flags);
      lcInfo.layerCluster_rechit_radius.push_back(rhInfo.rechit_radius);
      lcInfo.layerCluster_rechit_layerclusterid.push_back(rhInfo.rechit_layerclusterid);
      lcInfo.layerCluster_rechit_tracksterid.push_back(rhInfo.rechit_tracksterid);

    }//end of loop over layerClusters

  }// end of fillLayerClustersInfo


  //---------------------------------------------------------------------------------------------------
  //CaloParticles 
  //---------------------------------------------------------------------------------------------------
  struct caloParticlesInfo {

    // The PDG ID of the first associated gen particle. If there are no
    // gen particles associated then it returns type() from the first SimTrack.
    std::vector<int> caloparticle_pdgId;
    // EncodedEventId is taken from the first SimTrack only, but there shouldn't be any
    // SimTracks from different crossings in the SimCluster. 
    std::vector<EncodedEventId> caloparticle_eventId;
    // The index of the simCluster in the event collection
    std::vector<unsigned int> caloparticle_id;
    //particleId
    std::vector<uint64_t> caloparticle_particleId;
    // Electric charge. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_charge;
    // Magnitude of momentum vector. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_p;
    // Energy. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_energy;
    // Transverse energy. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_et;
    // Mass. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_mass;
    // Transverse mass. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_mt;
    // Transverse momentum. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_pt;
    // Momentum azimuthal angle. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_phi;
    // Momentum polar angle. Note this is taken from the first SimTrack only.
    std::vector<float> caloparticle_theta;
    // Momentum pseudorapidity. Note this is taken from the simtrack before the calorimeter
    std::vector<float> caloparticle_eta;
    // Rapidity. Note this is taken from the simtrack before the calorimeter
    std::vector<float> caloparticle_rapidity;
    // Returns status() from the first gen particle, or -99 if there are no gen particles attached.
    std::vector<int> caloparticle_status;
    // is long lived?
    std::vector<bool> caloparticle_longLived;
    // Gives the total number of SimHits, in the cluster 
    std::vector<int> caloparticle_numberOfSimHits;
    // returns the accumulated sim energy in the cluster from PCaloHit energy
    std::vector<float> caloparticle_simEnergy;
     // returns the reconstructable energy of the CaloParticle
    std::vector<float> caloparticle_sumEnergy;
   // simHit detids related to simCluster 
    std::vector<std::vector<uint32_t>> caloparticle_hits;
    // simCluster indices with at least one matched simhit with rechit via detid
    std::vector<std::vector<uint32_t>> caloparticle_matched_hits;
    // simHits/recHits cell thickness
    std::vector<std::vector<int>> caloparticle_hits_thickness;
    // fractions of the related matched recHits
    std::vector<std::vector<float>> caloparticle_fractions;
    // layers the simcluster expandes to, could go up to 100 depending of geometry
    // Will use a set because there will be many duplicates per simCluster coming 
    // from the corresponding simhits
    std::vector<std::vector<unsigned int>> caloparticle_layers;
    // wafer_u the matched rechit detid belongs too
    std::vector<std::vector<int>> caloparticle_wafers_u;
    // wafer_v the matched rechit detid belongs too
    std::vector<std::vector<int>> caloparticle_wafers_v;
    // cell_u the matched rechit detid belongs too
    std::vector<std::vector<int>> caloparticle_cells_u;
    // cell_v the matched rechit detid belongs too
    std::vector<std::vector<int>> caloparticle_cells_v;
  
  };

  void initCaloParticlesInfo(caloParticlesInfo &cpInfo) {
    cpInfo.caloparticle_pdgId.clear();
    cpInfo.caloparticle_eventId.clear();
    cpInfo.caloparticle_id.clear();
    cpInfo.caloparticle_particleId.clear();
    cpInfo.caloparticle_charge.clear();
    cpInfo.caloparticle_p.clear();
    cpInfo.caloparticle_energy.clear();
    cpInfo.caloparticle_et.clear();
    cpInfo.caloparticle_mass.clear();
    cpInfo.caloparticle_mt.clear();
    cpInfo.caloparticle_pt.clear();
    cpInfo.caloparticle_phi.clear();
    cpInfo.caloparticle_theta.clear();
    cpInfo.caloparticle_eta.clear();
    cpInfo.caloparticle_rapidity.clear();
    cpInfo.caloparticle_status.clear();
    cpInfo.caloparticle_longLived.clear();
    cpInfo.caloparticle_numberOfSimHits.clear();
    cpInfo.caloparticle_simEnergy.clear();
    cpInfo.caloparticle_sumEnergy.clear();
    cpInfo.caloparticle_hits.clear();
    cpInfo.caloparticle_matched_hits.clear();
    cpInfo.caloparticle_hits_thickness.clear();
    cpInfo.caloparticle_fractions.clear();
    cpInfo.caloparticle_layers.clear();
    cpInfo.caloparticle_wafers_u.clear();
    cpInfo.caloparticle_wafers_v.clear();
    cpInfo.caloparticle_cells_u.clear();
    cpInfo.caloparticle_cells_v.clear();

    }

  void createCaloParticlesBranches(TTree *fTree, caloParticlesInfo &cpInfo, simClustersInfo &scInfo) {
    fTree->Branch("caloparticle_pdgId",&cpInfo.caloparticle_pdgId);
    fTree->Branch("caloparticle_eventId",&cpInfo.caloparticle_eventId);
    fTree->Branch("caloparticle_id",&cpInfo.caloparticle_id);
    fTree->Branch("caloparticle_particleId",&cpInfo.caloparticle_particleId);
    fTree->Branch("caloparticle_charge",&cpInfo.caloparticle_charge);
    fTree->Branch("caloparticle_caloparticle_p",&cpInfo.caloparticle_p);
    fTree->Branch("caloparticle_energy",&cpInfo.caloparticle_energy);
    fTree->Branch("caloparticle_et",&cpInfo.caloparticle_et);
    fTree->Branch("caloparticle_mass",&cpInfo.caloparticle_mass);
    fTree->Branch("caloparticle_mt",&cpInfo.caloparticle_mt);
    fTree->Branch("caloparticle_pt",&cpInfo.caloparticle_pt);
    fTree->Branch("caloparticle_phi",&cpInfo.caloparticle_phi);
    fTree->Branch("caloparticle_theta",&cpInfo.caloparticle_theta);
    fTree->Branch("caloparticle_eta",&cpInfo.caloparticle_eta);
    fTree->Branch("caloparticle_rapidity",&cpInfo.caloparticle_rapidity);
    fTree->Branch("caloparticle_status",&cpInfo.caloparticle_status);
    fTree->Branch("caloparticle_longLived",&cpInfo.caloparticle_longLived);
    fTree->Branch("caloparticle_numberOfSimHits",&cpInfo.caloparticle_numberOfSimHits);
    fTree->Branch("caloparticle_simEnergy",&cpInfo.caloparticle_simEnergy);
    fTree->Branch("caloparticle_sumEnergy",&cpInfo.caloparticle_sumEnergy);
    fTree->Branch("caloparticle_hits",&cpInfo.caloparticle_hits);
    fTree->Branch("caloparticle_matched_hits",&cpInfo.caloparticle_matched_hits);
    fTree->Branch("caloparticle_hits_thickness",&cpInfo.caloparticle_hits_thickness);
    fTree->Branch("caloparticle_fractions",&cpInfo.caloparticle_fractions);
    fTree->Branch("caloparticle_layers",&cpInfo.caloparticle_layers);
    fTree->Branch("caloparticle_wafers_u",&cpInfo.caloparticle_wafers_u);
    fTree->Branch("caloparticle_wafers_v",&cpInfo.caloparticle_wafers_v);
    fTree->Branch("caloparticle_cells_u",&cpInfo.caloparticle_cells_u);
    fTree->Branch("caloparticle_cells_v",&cpInfo.caloparticle_cells_v);
    fTree->Branch("caloparticle_simcluster_pdgId",&scInfo.simcluster_pdgId);
    fTree->Branch("caloparticle_simcluster_id",&scInfo.simcluster_id);
    fTree->Branch("caloparticle_simcluster_particleId",&scInfo.simcluster_particleId);
    fTree->Branch("caloparticle_simcluster_charge",&scInfo.simcluster_charge);
    fTree->Branch("caloparticle_simcluster_simcluster_p",&scInfo.simcluster_p);
    fTree->Branch("caloparticle_simcluster_energy",&scInfo.simcluster_energy);
    fTree->Branch("caloparticle_simcluster_et",&scInfo.simcluster_et);
    fTree->Branch("caloparticle_simcluster_mass",&scInfo.simcluster_mass);
    fTree->Branch("caloparticle_simcluster_mt",&scInfo.simcluster_mt);
    fTree->Branch("caloparticle_simcluster_pt",&scInfo.simcluster_pt);
    fTree->Branch("caloparticle_simcluster_phi",&scInfo.simcluster_phi);
    fTree->Branch("caloparticle_simcluster_theta",&scInfo.simcluster_theta);
    fTree->Branch("caloparticle_simcluster_eta",&scInfo.simcluster_eta);
    fTree->Branch("caloparticle_simcluster_rapidity",&scInfo.simcluster_rapidity);
    fTree->Branch("caloparticle_simcluster_status",&scInfo.simcluster_status);
    fTree->Branch("caloparticle_simcluster_longLived",&scInfo.simcluster_longLived);
    fTree->Branch("caloparticle_simcluster_numberOfSimHits",&scInfo.simcluster_numberOfSimHits);
    fTree->Branch("caloparticle_simcluster_simEnergy",&scInfo.simcluster_simEnergy);
 
  }

   void fillCaloParticlesInfo(caloParticlesInfo &cpInfo,
			      simClustersInfo &scInfo, 
			      std::vector<CaloParticle> const& cP,
			      std::vector<size_t> const& cPSelectedIndices,
			      std::unordered_map<DetId, const HGCRecHit*> const& hitMap,
			      std::shared_ptr<hgcal::RecHitTools> recHitTools, 
			      unsigned int layers
			      ) {

    for (const auto& cpId : cPSelectedIndices) {

      std::map<unsigned int, float> cPEnergyOnLayer;
      for (unsigned int layerId = 0; layerId < layers * 2; ++layerId) {cPEnergyOnLayer[layerId] = 0;}

      std::vector<uint32_t> hits;
      std::vector<uint32_t> matched_hits;
      std::vector<int> hits_thickness;
      std::vector<float> fractions;
      std::vector<unsigned int> caloparticle_layers;
      std::vector<int> wafers_u;
      std::vector<int> wafers_v;
      std::vector<int> cells_u;
      std::vector<int> cells_v;

      //Clear here once for each caloParticle before looping in simClusters and simHits. 
      hits.clear();
      matched_hits.clear();
      hits_thickness.clear();
      fractions.clear();
      caloparticle_layers.clear();
      wafers_u.clear();
      wafers_v.clear();
      cells_u.clear();
      cells_v.clear();
      
      const SimClusterRefVector& simClusterRefVector = cP[cpId].simClusters();
      //auto nSimClusters = cP[cpId].simClusters().size();
      int countSimClusters = 0;
      for (const auto& it_sc : simClusterRefVector) {
	const SimCluster& simCluster = (*(it_sc));
	const auto& hits_and_fractions = simCluster.hits_and_fractions();
	//For the hits thickness of the current simCluster of CaloParticle
	double thickness = 0.;

  	for (const auto& it_haf : hits_and_fractions) {
	  const DetId hitid = (it_haf.first);
	  const int cpLayerId =
	    recHitTools->getLayerWithOffset(hitid) + layers * ((recHitTools->zside(hitid) + 1) >> 1) - 1;
	  std::unordered_map<DetId, const HGCRecHit*>::const_iterator itcheck = hitMap.find(hitid);
	  if (itcheck != hitMap.end()) {
	    const HGCRecHit* hit = itcheck->second;
	    cPEnergyOnLayer[cpLayerId] += it_haf.second * hit->energy();
	    //If there is a corresponding recHit save the detid 
	    matched_hits.push_back(hitid);
	  }

	  //We keep the detid of all relevant simhits of simClusters of the current caloParticle
	  hits.push_back(hitid);
	  //We keep the fractions of all relevant simhits of simClusters of the current caloParticle
	  fractions.push_back(it_haf.second);
	  // The layers the current caloParticle expands. We will save it per hit to compare and know where 
	  // each hit is later on the analysis of the tree. 
	  caloparticle_layers.push_back(cpLayerId);
 
	  if (hitid.det() == DetId::Forward || hitid.det() == DetId::HGCalEE || hitid.det() == DetId::HGCalHSi) {
	    thickness = recHitTools->getSiThickness(hitid);
	    std::pair<int, int> this_wafer = recHitTools->getWafer(hitid);
	    std::pair<int, int> this_cell = recHitTools->getCell(hitid);
	    wafers_u.push_back(this_wafer.first);
	    wafers_v.push_back(this_wafer.second);
	    cells_u.push_back(this_cell.first);
	    cells_v.push_back(this_cell.second);
	  } else {
	    thickness = std::numeric_limits<std::float_t>::max();
	    wafers_u.push_back(std::numeric_limits<unsigned int>::max());
	    wafers_v.push_back(std::numeric_limits<unsigned int>::max());
	    cells_u.push_back(std::numeric_limits<unsigned int>::max());
	    cells_v.push_back(std::numeric_limits<unsigned int>::max());
	  }

	  hits_thickness.push_back(thickness);
	} //end of loop over hits of the simCluster of the CaloParticle

	scInfo.simcluster_pdgId.push_back(simCluster.pdgId());
	scInfo.simcluster_id.push_back(countSimClusters);
	scInfo.simcluster_particleId.push_back(simCluster.particleId());
	scInfo.simcluster_charge.push_back(simCluster.charge());
	scInfo.simcluster_p.push_back(simCluster.p());
	scInfo.simcluster_energy.push_back(simCluster.energy());
	scInfo.simcluster_et.push_back(simCluster.et());
	scInfo.simcluster_mass.push_back(simCluster.mass());
	scInfo.simcluster_mt.push_back(simCluster.mt());
	scInfo.simcluster_pt.push_back(simCluster.pt());
	scInfo.simcluster_phi.push_back(simCluster.phi());
	scInfo.simcluster_theta.push_back(simCluster.theta());
	scInfo.simcluster_eta.push_back(simCluster.eta());
	scInfo.simcluster_rapidity.push_back(simCluster.rapidity());
	scInfo.simcluster_status.push_back(simCluster.status());
	scInfo.simcluster_longLived.push_back(simCluster.longLived());
	scInfo.simcluster_numberOfSimHits.push_back(simCluster.numberOfSimHits());
	scInfo.simcluster_simEnergy.push_back(simCluster.simEnergy());

	++countSimClusters;
      }//end of loop over simClusters of the CaloParticle

      //Calculate sum energy per-layer
      auto i = cPEnergyOnLayer.begin();
      double sum_energy = 0.0;
      while (i != cPEnergyOnLayer.end()) {
	sum_energy += i->second;
	i++;
      }

      cpInfo.caloparticle_pdgId.push_back(cP[cpId].pdgId());
      cpInfo.caloparticle_eventId.push_back(cP[cpId].eventId());
      cpInfo.caloparticle_id.push_back(cpId);
      cpInfo.caloparticle_particleId.push_back(cP[cpId].particleId());
      cpInfo.caloparticle_charge.push_back(cP[cpId].charge());
      cpInfo.caloparticle_p.push_back(cP[cpId].p());
      cpInfo.caloparticle_energy.push_back(cP[cpId].energy());
      cpInfo.caloparticle_et.push_back(cP[cpId].et());
      cpInfo.caloparticle_mass.push_back(cP[cpId].mass());
      cpInfo.caloparticle_mt.push_back(cP[cpId].mt());
      cpInfo.caloparticle_pt.push_back(cP[cpId].pt());
      cpInfo.caloparticle_phi.push_back(cP[cpId].phi());
      cpInfo.caloparticle_theta.push_back(cP[cpId].theta());
      cpInfo.caloparticle_eta.push_back(cP[cpId].eta());
      cpInfo.caloparticle_rapidity.push_back(cP[cpId].rapidity());
      cpInfo.caloparticle_status.push_back(cP[cpId].status());
      cpInfo.caloparticle_longLived.push_back(cP[cpId].longLived());
      cpInfo.caloparticle_numberOfSimHits.push_back(cP[cpId].numberOfSimHits());
      cpInfo.caloparticle_simEnergy.push_back(cP[cpId].simEnergy());
      cpInfo.caloparticle_sumEnergy.push_back(sum_energy);
      cpInfo.caloparticle_hits.push_back(hits);
      cpInfo.caloparticle_hits_thickness.push_back(hits_thickness);
      cpInfo.caloparticle_matched_hits.push_back(matched_hits);
      cpInfo.caloparticle_fractions.push_back(fractions);
      cpInfo.caloparticle_layers.push_back(caloparticle_layers);
      cpInfo.caloparticle_wafers_u.push_back(wafers_u);
      cpInfo.caloparticle_wafers_v.push_back(wafers_v);
      cpInfo.caloparticle_cells_u.push_back(cells_u);
      cpInfo.caloparticle_cells_v.push_back(cells_v);


    }//end of loop over caloParticles

    
   }

  //---------------------------------------------------------------------------------------------------
  //Tracksters
  //---------------------------------------------------------------------------------------------------
   struct trackstersInfo {
 
     std::vector<float> trackster_x;
     std::vector<float> trackster_y;
     std::vector<float> trackster_z;
     std::vector<float> trackster_eta;
     std::vector<float> trackster_phi;
     std::vector<unsigned int> trackster_firstlayer;
     std::vector<unsigned int> trackster_lastlayer;
     std::vector<unsigned int> trackster_layersnum;
     std::vector<int> trackster_id;
     std::vector<int> trackster_seedIndex;
     std::vector<float> trackster_time;
     std::vector<float> trackster_timeError;
     std::vector<float> trackster_regr_energy;
     std::vector<float> trackster_raw_energy;
     std::vector<float> trackster_raw_em_energy;
     std::vector<float> trackster_raw_pt;
     std::vector<float> trackster_raw_em_pt;
     std::vector<float> trackster_id_prob_1;
     std::vector<float> trackster_id_prob_2;
     std::vector<float> trackster_id_prob_3;
     std::vector<float> trackster_id_prob_4;
     std::vector<float> trackster_id_prob_5;
     std::vector<float> trackster_id_prob_6;
     std::vector<float> trackster_id_prob_7;
     std::vector<float> trackster_id_prob_8;
     std::vector<std::vector<unsigned int> > trackster_numedges;
     std::vector<std::vector<unsigned int> > trackster_edge_layerin_id;
     std::vector<std::vector<unsigned int> > trackster_edge_layerout_id;
     std::vector<std::vector<unsigned int> > trackster_edge_layerin;
     std::vector<std::vector<unsigned int> > trackster_edge_layerout;
     std::vector<std::vector<float> > trackster_delta_energy;
     std::vector<std::vector<float> > trackster_delta_energy_relative;
     std::vector<std::vector<float> > trackster_delta_layer;
     std::vector<std::vector<float> > trackster_angle_alpha;
     std::vector<std::vector<float> > trackster_angle_alpha_alternative;
     std::vector<std::vector<float> > trackster_angle_beta;

     std::vector<std::vector<float> > trackster_Raw_Energy;
     std::vector<std::vector<unsigned int> > trackster_numberOfHitsInTS;
     std::vector<std::vector<unsigned int> > trackster_numberOfNoiseHitsInTS;
     std::vector<std::vector<int> > trackster_maxCPId_byNumberOfHits;
     std::vector<std::vector<unsigned int> > trackster_maxCPNumberOfHitsInTS;
     std::vector<std::vector<int> > trackster_maxCPId_byEnergy;
     std::vector<std::vector<float> > trackster_maxEnergySharedTSandCP;
     std::vector<std::vector<float> > trackster_totalCPEnergyFromLayerCP;
     std::vector<std::vector<float> > trackster_energyFractionOfTSinCP;
     std::vector<std::vector<float> > trackster_energyFractionOfCPinTS;

     std::vector<std::vector<unsigned int> > trackster_cpId;
     std::vector<std::vector<unsigned int> > trackster_scId;
     std::vector<std::vector<unsigned int> > trackster_Id;
     std::vector<std::vector<unsigned int> > trackster_numofvertices;
     std::vector<std::vector<float> > trackster_score_trackster2caloparticle;
     std::vector<std::vector<float> > trackster_sharedenergy_trackster2caloparticle;
     std::vector<std::vector<float> > trackster_score_trackster2bestCaloparticle;
     std::vector<std::vector<float> > trackster_sharedenergy_trackster2bestCaloparticle;
     std::vector<std::vector<float> > trackster_trackster2bestCaloparticle_eta;
     std::vector<std::vector<float> > trackster_trackster2bestCaloparticle_phi;
     std::vector<std::vector<float> > trackster_score_trackster2bestCaloparticle2;
     std::vector<std::vector<float> > trackster_sharedenergy_trackster2bestCaloparticle2;
     
     std::vector<std::vector<unsigned int> > trackster_sts_cpId;
     std::vector<std::vector<unsigned int> > trackster_sts_id;
     std::vector<std::vector<unsigned int> > trackster_sts_ts_id;
     
     std::vector<std::vector<float> > trackster_sts_SimEnergy;
     std::vector<std::vector<float> > trackster_sts_SimEnergyWeight;
     std::vector<std::vector<float> > trackster_sts_trackster_raw_energy;
     std::vector<std::vector<float> > trackster_sts_score_caloparticle2trackster;
     std::vector<std::vector<float> > trackster_sts_sharedenergy_caloparticle2trackster;
     std::vector<std::vector<float> > trackster_sts_eta;
     std::vector<std::vector<float> > trackster_sts_phi;
     std::vector<std::vector<float> > trackster_sts_pt;
     
     std::vector<std::vector<float> > trackster_sts_raw_energy;
     std::vector<std::vector<float> > trackster_sts_scorePur_caloparticle2trackster;
     std::vector<std::vector<float> > trackster_sts_sharedenergy_caloparticle2trackster_assoc;
     std::vector<std::vector<float> > trackster_sts_besttrackster_raw_energy;
     std::vector<std::vector<float> > trackster_sts_scoreDupl_caloparticle2trackster;
     std::vector<std::vector<float> > trackster_sts_sharedenergy_caloparticle2trackster_assoc2;

     
   };

   void initTrackstersInfo(trackstersInfo &trInfo) {
     trInfo.trackster_x.clear();
     trInfo.trackster_y.clear();
     trInfo.trackster_z.clear();
     trInfo.trackster_eta.clear();
     trInfo.trackster_phi.clear();
     trInfo.trackster_firstlayer.clear();
     trInfo.trackster_lastlayer.clear();
     trInfo.trackster_layersnum.clear();
     trInfo.trackster_id.clear();
     trInfo.trackster_seedIndex.clear();
     trInfo.trackster_time.clear();
     trInfo.trackster_timeError.clear();
     trInfo.trackster_regr_energy.clear();
     trInfo.trackster_raw_energy.clear();
     trInfo.trackster_raw_em_energy.clear();
     trInfo.trackster_raw_pt.clear();
     trInfo.trackster_raw_em_pt.clear();
     trInfo.trackster_id_prob_1.clear();
     trInfo.trackster_id_prob_2.clear();
     trInfo.trackster_id_prob_3.clear();
     trInfo.trackster_id_prob_4.clear();
     trInfo.trackster_id_prob_5.clear();
     trInfo.trackster_id_prob_6.clear();
     trInfo.trackster_id_prob_7.clear();
     trInfo.trackster_id_prob_8.clear();
     trInfo.trackster_numedges.clear();
     trInfo.trackster_edge_layerin_id.clear();
     trInfo.trackster_edge_layerout_id.clear();
     trInfo.trackster_edge_layerin.clear();
     trInfo.trackster_edge_layerout.clear();
     trInfo.trackster_delta_energy.clear();
     trInfo.trackster_delta_energy_relative.clear();
     trInfo.trackster_delta_layer.clear();
     trInfo.trackster_angle_alpha.clear();
     trInfo.trackster_angle_alpha_alternative.clear();
     trInfo.trackster_angle_beta.clear();
     
     trInfo.trackster_Raw_Energy.clear();
     trInfo.trackster_numberOfHitsInTS.clear();
     trInfo.trackster_numberOfNoiseHitsInTS.clear();
     trInfo.trackster_maxCPId_byNumberOfHits.clear();
     trInfo.trackster_maxCPNumberOfHitsInTS.clear();
     trInfo.trackster_maxCPId_byEnergy.clear();
     trInfo.trackster_maxEnergySharedTSandCP.clear();
     trInfo.trackster_totalCPEnergyFromLayerCP.clear();
     trInfo.trackster_energyFractionOfTSinCP.clear();
     trInfo.trackster_energyFractionOfCPinTS.clear();

     trInfo.trackster_cpId.clear();
     trInfo.trackster_scId.clear();
     trInfo.trackster_Id.clear();
     trInfo.trackster_numofvertices.clear();
     trInfo.trackster_score_trackster2caloparticle.clear();
     trInfo.trackster_sharedenergy_trackster2caloparticle.clear();
     trInfo.trackster_score_trackster2bestCaloparticle.clear();
     trInfo.trackster_sharedenergy_trackster2bestCaloparticle.clear();
     trInfo.trackster_trackster2bestCaloparticle_eta.clear();
     trInfo.trackster_trackster2bestCaloparticle_phi.clear();
     trInfo.trackster_score_trackster2bestCaloparticle2.clear();
     trInfo.trackster_sharedenergy_trackster2bestCaloparticle2.clear();

     trInfo.trackster_sts_cpId.clear();
     trInfo.trackster_sts_id.clear();
     trInfo.trackster_sts_ts_id.clear();
     trInfo.trackster_sts_SimEnergy.clear();
     trInfo.trackster_sts_SimEnergyWeight.clear();
     trInfo.trackster_sts_trackster_raw_energy.clear();
     trInfo.trackster_sts_score_caloparticle2trackster.clear();
     trInfo.trackster_sts_sharedenergy_caloparticle2trackster.clear();
     trInfo.trackster_sts_eta.clear();
     trInfo.trackster_sts_phi.clear();
     trInfo.trackster_sts_pt.clear();
     trInfo.trackster_sts_raw_energy.clear();
     trInfo.trackster_sts_scorePur_caloparticle2trackster.clear();
     trInfo.trackster_sts_sharedenergy_caloparticle2trackster_assoc.clear();
     trInfo.trackster_sts_besttrackster_raw_energy.clear();
     trInfo.trackster_sts_scoreDupl_caloparticle2trackster.clear();
     trInfo.trackster_sts_sharedenergy_caloparticle2trackster_assoc2.clear();


     

   }

  void createTrackstersBranches(TTree *fTree, trackstersInfo &trInfo) {
    fTree->Branch("trackster_x", &trInfo.trackster_x);
    fTree->Branch("trackster_y", &trInfo.trackster_y);
    fTree->Branch("trackster_z", &trInfo.trackster_z);
    fTree->Branch("trackster_eta", &trInfo.trackster_eta);
    fTree->Branch("trackster_phi", &trInfo.trackster_phi);
    fTree->Branch("trackster_firstlayer", &trInfo.trackster_firstlayer);
    fTree->Branch("trackster_lastlayer", &trInfo.trackster_lastlayer);
    fTree->Branch("trackster_layersnum", &trInfo.trackster_layersnum);
    fTree->Branch("trackster_id", &trInfo.trackster_id);
    fTree->Branch("trackster_seedIndex", &trInfo.trackster_seedIndex);
    fTree->Branch("trackster_time", &trInfo.trackster_time);
    fTree->Branch("trackster_timeError", &trInfo.trackster_timeError);
    fTree->Branch("trackster_regr_energy", &trInfo.trackster_regr_energy);
    fTree->Branch("trackster_raw_energy", &trInfo.trackster_raw_energy);
    fTree->Branch("trackster_raw_em_energy", &trInfo.trackster_raw_em_energy);
    fTree->Branch("trackster_raw_pt", &trInfo.trackster_raw_pt);
    fTree->Branch("trackster_raw_em_pt", &trInfo.trackster_raw_em_pt);
    fTree->Branch("trackster_id_prob_1", &trInfo.trackster_id_prob_1);
    fTree->Branch("trackster_id_prob_2", &trInfo.trackster_id_prob_2);
    fTree->Branch("trackster_id_prob_3", &trInfo.trackster_id_prob_3);
    fTree->Branch("trackster_id_prob_4", &trInfo.trackster_id_prob_4);
    fTree->Branch("trackster_id_prob_5", &trInfo.trackster_id_prob_5);
    fTree->Branch("trackster_id_prob_6", &trInfo.trackster_id_prob_6);
    fTree->Branch("trackster_id_prob_7", &trInfo.trackster_id_prob_7);
    fTree->Branch("trackster_id_prob_8", &trInfo.trackster_id_prob_8);
    fTree->Branch("trackster_numedges", &trInfo.trackster_numedges);
    fTree->Branch("trackster_edge_layerin_id", &trInfo.trackster_edge_layerin_id);
    fTree->Branch("trackster_edge_layerout_id", &trInfo.trackster_edge_layerout_id);
    fTree->Branch("trackster_edge_layerin", &trInfo.trackster_edge_layerin);
    fTree->Branch("trackster_edge_layerout", &trInfo.trackster_edge_layerout);
    fTree->Branch("trackster_delta_energy", &trInfo.trackster_delta_energy);
    fTree->Branch("trackster_delta_energy_relative", &trInfo.trackster_delta_energy_relative);
    fTree->Branch("trackster_delta_layer", &trInfo.trackster_delta_layer);
    fTree->Branch("trackster_angle_alpha", &trInfo.trackster_angle_alpha);
    fTree->Branch("trackster_angle_alpha_alternative", &trInfo.trackster_angle_alpha_alternative);
    fTree->Branch("trackster_angle_beta", &trInfo.trackster_angle_beta);

    fTree->Branch("trackster_numberOfNoiseHitsInTS", &trInfo.trackster_numberOfNoiseHitsInTS);
    fTree->Branch("trackster_maxCPId_byNumberOfHits", &trInfo.trackster_maxCPId_byNumberOfHits);
    fTree->Branch("trackster_maxCPNumberOfHitsInTS", &trInfo.trackster_maxCPNumberOfHitsInTS);
    fTree->Branch("trackster_maxCPId_byEnergy", &trInfo.trackster_maxCPId_byEnergy);
    fTree->Branch("trackster_maxEnergySharedTSandCP", &trInfo.trackster_maxEnergySharedTSandCP);
    fTree->Branch("trackster_totalCPEnergyFromLayerCP", &trInfo.trackster_totalCPEnergyFromLayerCP);
    fTree->Branch("trackster_energyFractionOfTSinCP", &trInfo.trackster_energyFractionOfTSinCP);
    fTree->Branch("trackster_energyFractionOfCPinTS", &trInfo.trackster_energyFractionOfCPinTS);

    fTree->Branch("trackster_cpId", &trInfo.trackster_cpId);
    fTree->Branch("trackster_scId", &trInfo.trackster_scId);
    fTree->Branch("trackster_Id", &trInfo.trackster_Id);
    fTree->Branch("trackster_numofvertices", &trInfo.trackster_numofvertices);
    fTree->Branch("trackster_Raw_Energy", &trInfo.trackster_Raw_Energy);
    fTree->Branch("trackster_numberOfHitsInTS", &trInfo.trackster_numberOfHitsInTS);
    fTree->Branch("trackster_score_trackster2caloparticle", &trInfo.trackster_score_trackster2caloparticle);
    fTree->Branch("trackster_sharedenergy_trackster2caloparticle", &trInfo.trackster_sharedenergy_trackster2caloparticle);
    fTree->Branch("trackster_score_trackster2bestCaloparticle", &trInfo.trackster_score_trackster2bestCaloparticle);
    fTree->Branch("trackster_sharedenergy_trackster2bestCaloparticle", &trInfo.trackster_sharedenergy_trackster2bestCaloparticle);
    fTree->Branch("trackster_trackster2bestCaloparticle_eta", &trInfo.trackster_trackster2bestCaloparticle_eta);
    fTree->Branch("trackster_trackster2bestCaloparticle_phi", &trInfo.trackster_trackster2bestCaloparticle_phi);
    fTree->Branch("trackster_score_trackster2bestCaloparticle2", &trInfo.trackster_score_trackster2bestCaloparticle2);
    fTree->Branch("trackster_sharedenergy_trackster2bestCaloparticle2", &trInfo.trackster_sharedenergy_trackster2bestCaloparticle2);
    
    fTree->Branch("trackster_sts_cpId", &trInfo.trackster_sts_cpId);
    fTree->Branch("trackster_sts_id", &trInfo.trackster_sts_id);
    fTree->Branch("trackster_sts_ts_id", &trInfo.trackster_sts_ts_id);
    fTree->Branch("trackster_sts_SimEnergy", &trInfo.trackster_sts_SimEnergy);
    fTree->Branch("trackster_sts_SimEnergyWeight", &trInfo.trackster_sts_SimEnergyWeight);
    fTree->Branch("trackster_sts_trackster_raw_energy", &trInfo.trackster_sts_trackster_raw_energy);
    fTree->Branch("trackster_sts_score_caloparticle2trackster", &trInfo.trackster_sts_score_caloparticle2trackster);
    fTree->Branch("trackster_sts_sharedenergy_caloparticle2trackster", &trInfo.trackster_sts_sharedenergy_caloparticle2trackster);
    fTree->Branch("trackster_sts_eta", &trInfo.trackster_sts_eta);
    fTree->Branch("trackster_sts_phi", &trInfo.trackster_sts_phi);
    fTree->Branch("trackster_sts_pt", &trInfo.trackster_sts_pt);
    fTree->Branch("trackster_sts_raw_energy", &trInfo.trackster_sts_raw_energy);
    fTree->Branch("trackster_sts_scorePur_caloparticle2trackster", &trInfo.trackster_sts_scorePur_caloparticle2trackster);
    fTree->Branch("trackster_sts_sharedenergy_caloparticle2trackster_assoc", &trInfo.trackster_sts_sharedenergy_caloparticle2trackster_assoc);
    fTree->Branch("trackster_sts_besttrackster_raw_energy", &trInfo.trackster_sts_besttrackster_raw_energy);
    fTree->Branch("trackster_sts_scoreDupl_caloparticle2trackster", &trInfo.trackster_sts_scoreDupl_caloparticle2trackster);
    fTree->Branch("trackster_sts_sharedenergy_caloparticle2trackster_assoc2", &trInfo.trackster_sts_sharedenergy_caloparticle2trackster_assoc2);

  }

  void fillTrackstersInfo(trackstersInfo &trInfo,
			  simClustersInfo &scInfo, 
  			  layerClustersInfo &lcInfo,
  			  recHitInfo &rhInfo_sc,
  			  recHitInfo &rhInfo_lc,
  			  const ticl::TracksterCollection& tracksters,
			  std::vector<SimCluster> const& simClusters,
  			  const reco::CaloClusterCollection& clusters,
  			  const std::vector<TICLSeedingRegion>& ticlSeedingGlobal,
  			  const std::vector<TICLSeedingRegion>& ticlSeedingTrk,
			  std::unordered_map<DetId, std::vector<const HGCUncalibratedRecHit *>> const& UnCalibHitMap,
  			  std::unordered_map<DetId, const HGCRecHit*> const& hitMap,
  			  std::shared_ptr<hgcal::RecHitTools> recHitTools,
			  const bool doEdges, 
  			  unsigned int layers,
			  std::map<double, double> cummatbudg){

    const auto nTracksters = tracksters.size();

    // Loop through Tracksters
    for (unsigned int tstId = 0; tstId < nTracksters; ++tstId) {
      const auto& trackster = tracksters[tstId];
      //Will check later also if there are tracksters with zero LCs. 
      if (trackster.vertices().empty())
  	continue;

      //For the layers the Trackster expands to. Will use a set because there would be many
      //duplicates and then go back to vector for random access, since they say it is faster.
      std::set<int> trackster_layers;
      
      //This here is problematic. It doesn't relate to tracksters but it works as is since 
      //in the calibration only one trackster is created. We should keep rhInfo_lc and 
      //tstId as in fillLayerClustersInfo below e.g. as map tstId -> rhInfo_lc, then 
      //feed that into fillSimClustersInfo and disregard rechits that aren't not only 
      //matched to simhits but also don't belong to any trackster tstId. 
      //So, first I will go to SimCluster object and then if we want to impose apart from 
      //the simhit-rechit matching the criterion to belong to a LC of a trackster, i will 
      //come back here and add this and then uncomment below. 
      //fillSimClustersInfo(scInfo, rhInfo_sc, simClusters, UnCalibHitMap, hitMap, recHitTools, layers, cummatbudg);

      //Loop through layer clusters of the trackster
      for (const auto lcId : tracksters[tstId].vertices()) {
  	const std::vector<std::pair<DetId, float>>& hits_and_fractions = clusters[lcId].hitsAndFractions();
  	const auto firstHitDetId = hits_and_fractions[0].first;
  	//The layer that the layer cluster belongs to
  	int layerid =
  	  recHitTools->getLayerWithOffset(firstHitDetId) + layers * ((recHitTools->zside(firstHitDetId) + 1) >> 1) - 1;

	fillLayerClustersInfo(lcInfo, rhInfo_lc, clusters, lcId, tstId, true, UnCalibHitMap, hitMap, recHitTools, layers, cummatbudg);

  	trackster_layers.insert(layerid);


      }  // end of loop through layerClusters

      //Per Trackster quantities
      trInfo.trackster_x.push_back(tracksters[tstId].barycenter().x());
      trInfo.trackster_y.push_back(tracksters[tstId].barycenter().y());
      trInfo.trackster_z.push_back(tracksters[tstId].barycenter().z());
      trInfo.trackster_eta.push_back(tracksters[tstId].barycenter().eta());
      trInfo.trackster_phi.push_back(tracksters[tstId].barycenter().phi());
      trInfo.trackster_firstlayer.push_back((float)*trackster_layers.begin());
      trInfo.trackster_lastlayer.push_back((float)*trackster_layers.rbegin());
      trInfo.trackster_layersnum.push_back((float)trackster_layers.size());
      trInfo.trackster_id.push_back(tstId);

      // -------------------------------------------------------------------------------------
      // Plots on edges
      if (doEdges){
	unsigned int numedges = 0;
	std::vector<unsigned int> tmp_trackster_numedges;
	std::vector<unsigned int> tmp_trackster_edge_layerin;
	std::vector<unsigned int> tmp_trackster_edge_layerout;
	std::vector<unsigned int> tmp_trackster_edge_layerin_id;
	std::vector<unsigned int> tmp_trackster_edge_layerout_id;
	std::vector<float> tmp_trackster_delta_energy;
	std::vector<float> tmp_trackster_delta_energy_relative;
	std::vector<float> tmp_trackster_delta_layer;
	std::vector<float> tmp_trackster_angle_alpha;
	std::vector<float> tmp_trackster_angle_alpha_alternative;
	std::vector<float> tmp_trackster_angle_beta;
	tmp_trackster_numedges.clear();
	tmp_trackster_edge_layerin.clear();
	tmp_trackster_edge_layerout.clear();
	tmp_trackster_edge_layerin_id.clear();
	tmp_trackster_edge_layerout_id.clear();
	tmp_trackster_delta_energy.clear();
	tmp_trackster_delta_energy_relative.clear();
	tmp_trackster_delta_layer.clear();
	tmp_trackster_angle_alpha.clear();
	tmp_trackster_angle_alpha_alternative.clear();
	tmp_trackster_angle_beta.clear();
      
	for (const auto& edge : tracksters[tstId].edges()) {
	  auto& ic = clusters[edge[0]];
	  auto& oc = clusters[edge[1]];
	  auto const& cl_in = ic.hitsAndFractions()[0].first;
	  auto const& cl_out = oc.hitsAndFractions()[0].first;
	  auto const layer_in =
	    recHitTools->getLayerWithOffset(cl_in) + layers * ((recHitTools->zside(cl_in) + 1) >> 1) - 1;
	  auto const layer_out = 
	    recHitTools->getLayerWithOffset(cl_out) + layers * ((recHitTools->zside(cl_out) + 1) >> 1) - 1;

	  tmp_trackster_edge_layerin_id.push_back(edge[0]);
	  tmp_trackster_edge_layerout_id.push_back(edge[1]);
	  tmp_trackster_edge_layerin.push_back(layer_in);
	  tmp_trackster_edge_layerout.push_back(layer_out);

	  tmp_trackster_delta_energy.push_back(oc.energy() - ic.energy());
	  tmp_trackster_delta_energy_relative.push_back((oc.energy() - ic.energy()) / ic.energy());
	  tmp_trackster_delta_layer.push_back(layer_out - layer_in);

	  // Alpha angles
	  const auto& outer_outer_pos = oc.position();
	  const auto& outer_inner_pos = ic.position();
	  const auto& seed = tracksters[tstId].seedIndex();
	  auto seedGlobalPos = math::XYZPoint(
					      ticlSeedingGlobal[0].origin.x(), ticlSeedingGlobal[0].origin.y(), ticlSeedingGlobal[0].origin.z());
	  auto seedDirectionPos = outer_inner_pos;
	  if (tracksters[tstId].seedID().id() != 0) {
	    // Seed to trackster association is, at present, rather convoluted.
	    for (auto const& s : ticlSeedingTrk) {
	      if (s.index == seed) {
		seedGlobalPos = math::XYZPoint(s.origin.x(), s.origin.y(), s.origin.z());
		seedDirectionPos =
		  math::XYZPoint(s.directionAtOrigin.x(), s.directionAtOrigin.y(), s.directionAtOrigin.z());
		break;
	      }
	    }
	  }

	  auto alpha = (outer_inner_pos - seedGlobalPos).Dot(outer_outer_pos - outer_inner_pos) /
	    sqrt((outer_inner_pos - seedGlobalPos).Mag2() * (outer_outer_pos - outer_inner_pos).Mag2());
	  auto alpha_alternative = (outer_outer_pos - seedGlobalPos).Dot(seedDirectionPos) /
	    sqrt((outer_outer_pos - seedGlobalPos).Mag2() * seedDirectionPos.Mag2());
	
	  tmp_trackster_angle_alpha.push_back(alpha);
	  tmp_trackster_angle_alpha_alternative.push_back(alpha_alternative);

	  // Beta angle is usually computed using 2 edges. Another inner loop
	  // is therefore needed.	
	  std::vector<std::array<unsigned int, 2>> innerDoublets;
	  std::vector<std::array<unsigned int, 2>> outerDoublets;

	  for (const auto& otherEdge : tracksters[tstId].edges()) {
	    if (otherEdge[1] == edge[0]) {
	      innerDoublets.push_back(otherEdge);
	    }
	    if (edge[1] == otherEdge[0]) {
	      outerDoublets.push_back(otherEdge);
	    }
	  } 

	  for (const auto& inner : innerDoublets) {
	    const auto& inner_ic = clusters[inner[0]];
	    const auto& inner_inner_pos = inner_ic.position();
	    auto beta = (outer_inner_pos - inner_inner_pos).Dot(outer_outer_pos - inner_inner_pos) /
	      sqrt((outer_inner_pos - inner_inner_pos).Mag2() * (outer_outer_pos - inner_inner_pos).Mag2());

	    tmp_trackster_angle_beta.push_back(beta);
	  }

	  tmp_trackster_numedges.push_back(numedges);
	  ++numedges;
	} //end of loop through edges

	trInfo.trackster_numedges.push_back(tmp_trackster_numedges);
	trInfo.trackster_edge_layerin_id.push_back(tmp_trackster_edge_layerin_id);
	trInfo.trackster_edge_layerout_id.push_back(tmp_trackster_edge_layerout_id);
	trInfo.trackster_edge_layerin.push_back(tmp_trackster_edge_layerin);
	trInfo.trackster_edge_layerout.push_back(tmp_trackster_edge_layerout);
	trInfo.trackster_delta_energy.push_back(tmp_trackster_delta_energy);
	trInfo.trackster_delta_energy_relative.push_back(tmp_trackster_delta_energy_relative);
	trInfo.trackster_delta_layer.push_back(tmp_trackster_delta_layer);
	trInfo.trackster_angle_alpha.push_back(tmp_trackster_angle_alpha);
	trInfo.trackster_angle_alpha_alternative.push_back(tmp_trackster_angle_alpha_alternative);
	trInfo.trackster_angle_beta.push_back(tmp_trackster_angle_beta);

	trInfo.trackster_time.push_back(tracksters[tstId].time());
	trInfo.trackster_timeError.push_back(tracksters[tstId].timeError());
	trInfo.trackster_regr_energy.push_back(tracksters[tstId].regressed_energy());
	trInfo.trackster_raw_energy.push_back(tracksters[tstId].raw_energy());
	trInfo.trackster_raw_em_energy.push_back(tracksters[tstId].raw_em_energy());
	trInfo.trackster_raw_pt.push_back(tracksters[tstId].raw_pt());
	trInfo.trackster_raw_em_pt.push_back(tracksters[tstId].raw_em_pt());
	const auto& probs = tracksters[tstId].id_probabilities();
	std::vector<int> sorted_probs_idx(probs.size());
	std::iota(begin(sorted_probs_idx), end(sorted_probs_idx), 0);
	std::sort(begin(sorted_probs_idx), end(sorted_probs_idx), [&probs](int i, int j) { return probs[i] > probs[j]; });

	trInfo.trackster_id_prob_1.push_back(sorted_probs_idx[1]);
	trInfo.trackster_id_prob_2.push_back(sorted_probs_idx[2]);
	trInfo.trackster_id_prob_3.push_back(sorted_probs_idx[3]);
	trInfo.trackster_id_prob_4.push_back(sorted_probs_idx[4]);
	trInfo.trackster_id_prob_5.push_back(sorted_probs_idx[5]);
	trInfo.trackster_id_prob_6.push_back(sorted_probs_idx[6]);
	trInfo.trackster_id_prob_7.push_back(sorted_probs_idx[7]);
	trInfo.trackster_id_prob_8.push_back(sorted_probs_idx[8]);

      } // end of plots on edges
      // -------------------------------------------------------------------------------------

    } //end of loop over tracksters

  }// end of fill tracksters info


  void prepare_tracksters_to_SimTracksters( trackstersInfo &trInfo,
					    const ticl::TracksterCollection& tracksters,
					    const reco::CaloClusterCollection& layerClusters,
					    const ticl::TracksterCollection& simTSs,
					    const int i,
					    const ticl::TracksterCollection& simTSs_fromCP,
					    const std::map<unsigned int, std::vector<unsigned int>>& cpToSc_SimTrackstersMap,
					    std::vector<SimCluster> const& sC,
					    const edm::ProductID& cPHandle_id,
					    std::vector<CaloParticle> const& cP,
					    std::vector<size_t> const& cPIndices,
					    std::vector<size_t> const& cPSelectedIndices,
					    std::unordered_map<DetId, const HGCRecHit*> const& hitMap,
					    unsigned int layers) {
    
    const auto nTracksters = tracksters.size();
    const auto nSimTracksters = simTSs.size();

    /* std::vector<std::vector<float> > tmp_trackster_raw_energy; */
    /* std::vector<std::vector<unsigned int> > tmp_trackster_numberOfHitsInTS; */
    /* std::vector<std::vector<unsigned int> > tmp_trackster_numberOfNoiseHitsInTS; */
    /* std::vector<std::vector<int> > tmp_trackster_maxCPId_byNumberOfHits; */
    /* std::vector<std::vector<unsigned int> > tmp_trackster_maxCPNumberOfHitsInTS; */
    /* std::vector<std::vector<int> > tmp_trackster_maxCPId_byEnergy; */
    /* std::vector<std::vector<float> > tmp_trackster_maxEnergySharedTSandCP; */
    /* std::vector<std::vector<float> > tmp_trackster_totalCPEnergyFromLayerCP; */
    /* std::vector<std::vector<float> > tmp_trackster_energyFractionOfTSinCP; */
    /* std::vector<std::vector<float> > tmp_trackster_energyFractionOfCPinTS; */
    std::vector<float> tmp_trackster_raw_energy;
    std::vector<unsigned int> tmp_trackster_numberOfHitsInTS;
    std::vector<unsigned int> tmp_trackster_numberOfNoiseHitsInTS;
    std::vector<int> tmp_trackster_maxCPId_byNumberOfHits;
    std::vector<unsigned int> tmp_trackster_maxCPNumberOfHitsInTS;
    std::vector<int> tmp_trackster_maxCPId_byEnergy;
    std::vector<float> tmp_trackster_maxEnergySharedTSandCP;
    std::vector<float> tmp_trackster_totalCPEnergyFromLayerCP;
    std::vector<float> tmp_trackster_energyFractionOfTSinCP;
    std::vector<float> tmp_trackster_energyFractionOfCPinTS;
    
    tmp_trackster_raw_energy.clear();
    tmp_trackster_numberOfHitsInTS.clear();
    tmp_trackster_numberOfNoiseHitsInTS.clear();
    tmp_trackster_maxCPId_byNumberOfHits.clear();
    tmp_trackster_maxCPNumberOfHitsInTS.clear();
    tmp_trackster_maxCPId_byEnergy.clear();
    tmp_trackster_maxEnergySharedTSandCP.clear();
    tmp_trackster_totalCPEnergyFromLayerCP.clear();
    tmp_trackster_energyFractionOfTSinCP.clear();
    tmp_trackster_energyFractionOfCPinTS.clear();

    tmp_trackster_raw_energy.resize(nTracksters);
    tmp_trackster_numberOfHitsInTS.resize(nTracksters);
    tmp_trackster_numberOfNoiseHitsInTS.resize(nTracksters);
    tmp_trackster_maxCPId_byNumberOfHits.resize(nTracksters);
    tmp_trackster_maxCPNumberOfHitsInTS.resize(nTracksters);
    tmp_trackster_maxCPId_byEnergy.resize(nTracksters);
    tmp_trackster_maxEnergySharedTSandCP.resize(nTracksters);
    tmp_trackster_totalCPEnergyFromLayerCP.resize(nTracksters);
    tmp_trackster_energyFractionOfTSinCP.resize(nTracksters);
    tmp_trackster_energyFractionOfCPinTS.resize(nTracksters);
  
    std::unordered_map<DetId, std::vector<detIdInfoInCluster>> detIdSimTSId_Map;
    std::unordered_map<DetId, std::vector<detIdInfoInTrackster>> detIdToTracksterId_Map;
    std::vector<int> tracksters_FakeMerge(nTracksters, 0);
    std::vector<int> tracksters_PurityDuplicate(nTracksters, 0);

    // This vector contains the ids of the SimTracksters contributing with at least one hit to the Trackster and the reconstruction error
    //stsInTrackster[trackster][STSids]
    //Connects a Trackster with all related SimTracksters.
    std::vector<std::vector<std::pair<unsigned int, float>>> stsInTrackster;
    stsInTrackster.resize(nTracksters);

    //cPOnLayer[caloparticle][layer]
    //This defines a "calo particle on layer" concept. It is only filled in case
    //that calo particle has a reconstructed hit related via detid. So, a cPOnLayer[i][j] connects a
    //specific calo particle i in layer j with:
    //1. the sum of all rechits energy times fraction of the relevant simhit in layer j related to that calo particle i.
    //2. the hits and fractions of that calo particle i in layer j.
    //3. the layer clusters with matched rechit id.
    std::unordered_map<int, caloParticleOnLayer> cPOnLayer;
    std::unordered_map<int, std::vector<caloParticleOnLayer>> sCOnLayer;
    //Consider CaloParticles coming from the hard scatterer, excluding the PU contribution.
    for (const auto cpIndex : cPIndices) {
      cPOnLayer[cpIndex].caloParticleId = cpIndex;
      cPOnLayer[cpIndex].energy = 0.f;
      cPOnLayer[cpIndex].hits_and_fractions.clear();
      //const auto nSC_inCP = cP[cpIndex].simClusters().size();
      const auto nSC_inCP = sC.size();
      sCOnLayer[cpIndex].resize(nSC_inCP);
      for (unsigned int iSC = 0; iSC < nSC_inCP; iSC++) {
	sCOnLayer[cpIndex][iSC].caloParticleId = cpIndex;
	sCOnLayer[cpIndex][iSC].energy = 0.f;
	sCOnLayer[cpIndex][iSC].hits_and_fractions.clear();
      }
    }

    auto getCPId = [](const ticl::Trackster& simTS,
		      const unsigned int iSTS,
		      const edm::ProductID& cPHandle_id,
		      const std::map<unsigned int, std::vector<unsigned int>>& cpToSc_SimTrackstersMap,
		      const ticl::TracksterCollection& simTSs_fromCP) {
								       unsigned int cpId = -1;

								       const auto productID = simTS.seedID();
								       if (productID == cPHandle_id) {
									 cpId = simTS.seedIndex();
								       } else {  // SimTrackster from SimCluster
									 const auto findSimTSFromCP = std::find_if(
														   std::begin(cpToSc_SimTrackstersMap),
														   std::end(cpToSc_SimTrackstersMap),
														   [&](const std::pair<unsigned int, std::vector<unsigned int>>& cpToScs) {
														     return std::find(std::begin(cpToScs.second), std::end(cpToScs.second), iSTS) != std::end(cpToScs.second);
														   });
									 if (findSimTSFromCP != std::end(cpToSc_SimTrackstersMap)) {
									   cpId = simTSs_fromCP[findSimTSFromCP->first].seedIndex();
									 }
								       }

								       return cpId;
    };

    auto getLCId = [](const std::vector<unsigned int>& tst_vertices,
		      const reco::CaloClusterCollection& layerClusters,
		      const DetId& hitid) {
					   unsigned int lcId = -1;
					   std::for_each(std::begin(tst_vertices), std::end(tst_vertices), [&](unsigned int idx) {
					       const auto& lc_haf = layerClusters[idx].hitsAndFractions();
					       const auto& hitFound = std::find_if(std::begin(lc_haf),
										   std::end(lc_haf),
										   [&hitid](const std::pair<DetId, float>& v) { return v.first == hitid; });
					       if (hitFound != lc_haf.end())  // not all hits may be clusterized
						 lcId = idx;
					     });
					   return lcId;
    };

    for (unsigned int iSTS = 0; iSTS < nSimTracksters; ++iSTS) {
      const auto cpId = getCPId(simTSs[iSTS], iSTS, cPHandle_id, cpToSc_SimTrackstersMap, simTSs_fromCP);
      if (std::find(cPIndices.begin(), cPIndices.end(), cpId) == cPIndices.end())
	continue;

      // Loop through SimClusters
      for (const auto& simCluster : cP[cpId].simClusters()) {
	auto iSim = simTSs[iSTS].seedIndex();
	if (simTSs[iSTS].seedID() != cPHandle_id) {  // SimTrackster from SimCluster
	  if (iSim != (&(*simCluster) - &(sC[0])))
	    continue;
	} else
	  iSim = 0;

	for (const auto& it_haf : simCluster->hits_and_fractions()) {
	  const auto hitid = (it_haf.first);
	  const auto lcId = getLCId(simTSs[iSTS].vertices(), layerClusters, hitid);
	  //V9:maps the layers in -z: 0->51 and in +z: 52->103
	  //V10:maps the layers in -z: 0->49 and in +z: 50->99
	  const auto itcheck = hitMap.find(hitid);
	  //Checks whether the current hit belonging to sim cluster has a reconstructed hit.
	  if ((i == 0 && itcheck != hitMap.end()) || (i > 0 && int(lcId) >= 0)) {
	    const auto elemId = (i == 0) ? hitid : lcId;
	    const auto iLC = std::find(simTSs[iSTS].vertices().begin(), simTSs[iSTS].vertices().end(), lcId);
	    const auto lcFraction =
	      1.f / simTSs[iSTS].vertex_multiplicity(std::distance(std::begin(simTSs[iSTS].vertices()), iLC));
	    const auto elemFr = (i == 0) ? it_haf.second : lcFraction;
	    //Since the current hit from sim cluster has a reconstructed hit with the same detid,
	    //make a map that will connect a detid with:
	    //1. the CaloParticles that have a SimCluster with sim hits in that cell via caloparticle id.
	    //2. the sum of all SimHits fractions that contributes to that detid.
	    //So, keep in mind that in case of multiple CaloParticles contributing in the same cell
	    //the fraction is the sum over all calo particles. So, something like:
	    //detid: (caloparticle 1, sum of hits fractions in that detid over all cp) , (caloparticle 2, sum of hits fractions in that detid over all cp), (caloparticle 3, sum of hits fractions in that detid over all cp) ...
	    if (detIdSimTSId_Map.find(elemId) == detIdSimTSId_Map.end()) {
	      detIdSimTSId_Map[elemId] = std::vector<detIdInfoInCluster>();
	      detIdSimTSId_Map[elemId].emplace_back(detIdInfoInCluster{iSTS, elemFr});
	    } else {
	      auto findSTSIt =
		std::find(detIdSimTSId_Map[elemId].begin(),
			  detIdSimTSId_Map[elemId].end(),
			  detIdInfoInCluster{
			    iSTS, 0});  // only the first element is used for the matching (overloaded operator==)
	      if (findSTSIt != detIdSimTSId_Map[elemId].end()) {
		if (i == 0)
		  findSTSIt->fraction += elemFr;
	      } else {
		detIdSimTSId_Map[elemId].emplace_back(detIdInfoInCluster{iSTS, elemFr});
	      }
	    }
	    const auto hitEn = itcheck->second->energy();
	    //Since the current hit from sim cluster has a reconstructed hit with the same detid,
	    //fill the cPOnLayer[caloparticle][layer] object with energy (sum of all rechits energy times fraction
	    //of the relevant simhit) and keep the hit (detid and fraction) that contributed.
	    cPOnLayer[cpId].energy += it_haf.second * hitEn;
	    sCOnLayer[cpId][iSim].energy += lcFraction * hitEn;
	    // Need to compress the hits and fractions in order to have a
	    // reasonable score between CP and LC. Imagine, for example, that a
	    // CP has detID X used by 2 SimClusters with different fractions. If
	    // a single LC uses X with fraction 1 and is compared to the 2
	    // contributions separately, it will be assigned a score != 0, which
	    // is wrong.
	    auto& haf = cPOnLayer[cpId].hits_and_fractions;
	    auto found = std::find_if(std::begin(haf),
				      std::end(haf),
				      [&hitid](const std::pair<DetId, float>& v) { return v.first == hitid; });
	    if (found != haf.end())
	      found->second += it_haf.second;
	    else
	      haf.emplace_back(hitid, it_haf.second);
	    // Same for sCOnLayer
	    auto& haf_sc = sCOnLayer[cpId][iSim].hits_and_fractions;
	    auto found_sc = std::find_if(std::begin(haf_sc),
					 std::end(haf_sc),
					 [&hitid](const std::pair<DetId, float>& v) { return v.first == hitid; });
	    if (found_sc != haf_sc.end())
	      found_sc->second += it_haf.second;
	    else
	      haf_sc.emplace_back(hitid, it_haf.second);
	  }
	}  // end of loop through SimHits
      }    // end of loop through SimClusters
    }      // end of loop through SimTracksters

    auto apply_LCMultiplicity = [](const ticl::Trackster& trackster, const reco::CaloClusterCollection& layerClusters) {
															std::vector<std::pair<DetId, float>> hits_and_fractions_norm;
															int lcInTst = 0;
															std::for_each(std::begin(trackster.vertices()), std::end(trackster.vertices()), [&](unsigned int idx) {
															    const auto fraction = 1.f / trackster.vertex_multiplicity(lcInTst++);
															    for (const auto& cell : layerClusters[idx].hitsAndFractions()) {
															      hits_and_fractions_norm.emplace_back(
																				   cell.first, cell.second * fraction);  // cell.second is the hit fraction in the layerCluster
															    }
															  });
															return hits_and_fractions_norm;
    };

    //auto ScoreCutSTStoTSPurDup = ScoreCutSTStoTSPurDup_[0];
    auto ScoreCutTStoSTSFakeMerge = ScoreCutTStoSTSFakeMerge_[0];

    // Loop through Tracksters
    for (unsigned int tstId = 0; tstId < nTracksters; ++tstId) {
      const auto& tst = tracksters[tstId];
      if (tstId == 0)
	if ((i > 0) && (tst.ticlIteration() == ticl::Trackster::SIM)) {
	  //ScoreCutSTStoTSPurDup = ScoreCutSTStoTSPurDup_[i];
	  ScoreCutTStoSTSFakeMerge = ScoreCutTStoSTSFakeMerge_[i];
	}

      if (tst.vertices().empty())
	continue;

      std::unordered_map<unsigned, float> CPEnergyInTS;
      int maxCPId_byNumberOfHits = -1;
      unsigned int maxCPNumberOfHitsInTS = 0;
      int maxCPId_byEnergy = -1;
      float maxEnergySharedTSandCP = 0.f;
      float energyFractionOfTSinCP = 0.f;
      float energyFractionOfCPinTS = 0.f;

      //In case of matched rechit-simhit, so matched
      //CaloParticle-LayerCluster-Trackster, he counts and saves the number of
      //rechits related to the maximum energy CaloParticle out of all
      //CaloParticles related to that layer cluster and Trackster.

      std::unordered_map<unsigned, unsigned> occurrencesCPinTS;
      unsigned int numberOfNoiseHitsInTS = 0;
      unsigned int numberOfHaloHitsInTS = 0;

      const auto tst_hitsAndFractions = apply_LCMultiplicity(tst, layerClusters);
      const auto numberOfHitsInTS = tst_hitsAndFractions.size();

      //hitsToCaloParticleId is a vector of ints, one for each rechit of the
      //layer cluster under study. If negative, there is no simhit from any CaloParticle related.
      //If positive, at least one CaloParticle has been found with matched simhit.
      //In more detail:
      // 1. hitsToCaloParticleId[iHit] = -3
      //    TN:  These represent Halo Cells(N) that have not been
      //    assigned to any CaloParticle (hence the T).
      // 2. hitsToCaloParticleId[iHit] = -2
      //    FN: There represent Halo Cells(N) that have been assigned
      //    to a CaloParticle (hence the F, since those should have not been marked as halo)
      // 3. hitsToCaloParticleId[iHit] = -1
      //    FP: These represent Real Cells(P) that have not been
      //    assigned to any CaloParticle (hence the F, since these are fakes)
      // 4. hitsToCaloParticleId[iHit] >= 0
      //    TP There represent Real Cells(P) that have been assigned
      //    to a CaloParticle (hence the T)
      std::vector<int> hitsToCaloParticleId(numberOfHitsInTS);

      //Loop through the hits of the trackster under study
      for (unsigned int iHit = 0; iHit < numberOfHitsInTS; iHit++) {
	const auto rh_detid = tst_hitsAndFractions[iHit].first;
	const auto rhFraction = tst_hitsAndFractions[iHit].second;

	const auto lcId_r = getLCId(tst.vertices(), layerClusters, rh_detid);
	const auto iLC_r = std::find(tst.vertices().begin(), tst.vertices().end(), lcId_r);
	const auto lcFraction_r = 1.f / tst.vertex_multiplicity(std::distance(std::begin(tst.vertices()), iLC_r));

	//Make a map that will connect a detid (that belongs to a rechit of the layer cluster under study,
	//no need to save others) with:
	//1. the layer clusters that have rechits in that detid
	//2. the fraction of the rechit of each layer cluster that contributes to that detid.
	//So, something like:
	//detid: (layer cluster 1, hit fraction) , (layer cluster 2, hit fraction), (layer cluster 3, hit fraction) ...
	//here comparing with the calo particle map above the
	if (detIdToTracksterId_Map.find(rh_detid) == detIdToTracksterId_Map.end()) {
	  detIdToTracksterId_Map[rh_detid] = std::vector<detIdInfoInTrackster>();
	  detIdToTracksterId_Map[rh_detid].emplace_back(
							detIdInfoInTrackster{tstId, lcId_r, rhFraction});
	} else {
	  auto findTSIt =
	    std::find(detIdToTracksterId_Map[rh_detid].begin(),
		      detIdToTracksterId_Map[rh_detid].end(),
		      detIdInfoInTrackster{
			tstId, 0, 0});  // only the first element is used for the matching (overloaded operator==)
	  if (findTSIt != detIdToTracksterId_Map[rh_detid].end()) {
	    if (i == 0)
	      findTSIt->fraction += rhFraction;
	  } else {
	    detIdToTracksterId_Map[rh_detid].emplace_back(detIdInfoInTrackster{tstId, lcId_r, rhFraction});
	  }
	}

	// if the fraction is zero or the hit does not belong to any calo
	// particle, set the caloparticleId for the hit to -1 this will
	// contribute to the number of noise hits
	// MR Remove the case in which the fraction is 0, since this could be a
	// real hit that has been marked as halo.
	if (rhFraction == 0.) {
	  hitsToCaloParticleId[iHit] = -2;
	  numberOfHaloHitsInTS++;
	}

	// Check whether the RecHit of the trackster under study has a SimHit in the same cell
	const auto elemId = (i == 0) ? rh_detid.rawId() : lcId_r;
	const auto recoFr = (i == 0) ? rhFraction : lcFraction_r;
	const auto& hit_find_in_STS = detIdSimTSId_Map.find(elemId);
	if (hit_find_in_STS == detIdSimTSId_Map.end()) {
	  hitsToCaloParticleId[iHit] -= 1;
	} else {
	  // Since the hit is belonging to the layer cluster, it must be also in the rechits map
	  const auto hitEn = hitMap.find(rh_detid)->second->energy();
	  //const auto layerId =
	  //recHitTools_->getLayerWithOffset(rh_detid) + layers * ((recHitTools_->zside(rh_detid) + 1) >> 1) - 1;
	  //0;

	  auto maxCPEnergyInTS = 0.f;
	  auto maxCPId = -1;
	  for (const auto& h : hit_find_in_STS->second) {
	    const auto shared_fraction = std::min(recoFr, h.fraction);
	    const auto iSTS = h.clusterId;
	    const auto& simTS = simTSs[iSTS];
	    auto iSim = simTS.seedIndex();
	    if (simTSs[iSTS].seedID() == cPHandle_id)  // SimTrackster from CaloParticle
	      iSim = 0;

	    // SimTrackster with simHits connected via detid with the rechit under study
	    //So, from all layers clusters, find the rechits that are connected with a calo particle and save/calculate the
	    //energy of that calo particle as the sum over all rechits of the rechits energy weighted
	    //by the caloparticle's fraction related to that rechit.
	    const auto cpId = getCPId(simTS, iSTS, cPHandle_id, cpToSc_SimTrackstersMap, simTSs_fromCP);
	    if (std::find(cPIndices.begin(), cPIndices.end(), cpId) == cPIndices.end())
	      continue;

	    CPEnergyInTS[cpId] += shared_fraction * hitEn;
	    //Here cPOnLayer[caloparticle][layer] describe above is set.
	    //Here for Tracksters with matched rechit the CP fraction times hit energy is added and saved .
	    cPOnLayer[cpId].layerClusterIdToEnergyAndScore[tstId].first += shared_fraction * hitEn;
	    sCOnLayer[cpId][iSim].layerClusterIdToEnergyAndScore[tstId].first += shared_fraction * hitEn;
	    cPOnLayer[cpId].layerClusterIdToEnergyAndScore[tstId].second = FLT_MAX;
	    sCOnLayer[cpId][iSim].layerClusterIdToEnergyAndScore[tstId].second = FLT_MAX;
	    //stsInTrackster[trackster][STSids]
	    //Connects a Trackster with all related SimTracksters.
	    stsInTrackster[tstId].emplace_back(iSTS, FLT_MAX);
	    //From all CaloParticles related to a layer cluster, it saves id and energy of the calo particle
	    //that after simhit-rechit matching in layer has the maximum energy.
	    if (shared_fraction > maxCPEnergyInTS) {
	      //energy is used only here. cpid is saved for Tracksters
	      maxCPEnergyInTS = CPEnergyInTS[cpId];
	      maxCPId = cpId;
	    }
	  }
	  //Keep in mind here maxCPId could be zero. So, below ask for negative not including zero to count noise.
	  hitsToCaloParticleId[iHit] = maxCPId;
	}

      }  //end of loop through rechits of the layer cluster.

      //Loop through all rechits to count how many of them are noise and how many are matched.
      //In case of matched rechit-simhit, he counts and saves the number of rechits related to the maximum energy CaloParticle.
      for (auto c : hitsToCaloParticleId) {
	if (c < 0)
	  numberOfNoiseHitsInTS++;
	else
	  occurrencesCPinTS[c]++;
      }

      //Below from all maximum energy CaloParticles, he saves the one with the largest amount
      //of related rechits.
      for (auto& c : occurrencesCPinTS) {
	if (c.second > maxCPNumberOfHitsInTS) {
	  maxCPId_byNumberOfHits = c.first;
	  maxCPNumberOfHitsInTS = c.second;
	}
      }

      //Find the CaloParticle that has the maximum energy shared with the Trackster under study.
      for (auto& c : CPEnergyInTS) {
	if (c.second > maxEnergySharedTSandCP) {
	  maxCPId_byEnergy = c.first;
	  maxEnergySharedTSandCP = c.second;
	}
      }
      //The energy of the CaloParticle that found to have the maximum energy shared with the Trackster under study.
      float totalCPEnergyFromLayerCP = 0.f;
      if (maxCPId_byEnergy >= 0) {
	//Loop through all layers
	//for (unsigned int j = 0; j < layers * 2; ++j) {
	totalCPEnergyFromLayerCP += cPOnLayer[maxCPId_byEnergy].energy;
	//}
	energyFractionOfCPinTS = maxEnergySharedTSandCP / totalCPEnergyFromLayerCP;
	if (tst.raw_energy() > 0.f) {
	  energyFractionOfTSinCP = maxEnergySharedTSandCP / tst.raw_energy();
	}
      }

      tmp_trackster_raw_energy[tstId] = tst.raw_energy();
      tmp_trackster_numberOfHitsInTS[tstId] = numberOfHitsInTS;
      tmp_trackster_numberOfNoiseHitsInTS[tstId] = numberOfNoiseHitsInTS;
      tmp_trackster_maxCPId_byNumberOfHits[tstId] = maxCPId_byNumberOfHits;
      tmp_trackster_maxCPNumberOfHitsInTS[tstId] = maxCPNumberOfHitsInTS;
      tmp_trackster_maxCPId_byEnergy[tstId] = maxCPId_byEnergy;
      tmp_trackster_maxEnergySharedTSandCP[tstId] = maxEnergySharedTSandCP;
      tmp_trackster_totalCPEnergyFromLayerCP[tstId] = totalCPEnergyFromLayerCP;
      tmp_trackster_energyFractionOfTSinCP[tstId] = energyFractionOfTSinCP;
      tmp_trackster_energyFractionOfCPinTS[tstId] = energyFractionOfCPinTS;

      /* Std::Cout << std::setw(12) << "Trackster\t" << std::setw(10) << "energy\t" << std::setw(5) */
      /*                            << "nhits\t" << std::setw(12) << "noise hits\t" << std::setw(22) */
      /*                            << "maxCPId_byNumberOfHits\t" << std::setw(8) << "nhitsCP\t" << std::setw(16) */
      /*                            << "maxCPId_byEnergy\t" << std::setw(23) << "maxEnergySharedTSandCP\t" << std::setw(22) */
      /*                            << "totalCPEnergyFromAllLayerCP\t" << std::setw(22) << "energyFractionOfTSinCP\t" */
      /*                            << std::setw(25) << "energyFractionOfCPinTS\t" << std::endl; */
      /* std::cout << std::setw(12) << tstId << "\t"   */
      /*                            << std::setw(10) << tst.raw_energy() << "\t" << std::setw(5) << numberOfHitsInTS << "\t" */
      /*                            << std::setw(12) << numberOfNoiseHitsInTS << "\t" << std::setw(22) */
      /*                            << maxCPId_byNumberOfHits << "\t" << std::setw(8) << maxCPNumberOfHitsInTS << "\t" */
      /*                            << std::setw(16) << maxCPId_byEnergy << "\t" << std::setw(23) << maxEnergySharedTSandCP */
      /*                            << "\t" << std::setw(22) << totalCPEnergyFromLayerCP << "\t" << std::setw(22) */
      /*                            << energyFractionOfTSinCP << "\t" << std::setw(25) << energyFractionOfCPinTS */
      /*                            << std::endl; */

    }  // end of loop through Tracksters

  
    // Loop through Tracksters
    for (unsigned int tstId = 0; tstId < nTracksters; ++tstId) {
      const auto& tst = tracksters[tstId];
      if (tst.vertices().empty())
	continue;

      // find the unique SimTrackster ids contributing to the Trackster
      //stsInTrackster[trackster][STSids]
      std::sort(stsInTrackster[tstId].begin(), stsInTrackster[tstId].end());
      const auto last = std::unique(stsInTrackster[tstId].begin(), stsInTrackster[tstId].end());
      stsInTrackster[tstId].erase(last, stsInTrackster[tstId].end());

      if (tst.raw_energy() == 0. && !stsInTrackster[tstId].empty()) {
	//Loop through all SimTracksters contributing to Trackster tstId
	for (auto& stsPair : stsInTrackster[tstId]) {
	  // In case of a Trackster with zero energy but related SimTracksters the score is set to 1
	  stsPair.second = 1.;
	  /* std::cout << "Trackster Id:\t" << tstId << "\tSimTrackster id:\t" << stsPair.first */
	  /*                            << "\tscore\t" << stsPair.second << std::endl; */
	}
	continue;
      }

      const auto tst_hitsAndFractions = apply_LCMultiplicity(tst, layerClusters);

      // Compute the correct normalization
      float tracksterEnergy = 0.f, invTracksterEnergyWeight = 0.f;
      for (const auto& haf : tst_hitsAndFractions) {
	float hitFr = 0.f;
	if (i == 0) {
	  hitFr = haf.second;
	} else {
	  const auto lcId = getLCId(tst.vertices(), layerClusters, haf.first);
	  const auto iLC = std::find(tst.vertices().begin(), tst.vertices().end(), lcId);
	  hitFr = 1.f / tst.vertex_multiplicity(std::distance(std::begin(tst.vertices()), iLC));
	}
	tracksterEnergy += hitFr * hitMap.at(haf.first)->energy();
	invTracksterEnergyWeight += pow(hitFr * hitMap.at(haf.first)->energy(), 2);
      }
      if (invTracksterEnergyWeight)
	invTracksterEnergyWeight = 1.f / invTracksterEnergyWeight;

      for (const auto& haf : tst_hitsAndFractions) {
	const auto rh_detid = haf.first;
	unsigned int elemId = 0;
	float rhFraction = 0.f;
	if (i == 0) {
	  elemId = rh_detid.rawId();
	  rhFraction = haf.second;
	} else {
	  const auto lcId = getLCId(tst.vertices(), layerClusters, rh_detid);
	  elemId = lcId;
	  const auto iLC = std::find(tst.vertices().begin(), tst.vertices().end(), lcId);
	  rhFraction = 1.f / tst.vertex_multiplicity(std::distance(std::begin(tst.vertices()), iLC));
	}

	bool hitWithNoSTS = false;
	if (detIdSimTSId_Map.find(elemId) == detIdSimTSId_Map.end())
	  hitWithNoSTS = true;
	const HGCRecHit* hit = hitMap.find(rh_detid)->second;
	const auto hitEnergyWeight = pow(hit->energy(), 2);

	for (auto& stsPair : stsInTrackster[tstId]) {
	  float cpFraction = 0.f;
	  if (!hitWithNoSTS) {
	    const auto& findSTSIt = std::find(
					      detIdSimTSId_Map[elemId].begin(),
					      detIdSimTSId_Map[elemId].end(),
					      detIdInfoInCluster{
						stsPair.first, 0.f});  // only the first element is used for the matching (overloaded operator==)
	    if (findSTSIt != detIdSimTSId_Map[elemId].end())
	      cpFraction = findSTSIt->fraction;
	  }
	  if (stsPair.second == FLT_MAX) {
	    stsPair.second = 0.f;
	  }
	  stsPair.second +=
	    std::min(pow(rhFraction - cpFraction, 2), pow(rhFraction, 2)) * hitEnergyWeight * invTracksterEnergyWeight;
	}
      }  // end of loop through trackster rechits

      //In case of a Trackster with some energy but none related CaloParticles print some info.
      /* if (stsInTrackster[tstId].empty()) */
      /* 	std::cout << "Trackster Id: " << tstId << "\tSimTrackster id: -1" */
      /* 		  << "\tscore: -1\n"; */

      tracksters_FakeMerge[tstId] =
        std::count_if(std::begin(stsInTrackster[tstId]),
                      std::end(stsInTrackster[tstId]),
                      [ScoreCutTStoSTSFakeMerge](const auto& obj) { return obj.second < ScoreCutTStoSTSFakeMerge; });

      const auto score = std::min_element(std::begin(stsInTrackster[tstId]),
					  std::end(stsInTrackster[tstId]),
					  [](const auto& obj1, const auto& obj2) { return obj1.second < obj2.second; });
      float score2 = -1;
      float sharedEneFrac2 = 0;

      //Getting ready to save info for all related simtrackster. 
      std::vector<unsigned int> tmp_trackster_cpId;
      std::vector<unsigned int> tmp_trackster_scId;
      std::vector<unsigned int> tmp_trackster_id;
      std::vector<unsigned int> tmp_trackster_numofvertices;
      std::vector<float> tmp_trackster_score_trackster2caloparticle;
      std::vector<float> tmp_trackster_sharedenergy_trackster2caloparticle;
      std::vector<float> tmp_trackster_score_trackster2bestCaloparticle;
      std::vector<float> tmp_trackster_sharedenergy_trackster2bestCaloparticle;
      std::vector<float> tmp_trackster_trackster2bestCaloparticle_eta;
      std::vector<float> tmp_trackster_trackster2bestCaloparticle_phi;
      std::vector<float> tmp_trackster_score_trackster2bestCaloparticle2;
      std::vector<float> tmp_trackster_sharedenergy_trackster2bestCaloparticle2;
      tmp_trackster_cpId.clear();
      tmp_trackster_scId.clear();
      tmp_trackster_id.clear();
      tmp_trackster_numofvertices.clear();
      tmp_trackster_score_trackster2caloparticle.clear();
      tmp_trackster_sharedenergy_trackster2caloparticle.clear();
      tmp_trackster_score_trackster2bestCaloparticle.clear();
      tmp_trackster_sharedenergy_trackster2bestCaloparticle.clear();
      tmp_trackster_trackster2bestCaloparticle_eta.clear();
      tmp_trackster_trackster2bestCaloparticle_phi.clear();
      tmp_trackster_score_trackster2bestCaloparticle2.clear();
      tmp_trackster_sharedenergy_trackster2bestCaloparticle2.clear();
      //Also the per Trackster output above we want to save it per simtrackster
      //to end up with a single dataframe.
      std::vector<float> tmp_tmp_trackster_raw_energy;
      std::vector<unsigned int> tmp_tmp_trackster_numberOfHitsInTS;
      std::vector<unsigned int> tmp_tmp_trackster_numberOfNoiseHitsInTS;
      std::vector<int> tmp_tmp_trackster_maxCPId_byNumberOfHits;
      std::vector<unsigned int> tmp_tmp_trackster_maxCPNumberOfHitsInTS;
      std::vector<int> tmp_tmp_trackster_maxCPId_byEnergy;
      std::vector<float> tmp_tmp_trackster_maxEnergySharedTSandCP;
      std::vector<float> tmp_tmp_trackster_totalCPEnergyFromLayerCP;
      std::vector<float> tmp_tmp_trackster_energyFractionOfTSinCP;
      std::vector<float> tmp_tmp_trackster_energyFractionOfCPinTS;
      tmp_tmp_trackster_raw_energy.clear();
      tmp_tmp_trackster_numberOfHitsInTS.clear();
      tmp_tmp_trackster_numberOfNoiseHitsInTS.clear();
      tmp_tmp_trackster_maxCPId_byNumberOfHits.clear();
      tmp_tmp_trackster_maxCPNumberOfHitsInTS.clear();
      tmp_tmp_trackster_maxCPId_byEnergy.clear();
      tmp_tmp_trackster_maxEnergySharedTSandCP.clear();
      tmp_tmp_trackster_totalCPEnergyFromLayerCP.clear();
      tmp_tmp_trackster_energyFractionOfTSinCP.clear();
      tmp_tmp_trackster_energyFractionOfCPinTS.clear();

    
      for (const auto& stsPair : stsInTrackster[tstId]) {
	const auto iSTS = stsPair.first;
	const auto iScore = stsPair.second;
	const auto cpId = getCPId(simTSs[iSTS], iSTS, cPHandle_id, cpToSc_SimTrackstersMap, simTSs_fromCP);
	//if (std::find(cPIndices.begin(), cPIndices.end(), cpId) == cPIndices.end())
	//  continue;
	auto iSim = simTSs[iSTS].seedIndex();
	if (simTSs[iSTS].seedID() == cPHandle_id)  // SimTrackster from CaloParticle
	  iSim = 0;
	const auto& simOnLayer = (i == 0) ? cPOnLayer[cpId] : sCOnLayer[cpId][iSim];

	float sharedeneCPallLayers = 0.;
	//for (unsigned int j = 0; j < layers * 2; ++j)
	sharedeneCPallLayers += simOnLayer.layerClusterIdToEnergyAndScore.count(tstId)
	  ? simOnLayer.layerClusterIdToEnergyAndScore.at(tstId).first
	  : 0;
	if (tracksterEnergy == 0)
	  continue;
	const auto sharedEneFrac = sharedeneCPallLayers / tracksterEnergy;
	/* std::cout << "\nTrackster id: " << tstId << " (" << tst.vertices().size() << " vertices)" */
	/*                            << "\tSimTrackster Id: " << iSTS << " (" << simTSs[iSTS].vertices().size() */
	/*                            << " vertices)" */
	/*                            << " (CP id: " << cpId << ")\tscore: " << iScore */
	/*                            << "\tsharedeneCPallLayers: " << sharedeneCPallLayers << std::endl; */
	tmp_trackster_cpId.push_back(cpId);
	tmp_trackster_scId.push_back(iSim);
	tmp_trackster_id.push_back(tstId);
	tmp_trackster_numofvertices.push_back(tst.vertices().size());

	tmp_tmp_trackster_raw_energy.push_back(tmp_trackster_raw_energy[tstId]);
	tmp_tmp_trackster_numberOfHitsInTS.push_back(tmp_trackster_numberOfHitsInTS[tstId]);
	tmp_tmp_trackster_numberOfNoiseHitsInTS.push_back(tmp_trackster_numberOfNoiseHitsInTS[tstId]);
	tmp_tmp_trackster_maxCPId_byNumberOfHits.push_back(tmp_trackster_maxCPId_byNumberOfHits[tstId]);
	tmp_tmp_trackster_maxCPNumberOfHitsInTS.push_back(tmp_trackster_maxCPNumberOfHitsInTS[tstId]);
	tmp_tmp_trackster_maxCPId_byEnergy.push_back(tmp_trackster_maxCPId_byEnergy[tstId]);
	tmp_tmp_trackster_maxEnergySharedTSandCP.push_back(tmp_trackster_maxEnergySharedTSandCP[tstId]);
	tmp_tmp_trackster_totalCPEnergyFromLayerCP.push_back(tmp_trackster_totalCPEnergyFromLayerCP[tstId]);
	tmp_tmp_trackster_energyFractionOfTSinCP.push_back(tmp_trackster_energyFractionOfTSinCP[tstId]);
	tmp_tmp_trackster_energyFractionOfCPinTS.push_back(tmp_trackster_energyFractionOfCPinTS[tstId]);

	tmp_trackster_score_trackster2caloparticle.push_back(iScore);
	tmp_trackster_sharedenergy_trackster2caloparticle.push_back(sharedEneFrac);
 
	if (iSTS == score->first) {

	  tmp_trackster_score_trackster2bestCaloparticle.push_back(iScore);
	  tmp_trackster_sharedenergy_trackster2bestCaloparticle.push_back(sharedEneFrac);
	  tmp_trackster_trackster2bestCaloparticle_eta.push_back(tst.barycenter().eta());
	  tmp_trackster_trackster2bestCaloparticle_phi.push_back(tst.barycenter().phi());

	} else if (score2 < 0 || iScore < score2) {
	  score2 = iScore;
	  sharedEneFrac2 = sharedEneFrac;
	}
      }  // end of loop through SimTracksters associated to Trackster
      if (score2 > -1) {

	tmp_trackster_score_trackster2bestCaloparticle2.push_back(score2);
	tmp_trackster_sharedenergy_trackster2bestCaloparticle2.push_back(sharedEneFrac2);
       
      }

      trInfo.trackster_cpId.push_back(tmp_trackster_cpId);
      trInfo.trackster_scId.push_back(tmp_trackster_scId);
      trInfo.trackster_Id.push_back(tmp_trackster_id);
      trInfo.trackster_numofvertices.push_back(tmp_trackster_numofvertices);
    
      trInfo.trackster_score_trackster2caloparticle.push_back(tmp_trackster_score_trackster2caloparticle);
      trInfo.trackster_sharedenergy_trackster2caloparticle.push_back(tmp_trackster_sharedenergy_trackster2caloparticle);
      trInfo.trackster_score_trackster2bestCaloparticle.push_back(tmp_trackster_score_trackster2bestCaloparticle);
      trInfo.trackster_sharedenergy_trackster2bestCaloparticle.push_back(tmp_trackster_sharedenergy_trackster2bestCaloparticle);
      trInfo.trackster_trackster2bestCaloparticle_eta.push_back(tmp_trackster_trackster2bestCaloparticle_eta);
      trInfo.trackster_trackster2bestCaloparticle_phi.push_back(tmp_trackster_trackster2bestCaloparticle_phi);
      trInfo.trackster_score_trackster2bestCaloparticle2.push_back(tmp_trackster_score_trackster2bestCaloparticle2);
      trInfo.trackster_sharedenergy_trackster2bestCaloparticle2.push_back(tmp_trackster_sharedenergy_trackster2bestCaloparticle2);
      trInfo.trackster_Raw_Energy.push_back(tmp_tmp_trackster_raw_energy);
      trInfo.trackster_numberOfHitsInTS.push_back(tmp_tmp_trackster_numberOfHitsInTS);
      trInfo.trackster_numberOfNoiseHitsInTS.push_back(tmp_tmp_trackster_numberOfNoiseHitsInTS);
      trInfo.trackster_maxCPId_byNumberOfHits.push_back(tmp_tmp_trackster_maxCPId_byNumberOfHits);
      trInfo.trackster_maxCPNumberOfHitsInTS.push_back(tmp_tmp_trackster_maxCPNumberOfHitsInTS);
      trInfo.trackster_maxCPId_byEnergy.push_back(tmp_tmp_trackster_maxCPId_byEnergy);
      trInfo.trackster_maxEnergySharedTSandCP.push_back(tmp_tmp_trackster_maxEnergySharedTSandCP);
      trInfo.trackster_totalCPEnergyFromLayerCP.push_back(tmp_tmp_trackster_totalCPEnergyFromLayerCP);
      trInfo.trackster_energyFractionOfTSinCP.push_back(tmp_tmp_trackster_energyFractionOfTSinCP);
      trInfo.trackster_energyFractionOfCPinTS.push_back(tmp_tmp_trackster_energyFractionOfCPinTS);

    
    }  // end of loop through Tracksters

    std::unordered_map<unsigned int, std::vector<float>> score3d;
    std::unordered_map<unsigned int, std::vector<float>> tstSharedEnergy;

    for (unsigned int iSTS = 0; iSTS < nSimTracksters; ++iSTS) {
      score3d[iSTS].resize(nTracksters);
      tstSharedEnergy[iSTS].resize(nTracksters);
      for (unsigned int j = 0; j < nTracksters; ++j) {
	score3d[iSTS][j] = FLT_MAX;
	tstSharedEnergy[iSTS][j] = 0.f;
      }
    }

    // Fill the plots to compute the different metrics linked to
    // gen-level, namely efficiency, purity and duplicate. In this loop should restrict
    // only to the selected caloParaticles.
    for (unsigned int iSTS = 0; iSTS < nSimTracksters; ++iSTS) {
      const auto& sts = simTSs[iSTS];
      const auto& cpId = getCPId(sts, iSTS, cPHandle_id, cpToSc_SimTrackstersMap, simTSs_fromCP);
      if (i == 0 && std::find(cPSelectedIndices.begin(), cPSelectedIndices.end(), cpId) == cPSelectedIndices.end())
	continue;

      const auto& hafLC = apply_LCMultiplicity(sts, layerClusters);
      float SimEnergy_LC = 0.f;
      for (const auto& haf : hafLC) {
	const auto lcId = getLCId(sts.vertices(), layerClusters, haf.first);
	const auto iLC = std::find(sts.vertices().begin(), sts.vertices().end(), lcId);
	SimEnergy_LC +=
          hitMap.at(haf.first)->energy() / sts.vertex_multiplicity(std::distance(std::begin(sts.vertices()), iLC));
      }

      auto iSim = sts.seedIndex();
      if (sts.seedID() == cPHandle_id)  // SimTrackster from CaloParticle
	iSim = 0;
      auto& simOnLayer = (i == 0) ? cPOnLayer[cpId] : sCOnLayer[cpId][iSim];

      // Keep the Trackster ids that are related to
      // SimTrackster under study for the final filling of the score
      std::set<unsigned int> stsId_tstId_related;
      auto& score3d_iSTS = score3d[iSTS];

      float SimEnergy = 0.f;
      float SimEnergyWeight = 0.f, hitsEnergyWeight = 0.f;
      //for (unsigned int layerId = 0; layerId < 1/*layers * 2*/; ++layerId) {
      const auto SimNumberOfHits = simOnLayer.hits_and_fractions.size();
      if (SimNumberOfHits == 0)
	continue;
      SimEnergy += simOnLayer.energy;
      //int tstWithMaxEnergyInCP = -1;
      //This is the maximum energy related to Trackster per layer.
      //float maxEnergyTSperlayerinSim = 0.f;
      //float SimEnergyFractionInTSperlayer = 0.f;
      //Remember and not confused by name. layerClusterIdToEnergyAndScore contains the Trackster id.
      //for (const auto& tst : simOnLayer.layerClusterIdToEnergyAndScore) {
	//if (tst.second.first > maxEnergyTSperlayerinSim) {
	//  maxEnergyTSperlayerinSim = tst.second.first;
	  //tstWithMaxEnergyInCP = tst.first;
	//}
      //}
      //if (SimEnergy > 0.f)
	//SimEnergyFractionInTSperlayer = maxEnergyTSperlayerinSim / SimEnergy;

      /* std::cout << std::setw(12) << "caloparticle\t" << std::setw(15) << "cp total energy\t" */
      /* 				 << std::setw(15) << "cpEnergyOnLayer\t" << std::setw(14) << "CPNhitsOnLayer\t" */
      /* 				 << std::setw(18) << "tstWithMaxEnergyInCP\t" << std::setw(15) << "maxEnergyTSinCP\t" */
      /* 				 << std::setw(20) << "CPEnergyFractionInTS" */
      /* 				 << "\n"; */
      /* std::cout << std::setw(12) << cpId << "\t" << std::setw(15) << sts.raw_energy() << "\t" */
      /* 				 << std::setw(15) << SimEnergy << "\t" << std::setw(14) << SimNumberOfHits << "\t" */
      /* 				 << std::setw(18) << tstWithMaxEnergyInCP << "\t" << std::setw(15) */
      /* 				 << maxEnergyTSperlayerinSim << "\t" << std::setw(20) << SimEnergyFractionInTSperlayer */
      /* 				 << "\n"; */

      for (const auto& haf : ((i == 0) ? simOnLayer.hits_and_fractions : hafLC)) {
	const auto& hitDetId = haf.first;
	// Compute the correct normalization
	// Need to loop on the simOnLayer data structure since this is the
	// only one that has the compressed information for multiple usage
	// of the same DetId by different SimClusters by a single CaloParticle.
	SimEnergyWeight += pow(haf.second * hitMap.at(hitDetId)->energy(), 2);

	const auto lcId = getLCId(sts.vertices(), layerClusters, hitDetId);
	float cpFraction = 0.f;
	if (i == 0) {
	  cpFraction = haf.second;
	} else {
	  const auto iLC = std::find(sts.vertices().begin(), sts.vertices().end(), lcId);
	  cpFraction = 1.f / sts.vertex_multiplicity(std::distance(std::begin(sts.vertices()), iLC));
	}
	if (cpFraction == 0.f)
	  continue;  // hopefully this should never happen

	bool hitWithNoTS = false;
	if (detIdToTracksterId_Map.find(hitDetId) == detIdToTracksterId_Map.end())
	  hitWithNoTS = true;
	const HGCRecHit* hit = hitMap.find(hitDetId)->second;
	const auto hitEnergyWeight = pow(hit->energy(), 2);
	hitsEnergyWeight += pow(cpFraction, 2) * hitEnergyWeight;

	for (auto& tsPair : simOnLayer.layerClusterIdToEnergyAndScore) {
	  const auto tstId = tsPair.first;
	  stsId_tstId_related.insert(tstId);

	  float tstFraction = 0.f;
	  if (!hitWithNoTS) {
	    const auto findTSIt =
              std::find(detIdToTracksterId_Map[hitDetId].begin(),
                        detIdToTracksterId_Map[hitDetId].end(),
                        detIdInfoInTrackster{
			  tstId, 0, 0.f});  // only the first element is used for the matching (overloaded operator==)
	    if (findTSIt != detIdToTracksterId_Map[hitDetId].end()) {
	      if (i == 0) {
		tstFraction = findTSIt->fraction;
	      } else {
		const auto iLC = std::find(
					   tracksters[tstId].vertices().begin(), tracksters[tstId].vertices().end(), findTSIt->clusterId);
		if (iLC != tracksters[tstId].vertices().end()) {
		  tstFraction = 1.f / tracksters[tstId].vertex_multiplicity(
									    std::distance(std::begin(tracksters[tstId].vertices()), iLC));
		}
	      }
	    }
	  }
	  // Here do not divide as before by the trackster energy weight. Should sum first
	  // over all layers and divide with the total CP energy over all layers.
	  if (tsPair.second.second == FLT_MAX) {
	    tsPair.second.second = 0.f;
	  }
	  tsPair.second.second += std::min(pow(tstFraction - cpFraction, 2), pow(cpFraction, 2)) * hitEnergyWeight;

	  /* std::cout << "\nTracksterId:\t" << tstId << "\tSimTracksterId:\t" << iSTS << "\tcpId:\t" */
	  /* 			     << cpId << "\ttstfraction, cpfraction:\t" << tstFraction << ", " << cpFraction */
	  /* 			     << "\thitEnergyWeight:\t" << hitEnergyWeight << "\tadded delta:\t" */
	  /* 			     << pow((tstFraction - cpFraction), 2) * hitEnergyWeight */
	  /* 			     << "\tcurrent Sim-score numerator:\t" << tsPair.second.second */
	  /* 			     << "\tshared Sim energy:\t" << tsPair.second.first << '\n'; */
	}
      }  // end of loop through SimCluster SimHits on current layer

      /* if (simOnLayer.layerClusterIdToEnergyAndScore.empty()) */
      /* 	std::cout << "CP Id:\t" << cpId << "\tTS id:\t-1" */
      /* 				   << " Sub score in \t -1\n"; */

      for (const auto& tsPair : simOnLayer.layerClusterIdToEnergyAndScore) {
	const auto tstId = tsPair.first;
	// 3D score here without the denominator at this point
	if (score3d_iSTS[tstId] == FLT_MAX) {
	  score3d_iSTS[tstId] = 0.f;
	}
	score3d_iSTS[tstId] += tsPair.second.second;
	tstSharedEnergy[iSTS][tstId] += tsPair.second.first;
      }
      //} // end of loop through layers

      const auto scoreDenom = (i == 0) ? SimEnergyWeight : hitsEnergyWeight;
      const auto energyDenom = (i == 0) ? SimEnergy : SimEnergy_LC;

      const auto sts_eta = sts.barycenter().eta();
      const auto sts_phi = sts.barycenter().phi();
      const auto sts_en = sts.raw_energy();
      const auto sts_pt = sts.raw_pt();

      //Loop through related Tracksters here
      //Getting ready to save info for all related tracksters. 
      std::vector<unsigned int> tmp_trackster_sts_cpId;
      std::vector<unsigned int> tmp_trackster_sts_id;
      std::vector<unsigned int> tmp_trackster_sts_ts_id;
      std::vector<float> tmp_trackster_sts_SimEnergy;
      std::vector<float> tmp_trackster_sts_SimEnergyWeight;
      std::vector<float> tmp_trackster_sts_trackster_raw_energy;
      std::vector<float> tmp_trackster_sts_score_caloparticle2trackster;
      std::vector<float> tmp_trackster_sts_sharedenergy_caloparticle2trackster;
      std::vector<float> tmp_trackster_sts_eta;
      std::vector<float> tmp_trackster_sts_phi;
      std::vector<float> tmp_trackster_sts_pt;
      std::vector<float> tmp_trackster_sts_raw_energy;
      std::vector<float> tmp_trackster_sts_scorePur_caloparticle2trackster;
      std::vector<float> tmp_trackster_sts_sharedenergy_caloparticle2trackster_assoc;
      std::vector<float> tmp_trackster_sts_besttrackster_raw_energy;
      std::vector<float> tmp_trackster_sts_scoreDupl_caloparticle2trackster;
      std::vector<float> tmp_trackster_sts_sharedenergy_caloparticle2trackster_assoc2;

      tmp_trackster_sts_cpId.clear();
      tmp_trackster_sts_id.clear();
      tmp_trackster_sts_ts_id.clear();
      tmp_trackster_sts_SimEnergy.clear();
      tmp_trackster_sts_SimEnergyWeight.clear();
      tmp_trackster_sts_trackster_raw_energy.clear();
      tmp_trackster_sts_score_caloparticle2trackster.clear();
      tmp_trackster_sts_sharedenergy_caloparticle2trackster.clear();
      tmp_trackster_sts_eta.clear();
      tmp_trackster_sts_phi.clear();
      tmp_trackster_sts_pt.clear();
      tmp_trackster_sts_raw_energy.clear();
      tmp_trackster_sts_scorePur_caloparticle2trackster.clear();
      tmp_trackster_sts_sharedenergy_caloparticle2trackster_assoc.clear();
      tmp_trackster_sts_besttrackster_raw_energy.clear();
      tmp_trackster_sts_scoreDupl_caloparticle2trackster.clear();
      tmp_trackster_sts_sharedenergy_caloparticle2trackster_assoc2.clear();

      // In case the threshold to associate a CaloParticle to a Trackster is
      // below 50%, there could be cases in which the CP is linked to more than
      // one tracksters, leading to efficiencies >1. This boolean is used to
      // avoid "over counting".
      bool cp_considered_efficient = false;
      bool cp_considered_pure = false;
      for (const auto tstId : stsId_tstId_related) {
	// Now time for the denominator
	score3d_iSTS[tstId] /= scoreDenom;
	const auto tstSharedEnergyFrac = tstSharedEnergy[iSTS][tstId] / energyDenom;
	/* std::cout << "STS id: " << iSTS << "\t(CP id: " << cpId << ")\tTS id: " << tstId */
	/* 			   << "\nSimEnergy: " << energyDenom << "\tSimEnergyWeight: " << SimEnergyWeight */
	/* 			   << "\tTrackster energy: " << tracksters[tstId].raw_energy() */
	/* 			   << "\nscore: " << score3d_iSTS[tstId] */
	/* 			   << "\tshared energy: " << tstSharedEnergy[iSTS][tstId] */
	/* 			   << "\tshared energy fraction: " << tstSharedEnergyFrac << "\n"; */

	tmp_trackster_sts_cpId.push_back(cpId);
	tmp_trackster_sts_id.push_back(iSTS);
	tmp_trackster_sts_ts_id.push_back(tstId);
	tmp_trackster_sts_SimEnergy.push_back(energyDenom);
	tmp_trackster_sts_SimEnergyWeight.push_back(SimEnergyWeight);
	tmp_trackster_sts_trackster_raw_energy.push_back(tracksters[tstId].raw_energy());
	tmp_trackster_sts_score_caloparticle2trackster.push_back(score3d_iSTS[tstId]);
	tmp_trackster_sts_sharedenergy_caloparticle2trackster.push_back(tstSharedEnergyFrac);
	tmp_trackster_sts_eta.push_back(sts_eta);
	tmp_trackster_sts_phi.push_back(sts_phi);
	tmp_trackster_sts_pt.push_back(sts_pt);
	tmp_trackster_sts_raw_energy.push_back(sts_en);

      }  // end of loop through Tracksters related to SimTrackster

      const auto best = std::min_element(std::begin(score3d_iSTS), std::end(score3d_iSTS));
      if (best != score3d_iSTS.end()) {
	const auto bestTstId = std::distance(std::begin(score3d_iSTS), best);
	const auto bestTstSharedEnergyFrac = tstSharedEnergy[iSTS][bestTstId] / energyDenom;

	tmp_trackster_sts_scorePur_caloparticle2trackster.push_back(*best);
	tmp_trackster_sts_sharedenergy_caloparticle2trackster_assoc.push_back(bestTstSharedEnergyFrac);
	tmp_trackster_sts_besttrackster_raw_energy.push_back(tracksters[bestTstId].raw_energy());

	/* std::cout << count << " " << sts_eta << " " << sts_phi << " " */
	/* 			   << tracksters[bestTstId].raw_energy() << " " << sts.raw_energy() << " " */
	/* 			   << bestTstSharedEnergyFrac << "\n"; */

	if (score3d_iSTS.size() > 1) {
	  auto best2 = (best == score3d_iSTS.begin()) ? std::next(best, 1) : score3d_iSTS.begin();
	  for (auto tstId = score3d_iSTS.begin(); tstId != score3d_iSTS.end() && tstId != best; tstId++)
	    if (*tstId < *best2)
	      best2 = tstId;
	  const auto best2TstId = std::distance(std::begin(score3d_iSTS), best2);
	  const auto best2TstSharedEnergyFrac = tstSharedEnergy[iSTS][best2TstId] / energyDenom;

	  tmp_trackster_sts_scoreDupl_caloparticle2trackster.push_back(*best2);
	  tmp_trackster_sts_sharedenergy_caloparticle2trackster_assoc2.push_back(best2TstSharedEnergyFrac);

	}
      }

      trInfo.trackster_sts_cpId.push_back(tmp_trackster_sts_cpId);
      trInfo.trackster_sts_id.push_back(tmp_trackster_sts_id);
      trInfo.trackster_sts_ts_id.push_back(tmp_trackster_sts_ts_id);
      trInfo.trackster_sts_SimEnergy.push_back(tmp_trackster_sts_SimEnergy);
      trInfo.trackster_sts_SimEnergyWeight.push_back(tmp_trackster_sts_SimEnergyWeight);
      trInfo.trackster_sts_trackster_raw_energy.push_back(tmp_trackster_sts_trackster_raw_energy);
      trInfo.trackster_sts_score_caloparticle2trackster.push_back(tmp_trackster_sts_score_caloparticle2trackster);
      trInfo.trackster_sts_sharedenergy_caloparticle2trackster.push_back(tmp_trackster_sts_sharedenergy_caloparticle2trackster);
      trInfo.trackster_sts_eta.push_back(tmp_trackster_sts_eta);
      trInfo.trackster_sts_phi.push_back(tmp_trackster_sts_phi);
      trInfo.trackster_sts_pt.push_back(tmp_trackster_sts_pt);
      trInfo.trackster_sts_raw_energy.push_back(tmp_trackster_sts_raw_energy);
      trInfo.trackster_sts_scorePur_caloparticle2trackster.push_back(tmp_trackster_sts_scorePur_caloparticle2trackster);
      trInfo.trackster_sts_sharedenergy_caloparticle2trackster_assoc.push_back(tmp_trackster_sts_sharedenergy_caloparticle2trackster_assoc);
      trInfo.trackster_sts_besttrackster_raw_energy.push_back(tmp_trackster_sts_besttrackster_raw_energy);
      trInfo.trackster_sts_scoreDupl_caloparticle2trackster.push_back(tmp_trackster_sts_scoreDupl_caloparticle2trackster);
      trInfo.trackster_sts_sharedenergy_caloparticle2trackster_assoc2.push_back(tmp_trackster_sts_sharedenergy_caloparticle2trackster_assoc2);

      
    }  // end of loop through SimTracksters


  
  } // end of prepare_tracksters_to_SimTracksters


  



} // end namespace

#endif
