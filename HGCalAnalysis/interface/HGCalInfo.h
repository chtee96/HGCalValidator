#ifndef Validation_HGCalValidation_HGCalInfo_h
#define Validation_HGCalValidation_HGCalInfo_h

#include <vector>

#include <TTree.h>

// Forward declarations
class EncodedEventId;

namespace hgcal_validation
{

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
			   std::unordered_map<DetId, const HGCRecHit*> const& hitMap,
			   std::shared_ptr<hgcal::RecHitTools> recHitTools, 
			   unsigned int layers){


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

	  matched_hits.push_back(sh_detid.rawId());
	  //std::cout << "1 " << sh_detid.rawId() << std::endl;
	  /* std::cout << recHitTools->getLayerWithOffset(itcheck->first) + layers * ((recHitTools->zside(itcheck->first) + 1) >> 1) - 1 << std::endl; */
	  /* std::cout << "lastLayerEE() " << recHitTools->lastLayerEE() << " lastLayerFH() " << recHitTools->lastLayerFH() << " firstLayerBH() " << recHitTools->firstLayerBH() << " lastLayerBH() " << recHitTools->lastLayerBH() << " recHitTools->zside " << recHitTools->zside(itcheck->first) << std::endl; */

	  rhInfo.rechit_eta.push_back(recHitTools->getEta(position));
	  rhInfo.rechit_phi.push_back(recHitTools->getPhi(position));
	  rhInfo.rechit_pt.push_back(recHitTools->getPt(position,energy));
	  rhInfo.rechit_energy.push_back(energy);
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
			     unsigned int layers){

    
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
  			  unsigned int layers){

    const auto nTracksters = tracksters.size();

    /* auto apply_LCMultiplicity = [](const ticl::Trackster& trackster, const reco::CaloClusterCollection& layerClusters) { */
    /*   std::vector<std::pair<DetId, float>> hits_and_fractions_norm; */
    /*   int lcInTst = 0; */
    /*   std::for_each(std::begin(trackster.vertices()), std::end(trackster.vertices()), [&](unsigned int idx) { */
    /* 	  const auto fraction = 1.f / trackster.vertex_multiplicity(lcInTst++); */
    /* 	  for (const auto& cell : layerClusters[idx].hitsAndFractions()) { */
    /* 	    hits_and_fractions_norm.emplace_back(cell.first, cell.second * fraction); // cell.second is the hit fraction in the layerCluster */
    /* 	  } */
    /* 	}); */
    /*   return hits_and_fractions_norm; */
    /* }; */

    /* auto getLCId = [](const std::vector<unsigned int>& tst_vertices, const reco::CaloClusterCollection& layerClusters, const DetId& hitid) { */
    /*   unsigned int lcId = -1; */
    /*   std::for_each(std::begin(tst_vertices), std::end(tst_vertices), [&](unsigned int idx) { */
    /* 	  const auto& lc_haf = layerClusters[idx].hitsAndFractions(); */
    /* 	  const auto& hitFound = std::find_if(std::begin(lc_haf), std::end(lc_haf), [&hitid](const std::pair<DetId, float>& v) { */
    /* 	      return v.first == hitid; }); */
    /* 	  if (hitFound != lc_haf.end()) // not all hits may be clusterized */
    /* 	    lcId = idx; */
    /* 	}); */
    /*   return lcId; */
    /* }; */

    // Loop through Tracksters
    for (unsigned int tstId = 0; tstId < nTracksters; ++tstId) {
      const auto& trackster = tracksters[tstId];
      //Will check later also if there are tracksters with zero LCs. 
      if (trackster.vertices().empty())
  	continue;

      //For the layers the Trackster expands to. Will use a set because there would be many
      //duplicates and then go back to vector for random access, since they say it is faster.
      std::set<int> trackster_layers;
      
      //Loop through layer clusters of the trackster
      for (const auto lcId : tracksters[tstId].vertices()) {
  	const std::vector<std::pair<DetId, float>>& hits_and_fractions = clusters[lcId].hitsAndFractions();
  	const auto firstHitDetId = hits_and_fractions[0].first;
  	//The layer that the layer cluster belongs to
  	int layerid =
  	  recHitTools->getLayerWithOffset(firstHitDetId) + layers * ((recHitTools->zside(firstHitDetId) + 1) >> 1) - 1;

  	trackster_layers.insert(layerid);

	fillSimClustersInfo(scInfo, rhInfo_sc, simClusters, hitMap, recHitTools, layers);
	fillLayerClustersInfo(lcInfo, rhInfo_lc, clusters, lcId, tstId, true, UnCalibHitMap, hitMap, recHitTools, layers);

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

      }
 
      /* const auto tst_hitsAndFractions = apply_LCMultiplicity(trackster, clusters); */
      /* const auto numberOfHitsInTS = tst_hitsAndFractions.size(); */

      //Loop through the hits of the trackster under study
      /* for (unsigned int hitId = 0; hitId < numberOfHitsInTS; hitId++) { */
      /* 	const auto rh_detid = tst_hitsAndFractions[hitId].first; */
      /* 	const auto rhFraction = tst_hitsAndFractions[hitId].second; */

      /* 	/\* const auto lcId_r = getLCId(trackster.vertices(), clusters, rh_detid); *\/ */
      /* 	/\* const auto iLC_r = std::find(trackster.vertices().begin(), trackster.vertices().end(), lcId_r); *\/ */
      /* 	/\* const auto lcFraction_r = 1.f / trackster.vertex_multiplicity(std::distance(std::begin(trackster.vertices()), iLC_r)); *\/ */

      /* } //end of loop over hits of the trackster */



    } //end of loop over tracksters


  }



} // end namespace

#endif
