#include "NtupleMakerAlg.h"
#include "TrkTrack/Track.h"
#include "TrackerRIO_OnTrack/FaserSCT_ClusterOnTrack.h"
#include "TrackerIdentifier/FaserSCT_ID.h"
#include "ScintIdentifier/VetoNuID.h"
#include "ScintIdentifier/VetoID.h"
#include "TrackerPrepRawData/FaserSCT_Cluster.h"
#include "Identifier/Identifier.h"
#include "TrackerReadoutGeometry/SCT_DetectorManager.h"
#include "TrackerReadoutGeometry/SiDetectorElement.h"
#include "TrackerPrepRawData/FaserSCT_Cluster.h"
#include "xAODTruth/TruthParticle.h"
#include <cmath>



NtupleMakerAlg::NtupleMakerAlg(const std::string &name, 
    ISvcLocator *pSvcLocator)
  : AthReentrantAlgorithm(name, pSvcLocator), 
  AthHistogramming(name),
  m_histSvc("THistSvc/THistSvc", name) {}


StatusCode NtupleMakerAlg::initialize() 
{
  ATH_CHECK(m_truthEventContainer.initialize());
  ATH_CHECK(m_truthParticleContainer.initialize());
  ATH_CHECK(m_trackCollection.initialize());
  ATH_CHECK(m_vetoNuContainer.initialize());
  ATH_CHECK(m_vetoContainer.initialize());
  ATH_CHECK(m_clusterContainer.initialize());
  ATH_CHECK(m_simDataCollection.initialize());

  ATH_CHECK(detStore()->retrieve(m_sctHelper,       "FaserSCT_ID"));
  ATH_CHECK(detStore()->retrieve(m_vetoNuHelper,    "VetoNuID"));
  ATH_CHECK(detStore()->retrieve(m_vetoHelper,      "VetoID"));

  ATH_CHECK(detStore()->retrieve(m_detMgr, "SCT"));

  m_tree = new TTree("tree", "tree");
  m_tree->Branch("run_number", &m_run_number, "run_number/I");
  m_tree->Branch("event_number", &m_event_number, "event_number/I");

  m_tree->Branch("VetoNuPmt0", &m_vetoNu0, "vetoNu0/D");
  m_tree->Branch("VetoNuPmt1", &m_vetoNu1, "vetoNu1/D");

  m_tree->Branch("VetoPmt00",  &m_veto00,  "veto00/D");
  m_tree->Branch("VetoPmt01",  &m_veto01,  "veto01/D");
  m_tree->Branch("VetoUpstream", &m_vetoUpstream, "vetoUpstream/D");
  m_tree->Branch("VetoPmt10",  &m_veto10,  "veto10/D");
  m_tree->Branch("VetoPmt11",  &m_veto11,  "veto11/D");
  m_tree->Branch("VetoDownstream", &m_vetoDownstream, "vetoDownstream/D");

  m_tree->Branch("sta0nclus", &m_sta0nclus, "sta0nclus/I");
  m_tree->Branch("sta1nclus", &m_sta1nclus, "sta1nclus/I");
  m_tree->Branch("sta2nclus", &m_sta2nclus, "sta2nclus/I");
  m_tree->Branch("sta3nclus", &m_sta3nclus, "sta3nclus/I");

  m_tree->Branch("x", &m_x);
  m_tree->Branch("y", &m_y);
  m_tree->Branch("z", &m_z);
  m_tree->Branch("px", &m_px);
  m_tree->Branch("py", &m_py);
  m_tree->Branch("pz", &m_pz);
  m_tree->Branch("p", &m_p);
  m_tree->Branch("charge", &m_charge);
  m_tree->Branch("chi2", &m_chi2);
  m_tree->Branch("ndof", &m_ndof);
  m_tree->Branch("nclus0", &m_nclus0);
  m_tree->Branch("nclus1", &m_nclus1);
  m_tree->Branch("nclus2", &m_nclus2);
  m_tree->Branch("nclus3", &m_nclus3);
  m_tree->Branch("nclus", &m_nclus);
  m_tree->Branch("nlayer", &m_nlayer);
  m_tree->Branch("nlayer_ift", &m_nlayer_ift);
  m_tree->Branch("nsta", &m_nsta);
  m_tree->Branch("firstlayerId", &m_firstlayer);
  m_tree->Branch("lastlayerId", &m_lastlayer);
  m_tree->Branch("firstclus_x", &m_first_x);
  m_tree->Branch("firstclus_y", &m_first_y);
  m_tree->Branch("firstclus_z", &m_first_z);
  m_tree->Branch("lastclus_x", &m_last_x);
  m_tree->Branch("lastclus_y", &m_last_y);
  m_tree->Branch("lastclus_z", &m_last_z);
  m_tree->Branch("ntruthmatchedclus", &m_ntruthmatchedclus);
  m_tree->Branch("ntruthmatchedclus0", &m_ntruthmatchedclus0);
  m_tree->Branch("ntruthmatchedclus1", &m_ntruthmatchedclus1);
  m_tree->Branch("ntruthmatchedclus2", &m_ntruthmatchedclus2);
  m_tree->Branch("ntruthmatchedclus3", &m_ntruthmatchedclus3);
  m_tree->Branch("goodStations", &m_ngoodsta, "ngoodsta/I");
  m_tree->Branch("longTracks", &m_longTracks, "longTracks/I");
  m_tree->Branch("pTruthLepton", &m_truthLeptonMomentum, "pTruthLepton/D");
  m_tree->Branch("truthBarcode", &m_truthBarcode, "truthBarcode/I");
  m_tree->Branch("truthPdg", &m_truthPdg, "truthPdg/I");
  m_tree->Branch("ntruthmatchedclus_layer", &m_ntruthmatchedclus_layer);
  m_tree->Branch("ntruthmatchedclus_event", &m_ntruthmatchedclus_event, "truthmatchedclus/I");

  ATH_CHECK(histSvc()->regTree("/HIST2/tree", m_tree));

  return StatusCode::SUCCESS;
}


StatusCode NtupleMakerAlg::execute(const EventContext &ctx) const 
{

  clearTree();
    bool ifsave=false;

  m_run_number = ctx.eventID().run_number();
  m_event_number = ctx.eventID().event_number();

  bool realData = true;

  // Work out effective cross section for MC
  SG::ReadHandle<xAOD::TruthEventContainer> truthEventContainer { m_truthEventContainer, ctx };
  if (truthEventContainer.isValid() && truthEventContainer->size() > 0)
  {
    realData = false;
  }

  // Find the primary lepton (if any)
  SG::ReadHandle<xAOD::TruthParticleContainer> truthParticleContainer { m_truthParticleContainer, ctx };
  if (truthParticleContainer.isValid() && truthParticleContainer->size() > 0)
  {
    for (auto particle : *truthParticleContainer)
    {
      if ( particle->absPdgId() == 11 || particle->absPdgId() == 13 || particle->absPdgId() == 15 )
      {
	if (particle->status() == 1 && (particle->nParents() == 0 || particle->nParents() == 2) )
	  m_truthLeptonMomentum = particle->p4().P();
	  m_truthLeptonPx = particle->px().P();
	  m_truthLeptonPy = particle->py().P();
	  m_truthLeptonPz = particle->pz().P();
	break;
      }
    }
  }

  SG::ReadHandle<xAOD::WaveformHitContainer> vetoNuContainer { m_vetoNuContainer, ctx };
  ATH_CHECK(vetoNuContainer.isValid());

  // If real data, check for blinding before we do anything else

  for (auto hit : *vetoNuContainer)
  {
    if (!waveformHitOK(hit)) continue;
    auto id = hit->identify();
    if (m_vetoNuHelper->plate(id) == 0)
    {
      m_vetoNu0 += hit->integral();
    }
    else if (m_vetoNuHelper->plate(id) == 1)
    {
      m_vetoNu1 += hit->integral();
    }
    else
    {
      ATH_MSG_FATAL("Invalid VetoNu plate number: " << m_vetoNuHelper->plate(id));
      return StatusCode::FAILURE;
    }
  }

  SG::ReadHandle<xAOD::WaveformHitContainer> vetoContainer { m_vetoContainer, ctx };
  ATH_CHECK(vetoContainer.isValid());

  for (auto hit : *vetoContainer)
  {
    if (!waveformHitOK(hit)) continue;
    auto id = hit->identify();
    auto station = m_vetoHelper->station(id);
    auto plate   = m_vetoHelper->plate(id);
    if (station == 0)
    {
      if (plate == 0)
      {
	m_veto00 += hit->integral();
	m_vetoUpstream += hit->integral();
      }
      else if (plate == 1)
      {
	m_veto01 += hit->integral();
	m_vetoUpstream += hit->integral();
      }
      else
      {
	ATH_MSG_FATAL("Invalid Veto plate number: " << plate);
      }
    }
    else if (station == 1)
    {
      if (plate == 0)
      {
	m_veto10 += hit->integral();
	m_vetoDownstream += hit->integral();
      }
      else if (plate == 1)
      {
	m_veto11 += hit->integral();
	m_vetoDownstream += hit->integral();
      }
      else
      {
	ATH_MSG_FATAL("Invalid Veto plate number: " << plate);
      }
    }
    else
    {
      ATH_MSG_FATAL("Invalid Veto station number: " << station);
      return StatusCode::FAILURE;
    }
  }



  SG::ReadHandle<Tracker::FaserSCT_ClusterContainer> clusterContainer { m_clusterContainer, ctx };
  ATH_CHECK(clusterContainer.isValid());
  SG::ReadHandle<TrackerSimDataCollection> simDataCollection {m_simDataCollection, ctx};
  m_ntruthmatchedclus_event = 0;
  for(int i=0;i<12;i++){
    m_ntruthmatchedclus_layer.push_back(0);
  }

  std::map<Identifier,std::pair<int, int>> identmap;
  for (auto collection : *clusterContainer)
  {
    Identifier id = collection->identify();
    int station = m_sctHelper->station(id);
    int clusters = (int) collection->size();
    int layerid = (int)(m_sctHelper->layer(id)) + 3*station;
    switch (station)
    {
      case 0:
	m_sta0nclus += clusters;
	break;
      case 1:
	m_sta1nclus += clusters;
	break;
      case 2:
	m_sta2nclus += clusters;
	break;
      case 3:
	m_sta3nclus += clusters;
	break;
      default:
	ATH_MSG_FATAL("Unknown tracker station number " << station);
	break;
    }
    //check the truthmatched-cluster
    if (simDataCollection.isValid()){
      for(auto cluster : *collection){
	if (cluster != nullptr)
	{
	  auto idRDO = cluster->identify();

	  if (simDataCollection->count(idRDO) > 0){
	    identmap.emplace(idRDO,std::make_pair(layerid,1));
	    m_ntruthmatchedclus_layer[layerid] +=1;
	    m_ntruthmatchedclus_event +=1;
	  }
	  else
	    identmap.emplace(idRDO,std::make_pair(layerid,0));
	}
      }
    }
  }

  //good station  ihas n cluster larger than 2
  m_ngoodsta = (m_sta0nclus>2) + (m_sta1nclus>2) + (m_sta2nclus>2) + (m_sta3nclus>2);


  SG::ReadHandle<TrackCollection> trackCollection {m_trackCollection, ctx};
  ATH_CHECK(trackCollection.isValid());

  const Trk::TrackParameters* candidateParameters {nullptr};
  const Trk::Track* candidateTrack {nullptr};

  for (const Trk::Track* track : *trackCollection)
  {
    if (track == nullptr) continue;
    std::vector<int> stationMap;
    std::set<std::pair<int, int>> layerMap;

    int firstlayer=99;
    int lastlayer=-1;
    double first_x(999.);
    double first_y(999.);
    double first_z(999.);
    double last_x(999.);
    double last_y(999.);
    double last_z(999.);
    int ntruthmatchedclus[4]={0,0,0,0};
    // Check for hit in the three downstream stations
    for (auto measurement : *(track->measurementsOnTrack()))
    {
      const Tracker::FaserSCT_ClusterOnTrack* cluster = dynamic_cast<const Tracker::FaserSCT_ClusterOnTrack*>(measurement);
      if (cluster != nullptr)
      {
	Identifier id = cluster->identify();
	int station = m_sctHelper->station(id);
	int layer = m_sctHelper->layer(id);
	stationMap.push_back(station);
	layerMap.emplace(station, layer);
	int layerid = 3*station+layer;
	const Tracker::FaserSCT_Cluster* clu = cluster->prepRawData();
	auto match = identmap.find(clu->identify());
	if(match != identmap.end()){
	  ntruthmatchedclus[station]+=1;
	}
	if(layerid<firstlayer){
	  firstlayer=layerid;
	  auto glopos = clu->globalPosition();
	  first_x=glopos.x();
	  first_y=glopos.y();
	  first_z=glopos.z();
	}
	if(layerid>lastlayer){
	  lastlayer=layerid;
	  auto glopos = clu->globalPosition();
	  last_x=glopos.x();
	  last_y=glopos.y();
	  last_z=glopos.z();
	}
      }
    }

    auto countsta = [&](int sta){int n =0;for(auto iter: layerMap){if(iter.first == sta) n = 1;}return n;};
    int nsta = countsta(0)+countsta(1)+countsta(2)+countsta(3);
    int nLayers = std::count_if(layerMap.begin(), layerMap.end(), [](std::pair<int,int> p){return p.first != 0;});
    int nLayers_ift = std::count_if(layerMap.begin(), layerMap.end(), [](std::pair<int,int> p){return p.first == 0;});
    //    if (nLayers < m_minTrackerLayers) continue;
    if(nsta>2)
    m_longTracks++;
    const Trk::TrackParameters* upstreamParameters {nullptr};
    double chi2=9999.;
    double ndof=9999.;
    for (auto params : *(track->trackParameters()))
    {
      if (upstreamParameters == nullptr || params->position().z() < upstreamParameters->position().z()) upstreamParameters = params;
    }
    if (candidateParameters == nullptr || upstreamParameters->momentum().mag() > candidateParameters->momentum().mag())
    {
      candidateParameters = upstreamParameters;
      candidateTrack = track;
      chi2 = track->fitQuality()->chiSquared();
      ndof = track->fitQuality()->numberDoF();
    }
    if(chi2<1000){
    m_nclus0.push_back(std::count_if(stationMap.begin(), stationMap.end(), [](int p){return p==0;}));
    m_nclus1.push_back(std::count_if(stationMap.begin(), stationMap.end(), [](int p){return p==1;}));
    m_nclus2.push_back(std::count_if(stationMap.begin(), stationMap.end(), [](int p){return p==2;}));
    m_nclus3.push_back(std::count_if(stationMap.begin(), stationMap.end(), [](int p){return p==3;}));
    m_nclus.push_back(stationMap.size());
    m_nlayer.push_back(nLayers);
    m_nlayer_ift.push_back(nLayers_ift);
    m_chi2.push_back(chi2);
    m_ndof.push_back(ndof);
    m_nsta.push_back(nsta);
    m_first_x.push_back(first_x);
    m_first_y.push_back(first_y);
    m_first_z.push_back(first_z);
    m_last_x.push_back(last_x);
    m_last_y.push_back(last_y);
    m_last_z.push_back(last_z);
    m_firstlayer.push_back(firstlayer);
    m_lastlayer.push_back(lastlayer);
    m_x.push_back(candidateParameters->position().x());
    m_y.push_back(candidateParameters->position().y());
    m_z.push_back(candidateParameters->position().z());
    m_px.push_back(candidateParameters->momentum().x());
    m_py.push_back(candidateParameters->momentum().y());
    m_pz.push_back(candidateParameters->momentum().z());
    m_p.push_back(sqrt(m_px.back() * m_px.back() + m_py.back() * m_py.back() + m_pz.back() * m_pz.back()));
    m_charge.push_back((int) candidateParameters->charge());
    m_ntruthmatchedclus.push_back(ntruthmatchedclus[0]+ntruthmatchedclus[1]+ntruthmatchedclus[2]+ntruthmatchedclus[3]);
    m_ntruthmatchedclus0.push_back(ntruthmatchedclus[0]);
    m_ntruthmatchedclus1.push_back(ntruthmatchedclus[1]);
    m_ntruthmatchedclus2.push_back(ntruthmatchedclus[2]);
    m_ntruthmatchedclus3.push_back(ntruthmatchedclus[3]);
    ifsave=true;
    }
  }

  //  SG::ReadHandle<TrackerSimDataCollection> simDataCollection {m_simDataCollection, ctx};
  //   ATH_MSG_INFO("SimData valid? " << simDataCollection.isValid());
  if (candidateTrack != nullptr && simDataCollection.isValid())
  {
    std::map<int, float> truthMap;
    for (auto measurement : *(candidateTrack->measurementsOnTrack()))
    {
      const Tracker::FaserSCT_ClusterOnTrack* cluster = dynamic_cast<const Tracker::FaserSCT_ClusterOnTrack*>(measurement);
      if (cluster != nullptr)
      {
	// ATH_MSG_INFO("ClusterOnTrack is OK");
	cluster->dump(msg());

	// Hack to work around issue with cluster navigation

	auto idRDO = cluster->identify();

	if (simDataCollection->count(idRDO) > 0)
	{
	  // ATH_MSG_INFO("rdo entry found");
	  const auto& simdata = simDataCollection->find(idRDO)->second;
	  const auto& deposits = simdata.getdeposits();
	  //loop through deposits and record contributions
	  HepMcParticleLink primary{};
	  for( const auto& depositPair : deposits)
	  {
	    // ATH_MSG_INFO("Deposit found");
	    float eDep = depositPair.second;
	    int barcode = depositPair.first->barcode();
	    if (truthMap.count(barcode) > 0)
	    {
	      truthMap[barcode] += eDep;
	    }
	    else
	    {
	      truthMap[barcode] = eDep;
	    }
	  }
	}


      }
    }
    std::vector<std::pair<int, float>> truth(truthMap.begin(), truthMap.end());
    std::sort(truth.begin(), truth.end(), [](auto v1, auto v2) { return v1.second > v2.second; });
    if (truth.size()>0) ATH_MSG_ALWAYS("Selected track truth info:");
    for (auto v : truth)
    {
      auto truthParticle = (*(std::find_if(truthParticleContainer->cbegin(), truthParticleContainer->cend(), [v](const xAOD::TruthParticle* p){ return p->barcode() == v.first; })));
      if (m_truthPdg == 0) m_truthPdg = truthParticle->pdgId();
      if (m_truthBarcode == 0) m_truthBarcode = v.first;
      ATH_MSG_ALWAYS("truth info: barcode = " << v.first << " ( " << truthParticle->p4().P()/1000 << " GeV/c, Id code = " << truthParticle->pdgId() << ") -> deposited energy: " << v.second/1000);
    }
  }

  //  if (candidateParameters != nullptr)
  //  {
  // }

  // Here we apply the signal selection
  // Very simple/unrealistic to start
  //  if (m_vetoUpstream == 0 || m_vetoDownstream == 0 ||
  //        candidateParameters == nullptr)
  //      return StatusCode::SUCCESS;

//  if(ifsave)
  m_tree->Fill();

  return StatusCode::SUCCESS;
}


StatusCode NtupleMakerAlg::finalize() 
{
  return StatusCode::SUCCESS;
}

bool NtupleMakerAlg::waveformHitOK(const xAOD::WaveformHit* hit) const
{
  if (hit->status_bit(xAOD::WaveformStatus::THRESHOLD_FAILED) || hit->status_bit(xAOD::WaveformStatus::SECONDARY)) return false;
  return true;
}

void
NtupleMakerAlg::clearTree() const
{
  m_run_number = 0;
  m_event_number = 0;
  m_vetoNu0 = 0;
  m_vetoNu1 = 0;
  m_veto00 = 0;
  m_veto01 = 0;
  m_veto10 = 0;
  m_veto11 = 0;
  m_vetoUpstream = 0;
  m_vetoDownstream = 0;
  m_sta0nclus = 0;
  m_sta1nclus = 0;
  m_sta2nclus = 0;
  m_sta3nclus = 0;
  m_chi2.clear();
  m_nclus0.clear();
  m_nclus1.clear();
  m_nclus2.clear();
  m_nclus3.clear();
  m_nclus.clear();
  m_nlayer.clear();
  m_nlayer_ift.clear();
  m_nsta.clear();
  m_firstlayer.clear();
  m_lastlayer.clear();
  m_first_x.clear();
  m_first_y.clear();
  m_first_z.clear();
  m_last_x.clear();
  m_last_y.clear();
  m_last_z.clear();
  m_ndof.clear();
  m_px.clear();
  m_py.clear();
  m_pz.clear();
  m_p.clear();
  m_charge.clear();
  m_x.clear();
  m_y.clear();
  m_z.clear();
  m_ntruthmatchedclus.clear();
  m_ntruthmatchedclus0.clear();
  m_ntruthmatchedclus1.clear();
  m_ntruthmatchedclus2.clear();
  m_ntruthmatchedclus3.clear();
  m_longTracks = 0;
  m_truthLeptonMomentum = 0;
  m_truthBarcode = 0;
  m_truthPdg = 0;
  m_ntruthmatchedclus_layer.clear();
  m_ntruthmatchedclus_event = 0;
}
