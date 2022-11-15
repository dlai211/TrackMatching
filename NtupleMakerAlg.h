#ifndef SKIMXAOD_NTUPLEMAKERALG_H
#define SKIMXAOD_NTUPLEMAKERALG_H

#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "AthenaBaseComps/AthHistogramming.h"
#include "TrkTrack/TrackCollection.h"
#include "xAODFaserWaveform/WaveformHitContainer.h"
#include "xAODFaserWaveform/WaveformHit.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "TrackerPrepRawData/FaserSCT_ClusterContainer.h"
#include "TrackerSimData/TrackerSimDataCollection.h"


class TTree;
class FaserSCT_ID;
class VetoNuID;
class VetoID;
namespace  TrackerDD
{
    class SCT_DetectorManager;
}

class NtupleMakerAlg : public AthReentrantAlgorithm, AthHistogramming {
public:
  NtupleMakerAlg(const std::string &name, ISvcLocator *pSvcLocator);
  virtual ~NtupleMakerAlg() = default;
  virtual StatusCode initialize() override;
  virtual StatusCode execute(const EventContext &ctx) const override;
  virtual StatusCode finalize() override;
  const ServiceHandle <ITHistSvc> &histSvc() const;

private:

  bool waveformHitOK(const xAOD::WaveformHit* hit) const;
  void clearTree() const;

  ServiceHandle <ITHistSvc> m_histSvc;

  SG::ReadHandleKey<xAOD::TruthEventContainer> m_truthEventContainer { this, "EventContainer", "TruthEvents", "Truth event container name." };
  SG::ReadHandleKey<xAOD::TruthParticleContainer> m_truthParticleContainer { this, "ParticleContainer", "TruthParticles", "Truth particle container name." };
  SG::ReadHandleKey<TrackerSimDataCollection> m_simDataCollection {this, "TrackerSimDataCollection", "SCT_SDO_Map"};

  SG::ReadHandleKey<TrackCollection> m_trackCollection { this, "TrackCollection", "CKFTrackCollection", "Input track collection name" };
  SG::ReadHandleKey<xAOD::WaveformHitContainer> m_vetoNuContainer { this, "VetoNuContainer", "VetoNuWaveformHits", "VetoNu hit container name" };
  SG::ReadHandleKey<xAOD::WaveformHitContainer> m_vetoContainer { this, "VetoContainer", "VetoWaveformHits", "Veto hit container name" };
  SG::ReadHandleKey<Tracker::FaserSCT_ClusterContainer> m_clusterContainer { this, "ClusterContainer", "SCT_ClusterContainer", "Tracker cluster container name" };

  const TrackerDD::SCT_DetectorManager* m_detMgr {nullptr};

  const FaserSCT_ID* m_sctHelper;
  const VetoNuID*    m_vetoNuHelper;
  const VetoID*      m_vetoHelper;

  IntegerProperty m_minTrackerLayers     { this, "MinTrackerLayers", 7, "Minimum number of layers with hits on track" };

  //DoubleProperty  m_genieLuminosity   { this, "GenieLuminosity", 150.0, "Genie luminosity in inverse fb." };

//   BooleanProperty m_enforceBlinding   { this, "EnforceBlinding", true, "Ignore data events with no VetoNu signals." };

  mutable TTree* m_tree;
  mutable unsigned int m_run_number;
  mutable unsigned int m_event_number;
  mutable double m_vetoNu0;
  mutable double m_vetoNu1;
  mutable double m_veto00;
  mutable double m_veto01;
  mutable double m_vetoUpstream;
  mutable double m_veto10;
  mutable double m_veto11;
  mutable double m_vetoDownstream;
  mutable int m_sta0nclus;
  mutable int m_sta1nclus;
  mutable int m_sta2nclus;
  mutable int m_sta3nclus;
  mutable int m_ngoodsta;

  mutable std::vector<double> m_x;
  mutable std::vector<double> m_y;
  mutable std::vector<double> m_z;
  mutable std::vector<double> m_px;
  mutable std::vector<double> m_py;
  mutable std::vector<double> m_pz;
  mutable std::vector<double> m_p;
  mutable std::vector<int>    m_charge;
  mutable std::vector<double> m_chi2;
  mutable std::vector<int>    m_ndof;
  mutable std::vector<int>    m_nclus0;
  mutable std::vector<int>    m_nclus1;
  mutable std::vector<int>    m_nclus2;
  mutable std::vector<int>    m_nclus3;
  mutable std::vector<int>    m_nclus;
  mutable std::vector<int>    m_nlayer;
  mutable std::vector<int>    m_nlayer_ift;
  mutable std::vector<int>    m_nsta;
  mutable std::vector<int>    m_firstlayer;
  mutable std::vector<int>    m_lastlayer;
  mutable std::vector<int>    m_first_x;
  mutable std::vector<int>    m_first_y;
  mutable std::vector<int>    m_first_z;
  mutable std::vector<int>    m_last_x;
  mutable std::vector<int>    m_last_y;
  mutable std::vector<int>    m_last_z;
  mutable std::vector<int>    m_ntruthmatchedclus;
  mutable std::vector<int>    m_ntruthmatchedclus0;
  mutable std::vector<int>    m_ntruthmatchedclus1;
  mutable std::vector<int>    m_ntruthmatchedclus2;
  mutable std::vector<int>    m_ntruthmatchedclus3;
  mutable int    m_longTracks;
  mutable double m_truthLeptonMomentum;
  mutable int    m_truthBarcode;
  mutable int    m_truthPdg;
  mutable std::vector<int>    m_ntruthmatchedclus_layer;
  mutable int    m_ntruthmatchedclus_event;

};

inline const ServiceHandle <ITHistSvc> &NtupleMakerAlg::histSvc() const {
  return m_histSvc;
}

#endif  // SKIMXAOD_NTUPLEMAKERALG_H
