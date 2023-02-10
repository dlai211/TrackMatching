#include "NtupleDumperAlg.h"
#include "TrkTrack/Track.h"
#include "TrackerRIO_OnTrack/FaserSCT_ClusterOnTrack.h"
#include "TrackerIdentifier/FaserSCT_ID.h"
#include "ScintIdentifier/VetoNuID.h"
#include "ScintIdentifier/VetoID.h"
#include "ScintIdentifier/TriggerID.h"
#include "ScintIdentifier/PreshowerID.h"
#include "FaserCaloIdentifier/EcalID.h"
#include "TrackerPrepRawData/FaserSCT_Cluster.h"
#include "TrackerSpacePoint/FaserSCT_SpacePoint.h"
#include "Identifier/Identifier.h"
#include "TrackerReadoutGeometry/SCT_DetectorManager.h"
#include "TrackerReadoutGeometry/SiDetectorElement.h"
#include "TrackerPrepRawData/FaserSCT_Cluster.h"
#include "xAODTruth/TruthParticle.h"
#include <cmath>
#include <TH1F.h>
#include <numeric>

constexpr float NaN = std::numeric_limits<double>::quiet_NaN();


NtupleDumperAlg::NtupleDumperAlg(const std::string &name, 
                                    ISvcLocator *pSvcLocator)
    : AthReentrantAlgorithm(name, pSvcLocator), 
      AthHistogramming(name),
      m_histSvc("THistSvc/THistSvc", name) {}


void NtupleDumperAlg::addBranch(const std::string &name,
				float* var) {
  m_tree->Branch(name.c_str(),var,(name+"/F").c_str());
}
void NtupleDumperAlg::addBranch(const std::string &name,
				unsigned int* var) {
  m_tree->Branch(name.c_str(),var,(name+"/I").c_str());
}

void NtupleDumperAlg::addWaveBranches(const std::string &name,
				      int nchannels,
				      int first) {
  for(int ch=0;ch<nchannels;ch++) {
    std::string base=name+std::to_string(ch)+"_";
    addBranch(base+"time",&m_wave_localtime[first]);
    addBranch(base+"peak",&m_wave_peak[first]);
    addBranch(base+"width",&m_wave_width[first]);
    addBranch(base+"charge",&m_wave_charge[first]);
    addBranch(base+"raw_peak",&m_wave_raw_peak[first]);
    addBranch(base+"raw_charge",&m_wave_raw_charge[first]);
    addBranch(base+"baseline",&m_wave_baseline_mean[first]);
    addBranch(base+"baseline_rms",&m_wave_baseline_rms[first]);
    addBranch(base+"status",&m_wave_status[first]);
    first++;
  }
}

void NtupleDumperAlg::FillWaveBranches(const xAOD::WaveformHitContainer &wave) const {
  for (auto hit : wave) {
    if ((hit->hit_status()&2)==0) { // dont store secoondary hits as they can overwrite the primary hit
      int ch=hit->channel();
      m_wave_localtime[ch]=hit->localtime()+m_clock_phase;
      m_wave_peak[ch]=hit->peak();
      m_wave_width[ch]=hit->width();
      m_wave_charge[ch]=hit->integral()/50;

      m_wave_raw_peak[ch]=hit->raw_peak();
      m_wave_raw_charge[ch]=hit->raw_integral()/50;
      m_wave_baseline_mean[ch]=hit->baseline_mean();
      m_wave_baseline_rms[ch]=hit->baseline_rms();
      m_wave_status[ch]=hit->hit_status();  
    }
  }
}

void NtupleDumperAlg::addCalibratedBranches(const std::string &name,
                      int nchannels,
                      int first) {
  for(int ch=0;ch<nchannels;ch++) {
    std::string base=name+std::to_string(ch)+"_";
    addBranch(base+"nMIP",&m_calibrated_nMIP[first]);
    addBranch(base+"E_dep",&m_calibrated_E_dep[first]);
    addBranch(base+"E_EM",&m_calibrated_E_EM[first]);
    first++;
  }
}

StatusCode NtupleDumperAlg::initialize() 
{
  ATH_CHECK(m_truthEventContainer.initialize());
  ATH_CHECK(m_truthParticleContainer.initialize());
  ATH_CHECK(m_lhcData.initialize());
  ATH_CHECK(m_trackCollection.initialize());

  if (!m_isMC) { // disable for MC for now
    ATH_CHECK(m_trackCollectionWithoutIFT.initialize()); 
  }

  ATH_CHECK(m_trackSegmentCollection.initialize());
  ATH_CHECK(m_vetoNuContainer.initialize());
  ATH_CHECK(m_vetoContainer.initialize());
  ATH_CHECK(m_triggerContainer.initialize());
  ATH_CHECK(m_preshowerContainer.initialize());
  ATH_CHECK(m_ecalContainer.initialize());
  ATH_CHECK(m_clusterContainer.initialize());
  ATH_CHECK(m_simDataCollection.initialize());
  ATH_CHECK(m_FaserTriggerData.initialize());
  ATH_CHECK(m_ClockWaveformContainer.initialize());
  ATH_CHECK(m_siHitCollectionKey.initialize());

  if (!m_isMC) { // disable for MC for now
    ATH_CHECK(m_preshowerCalibratedContainer.initialize());
    ATH_CHECK(m_ecalCalibratedContainer.initialize());
  }

  ATH_CHECK(detStore()->retrieve(m_sctHelper,       "FaserSCT_ID"));
  ATH_CHECK(detStore()->retrieve(m_vetoNuHelper,    "VetoNuID"));
  ATH_CHECK(detStore()->retrieve(m_vetoHelper,      "VetoID"));
  ATH_CHECK(detStore()->retrieve(m_triggerHelper,   "TriggerID"));
  ATH_CHECK(detStore()->retrieve(m_preshowerHelper, "PreshowerID"));
  ATH_CHECK(detStore()->retrieve(m_ecalHelper,      "EcalID"));

  ATH_CHECK(detStore()->retrieve(m_detMgr, "SCT"));
  ATH_CHECK(m_extrapolationTool.retrieve());
  ATH_CHECK(m_trackingGeometryTool.retrieve());
  ATH_CHECK(m_trackTruthMatchingTool.retrieve());
  ATH_CHECK(m_fiducialParticleTool.retrieve());

  ATH_CHECK(m_spacePointContainerKey.initialize());

  if (m_useFlukaWeights)
  {
    m_baseEventCrossSection = (m_flukaCrossSection * kfemtoBarnsPerMilliBarn)/m_flukaCollisions;
  }
  else if (m_useGenieWeights)
  {
    m_baseEventCrossSection = 1.0/m_genieLuminosity;
  }
  else
  {
    m_baseEventCrossSection = 1.0;
  }

  m_tree = new TTree("nt", "NtupleDumper tree");
  m_tree->Branch("run", &m_run_number, "run/I");
  m_tree->Branch("eventID", &m_event_number, "eventID/I");
  m_tree->Branch("eventTime", &m_event_time, "eventTime/I");
  m_tree->Branch("BCID", &m_bcid, "BCID/I");

  m_tree->Branch("fillNumber", &m_fillNumber, "fillNumber/I");
  m_tree->Branch("betaStar", &m_betaStar, "betaStar/F");
  m_tree->Branch("crossingAngle", &m_crossingAngle, "crossingAngle/F");
  m_tree->Branch("distanceToCollidingBCID", &m_distanceToCollidingBCID, "distanceToCollidingBCID/I");
  m_tree->Branch("distanceToUnpairedB1", &m_distanceToUnpairedB1, "distanceToUnpairedB1/I");
  m_tree->Branch("distanceToUnpairedB2", &m_distanceToUnpairedB2, "distanceToUnpairedB2/I");
  m_tree->Branch("distanceToInboundB1", &m_distanceToInboundB1, "distanceToInboundB1/I");
  m_tree->Branch("distanceToTrainStart", &m_distanceToTrainStart, "distanceToTrainStart/I");
  m_tree->Branch("distanceToPreviousColliding", &m_distanceToPreviousColliding, "distanceToPreviousColliding/I");

  m_tree->Branch("TBP", &m_tbp, "TBP/I");
  m_tree->Branch("TAP", &m_tap, "TAP/I");
  m_tree->Branch("inputBits", &m_inputBits, "inputBits/I");
  m_tree->Branch("inputBitsNext", &m_inputBitsNext, "inputBitsNext/I");

  addWaveBranches("VetoNu",2,4);
  addWaveBranches("VetoSt1",2,6);
  addWaveBranches("VetoSt2",1,14);
  addWaveBranches("Timing",4,8);
  addWaveBranches("Preshower",2,12);
  addWaveBranches("Calo",4,0);

  addCalibratedBranches("Calo",4,0);
  addBranch("Calo_total_nMIP", &m_calo_total_nMIP);
  addBranch("Calo_total_E_dep", &m_calo_total_E_dep);
  addBranch("Calo_total_E_EM", &m_calo_total_E_EM);

  addCalibratedBranches("Preshower",2,12);
  addBranch("Preshower_total_nMIP", &m_preshower_total_nMIP);
  addBranch("Preshower_total_E_dep", &m_preshower_total_E_dep);

  addBranch("nClusters0",&m_station0Clusters);
  addBranch("nClusters1",&m_station1Clusters);
  addBranch("nClusters2",&m_station2Clusters);
  addBranch("nClusters3",&m_station3Clusters);

  addBranch("SpacePoints",&m_nspacepoints);
  m_tree->Branch("SpacePoint_x", &m_spacepointX);
  m_tree->Branch("SpacePoint_y", &m_spacepointY);
  m_tree->Branch("SpacePoint_z", &m_spacepointZ);

  addBranch("TrackSegments",&m_ntracksegs);
  m_tree->Branch("TrackSegment_Chi2", &m_trackseg_Chi2);
  m_tree->Branch("TrackSegment_nDoF", &m_trackseg_DoF);
  m_tree->Branch("TrackSegment_x", &m_trackseg_x);
  m_tree->Branch("TrackSegment_y", &m_trackseg_y);
  m_tree->Branch("TrackSegment_z", &m_trackseg_z);
  m_tree->Branch("TrackSegment_px", &m_trackseg_px);
  m_tree->Branch("TrackSegment_py", &m_trackseg_py);
  m_tree->Branch("TrackSegment_pz", &m_trackseg_pz);

  m_tree->Branch("longTracks", &m_longTracks, "longTracks/I");
  m_tree->Branch("Track_Chi2", &m_Chi2);
  m_tree->Branch("Track_nDoF", &m_DoF);
  m_tree->Branch("Track_resx", &m_resx);
  m_tree->Branch("Track_x0", &m_xup);
  m_tree->Branch("Track_y0", &m_yup);
  m_tree->Branch("Track_z0", &m_zup);
  m_tree->Branch("Track_px0", &m_pxup);
  m_tree->Branch("Track_py0", &m_pyup);
  m_tree->Branch("Track_pz0", &m_pzup);
  m_tree->Branch("Track_p0", &m_pup);
  m_tree->Branch("Track_x1", &m_xdown);
  m_tree->Branch("Track_y1", &m_ydown);
  m_tree->Branch("Track_z1", &m_zdown);
  m_tree->Branch("Track_px1", &m_pxdown);
  m_tree->Branch("Track_py1", &m_pydown);
  m_tree->Branch("Track_pz1", &m_pzdown);
  m_tree->Branch("Track_p1", &m_pdown);
  m_tree->Branch("Track_charge", &m_charge);
  m_tree->Branch("Track_nLayers", &m_nLayers);

  m_tree->Branch("Track_InStation0",&m_nHit0);
  m_tree->Branch("Track_InStation1",&m_nHit1);
  m_tree->Branch("Track_InStation2",&m_nHit2);
  m_tree->Branch("Track_InStation3",&m_nHit3);

  m_tree->Branch("Track_X_atVetoNu", &m_xVetoNu);
  m_tree->Branch("Track_Y_atVetoNu", &m_yVetoNu);
  m_tree->Branch("Track_ThetaX_atVetoNu", &m_thetaxVetoNu);
  m_tree->Branch("Track_ThetaY_atVetoNu", &m_thetayVetoNu);

  m_tree->Branch("Track_X_atVetoStation1", &m_xVetoStation1);
  m_tree->Branch("Track_Y_atVetoStation1", &m_yVetoStation1);
  m_tree->Branch("Track_ThetaX_atVetoStation1", &m_thetaxVetoStation1);
  m_tree->Branch("Track_ThetaY_atVetoStation1", &m_thetayVetoStation1);

  m_tree->Branch("Track_X_atVetoStation2", &m_xVetoStation2);
  m_tree->Branch("Track_Y_atVetoStation2", &m_yVetoStation2);
  m_tree->Branch("Track_ThetaX_atVetoStation2", &m_thetaxVetoStation2);
  m_tree->Branch("Track_ThetaY_atVetoStation2", &m_thetayVetoStation2);

  m_tree->Branch("Track_X_atTrig", &m_xTrig);
  m_tree->Branch("Track_Y_atTrig", &m_yTrig);
  m_tree->Branch("Track_ThetaX_atTrig", &m_thetaxTrig);
  m_tree->Branch("Track_ThetaY_atTrig", &m_thetayTrig);

  m_tree->Branch("Track_X_atPreshower1", &m_xPreshower1);
  m_tree->Branch("Track_Y_atPreshower1", &m_yPreshower1);
  m_tree->Branch("Track_ThetaX_atPreshower1", &m_thetaxPreshower1);
  m_tree->Branch("Track_ThetaY_atPreshower1", &m_thetayPreshower1);

  m_tree->Branch("Track_X_atPreshower2", &m_xPreshower2);
  m_tree->Branch("Track_Y_atPreshower2", &m_yPreshower2);
  m_tree->Branch("Track_ThetaX_atPreshower2", &m_thetaxPreshower2);
  m_tree->Branch("Track_ThetaY_atPreshower2", &m_thetayPreshower2);

  m_tree->Branch("Track_X_atCalo", &m_xCalo);
  m_tree->Branch("Track_Y_atCalo", &m_yCalo);
  m_tree->Branch("Track_ThetaX_atCalo", &m_thetaxCalo);
  m_tree->Branch("Track_ThetaY_atCalo", &m_thetayCalo);

  m_tree->Branch("t_pdg", &m_t_pdg);
  m_tree->Branch("t_barcode", &m_t_barcode);
  m_tree->Branch("t_truthHitRatio", &m_t_truthHitRatio);
  m_tree->Branch("t_prodVtx_x", &m_t_prodVtx_x);
  m_tree->Branch("t_prodVtx_y", &m_t_prodVtx_y);
  m_tree->Branch("t_prodVtx_z", &m_t_prodVtx_z);
  m_tree->Branch("t_decayVtx_x", &m_t_decayVtx_x);
  m_tree->Branch("t_decayVtx_y", &m_t_decayVtx_y);
  m_tree->Branch("t_decayVtx_z", &m_t_decayVtx_z);
  m_tree->Branch("t_px", &m_t_px);
  m_tree->Branch("t_py", &m_t_py);
  m_tree->Branch("t_pz", &m_t_pz);
  m_tree->Branch("t_theta", &m_t_theta);
  m_tree->Branch("t_phi", &m_t_phi);
  m_tree->Branch("t_p", &m_t_p);
  m_tree->Branch("t_pT", &m_t_pT);
  m_tree->Branch("t_eta", &m_t_eta);
  m_tree->Branch("t_st0_x", &m_t_st_x[0]);
  m_tree->Branch("t_st0_y", &m_t_st_y[0]);
  m_tree->Branch("t_st0_z", &m_t_st_z[0]);
  m_tree->Branch("t_st1_x", &m_t_st_x[1]);
  m_tree->Branch("t_st1_y", &m_t_st_y[1]);
  m_tree->Branch("t_st1_z", &m_t_st_z[1]);
  m_tree->Branch("t_st2_x", &m_t_st_x[2]);
  m_tree->Branch("t_st2_y", &m_t_st_y[2]);
  m_tree->Branch("t_st2_z", &m_t_st_z[2]);
  m_tree->Branch("t_st3_x", &m_t_st_x[3]);
  m_tree->Branch("t_st3_y", &m_t_st_y[3]);
  m_tree->Branch("t_st3_z", &m_t_st_z[3]);
  m_tree->Branch("isFiducial", &m_isFiducial);

  m_tree->Branch("truthParticleBarcode", &m_truthParticleBarcode);
  m_tree->Branch("truthParticleMatchedTracks", &m_truthParticleMatchedTracks);
  m_tree->Branch("truthParticleIsFiducial", &m_truthParticleIsFiducial);

  m_tree->Branch("pTruthLepton", &m_truthLeptonMomentum, "pTruthLepton/D");
  m_tree->Branch("truthBarcode", &m_truthBarcode, "truthBarcode/I");
  m_tree->Branch("truthPdg", &m_truthPdg, "truthPdg/I");
  m_tree->Branch("CrossSection", &m_crossSection, "crossSection/D");

  // for mother + daughter particle truth infomation 

  m_tree->Branch("truthM_P", &m_truthM_P);
  m_tree->Branch("truthM_px", &m_truthM_px);
  m_tree->Branch("truthM_py", &m_truthM_py);
  m_tree->Branch("truthM_pz", &m_truthM_pz);
  m_tree->Branch("truthM_x", &m_truthM_x);
  m_tree->Branch("truthM_y", &m_truthM_y);
  m_tree->Branch("truthM_z", &m_truthM_z);

  m_tree->Branch("truthd0_P", &m_truthd0_P);
  m_tree->Branch("truthd0_px", &m_truthd0_px);
  m_tree->Branch("truthd0_py", &m_truthd0_py);
  m_tree->Branch("truthd0_pz", &m_truthd0_pz);
  m_tree->Branch("truthd0_x", &m_truthd0_x);
  m_tree->Branch("truthd0_y", &m_truthd0_y);
  m_tree->Branch("truthd0_z", &m_truthd0_z);

  m_tree->Branch("truthd1_P", &m_truthd1_P);
  m_tree->Branch("truthd1_px", &m_truthd1_px);
  m_tree->Branch("truthd1_py", &m_truthd1_py);
  m_tree->Branch("truthd1_pz", &m_truthd1_pz);
  m_tree->Branch("truthd1_x", &m_truthd1_x);
  m_tree->Branch("truthd1_y", &m_truthd1_y);
  m_tree->Branch("truthd1_z", &m_truthd1_z);

  ATH_CHECK(histSvc()->regTree("/HIST2/tree", m_tree));

  // Register histograms
  m_HistRandomCharge[0] = new TH1F("hRandomCharge0", "Calo ch0 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[1] = new TH1F("hRandomCharge1", "Calo ch1 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[2] = new TH1F("hRandomCharge2", "Calo ch2 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[3] = new TH1F("hRandomCharge3", "Calo ch3 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[4] = new TH1F("hRandomCharge4", "VetoNu ch4 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[5] = new TH1F("hRandomCharge5", "VetoNu ch5 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[6] = new TH1F("hRandomCharge6", "Veto ch6 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[7] = new TH1F("hRandomCharge7", "Veto ch7 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[8] = new TH1F("hRandomCharge8", "Trig ch8 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[9] = new TH1F("hRandomCharge9", "Trig ch9 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[10] = new TH1F("hRandomCharge10", "Trig ch10 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[11] = new TH1F("hRandomCharge11", "Trig ch11 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[12] = new TH1F("hRandomCharge12", "Preshower ch12 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[13] = new TH1F("hRandomCharge13", "Preshower ch13 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);
  m_HistRandomCharge[14] = new TH1F("hRandomCharge14", "Veto ch14 Charge from Random Events;charge (pC);Events/bin", 100, -1.0, 1.0);

  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge0", m_HistRandomCharge[0]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge1", m_HistRandomCharge[1]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge2", m_HistRandomCharge[2]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge3", m_HistRandomCharge[3]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge4", m_HistRandomCharge[4]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge5", m_HistRandomCharge[5]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge6", m_HistRandomCharge[6]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge7", m_HistRandomCharge[7]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge8", m_HistRandomCharge[8]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge9", m_HistRandomCharge[9]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge10", m_HistRandomCharge[10]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge11", m_HistRandomCharge[11]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge12", m_HistRandomCharge[12]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge13", m_HistRandomCharge[13]));
  ATH_CHECK(histSvc()->regHist("/HIST2/RandomCharge14", m_HistRandomCharge[14]));

  if (m_doBlinding) {
    ATH_MSG_INFO("Blinding will be enforced for real data.");
  } else {
    ATH_MSG_INFO("Blinding will NOT be enforced for real data.");
  }

  return StatusCode::SUCCESS;
}


StatusCode NtupleDumperAlg::execute(const EventContext &ctx) const
{
  clearTree();

  // check if real data or simulation data
  bool isMC = false;
  SG::ReadHandle<xAOD::TruthEventContainer> truthEventContainer { m_truthEventContainer, ctx };
  if (truthEventContainer.isValid() && truthEventContainer->size() > 0)
  {
    isMC = true;
  }

  // if real data, store charge in histograms from random events and only fill ntuple from coincidence events
  if (!isMC) {
    SG::ReadHandle<xAOD::FaserTriggerData> triggerData(m_FaserTriggerData, ctx);
    m_tap=triggerData->tap();
    // unpack trigger word: 1=calo, 2=veotnu|veto1|preshower, 4=TimingLayer, 8=(VetoNu|Veto2)&Preshower, 16=random, 32=LED 
    bool trig_random = false;
    if ( m_tap == 16 ) trig_random = true;

    bool trig_coincidence_preshower_and_vetoes = false;
    if ( (m_tap&8) != 0 ) trig_coincidence_preshower_and_vetoes = true;

    bool trig_coincidence_timing_and_vetoesORpreshower = false;
    if ( ((m_tap&4)!=0) && ((m_tap&2)!=0) ) trig_coincidence_timing_and_vetoesORpreshower = true;

    bool trig_coincidence_timing_and_calo = false;
    if ( ((m_tap&4)!=0) && ((m_tap&1)!=0) ) trig_coincidence_timing_and_calo = true;

    bool trig_coincidence_vetoesORpreshower_and_calo = false;
    if ( ((m_tap&2)!=0) && ((m_tap&1)!=0) ) trig_coincidence_vetoesORpreshower_and_calo = true;

    // for random trigger, store charge of scintillators in histograms
    if (trig_random) {
      // Read in Waveform containers
      SG::ReadHandle<xAOD::WaveformHitContainer> vetoNuContainer { m_vetoNuContainer, ctx };
      ATH_CHECK(vetoNuContainer.isValid());

      SG::ReadHandle<xAOD::WaveformHitContainer> vetoContainer { m_vetoContainer, ctx };
      ATH_CHECK(vetoContainer.isValid());

      SG::ReadHandle<xAOD::WaveformHitContainer> triggerContainer { m_triggerContainer, ctx };
      ATH_CHECK(triggerContainer.isValid());

      SG::ReadHandle<xAOD::WaveformHitContainer> preshowerContainer { m_preshowerContainer, ctx };
      ATH_CHECK(preshowerContainer.isValid());

      SG::ReadHandle<xAOD::WaveformHitContainer> ecalContainer { m_ecalContainer, ctx };
      ATH_CHECK(ecalContainer.isValid());

      // Fill histograms
      if (vetoNuContainer.isValid()) {
        for (auto hit : *vetoNuContainer) {
          int ch=hit->channel();
          m_HistRandomCharge[ch]->Fill(hit->raw_integral()/50.0);
        }
      }
      if (vetoContainer.isValid()) {
        for (auto hit : *vetoContainer) {
          int ch=hit->channel();
          m_HistRandomCharge[ch]->Fill(hit->raw_integral()/50.0);
        }
      }
      if (triggerContainer.isValid()) {
        for (auto hit : *triggerContainer) {
          int ch=hit->channel();
          m_HistRandomCharge[ch]->Fill(hit->raw_integral()/50.0);
        }
      }
      if (preshowerContainer.isValid()) {
        for (auto hit : *preshowerContainer) {
          int ch=hit->channel();
          m_HistRandomCharge[ch]->Fill(hit->raw_integral()/50.0);
        }
      }
      if (ecalContainer.isValid()) {
        for (auto hit : *ecalContainer) {
          int ch=hit->channel();
          m_HistRandomCharge[ch]->Fill(hit->raw_integral()/50.0);
        }
      }
      if (!m_storeAllEvents) return StatusCode::SUCCESS;  // finished with this randomly triggered event

//     if ( !(trig_coincidence_preshower_and_vetoes || trig_coincidence_timing_and_vetoesORpreshower || trig_coincidence_timing_and_calo || trig_coincidence_vetoesORpreshower_and_calo) ) { 
//      // don't process events that fail to activate coincidence triggers
//      return StatusCode::SUCCESS;
    }

    // store trigger data in ntuple variables
    m_tbp=triggerData->tbp();
    m_tap=triggerData->tap();
    m_inputBits=triggerData->inputBits();
    m_inputBitsNext=triggerData->inputBitsNextClk();

    // load in LHC data
    SG::ReadHandle<xAOD::FaserLHCData> lhcData { m_lhcData, ctx };
    ATH_CHECK(lhcData.isValid());
    // don't process events that were not taken during "Stable Beams"
    if ( (!m_storeAllEvents) && !(lhcData->stableBeams()) ) return StatusCode::SUCCESS;
    // store interesting data in ntuple variables
    m_fillNumber = lhcData->fillNumber();
    m_betaStar = lhcData->betaStar();
    m_crossingAngle = lhcData->crossingAngle();
    m_distanceToCollidingBCID = lhcData->distanceToCollidingBCID();
    m_distanceToUnpairedB1 = lhcData->distanceToUnpairedB1();
    m_distanceToUnpairedB2 = lhcData->distanceToUnpairedB2();
    m_distanceToInboundB1 = lhcData->distanceToInboundB1();
    m_distanceToTrainStart = lhcData->distanceToTrainStart();
    m_distanceToPreviousColliding = lhcData->distanceToPreviousColliding();
    // debug print out all LHC data info available
    ATH_MSG_DEBUG("LHC data fillNumber = " << lhcData->fillNumber() );    
    ATH_MSG_DEBUG("LHC data machineMode = " << lhcData->machineMode() );
    ATH_MSG_DEBUG("LHC data beamMode = " << lhcData->beamMode() );
    ATH_MSG_DEBUG("LHC data beamType1 = " << lhcData->beamType1() );
    ATH_MSG_DEBUG("LHC data beamType2 = " << lhcData->beamType2() );
    ATH_MSG_DEBUG("LHC data betaStar = " << lhcData->betaStar() );
    ATH_MSG_DEBUG("LHC data crossingAngle = " << lhcData->crossingAngle() );
    ATH_MSG_DEBUG("LHC data stableBeams = " << lhcData->stableBeams() );
    ATH_MSG_DEBUG("LHC data injectionScheme = " << lhcData->injectionScheme() );
    ATH_MSG_DEBUG("LHC data numBunchBeam1 = " << lhcData->numBunchBeam1() );
    ATH_MSG_DEBUG("LHC data numBunchBeam2 = " << lhcData->numBunchBeam2() );
    ATH_MSG_DEBUG("LHC data numBunchColliding = " << lhcData->numBunchColliding() );
    ATH_MSG_DEBUG("LHC data distanceToCollidingBCID = " << lhcData->distanceToCollidingBCID() );
    ATH_MSG_DEBUG("LHC data distanceToUnpairedB1 = " << lhcData->distanceToUnpairedB1() );
    ATH_MSG_DEBUG("LHC data distanceToUnpairedB1 = " << lhcData->distanceToUnpairedB2() );
    ATH_MSG_DEBUG("LHC data distanceToInboundB1 = " << lhcData->distanceToInboundB1() );
    ATH_MSG_DEBUG("LHC data distanceToTrainStart = " << lhcData->distanceToTrainStart() );
    ATH_MSG_DEBUG("LHC data distanceToPreviousColliding = " << lhcData->distanceToPreviousColliding() );

    // correct waveform time with clock phase
    SG::ReadHandle<xAOD::WaveformClock> clockHandle(m_ClockWaveformContainer, ctx);
    ATH_CHECK(clockHandle.isValid());

    if (clockHandle->phase() < -2.0) { // wrap around clock pahse so -pi goes to pi
      m_clock_phase = ((clockHandle->phase() + 2*3.14159) / 3.14159) * 12.5;
    } else {
      m_clock_phase = (clockHandle->phase() / 3.14159) * 12.5;
    }

  } // done with processing only on real data

  m_run_number = ctx.eventID().run_number();
  m_event_number = ctx.eventID().event_number();
  m_event_time = ctx.eventID().time_stamp();
  m_bcid = ctx.eventID().bunch_crossing_id();

  if (isMC) { // if simulation find MC cross section and primary lepton
    // Work out effective cross section for MC
    if (m_useFlukaWeights)
    {
        double flukaWeight = truthEventContainer->at(0)->weights()[0];
        ATH_MSG_ALWAYS("Found fluka weight = " << flukaWeight);
        m_crossSection = m_baseEventCrossSection * flukaWeight;
    }
    else if (m_useGenieWeights)
    {
        m_crossSection = m_baseEventCrossSection;
    }
    else
    {
      //ATH_MSG_WARNING("Monte carlo event with no weighting scheme specified.  Setting crossSection (weight) to " << m_baseEventCrossSection << " fb.");
        m_crossSection = m_baseEventCrossSection;
    }
    
    // Find the M d0 and d1 truth information 
    SG::ReadHandle<xAOD::TruthParticleContainer> truthParticleContainer { m_truthParticleContainer, ctx };
    if (truthParticleContainer.isValid() && truthParticleContainer->size() > 0)
    {
      for (auto particle : *truthParticleContainer)
      {
        if ( particle->barcode() <= 3) 
        {

	if ( particle->barcode() == 1) // mother particle (A')
	  {

	  m_truthM_P.push_back(particle->p4().P());
	  m_truthM_px.push_back(particle->p4().X());
	  m_truthM_py.push_back(particle->p4().Y());
	  m_truthM_pz.push_back(particle->p4().Z());

	  if ( particle->hasProdVtx()) {
	    m_truthM_x.push_back(particle->prodVtx()->x());
	    m_truthM_y.push_back(particle->prodVtx()->y());
	    m_truthM_z.push_back(particle->prodVtx()->z());
	  } else {
	    m_truthM_x.push_back(NaN);
	    m_truthM_y.push_back(NaN);
	    m_truthM_z.push_back(NaN);
	  }

	  }
	if ( particle->pdgId() == 11) // daughter particle (positron)
	  {
	    m_truthd0_P.push_back(particle->p4().P());
	    m_truthd0_px.push_back(particle->p4().X());
	    m_truthd0_py.push_back(particle->p4().Y());
	    m_truthd0_pz.push_back(particle->p4().Z());

	  if ( particle->hasProdVtx()) {
	    m_truthd0_x.push_back(particle->prodVtx()->x());
	    m_truthd0_y.push_back(particle->prodVtx()->y());
	    m_truthd0_z.push_back(particle->prodVtx()->z());
	  } else {
	    m_truthd0_x.push_back(NaN);
	    m_truthd0_y.push_back(NaN);
	    m_truthd0_z.push_back(NaN);
	  }
	  }
	if ( particle->pdgId() == -11) // daughter particle (electron)
	  {
	    m_truthd1_P.push_back(particle->p4().P());
	    m_truthd1_px.push_back(particle->p4().X());
	    m_truthd1_py.push_back(particle->p4().Y());
	    m_truthd1_pz.push_back(particle->p4().Z());

	  if ( particle->hasProdVtx()) {
	    m_truthd1_x.push_back(particle->prodVtx()->x());
	    m_truthd1_y.push_back(particle->prodVtx()->y());
	    m_truthd1_z.push_back(particle->prodVtx()->z());
	  } else {
	    m_truthd1_x.push_back(NaN);
	    m_truthd1_y.push_back(NaN);
	    m_truthd1_z.push_back(NaN);
	  }
	  }
	}
      }
    }
  }

  if (!m_isMC) { // disable for MC for now
      // load in calibrated calo container
      SG::ReadHandle<xAOD::CalorimeterHitContainer> ecalCalibratedContainer { m_ecalCalibratedContainer, ctx };
      ATH_CHECK(ecalCalibratedContainer.isValid());
      for (auto hit : *ecalCalibratedContainer) {
        int ch=hit->channel();
        m_calibrated_nMIP[ch] = hit->Nmip();
        m_calibrated_E_dep[ch] = hit->E_dep(); 
        m_calibrated_E_EM[ch] = hit->E_EM();

        m_calo_total_nMIP += hit->Nmip();
        m_calo_total_E_dep += hit->E_dep();
        m_calo_total_E_EM += hit->E_EM();

        ATH_MSG_DEBUG("Calibrated calo: ch is " << ch << ", edep is " << hit->E_dep() << ", E_EM is " << hit->E_EM() << ", Nmip is " << hit->Nmip() << ", fit_to_raw_ratio is " << hit->fit_to_raw_ratio());

        //// the following is an example of how to access the linked waveform data from the calibrated data
        //auto measurements = &(hit->WaveformLinks())[0];
        //auto link_collections = measurements->getDataPtr();
        //auto link_collection = link_collections[0];
        //auto link_index = measurements->index();
        //auto link = link_collection[link_index];
        //if (link_collection != nullptr) {
        //  ATH_MSG_INFO("DEION TEST: wavelink status is " << link->hit_status() );
        //  ATH_MSG_INFO("DEION TEST: wavelink integral is " << link->integral() );
        //}
      }

      // load in calibrated preshower container
      SG::ReadHandle<xAOD::CalorimeterHitContainer> preshowerCalibratedContainer { m_preshowerCalibratedContainer, ctx };
      ATH_CHECK(preshowerCalibratedContainer.isValid());
      for (auto hit : *preshowerCalibratedContainer) {
        int ch=hit->channel();
        m_calibrated_nMIP[ch] = hit->Nmip();
        m_calibrated_E_dep[ch] = hit->E_dep();

        m_preshower_total_nMIP += hit->Nmip();
        m_preshower_total_E_dep += hit->E_dep();

        ATH_MSG_DEBUG("Calibrated preshower: ch is " << ch << ", edep is " << hit->E_dep() << ", E_EM is " << hit->E_EM() << ", Nmip is " << hit->Nmip() << ", fit_to_raw_ratio is " << hit->fit_to_raw_ratio());
      }
  }

  // process all waveeform data from all scintillator and calorimeter channels
  SG::ReadHandle<xAOD::WaveformHitContainer> vetoNuContainer { m_vetoNuContainer, ctx };
  ATH_CHECK(vetoNuContainer.isValid());

  SG::ReadHandle<xAOD::WaveformHitContainer> vetoContainer { m_vetoContainer, ctx };
  ATH_CHECK(vetoContainer.isValid());

  SG::ReadHandle<xAOD::WaveformHitContainer> triggerContainer { m_triggerContainer, ctx };
  ATH_CHECK(triggerContainer.isValid());

  SG::ReadHandle<xAOD::WaveformHitContainer> preshowerContainer { m_preshowerContainer, ctx };
  ATH_CHECK(preshowerContainer.isValid());

  SG::ReadHandle<xAOD::WaveformHitContainer> ecalContainer { m_ecalContainer, ctx };
  ATH_CHECK(ecalContainer.isValid());

  FillWaveBranches(*vetoNuContainer);
  FillWaveBranches(*vetoContainer);
  FillWaveBranches(*triggerContainer);
  FillWaveBranches(*preshowerContainer);
  FillWaveBranches(*ecalContainer);

  // enforce blinding such that events that miss the veto layers and have a large calo signal are skipped and not in the output root file
  if ( (!isMC) && m_doBlinding) {
    if ( m_calo_total_E_EM/1000.0 > 10.0 ) { // energy is in MeV so divide by 1000 to compare to 10 GeV
      if (m_wave_status[4] == 1 and m_wave_status[5] == 1 and m_wave_status[6] == 1 and m_wave_status[7] == 1 and m_wave_status[14] == 1) {  // hit status == 1 means it is below threshold. channles 4 and 5 are vetoNu, channels 6,7, and 14 are veto
        return StatusCode::SUCCESS;
      }
    }
  }

  // get geometry context
  FaserActsGeometryContext faserGeometryContext = m_trackingGeometryTool->getNominalGeometryContext();
  auto gctx = faserGeometryContext.context();

  // loop over clusters and store how many clusters are in each tracking station
  SG::ReadHandle<Tracker::FaserSCT_ClusterContainer> clusterContainer { m_clusterContainer, ctx };
  ATH_CHECK(clusterContainer.isValid());

  for (auto collection : *clusterContainer)
  {
    Identifier id = collection->identify();
    int station = m_sctHelper->station(id);
    int clusters = (int) collection->size();
    switch (station)
    {
      case 0:
        m_station0Clusters += clusters;
        // following lines commented out depict how to access cluster position
        //for (auto cluster : *collection) {
        //  if (cluster == nullptr) continue;
        //  auto pos = cluster->globalPosition();
        //  m_station0ClusterX.push_back(pos.x());
        //}
        break;
      case 1:
        m_station1Clusters += clusters;
        break;
      case 2:
        m_station2Clusters += clusters;
        break;
      case 3:
        m_station3Clusters += clusters;
        break;
      default:
        ATH_MSG_FATAL("Unknown tracker station number " << station);
        break;
    }
  }

  // loop over spacepoints and store each space point position
  SG::ReadHandle<FaserSCT_SpacePointContainer> spacePointContainer {m_spacePointContainerKey, ctx};
  ATH_CHECK(spacePointContainer.isValid());
  for (const FaserSCT_SpacePointCollection* spacePointCollection : *spacePointContainer) {
    m_nspacepoints += spacePointCollection->size();
    for (const Tracker::FaserSCT_SpacePoint *spacePoint: *spacePointCollection) {
      auto pos = spacePoint->globalPosition();
      m_spacepointX.push_back(pos.x());
      m_spacepointY.push_back(pos.y());
      m_spacepointZ.push_back(pos.z());
    }
  }

  // loop over track segments and store position, momentum, chi2, and nDOF for each segment
  SG::ReadHandle<TrackCollection> trackSegmentCollection {m_trackSegmentCollection, ctx};
  ATH_CHECK(trackSegmentCollection.isValid());
  for (const Trk::Track* trackSeg : *trackSegmentCollection) {
    if (trackSeg == nullptr) continue;
    m_ntracksegs += 1;
    m_trackseg_Chi2.push_back(trackSeg->fitQuality()->chiSquared());
    m_trackseg_DoF.push_back(trackSeg->fitQuality()->numberDoF());
    auto SegParameters = trackSeg->trackParameters()->front();
    const Amg::Vector3D SegPosition = SegParameters->position();
    const Amg::Vector3D SegMomentum = SegParameters->momentum();
    m_trackseg_x.push_back(SegPosition.x());
    m_trackseg_y.push_back(SegPosition.y());
    m_trackseg_z.push_back(SegPosition.z());
    m_trackseg_px.push_back(SegMomentum.x());
    m_trackseg_py.push_back(SegMomentum.y());
    m_trackseg_pz.push_back(SegMomentum.z());
  }

  // Write out all truth particle barcodes that have a momentum larger than MinMomentum (default is 50 GeV)
  std::map<int, size_t> truthParticleCount {};
  if (isMC) {
    SG::ReadHandle<xAOD::TruthParticleContainer> truthParticleContainer { m_truthParticleContainer, ctx };
    ATH_CHECK(truthParticleContainer.isValid() && truthParticleContainer->size() > 0);
    for (const xAOD::TruthParticle *tp : *truthParticleContainer) {
      if (tp->p4().P() > m_minMomentum)
        truthParticleCount[tp->barcode()] = 0;
    }
  }
  
    //loop over all the measurements and get the residual
  for (auto tsos : *(track->trackStateOnSurfaces())) {
    const Trk::TrackParameters* tsos_params = tsos->trackParameters();
    const Acts::BoundTrackParameters parameter = ATLASTrackParameterToActs(tsos_params);
    //only use the measurement
    if(tsos_params!=nullptr&&tsos->type(Trk::TrackStateOnSurface::Measurement)){
      const Tracker::FaserSCT_ClusterOnTrack* cluster = dynamic_cast<const Tracker::FaserSCT_ClusterOnTrack*>(tsos->measurementOnTrack());
      if (cluster != nullptr) {
        ATH_MSG_DEBUG("Track Parameter at the surface "<<tsos_params->parameters());
	      
	Amg::Vector3D global_fit(tsos_params->position().x(), tsos_params->position().y() ,tsos_params->position().z());
	auto local_fit=trans2.inverse()*global_fit;
	  //fill the positions
	m_resx.push_back(cluster->localParameters()[Trk::loc1] - local_fit.x());
	
  // loop over all reconstructed tracks and use only the tracks that have hits in all three tracking stations (excludes IFT)
  // store track parameters at most upstream measurement and at most downstream measurement
  // extrapolate track to all scintillator positions and store extrapolated position and angle
  SG::ReadHandle<TrackCollection> trackCollection {m_trackCollection, ctx}; // hack to look at old data that did not have track collection without ift
  if (!m_isMC) { // disable for MC for now
    SG::ReadHandle<TrackCollection> trackCollection {m_trackCollectionWithoutIFT, ctx}; // use track collection that excludes IFT
  }
  ATH_CHECK(trackCollection.isValid());
  for (const Trk::Track* track : *trackCollection)
  {
    if (track == nullptr) continue;

    std::set<std::pair<int, int>> layerMap;
    std::set<int> stationMap;
    // Check for hit in the three downstream stations
    for (auto measurement : *(track->measurementsOnTrack())) {
        const Tracker::FaserSCT_ClusterOnTrack* cluster = dynamic_cast<const Tracker::FaserSCT_ClusterOnTrack*>(measurement);
        if (cluster != nullptr) {
            Identifier id = cluster->identify();
            int station = m_sctHelper->station(id);
            int layer = m_sctHelper->layer(id);
            stationMap.emplace(station);
            layerMap.emplace(station, layer);
        }
    }
    if (stationMap.count(1) == 0 || stationMap.count(2) == 0 || stationMap.count(3) == 0) continue;

    const Trk::TrackParameters* upstreamParameters = track->trackParameters()->front();
    const Trk::TrackParameters* downstreamParameters = track->trackParameters()->back();

    if ((upstreamParameters == nullptr) || (downstreamParameters == nullptr)) continue;

	m_nLayers.push_back(layerMap.size());

    m_Chi2.push_back(track->fitQuality()->chiSquared());
    m_DoF.push_back(track->fitQuality()->numberDoF());

    m_nHit0.push_back(stationMap.count(0));
    m_nHit1.push_back(stationMap.count(1));
    m_nHit2.push_back(stationMap.count(2));
    m_nHit3.push_back(stationMap.count(3));

    m_charge.push_back( (int) upstreamParameters->charge() );

    m_xup.push_back(upstreamParameters->position().x());
    m_yup.push_back(upstreamParameters->position().y());
    m_zup.push_back(upstreamParameters->position().z());
    m_pxup.push_back(upstreamParameters->momentum().x());
    m_pyup.push_back(upstreamParameters->momentum().y());
    m_pzup.push_back(upstreamParameters->momentum().z());
    m_pup.push_back(sqrt( pow(upstreamParameters->momentum().x(),2) + pow(upstreamParameters->momentum().y(),2) + pow(upstreamParameters->momentum().z(),2) ));

    m_xdown.push_back(downstreamParameters->position().x());
    m_ydown.push_back(downstreamParameters->position().y());
    m_zdown.push_back(downstreamParameters->position().z());
    m_pxdown.push_back(downstreamParameters->momentum().x());
    m_pydown.push_back(downstreamParameters->momentum().y());
    m_pzdown.push_back(downstreamParameters->momentum().z());
    m_pdown.push_back(sqrt( pow(downstreamParameters->momentum().x(),2) + pow(downstreamParameters->momentum().y(),2) + pow(downstreamParameters->momentum().z(),2) ));

    if (isMC) { // if simulation, store track truth info as well
      auto [truthParticle, hitCount] = m_trackTruthMatchingTool->getTruthParticle(track);
      if (truthParticle != nullptr) {
        if (truthParticleCount.count(truthParticle->barcode()) > 0)
          truthParticleCount[truthParticle->barcode()] += 1;
        m_t_pdg.push_back(truthParticle->pdgId());
        m_t_barcode.push_back(truthParticle->barcode());
        // the track fit eats up 5 degrees of freedom, thus the number of hits on track is m_DoF + 5
        m_t_truthHitRatio.push_back(hitCount / (m_DoF.back() + 5));
        m_isFiducial.push_back(m_fiducialParticleTool->isFiducial(truthParticle->barcode()));
        auto positions = m_fiducialParticleTool->getTruthPositions(truthParticle->barcode());
        for (int station = 0; station < 4; ++station) {
          m_t_st_x[station].push_back(positions[station].x());
          m_t_st_y[station].push_back(positions[station].y());
          m_t_st_z[station].push_back(positions[station].z());
        }
        if (truthParticle->hasProdVtx()) {
          m_t_prodVtx_x.push_back(truthParticle->prodVtx()->x());
          m_t_prodVtx_y.push_back(truthParticle->prodVtx()->y());
          m_t_prodVtx_z.push_back(truthParticle->prodVtx()->z());
        } else {
          m_t_prodVtx_x.push_back(NaN);
          m_t_prodVtx_y.push_back(NaN);
          m_t_prodVtx_z.push_back(NaN);
        }
        if (truthParticle->hasDecayVtx()) {
          m_t_decayVtx_x.push_back(truthParticle->decayVtx()->x());
          m_t_decayVtx_y.push_back(truthParticle->decayVtx()->y());
          m_t_decayVtx_z.push_back(truthParticle->decayVtx()->z());
        } else {
          m_t_decayVtx_x.push_back(NaN);
          m_t_decayVtx_y.push_back(NaN);
          m_t_decayVtx_z.push_back(NaN);
        }
        m_t_px.push_back(truthParticle->px());
        m_t_py.push_back(truthParticle->py());
        m_t_pz.push_back(truthParticle->pz());
        m_t_theta.push_back(truthParticle->p4().Theta());
        m_t_phi.push_back(truthParticle->p4().Phi());
        m_t_p.push_back(truthParticle->p4().P());
        m_t_pT.push_back(truthParticle->p4().Pt());
        m_t_eta.push_back(truthParticle->p4().Eta());
      } else {
        ATH_MSG_WARNING("Can not find truthParticle.");
        setNaN();
      }
    } else {
      setNaN();
    }

    // fill extrapolation vectors with NaN, will get set to real number if the track extrapolation succeeds
    m_xVetoNu.push_back(NaN);
    m_yVetoNu.push_back(NaN);
    m_thetaxVetoNu.push_back(NaN);
    m_thetayVetoNu.push_back(NaN);
    m_xVetoStation1.push_back(NaN);
    m_yVetoStation1.push_back(NaN);
    m_thetaxVetoStation1.push_back(NaN);
    m_thetayVetoStation1.push_back(NaN);
    m_xVetoStation2.push_back(NaN);
    m_yVetoStation2.push_back(NaN);
    m_thetaxVetoStation2.push_back(NaN);
    m_thetayVetoStation2.push_back(NaN);
    m_xTrig.push_back(NaN);
    m_yTrig.push_back(NaN);
    m_thetaxTrig.push_back(NaN);
    m_thetayTrig.push_back(NaN);
    m_xPreshower1.push_back(NaN);
    m_yPreshower1.push_back(NaN);
    m_thetaxPreshower1.push_back(NaN);
    m_thetayPreshower1.push_back(NaN);
    m_xPreshower2.push_back(NaN);
    m_yPreshower2.push_back(NaN);
    m_thetaxPreshower2.push_back(NaN);
    m_thetayPreshower2.push_back(NaN);
    m_xCalo.push_back(NaN);
    m_yCalo.push_back(NaN);
    m_thetaxCalo.push_back(NaN);
    m_thetayCalo.push_back(NaN);

    // extrapolate track from first station
    if (stationMap.count(1) > 0) { // extrapolation crashes if the track parameters are is too far away to extrapolate
      Amg::Vector3D position = upstreamParameters->position();
      Amg::Vector3D momentum = upstreamParameters->momentum();
      Acts::BoundVector params = Acts::BoundVector::Zero();
      params[Acts::eBoundLoc0] = -position.y();
      params[Acts::eBoundLoc1] = position.x();
      params[Acts::eBoundPhi] = momentum.phi();
      params[Acts::eBoundTheta] = momentum.theta();
      params[Acts::eBoundQOverP] = upstreamParameters->charge() / momentum.mag();
      params[Acts::eBoundTime] = 0;
      auto startSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, position.z()), Acts::Vector3(0, 0, 1));
      Acts::BoundTrackParameters startParameters(std::move(startSurface), params, upstreamParameters->charge());

      auto targetSurface_VetoNu = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, -3112.0), Acts::Vector3(0, 0, 1)); // -3112 mm is z position of VetoNu planes touching
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_VetoNu =m_extrapolationTool->propagate(ctx, startParameters, *targetSurface_VetoNu, Acts::backward);
      if (targetParameters_VetoNu != nullptr) {
        auto targetPosition_VetoNu = targetParameters_VetoNu->position(gctx);
        auto targetMomentum_VetoNu = targetParameters_VetoNu->momentum();
        m_xVetoNu[m_longTracks] = targetPosition_VetoNu.x();
        m_yVetoNu[m_longTracks] = targetPosition_VetoNu.y();
        m_thetaxVetoNu[m_longTracks] = atan(targetMomentum_VetoNu[0]/targetMomentum_VetoNu[2]);
        m_thetayVetoNu[m_longTracks] = atan(targetMomentum_VetoNu[1]/targetMomentum_VetoNu[2]);
      } else {
        ATH_MSG_INFO("vetoNu null targetParameters");
      }

      auto targetSurface_Veto1 = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, -1769.65), Acts::Vector3(0, 0, 1)); // -1769.65 mm is z position of center of operational layer in Veto station 1
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Veto1 =m_extrapolationTool->propagate(ctx, startParameters, *targetSurface_Veto1, Acts::backward);
      if (targetParameters_Veto1 != nullptr) {
        auto targetPosition_Veto1 = targetParameters_Veto1->position(gctx);
        auto targetMomentum_Veto1 = targetParameters_Veto1->momentum();
        m_xVetoStation1[m_longTracks] = targetPosition_Veto1.x();
        m_yVetoStation1[m_longTracks] = targetPosition_Veto1.y();
        m_thetaxVetoStation1[m_longTracks] = atan(targetMomentum_Veto1[0]/targetMomentum_Veto1[2]);
        m_thetayVetoStation1[m_longTracks] = atan(targetMomentum_Veto1[1]/targetMomentum_Veto1[2]);
      } else {
        ATH_MSG_INFO("veto1 null targetParameters");
      }

      auto targetSurface_Veto2 = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, -1609.65), Acts::Vector3(0, 0, 1)); // -1609.65 mm is z position of where planes touch in Veto station 2
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Veto2 =m_extrapolationTool->propagate(ctx, startParameters, *targetSurface_Veto2, Acts::backward);
      if (targetParameters_Veto2 != nullptr) {
        auto targetPosition_Veto2 = targetParameters_Veto2->position(gctx);
        auto targetMomentum_Veto2 = targetParameters_Veto2->momentum();
        m_xVetoStation2[m_longTracks] = targetPosition_Veto2.x();
        m_yVetoStation2[m_longTracks] = targetPosition_Veto2.y();
        m_thetaxVetoStation2[m_longTracks] = atan(targetMomentum_Veto2[0]/targetMomentum_Veto2[2]);
        m_thetayVetoStation2[m_longTracks] = atan(targetMomentum_Veto2[1]/targetMomentum_Veto2[2]);
      } else {
        ATH_MSG_INFO("veto2 null targetParameters");
      }

      auto targetSurface_Trig = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, 0.0), Acts::Vector3(0, 0, 1)); // 0 mm is z position of Trig planes overlapping
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Trig =m_extrapolationTool->propagate(ctx, startParameters, *targetSurface_Trig, Acts::backward); // must extrapolate backsards to trig plane if track starts in station 1
      if (targetParameters_Trig != nullptr) {
        auto targetPosition_Trig = targetParameters_Trig->position(gctx);
        auto targetMomentum_Trig = targetParameters_Trig->momentum();
        m_xTrig[m_longTracks] = targetPosition_Trig.x();
        m_yTrig[m_longTracks] = targetPosition_Trig.y();
        m_thetaxTrig[m_longTracks] = atan(targetMomentum_Trig[0]/targetMomentum_Trig[2]);
        m_thetayTrig[m_longTracks] = atan(targetMomentum_Trig[1]/targetMomentum_Trig[2]);
      } else {
        ATH_MSG_INFO("Trig null targetParameters");
      }

    }

    // extrapolate track from tracking station 3
    if (stationMap.count(3) > 0) { // extrapolation crashes if the track does not end in the Station 3, as it is too far away to extrapolate
      Amg::Vector3D positionDown = downstreamParameters->position();
      Amg::Vector3D momentumDown = downstreamParameters->momentum();
      Acts::BoundVector paramsDown = Acts::BoundVector::Zero();
      paramsDown[Acts::eBoundLoc0] = -positionDown.y();
      paramsDown[Acts::eBoundLoc1] = positionDown.x();
      paramsDown[Acts::eBoundPhi] = momentumDown.phi();
      paramsDown[Acts::eBoundTheta] = momentumDown.theta();
      paramsDown[Acts::eBoundQOverP] = downstreamParameters->charge() / momentumDown.mag();
      paramsDown[Acts::eBoundTime] = 0;
      auto startSurfaceDown = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, positionDown.z()), Acts::Vector3(0, 0, 1));
      Acts::BoundTrackParameters startParametersDown(std::move(startSurfaceDown), paramsDown, downstreamParameters->charge());

      auto targetSurface_Preshower1 = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, 2582.68), Acts::Vector3(0, 0, 1)); // 2582.68  mm is z position of center of upstream preshower layer
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Preshower1 =m_extrapolationTool->propagate(ctx, startParametersDown, *targetSurface_Preshower1, Acts::forward);
      if (targetParameters_Preshower1 != nullptr) {
        auto targetPosition_Preshower1 = targetParameters_Preshower1->position(gctx);
        auto targetMomentum_Preshower1 = targetParameters_Preshower1->momentum();
        m_xPreshower1[m_longTracks] = targetPosition_Preshower1.x();
        m_yPreshower1[m_longTracks] = targetPosition_Preshower1.y();
        m_thetaxPreshower1[m_longTracks] = atan(targetMomentum_Preshower1[0]/targetMomentum_Preshower1[2]);
        m_thetayPreshower1[m_longTracks] = atan(targetMomentum_Preshower1[1]/targetMomentum_Preshower1[2]);
      } else {
        ATH_MSG_INFO("Preshower1 null targetParameters");
      }

      auto targetSurface_Preshower2 = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, 2657.68), Acts::Vector3(0, 0, 1)); // 2657.68  mm is z position of center of downstream preshower layer
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Preshower2 =m_extrapolationTool->propagate(ctx, startParametersDown, *targetSurface_Preshower2, Acts::forward);
      if (targetParameters_Preshower2 != nullptr) {
        auto targetPosition_Preshower2 = targetParameters_Preshower2->position(gctx);
        auto targetMomentum_Preshower2 = targetParameters_Preshower2->momentum();
        m_xPreshower2[m_longTracks] = targetPosition_Preshower2.x();
        m_yPreshower2[m_longTracks] = targetPosition_Preshower2.y();
        m_thetaxPreshower2[m_longTracks] = atan(targetMomentum_Preshower2[0]/targetMomentum_Preshower2[2]);
        m_thetayPreshower2[m_longTracks] =  atan(targetMomentum_Preshower2[1]/targetMomentum_Preshower2[2]);
      } else {
        ATH_MSG_INFO("Preshower2 null targetParameters");
      }

      auto targetSurface_Calo = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, 2760.0), Acts::Vector3(0, 0, 1)); // 2760  mm is estimated z position of calorimeter face
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Calo =m_extrapolationTool->propagate(ctx, startParametersDown, *targetSurface_Calo, Acts::forward);
      if (targetParameters_Calo != nullptr) {
        auto targetPosition_Calo = targetParameters_Calo->position(gctx);
        auto targetMomentum_Calo = targetParameters_Calo->momentum();
        m_xCalo[m_longTracks] = targetPosition_Calo.x();
        m_yCalo[m_longTracks] = targetPosition_Calo.y();
        m_thetaxCalo[m_longTracks] = atan(targetMomentum_Calo[0]/targetMomentum_Calo[2]) ;
        m_thetayCalo[m_longTracks] = atan(targetMomentum_Calo[1]/targetMomentum_Calo[2]) ;
      } else {
        ATH_MSG_INFO("Calo null targetParameters");
      }
    }

    m_longTracks++;
  }
  if ((!m_storeAllEvents) && m_longTracks == 0) return StatusCode::SUCCESS;

  if (isMC) {
    for (auto &tp : truthParticleCount) {
      m_truthParticleBarcode.push_back(tp.first);
      m_truthParticleMatchedTracks.push_back(tp.second);
      m_truthParticleIsFiducial.push_back(m_fiducialParticleTool->isFiducial(tp.first));
    }
  }

  // finished processing event, now fill ntuple tree
  m_tree->Fill();
  m_eventsPassed += 1;
  return StatusCode::SUCCESS;
}


StatusCode NtupleDumperAlg::finalize()
{
  ATH_MSG_INFO("Number of events passed Ntuple selectioon = " << m_eventsPassed);
  return StatusCode::SUCCESS;
}

bool NtupleDumperAlg::waveformHitOK(const xAOD::WaveformHit* hit) const
{
    if (hit->status_bit(xAOD::WaveformStatus::THRESHOLD_FAILED) || hit->status_bit(xAOD::WaveformStatus::SECONDARY)) return false;
    return true;
}

void
NtupleDumperAlg::clearTree() const
{
  // set all float variables to NaN
  // set all int variables to 999999 (can't set int to NaN, because NaN is a double)
  // set all counter variables to zero
  // set all trigger words to zero
  m_run_number = 999999; 
  m_event_number = 999999;
  m_event_time = 999999;
  m_bcid = 999999;

  m_fillNumber = 999999;
  m_betaStar = NaN;
  m_crossingAngle = NaN;
  m_distanceToCollidingBCID = 999999;
  m_distanceToUnpairedB1 = 999999;
  m_distanceToUnpairedB2 = 999999;
  m_distanceToInboundB1 = 999999;
  m_distanceToTrainStart = 999999;
  m_distanceToPreviousColliding = 999999;

  m_tbp=0;
  m_tap=0;
  m_inputBits=0;
  m_inputBitsNext=0;

  for(int ii=0;ii<15;ii++) {
      m_wave_localtime[ii]=NaN;
      m_wave_peak[ii]=NaN;
      m_wave_width[ii]=NaN;
      m_wave_charge[ii]=NaN;

      m_wave_raw_peak[ii]=NaN;
      m_wave_raw_charge[ii]=NaN;
      m_wave_baseline_mean[ii]=NaN;
      m_wave_baseline_rms[ii]=NaN;
      m_wave_status[ii]=1; // default = 1 means below threshold

      m_calibrated_nMIP[ii]=NaN;
      m_calibrated_E_dep[ii]=NaN;
      m_calibrated_E_EM[ii]=NaN;
  }

  m_calo_total_nMIP=0;
  m_calo_total_E_dep=0;
  m_calo_total_E_EM=0;

  m_preshower_total_nMIP=0;
  m_preshower_total_E_dep=0;

  m_clock_phase=0;

  m_station0Clusters = 0;
  m_station1Clusters = 0;
  m_station2Clusters = 0;
  m_station3Clusters = 0;
  m_crossSection = 0;

  m_nspacepoints = 0;
  m_spacepointX.clear();
  m_spacepointY.clear();
  m_spacepointZ.clear();

  m_ntracksegs = 0;
  m_trackseg_Chi2.clear();
  m_trackseg_DoF.clear();
  m_trackseg_x.clear();
  m_trackseg_y.clear();
  m_trackseg_z.clear();
  m_trackseg_px.clear();
  m_trackseg_py.clear();
  m_trackseg_pz.clear();
  
  m_resx.clear();
  m_xup.clear();
  m_yup.clear();
  m_zup.clear();
  m_pxup.clear();
  m_pyup.clear();
  m_pzup.clear();
  m_pup.clear();

  m_xdown.clear();
  m_ydown.clear();
  m_zdown.clear();
  m_pxdown.clear();
  m_pydown.clear();
  m_pzdown.clear();
  m_pdown.clear();

  m_Chi2.clear();
  m_DoF.clear();
  m_charge.clear();
  m_nLayers.clear();
  m_longTracks = 0;

  m_nHit0.clear();
  m_nHit1.clear();
  m_nHit2.clear();
  m_nHit3.clear();

  m_xVetoNu.clear();
  m_yVetoNu.clear();
  m_thetaxVetoNu.clear();
  m_thetayVetoNu.clear();

  m_xVetoStation1.clear();
  m_yVetoStation1.clear();
  m_thetaxVetoStation1.clear();
  m_thetayVetoStation1.clear();

  m_xVetoStation2.clear();
  m_yVetoStation2.clear();
  m_thetaxVetoStation2.clear();
  m_thetayVetoStation2.clear();

  m_xTrig.clear();
  m_yTrig.clear();
  m_thetaxTrig.clear();
  m_thetayTrig.clear();

  m_xPreshower1.clear();
  m_yPreshower1.clear();
  m_thetaxPreshower1.clear();
  m_thetayPreshower1.clear();

  m_xPreshower2.clear();
  m_yPreshower2.clear();
  m_thetaxPreshower2.clear();
  m_thetayPreshower2.clear();

  m_xCalo.clear();
  m_yCalo.clear();
  m_thetaxCalo.clear();
  m_thetayCalo.clear();

  m_t_pdg.clear();
  m_t_barcode.clear();
  m_t_truthHitRatio.clear();
  m_t_prodVtx_x.clear();
  m_t_prodVtx_y.clear();
  m_t_prodVtx_z.clear();
  m_t_decayVtx_x.clear();
  m_t_decayVtx_y.clear();
  m_t_decayVtx_z.clear();
  m_t_px.clear();
  m_t_py.clear();
  m_t_pz.clear();
  m_t_theta.clear();
  m_t_phi.clear();
  m_t_p.clear();
  m_t_pT.clear();
  m_t_eta.clear();
  m_isFiducial.clear();
  for (int station = 0; station < 4; ++station) {
    m_t_st_x[station].clear();
    m_t_st_y[station].clear();
    m_t_st_z[station].clear();
  }
  m_truthParticleBarcode.clear();
  m_truthParticleMatchedTracks.clear();
  m_truthParticleIsFiducial.clear();

  m_truthLeptonMomentum = 0;
  m_truthBarcode = 0;
  m_truthPdg = 0;
}

void NtupleDumperAlg::setNaN() const {
  m_t_pdg.push_back(0);
  m_t_barcode.push_back(-1);
  m_t_truthHitRatio.push_back(NaN);
  m_t_prodVtx_x.push_back(NaN);
  m_t_prodVtx_y.push_back(NaN);
  m_t_prodVtx_z.push_back(NaN);
  m_t_decayVtx_x.push_back(NaN);
  m_t_decayVtx_y.push_back(NaN);
  m_t_decayVtx_z.push_back(NaN);
  m_t_px.push_back(NaN);
  m_t_py.push_back(NaN);
  m_t_pz.push_back(NaN);
  m_t_theta.push_back(NaN);
  m_t_phi.push_back(NaN);
  m_t_p.push_back(NaN);
  m_t_pT.push_back(NaN);
  m_t_eta.push_back(NaN);
  for (int station = 0; station < 4; ++station) {
    m_t_st_x[station].push_back(NaN);
    m_t_st_y[station].push_back(NaN);
    m_t_st_z[station].push_back(NaN);
  }
  m_isFiducial.push_back(false);
}
