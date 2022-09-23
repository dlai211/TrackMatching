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

using namespace std;

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

StatusCode NtupleDumperAlg::initialize() 
{
  ATH_CHECK(m_truthEventContainer.initialize());
  ATH_CHECK(m_truthParticleContainer.initialize());
  ATH_CHECK(m_trackCollection.initialize());
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

  ATH_CHECK(detStore()->retrieve(m_sctHelper,       "FaserSCT_ID"));
  ATH_CHECK(detStore()->retrieve(m_vetoNuHelper,    "VetoNuID"));
  ATH_CHECK(detStore()->retrieve(m_vetoHelper,      "VetoID"));
  ATH_CHECK(detStore()->retrieve(m_triggerHelper,   "TriggerID"));
  ATH_CHECK(detStore()->retrieve(m_preshowerHelper, "PreshowerID"));
  ATH_CHECK(detStore()->retrieve(m_ecalHelper,      "EcalID"));

  ATH_CHECK(detStore()->retrieve(m_detMgr, "SCT"));
  ATH_CHECK(m_extrapolationTool.retrieve());
  ATH_CHECK(m_trackingGeometryTool.retrieve());

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
  addBranch("Calo_total_charge", &m_calo_total);
  addBranch("Calo_total_raw_charge", &m_calo_rawtotal);

  addBranch("Calo0_Edep", &m_Calo0_Edep);
  addBranch("Calo1_Edep", &m_Calo1_Edep);
  addBranch("Calo2_Edep", &m_Calo2_Edep);
  addBranch("Calo3_Edep", &m_Calo3_Edep);
  addBranch("Calo_total_Edep", &m_Calo_Total_Edep);
  addBranch("Preshower12_Edep", &m_Preshower12_Edep);
  addBranch("Preshower13_Edep", &m_Preshower13_Edep);

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
  m_tree->Branch("nClus",&n_clu);
  m_tree->Branch("Track_Chi2", &m_Chi2);
  m_tree->Branch("Track_nDoF", &m_DoF);
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

  m_tree->Branch("pTruthLepton", &m_truthLeptonMomentum, "pTruthLepton/D");
  m_tree->Branch("truthBarcode", &m_truthBarcode, "truthBarcode/I");
  m_tree->Branch("truthPdg", &m_truthPdg, "truthPdg/I");
  m_tree->Branch("CrossSection", &m_crossSection, "crossSection/D");

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

  m_MIP_sim_Edep_calo = 0.0585; // MIP deposits 0.0585 GeV of energy in calo
  m_MIP_sim_Edep_preshower = 0.004894; // MIP deposits 0.004894 GeV of energy in a preshower layer

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
  bool realData = true;
  SG::ReadHandle<xAOD::TruthEventContainer> truthEventContainer { m_truthEventContainer, ctx };
  if (truthEventContainer.isValid() && truthEventContainer->size() > 0)
  {
    realData = false;
  }

  // if real data, store charge in histograms from random events and only fill ntuple from coincidence events
  if (realData) { //no trigger simulation yet
    SG::ReadHandle<xAOD::FaserTriggerData> triggerData(m_FaserTriggerData, ctx);
    m_tap=triggerData->tap();
    if (m_tap==16) { // random trigger, store charge of scintillators in histograms
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

      return StatusCode::SUCCESS; // finished with this event

    } else if ( ((m_tap&8)==0) && (((m_tap&4)==0)||((m_tap&2)==0)) && (((m_tap&4)==0)||((m_tap&1)==0)) && (((m_tap&2)==0)||((m_tap&1)==0)) ) { // don't process events that don't trigger coincidence triggers: 1=calo, 2=veotnu|neto1|preshower, 4=TimingLayer, 8=(VetoNu|Veto2)&Preshower 
      return StatusCode::SUCCESS;
    }
    m_tbp=triggerData->tbp();
    m_tap=triggerData->tap();
    m_inputBits=triggerData->inputBits();
    m_inputBitsNext=triggerData->inputBitsNextClk();
  }

  m_run_number = ctx.eventID().run_number();
  m_event_number = ctx.eventID().event_number();
  m_event_time = ctx.eventID().time_stamp();
  m_bcid = ctx.eventID().bunch_crossing_id();

  if (!realData) { // if simulation find MC cross section and primary lepton
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

    // Find the primary lepton (if any)
    SG::ReadHandle<xAOD::TruthParticleContainer> truthParticleContainer { m_truthParticleContainer, ctx };
    if (truthParticleContainer.isValid() && truthParticleContainer->size() > 0)
    {
      for (auto particle : *truthParticleContainer)
      {
        if ( particle->absPdgId() == 11 || particle->absPdgId() == 13 || particle->absPdgId() == 15 )
        {
          if (particle->status() == 1 && (particle->nParents() == 0 || particle->nParents() == 2) )
          {
            m_truthLeptonMomentum = particle->p4().P();
            break;
          }
        }
      }
    }
  }

  if (realData) { // correct waveform time with clock phase
    SG::ReadHandle<xAOD::WaveformClock> clockHandle(m_ClockWaveformContainer, ctx);
    ATH_CHECK(clockHandle.isValid());

    if (clockHandle->phase() < -2.0) { // wrap around clock pahse so -pi goes to pi
      m_clock_phase = ((clockHandle->phase() + 3.14159) / 3.14159) * 12.5;
    } else {
      m_clock_phase = (clockHandle->phase() / 3.14159) * 12.5;
    }
  }

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
  
  m_calo_total=m_wave_charge[0]+m_wave_charge[1]+m_wave_charge[2]+m_wave_charge[3];
  m_calo_rawtotal=m_wave_raw_charge[0]+m_wave_raw_charge[1]+m_wave_raw_charge[2]+m_wave_raw_charge[3];

  // do calibration of calo channels from pC to GeV deposited
  if (m_CaloConfig == "High_gain") {
    m_Calo0_Edep = (m_wave_charge[0] / 23.709) * m_MIP_sim_Edep_calo;
    m_Calo1_Edep = (m_wave_charge[1] / 24.333) * m_MIP_sim_Edep_calo;
    m_Calo2_Edep = (m_wave_charge[2] / 24.409) * m_MIP_sim_Edep_calo;
    m_Calo3_Edep = (m_wave_charge[3] / 25.555) * m_MIP_sim_Edep_calo;
  } else if (m_CaloConfig == "Low_gain") { // assume low gain calo 
    m_Calo0_Edep = (m_wave_charge[0] / 0.7909) * m_MIP_sim_Edep_calo;
    m_Calo1_Edep = (m_wave_charge[1] / 0.8197) * m_MIP_sim_Edep_calo;
    m_Calo2_Edep = (m_wave_charge[2] / 0.8256) * m_MIP_sim_Edep_calo;
    m_Calo3_Edep = (m_wave_charge[3] / 0.8821) * m_MIP_sim_Edep_calo;
  } else {
   ATH_MSG_WARNING("Run config is neither High_gain nor Low_gain, it is " << m_CaloConfig << ", calo calibration will be zero"); 
  }
  m_Calo_Total_Edep = m_Calo0_Edep + m_Calo1_Edep + m_Calo2_Edep + m_Calo3_Edep;

  // do calibration of preshower channels from pC to GeV deposited
  m_Preshower12_Edep = (m_wave_charge[12] / 5.0) * m_MIP_sim_Edep_preshower; // 5 pC per MIP is rough measurement
  m_Preshower13_Edep = (m_wave_charge[12] / 5.0) * m_MIP_sim_Edep_preshower;

  if (realData && m_doBlinding) { // enforce blinding such that events with large calo signals are skipped and not in the output root file
    if ((m_Calo_Total_Edep/0.155) > 10.0) { // only save events with a shower less than a 10 GeV e- (assume 10 GeV electron deposits 15.5% of their energy in calo)
      return StatusCode::SUCCESS;
    }
  }

  SG::ReadHandle<Tracker::FaserSCT_ClusterContainer> clusterContainer { m_clusterContainer, ctx };
  ATH_CHECK(clusterContainer.isValid());

  FaserActsGeometryContext faserGeometryContext = m_trackingGeometryTool->getNominalGeometryContext();
  auto gctx = faserGeometryContext.context();

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

  SG::ReadHandle<TrackCollection> trackCollection {m_trackCollection, ctx};
  ATH_CHECK(trackCollection.isValid());
  const Trk::TrackParameters* candidateParameters {nullptr};
  const Trk::TrackParameters* candidateDownParameters {nullptr};
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

    int nLayers = std::count_if(layerMap.begin(), layerMap.end(), [](std::pair<int,int> p){return p.first != 0;});
    const Trk::TrackParameters* upstreamParameters = track->trackParameters()->front();
    const Trk::TrackParameters* downstreamParameters = track->trackParameters()->back();
    double chi2=9999.;
    double ndof=9999.;

    if (candidateParameters == nullptr || upstreamParameters->momentum().mag() > candidateParameters->momentum().mag())
    {
        candidateParameters = upstreamParameters;
        candidateDownParameters = downstreamParameters;
        chi2 = track->fitQuality()->chiSquared();
        ndof = track->fitQuality()->numberDoF();
    }

    if ((candidateParameters == nullptr) || (candidateDownParameters == nullptr)) continue;

	m_nLayers.push_back(nLayers);

    std:cout << " chi2 = " << chi2 << std::endl;
    m_Chi2.push_back(chi2);
    m_DoF.push_back(ndof);

    m_nHit0.push_back(stationMap.count(0));
    m_nHit1.push_back(stationMap.count(1));
    m_nHit2.push_back(stationMap.count(2));
    m_nHit3.push_back(stationMap.count(3));
    n_clu.push_back(stationMap.size());

    m_charge.push_back( (int) candidateParameters->charge() );

    m_xup.push_back(candidateParameters->position().x());
    m_yup.push_back(candidateParameters->position().y());
    m_zup.push_back(candidateParameters->position().z());
    m_pxup.push_back(candidateParameters->momentum().x());
    m_pyup.push_back(candidateParameters->momentum().y());
    m_pzup.push_back(candidateParameters->momentum().z());
    m_pup.push_back(sqrt( pow(candidateParameters->momentum().x(),2) + pow(candidateParameters->momentum().y(),2) + pow(candidateParameters->momentum().z(),2) ));

    m_xdown.push_back(candidateDownParameters->position().x());
    m_ydown.push_back(candidateDownParameters->position().y());
    m_zdown.push_back(candidateDownParameters->position().z());
    m_pxdown.push_back(candidateDownParameters->momentum().x());
    m_pydown.push_back(candidateDownParameters->momentum().y());
    m_pzdown.push_back(candidateDownParameters->momentum().z());
    m_pdown.push_back(sqrt( pow(candidateDownParameters->momentum().x(),2) + pow(candidateDownParameters->momentum().y(),2) + pow(candidateDownParameters->momentum().z(),2) ));

    // fill extrapolation vectors with filler values that get changed iif the track extrapolation succeeds
    m_xVetoNu.push_back(-10000);
    m_yVetoNu.push_back(-10000);
    m_thetaxVetoNu.push_back(-10000);
    m_thetayVetoNu.push_back(-10000);
    m_xVetoStation1.push_back(-10000);
    m_yVetoStation1.push_back(-10000);
    m_thetaxVetoStation1.push_back(-10000);
    m_thetayVetoStation1.push_back(-10000);
    m_xVetoStation2.push_back(-10000);
    m_yVetoStation2.push_back(-10000);
    m_thetaxVetoStation2.push_back(-10000);
    m_thetayVetoStation2.push_back(-10000);
    m_xTrig.push_back(-10000);
    m_yTrig.push_back(-10000);
    m_thetaxTrig.push_back(-10000);
    m_thetayTrig.push_back(-10000);
    m_xPreshower1.push_back(-10000);
    m_yPreshower1.push_back(-10000);
    m_thetaxPreshower1.push_back(-10000);
    m_thetayPreshower1.push_back(-10000);
    m_xPreshower2.push_back(-10000);
    m_yPreshower2.push_back(-10000);
    m_thetaxPreshower2.push_back(-10000);
    m_thetayPreshower2.push_back(-10000);
    m_xCalo.push_back(-10000);
    m_yCalo.push_back(-10000);
    m_thetaxCalo.push_back(-10000);
    m_thetayCalo.push_back(-10000);

    // extrapolate track from IFT
    if (stationMap.count(0) > 0) { // extrapolation crashes if the track does not start in the IFT, as it is too far away to extrapolate
      Amg::Vector3D position = candidateParameters->position();
      Amg::Vector3D momentum = candidateParameters->momentum();
      Acts::BoundVector params = Acts::BoundVector::Zero();
      params[Acts::eBoundLoc0] = -position.y();
      params[Acts::eBoundLoc1] = position.x();
      params[Acts::eBoundPhi] = momentum.phi();
      params[Acts::eBoundTheta] = momentum.theta();
      params[Acts::eBoundQOverP] = candidateParameters->charge() / momentum.mag();
      params[Acts::eBoundTime] = 0;
      auto startSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, position.z()), Acts::Vector3(0, 0, 1));
      Acts::BoundTrackParameters startParameters(std::move(startSurface), params, candidateParameters->charge());

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
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Veto1 =m_extrapolationTool->propagate(ctx, startParameters, *targetSurface_Veto1, Acts::forward);
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
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Veto2 =m_extrapolationTool->propagate(ctx, startParameters, *targetSurface_Veto2, Acts::forward);
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
      std::unique_ptr<const Acts::BoundTrackParameters> targetParameters_Trig =m_extrapolationTool->propagate(ctx, startParameters, *targetSurface_Trig, Acts::forward); // must extrapolate forward to trig plane if track starts in IFT
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
      Amg::Vector3D positionDown = candidateDownParameters->position();
      Amg::Vector3D momentumDown = candidateDownParameters->momentum();
      Acts::BoundVector paramsDown = Acts::BoundVector::Zero();
      paramsDown[Acts::eBoundLoc0] = -positionDown.y();
      paramsDown[Acts::eBoundLoc1] = positionDown.x();
      paramsDown[Acts::eBoundPhi] = momentumDown.phi();
      paramsDown[Acts::eBoundTheta] = momentumDown.theta();
      paramsDown[Acts::eBoundQOverP] = candidateDownParameters->charge() / momentumDown.mag();
      paramsDown[Acts::eBoundTime] = 0;
      auto startSurfaceDown = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3(0, 0, positionDown.z()), Acts::Vector3(0, 0, 1));
      Acts::BoundTrackParameters startParametersDown(std::move(startSurfaceDown), paramsDown, candidateDownParameters->charge());

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

  /*
  // Here we apply the signal selection
  // Very simple/unrealistic to start
  if (m_vetoUpstream == 0 || m_vetoDownstream == 0 ||
        m_triggerTotal == 0 ||
        m_preshower0 == 0 || m_preshower1 == 0 ||
        // m_ecalTotal == 0 ||
        candidateParameters == nullptr)
      return StatusCode::SUCCESS;
  */
  m_tree->Fill();

  return StatusCode::SUCCESS;
}


StatusCode NtupleDumperAlg::finalize() 
{
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
  m_run_number = 0;
  m_event_number = 0;
  m_event_time = 0;
  m_bcid = 0;

  m_tbp=0;
  m_tap=0;
  m_inputBits=0;
  m_inputBitsNext=0;

  for(int ii=0;ii<15;ii++) {
      m_wave_localtime[ii]=0;
      m_wave_peak[ii]=0;
      m_wave_width[ii]=0;
      m_wave_charge[ii]=0;

      m_wave_raw_peak[ii]=0;
      m_wave_raw_charge[ii]=0;
      m_wave_baseline_mean[ii]=0;
      m_wave_baseline_rms[ii]=0;
      m_wave_status[ii]=0;
  }

  m_calo_total=0;
  m_calo_rawtotal=0;

  m_Calo0_Edep=0;
  m_Calo1_Edep=0;
  m_Calo2_Edep=0;
  m_Calo3_Edep=0;
  m_Calo_Total_Edep=0;
  m_Preshower12_Edep=0;
  m_Preshower13_Edep=0;

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
  
  n_clu.clear();
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

  m_truthLeptonMomentum = 0;
  m_truthBarcode = 0;
  m_truthPdg = 0;
}
