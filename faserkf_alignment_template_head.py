#!/usr/bin/env python

import sys
from AthenaCommon.Logging import log, logging
from AthenaCommon.Constants import INFO
from AthenaCommon.Configurable import Configurable
from CalypsoConfiguration.AllConfigFlags import ConfigFlags
from AthenaConfiguration.TestDefaults import defaultTestFiles
from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
from AthenaPoolCnvSvc.PoolWriteConfig import PoolWriteCfg
from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
from TrackerPrepRawDataFormation.TrackerPrepRawDataFormationConfig import FaserSCT_ClusterizationCfg
from TrackerSpacePointFormation.TrackerSpacePointFormationConfig import TrackerSpacePointFinderCfg
from TrackerSegmentFit.TrackerSegmentFitConfig import SegmentFitAlgCfg
from FaserActsKalmanFilter.GhostBustersConfig import GhostBustersCfg
from TruthSeededTrackFinder.TruthSeededTrackFinderConfig import TruthSeededTrackFinderCfg
from FaserSCT_GeoModel.FaserSCT_GeoModelConfig import FaserSCT_GeometryCfg
from FaserActsGeometry.ActsGeometryConfig import ActsTrackingGeometryToolCfg
from FaserActsGeometry.ActsGeometryConfig import ActsExtrapolationToolCfg
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from MagFieldServices.MagFieldServicesConfig import MagneticFieldSvcCfg
from AthenaConfiguration.ComponentFactory import CompFactory

THistSvc=CompFactory.getComps("THistSvc")

CKF2Alignment, FaserActsExtrapolationTool, FaserActsTrackingGeometryTool=CompFactory.getComps("CKF2Alignment", "FaserActsExtrapolationTool","FaserActsTrackingGeometryTool")
from FaserActsGeometry.ActsGeometryConfig import ActsTrackingGeometrySvcCfg

def FaserActsAlignmentCondAlgCfg(flags, **kwargs):
   acc = ComponentAccumulator()
   acc.addCondAlgo(CompFactory.FaserActsAlignmentCondAlg(name="FaserActsAlignmentCondAlg", **kwargs))
   return acc

def CKF2Alignment_OutputESDCfg(flags):
    acc = ComponentAccumulator()
    itemList = ["xAOD::EventInfo#*",
                "xAOD::EventAuxInfo#*",
		"TrackCollection#*"
                ]
    acc.merge(OutputStreamCfg(flags, "ESD",itemList))
    ostream = acc.getEventAlgo("OutputStreamESD")
    ostream.TakeItemsFromInput = True
    return acc

def CKF2AlignmentCfg(flags, **kwargs):
    acc = FaserSCT_GeometryCfg(flags)
    acc.merge(MagneticFieldSvcCfg(flags))
    acts_tracking_geometry_svc = ActsTrackingGeometrySvcCfg(flags)
    acc.merge(acts_tracking_geometry_svc )
    track_seed_tool = CompFactory.CircleFitTrackSeedTool()
    sigma_loc0 = 1.9e-2
    sigma_loc1 = 9e-1
    sigma_phi = 3.3e-2
    sigma_theta = 2.9e-4
    p = 1000
    sigma_p = 0.1 * p
    sigma_qop = sigma_p / (p * p)
    initial_variance_inflation = [100, 100, 100, 100, 1000]
    track_seed_tool.covLoc0 = initial_variance_inflation[0] * sigma_loc1 * sigma_loc1
    track_seed_tool.covLoc1 = initial_variance_inflation[1] * sigma_loc0 * sigma_loc0
    track_seed_tool.covPhi = initial_variance_inflation[2] * sigma_phi * sigma_phi
    track_seed_tool.covTheta = initial_variance_inflation[3] * sigma_theta * sigma_theta
    track_seed_tool.covQOverP = initial_variance_inflation[4] * sigma_qop * sigma_qop
    track_seed_tool.std_cluster = 0.0231
    track_seed_tool.removeIFT = False
    track_seed_tool.TrackCollection = "Segments"
    trajectory_states_writer_tool = CompFactory.RootTrajectoryStatesWriterTool()
    trajectory_states_writer_tool.noDiagnostics = kwargs.pop("noDiagnostics", True)
    trajectory_states_writer_tool1 = CompFactory.RootTrajectoryStatesWriterTool()
    trajectory_states_writer_tool1.noDiagnostics = kwargs.pop("noDiagnostics", True)
    trajectory_states_writer_tool1.FilePath = "track_states_ckf1_666NAME666.root"
    
    trajectory_summary_writer_tool = CompFactory.RootTrajectorySummaryWriterTool()
    trajectory_summary_writer_tool.noDiagnostics = kwargs.pop("noDiagnostics", True)
    trajectory_summary_writer_tool1 = CompFactory.RootTrajectorySummaryWriterTool()
    trajectory_summary_writer_tool1.FilePath = "track_summary_ckf_666NAME666.root"
    trajectory_summary_writer_tool1.noDiagnostics = kwargs.pop("noDiagnostics", True)
    actsExtrapolationTool = CompFactory.FaserActsExtrapolationTool("FaserActsExtrapolationTool")
    actsExtrapolationTool.MaxSteps = 1000
    result, actsTrackingGeometryTool = ActsTrackingGeometryToolCfg(flags)
    actsExtrapolationTool.TrackingGeometryTool = actsTrackingGeometryTool
    acc.merge(result)
    #actsExtrapolationTool.TrackingGeometryTool = CompFactory.FaserActsTrackingGeometryTool("TrackingGeometryTool")
    
    trajectory_performance_writer_tool = CompFactory.PerformanceWriterTool("PerformanceWriterTool")
    trajectory_performance_writer_tool.ExtrapolationTool = actsExtrapolationTool
    trajectory_performance_writer_tool.noDiagnostics = kwargs.pop("noDiagnostics", True)
    #track_finder_tool = CompFactory.SPSimpleTrackFinderTool()
    #track_finder_tool = CompFactory.SegmentFitClusterTrackFinderTool()
    globalchi2fit = CompFactory.CKF2Alignment(**kwargs)
    globalchi2fit.BiasedResidual = False
    globalchi2fit.ExtrapolationTool = actsExtrapolationTool
    globalchi2fit.TrackingGeometryTool=actsTrackingGeometryTool
    kalman_fitter1 = CompFactory.KalmanFitterTool(name="fitterTool1")
    kalman_fitter1.noDiagnostics = True
    kalman_fitter1.ActsLogging = "INFO"
    kalman_fitter1.SummaryWriter = False
    kalman_fitter1.StatesWriter = False
    kalman_fitter1.SeedCovarianceScale = 10
    kalman_fitter1.isMC = False
    kalman_fitter1.RootTrajectoryStatesWriterTool = trajectory_states_writer_tool1
    kalman_fitter1.RootTrajectorySummaryWriterTool = trajectory_summary_writer_tool1
    globalchi2fit.KalmanFitterTool1 = kalman_fitter1
    globalchi2fit.TrackSeed = track_seed_tool
    globalchi2fit.ActsLogging = "INFO"
    globalchi2fit.isMC = False
    globalchi2fit.nMax = 10
    globalchi2fit.chi2Max = 100000
    globalchi2fit.maxSteps = 5000
    histSvc= CompFactory.THistSvc()
    histSvc.Output +=  [ "CKF2Alignment DATAFILE='kfalignment_666RUN666_666NAME666.root' OPT='RECREATE'"]
    acc.addService(histSvc)
    acc.addEventAlgo(globalchi2fit)
#    acc.merge(CKF2Alignment_OutputESDCfg(flags))
    return acc

if __name__ == "__main__":

   log.setLevel(INFO)
   Configurable.configurableRun3Behavior = True
   
   # Configure
#   ConfigFlags.Input.Files = [ '/eos/user/k/keli/Faser/alignment/misalign.RDO.pool.root']
#   ConfigFlags.Input.Files = [ '/eos/experiment/faser/rec/2022/p0008/008023/Faser-Physics-008023-00000-p0008-xAOD.root']
#   ConfigFlags.Input.Files = [ '/eos/home-k/keli/Faser/skim/skim666RUN666.DAOD.pool.root']
   ConfigFlags.Input.Files = [
