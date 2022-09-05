#!/usr/bin/env python

import sys
from AthenaCommon.Logging import log, logging
from AthenaCommon.Constants import DEBUG, VERBOSE, INFO
from AthenaCommon.Configurable import Configurable
from CalypsoConfiguration.AllConfigFlags import ConfigFlags
from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
from AthenaPoolCnvSvc.PoolWriteConfig import PoolWriteCfg
from FaserByteStreamCnvSvc.FaserByteStreamCnvSvcConfig import FaserByteStreamCnvSvcCfg
from TrackerPrepRawDataFormation.TrackerPrepRawDataFormationConfig import FaserSCT_ClusterizationCfg
from FaserActsKalmanFilter.CombinatorialKalmanFilterConfig import CombinatorialKalmanFilterCfg
from TrackerSegmentFit.TrackerSegmentFitConfig import SegmentFitAlgCfg

log.setLevel(DEBUG)
Configurable.configurableRun3Behavior = True

ConfigFlags.Input.Files = ['/eos/experiment/faser/raw/2022/008362/Faser-Physics-008362-00000.raw']
ConfigFlags.Output.ESDFileName = "CKF.ESD.pool.root"
ConfigFlags.addFlag("Output.xAODFileName", f"CKF.xAOD.root")
ConfigFlags.IOVDb.GlobalTag = "OFLCOND-FASER-02"
ConfigFlags.IOVDb.DatabaseInstance = "OFLP200"
ConfigFlags.Input.ProjectName = "data21"
ConfigFlags.Input.isMC = True
ConfigFlags.GeoModel.FaserVersion = "FASER-03"
ConfigFlags.Common.isOnline = False
ConfigFlags.GeoModel.Align.Dynamic = False
ConfigFlags.Beam.NumberOfCollisions = 0.
ConfigFlags.Detector.GeometryFaserSCT = True
ConfigFlags.lock()

acc = MainServicesCfg(ConfigFlags)
acc.merge(PoolWriteCfg(ConfigFlags))
acc.merge(FaserByteStreamCnvSvcCfg(ConfigFlags))
acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags, DataObjectName="SCT_LEVELMODE_RDOs", ClusterToolTimingPattern="X1X"))
acc.merge(SegmentFitAlgCfg(ConfigFlags))
acc.merge(CombinatorialKalmanFilterCfg(ConfigFlags, SummaryWriter=True, StatesWriter=True, PerformanceWriter=True))
acc.getEventAlgo("CombinatorialKalmanFilterAlg").OutputLevel = VERBOSE

from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
itemList = ["xAOD::EventInfo#*",
            "xAOD::EventAuxInfo#*",
            "FaserSCT_RDO_Container#*",
            "Tracker::FaserSCT_ClusterContainer#*",
            "TrackCollection#*"]
acc.merge(OutputStreamCfg(ConfigFlags, "xAOD", itemList))

print( "Writing out xAOD objects:" )
print( acc.getEventAlgo("OutputStreamxAOD").ItemList )

# logging.getLogger('forcomps').setLevel(VERBOSE)
# acc.foreach_component("*").OutputLevel = VERBOSE
# acc.foreach_component("*ClassID*").OutputLevel = INFO
# acc.getService("StoreGateSvc").Dump = True
# acc.getService("ConditionStore").Dump = True
# acc.printConfig(withDetails=True)
# ConfigFlags.dump()

# Hack to avoid problem with our use of MC databases when isMC = False
replicaSvc = acc.getService("DBReplicaSvc")
replicaSvc.COOLSQLiteVetoPattern = ""
replicaSvc.UseCOOLSQLite = True
replicaSvc.UseCOOLFrontier = False
replicaSvc.UseGeomSQLite = True

sc = acc.run(maxEvents=20)
sys.exit(not sc.isSuccess())
