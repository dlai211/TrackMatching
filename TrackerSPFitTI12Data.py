# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

import sys
from AthenaCommon.Logging import log, logging
from AthenaCommon.Constants import DEBUG, INFO, VERBOSE
from AthenaCommon.Configurable import Configurable
from CalypsoConfiguration.AllConfigFlags import ConfigFlags
from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
from AthenaPoolCnvSvc.PoolWriteConfig import PoolWriteCfg
from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
from WaveRecAlgs.WaveRecAlgsConfig import WaveformReconstructionCfg
from TrackerPrepRawDataFormation.TrackerPrepRawDataFormationConfig import FaserSCT_ClusterizationCfg
from TrackerSpacePointFormation.TrackerSpacePointFormationConfig import TrackerSpacePointFinderCfg
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from FaserSCT_GeoModel.FaserSCT_GeoModelConfig import FaserSCT_GeometryCfg

Tracker__TrackerSPFit, THistSvc=CompFactory.getComps("Tracker::TrackerSPFit", "THistSvc")


def TrackerSPFitBasicCfg(flags, **kwargs):
    """Return ComponentAccumulator for TrackerSPFit"""
    acc = FaserSCT_GeometryCfg(flags)
    kwargs.setdefault("SpacePointsSCTName", "SCT_SpacePointContainer")
    kwargs.setdefault("SCT_ClustersName", "SCT_ClusterContainer")
    kwargs.setdefault("MaxChi2", 100)
    kwargs.setdefault("UseBiasedResidual", True)
    acc.addEventAlgo(Tracker__TrackerSPFit(**kwargs))
   # attach ToolHandles
    return acc

def TrackerSPFit_OutputCfg(flags):
    """Return ComponentAccumulator with Output for SCT. Not standalone."""
    acc = ComponentAccumulator()
    acc.merge(OutputStreamCfg(flags, "ESD"))
    ostream = acc.getEventAlgo("OutputStreamESD")
    ostream.TakeItemsFromInput = True
    return acc

def TrackerSPFitCfg(flags, **kwargs):
    acc=TrackerSPFitBasicCfg(flags, **kwargs)
    histSvc= THistSvc()
    histSvc.Output += [ "TrackerSPFit DATAFILE='/eos/user/j/jlai/TestBeam/trackerspfit_8362_02.root' OPT='RECREATE'" ]
    acc.addService(histSvc)
    acc.merge(TrackerSPFit_OutputCfg(flags))
    return acc

if __name__ == "__main__":
  log.setLevel(DEBUG)
  Configurable.configurableRun3Behavior = True

  # Configure
  ConfigFlags.Input.Files = ['/eos/experiment/faser/rec/2022/p0008/008362/Faser-Physics-008362-00002-p0008-xAOD.root']
  #ConfigFlags.Input.Files = ['/eos/project/f/faser-commissioning/reco/Commissioning2021/r0005/004979/Faser-Physics-004979-00000-r0005-xAOD.root']
  ConfigFlags.Output.ESDFileName = "mySeeds_8362_02.ESD.pool.root"
  ConfigFlags.IOVDb.GlobalTag = "OFLCOND-FASER-02"             # Always needed; must match FaserVersion
  ConfigFlags.IOVDb.DatabaseInstance = "OFLP200"               # Use MC conditions for now
  ConfigFlags.Input.ProjectName = "data21"                     # Needed to bypass autoconfig
  ConfigFlags.Input.isMC = False
  # Needed to bypass autoconfig
  ConfigFlags.GeoModel.FaserVersion     = "FASERNU-02"           # FASER cosmic ray geometry (station 2 only)
#  ConfigFlags.Detector.GeometryFaserSCT = True
  ConfigFlags.Common.isOnline = False
  ConfigFlags.GeoModel.Align.Dynamic = False
  ConfigFlags.Beam.NumberOfCollisions = 0.

  ConfigFlags.lock()

  # Core components
  acc = MainServicesCfg(ConfigFlags)
  acc.merge(PoolReadCfg(ConfigFlags))
  acc.merge(PoolWriteCfg(ConfigFlags))

#  from FaserByteStreamCnvSvc.FaserByteStreamCnvSvcConfig import FaserByteStreamCnvSvcCfg
#  acc.merge(FaserByteStreamCnvSvcCfg(ConfigFlags))
#  acc.merge(WaveformReconstructionCfg(ConfigFlags))
#  acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags, DataObjectName="SCT_RDOs"))
#  acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags))
#  acc.merge(TrackerSpacePointFinderCfg(ConfigFlags))
  acc.merge(TrackerSPFitCfg(ConfigFlags))
#  from AthenaConfiguration.ComponentFactory import CompFactory
#  decoderTool = CompFactory.RawWaveformDecoderTool(name = "RawWaveformDecoderTool", 
##      CaloChannels = [0, 1, 2,3, 4, 5], 
##    PreshowerChannels = [6,7], 
##    TriggerChannels = [8, 9],
#      VetoChannels=[])
#  acc.addPublicTool(decoderTool)

  # Timing
  #acc.merge(MergeRecoTimingObjCfg(ConfigFlags))

#  replicaSvc = acc.getService("DBReplicaSvc")
#  replicaSvc.COOLSQLiteVetoPattern = ""
#  replicaSvc.UseCOOLSQLite = True
#  replicaSvc.UseCOOLFrontier = False
#  replicaSvc.UseGeomSQLite = True
  # Dump config
  logging.getLogger('forcomps').setLevel(VERBOSE)
  acc.foreach_component("*").OutputLevel = VERBOSE
  acc.foreach_component("*ClassID*").OutputLevel = INFO
  # acc.getCondAlgo("FaserSCT_AlignCondAlg").OutputLevel = VERBOSE
  # acc.getCondAlgo("FaserSCT_DetectorElementCondAlg").OutputLevel = VERBOSE
  acc.getService("StoreGateSvc").Dump = True
  acc.getService("ConditionStore").Dump = True
  acc.printConfig(withDetails=True)
  ConfigFlags.dump()

  # Execute and finish
  sc = acc.run(maxEvents=-1)

  # Success should be 0
  sys.exit(not sc.isSuccess())
