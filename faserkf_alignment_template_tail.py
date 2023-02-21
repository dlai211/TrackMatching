]
   ConfigFlags.IOVDb.GlobalTag = "OFLCOND-FASER-03"
   ConfigFlags.GeoModel.FaserVersion     = "FASERNU-03"           # FASER cosmic ray geometry (station 2 only)
   ConfigFlags.Input.isMC = True
   ConfigFlags.GeoModel.Align.Dynamic = False
   ConfigFlags.Output.doWriteESD = False
   ConfigFlags.Input.ProjectName = "data21"                     # Needed to bypass autoconfig
   ConfigFlags.Beam.NumberOfCollisions = 0.
   ConfigFlags.addFlag("Input.InitialTimeStamp", 0)
   ConfigFlags.fillFromArgs(sys.argv[1:])
#   ConfigFlags.Exec.SkipEvents = 666SKIP666
   ConfigFlags.Exec.MaxEvents = -1
   #ConfigFlags.Concurrency.NumThreads = 1
   ConfigFlags.lock()
   
   from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
   # Core components
   acc = MainServicesCfg(ConfigFlags)
   acc.merge(FaserGeometryCfg(ConfigFlags))
   acc.merge(PoolReadCfg(ConfigFlags))
   acc.merge(PoolWriteCfg(ConfigFlags))
   from FaserByteStreamCnvSvc.FaserByteStreamCnvSvcConfig import FaserByteStreamCnvSvcCfg
#   acc.merge(FaserByteStreamCnvSvcCfg(ConfigFlags))
#   acc.merge(FaserGeometryCfg(ConfigFlags))
#   from WaveRecAlgs.WaveRecAlgsConfig import WaveformReconstructionCfg
#   acc.merge(WaveformReconstructionCfg(ConfigFlags))
#   from WaveRecAlgs.WaveRecAlgsConfig import WaveformReconstructionOutputCfg    
#   acc.merge(WaveformReconstructionOutputCfg(ConfigFlags))
#   acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags, DataObjectName="SCT_RDOs"))
#   acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags, DataObjectName="SCT_LEVELMODE_RDOs", ClusterToolTimingPattern="X1X"))
   # Inner Detector
   acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags))
  # acc.merge(SegmentFitAlgCfg(ConfigFlags))
   acc.merge(TrackerSpacePointFinderCfg(ConfigFlags))
#   acc.merge(TrackerSPFitCfg(ConfigFlags))
#   acc.merge(TruthSeededTrackFinderCfg(ConfigFlags))
   acc.merge(SegmentFitAlgCfg(ConfigFlags, SharedHitFraction=0.61, MinClustersPerFit=5, TanThetaXZCut=0.083))
   acc.merge(GhostBustersCfg(ConfigFlags))
   acc.merge(CKF2AlignmentCfg(ConfigFlags))
   #from FaserActsKalmanFilter.CKF2Config import CKF2Cfg
   #acc.merge(CKF2Cfg(ConfigFlags))
#   replicaSvc = acc.getService("DBReplicaSvc")
#   replicaSvc.COOLSQLiteVetoPattern = ""
#   replicaSvc.UseCOOLSQLite = True
#   replicaSvc.UseCOOLFrontier = False
#   replicaSvc.UseGeomSQLite = True
#   acc.getEventAlgo("CKF2Alignment").AlignmentConstants = ConfigFlags.CKF2Alignment.AlignmentConstants
   logging.getLogger('forcomps').setLevel(INFO)
   acc.foreach_component("*").OutputLevel = INFO
   acc.foreach_component("*ClassID*").OutputLevel =INFO
   #acc.getService("StoreGateSvc").Dump = True
   #acc.getService("ConditionStore").Dump = True
   #acc.printConfig(withDetails=True)
   #ConfigFlags.dump()
   
   # Execute and finish
   sc = acc.run()
   #sc = acc.run(maxEvents=-1000000,skipEvents=0)
   sys.exit(not sc.isSuccess())
