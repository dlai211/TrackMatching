#!/usr/bin/env python3

"""
    Copyright (C) 2002-2022 CERN for the benefit of the FASER collaboration
"""

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from MagFieldServices.MagFieldServicesConfig import MagneticFieldSvcCfg


def NtupleDumperAlgCfg(flags, OutName, **kwargs):
    # Initialize GeoModel
    from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
    acc = FaserGeometryCfg(flags)

    acc.merge(MagneticFieldSvcCfg(flags))
    # acc.merge(FaserActsTrackingGeometrySvcCfg(flags))
    # acc.merge(FaserActsAlignmentCondAlgCfg(flags))

    actsExtrapolationTool = CompFactory.FaserActsExtrapolationTool("FaserActsExtrapolationTool")
    actsExtrapolationTool.MaxSteps = 10000
    actsExtrapolationTool.TrackingGeometryTool = CompFactory.FaserActsTrackingGeometryTool("TrackingGeometryTool")

    NtupleDumperAlg = CompFactory.NtupleDumperAlg("NtupleDumperAlg",**kwargs)
    NtupleDumperAlg.ExtrapolationTool = actsExtrapolationTool
    acc.addEventAlgo(NtupleDumperAlg)

    thistSvc = CompFactory.THistSvc()
    thistSvc.Output += [f"HIST2 DATAFILE='{OutName}' OPT='RECREATE'"]
    acc.addService(thistSvc)

    return acc

if __name__ == "__main__":

    import glob
    import sys
    import ROOT

    runno=int(sys.argv[1])
    num=int(sys.argv[2])
    filesPerJob=int(sys.argv[3])
    run_config=str(sys.argv[4]) 

    ptag="p0008"

    from AthenaCommon.Logging import log, logging
    from AthenaCommon.Constants import DEBUG, VERBOSE, INFO
    from AthenaCommon.Configurable import Configurable
    from CalypsoConfiguration.AllConfigFlags import ConfigFlags
    from AthenaConfiguration.TestDefaults import defaultTestFiles
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    # from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
    # Set up logging and new style config
    log.setLevel(INFO)
    Configurable.configurableRun3Behavior = True

    dataDir=f"/eos/experiment/faser/rec/2022/{ptag}/{runno:06d}"
    files=sorted(glob.glob(f"{dataDir}/Faser-Physics*"))
    fileListInitial=files[num*filesPerJob:(num+1)*filesPerJob]
    fileList=[]
    for fName in fileListInitial:
        try:
            fh=ROOT.TFile(fName)
            fileList.append(fName)
        except OSError:
            print("Warning bad file: ",fName)

    log.info(f"Analyzing Run {runno} files {num*filesPerJob} to {(num+1)*filesPerJob} (num={num})")
    log.info(f"Got {len(fileList)} files out of {len(fileListInitial)}")

    outName=f"Data-tuple-{run_config}-{runno:06d}-{num:05d}-{filesPerJob}.root"

    # Configure
    ConfigFlags.Input.Files = fileList
    ConfigFlags.IOVDb.GlobalTag = "OFLCOND-FASER-02"             # Always needed; must match FaserVersionS
    ConfigFlags.IOVDb.DatabaseInstance = "OFLP200"               # Use MC conditions for now
    ConfigFlags.Input.ProjectName = "data21"                     # Needed to bypass autoconfig
    ConfigFlags.Input.isMC = False                                # Needed to bypass autoconfig
    ConfigFlags.GeoModel.FaserVersion     = "FASERNU-03"           # FASER geometry
    ConfigFlags.Common.isOnline = False
    ConfigFlags.GeoModel.Align.Dynamic = False
    ConfigFlags.Beam.NumberOfCollisions = 0.

    ConfigFlags.Detector.GeometryFaserSCT = True

    ConfigFlags.lock()

    # Core components
    acc = MainServicesCfg(ConfigFlags)
    acc.merge(PoolReadCfg(ConfigFlags))

    # algorithm
    acc.merge(NtupleDumperAlgCfg(ConfigFlags, outName, UseFlukaWeights=True, CaloConfig=run_config))

    AthenaEventLoopMgr = CompFactory.AthenaEventLoopMgr()
    AthenaEventLoopMgr.EventPrintoutInterval=1000
    acc.addService(AthenaEventLoopMgr)

    # # Hack to avoid problem with our use of MC databases when isMC = False
    replicaSvc = acc.getService("DBReplicaSvc")
    replicaSvc.COOLSQLiteVetoPattern = ""
    replicaSvc.UseCOOLSQLite = True
    replicaSvc.UseCOOLFrontier = False
    replicaSvc.UseGeomSQLite = True

    # Timing
    #acc.merge(MergeRecoTimingObjCfg(ConfigFlags))

    # Dump config
    # logging.getLogger('forcomps').setLevel(VERBOSE)
    # acc.foreach_component("*").OutputLevel = VERBOSE
    # acc.foreach_component("*ClassID*").OutputLevel = INFO
    # acc.getCondAlgo("FaserSCT_AlignCondAlg").OutputLevel = VERBOSE
    # acc.getCondAlgo("FaserSCT_DetectorElementCondAlg").OutputLevel = VERBOSE
    # acc.getService("StoreGateSvc").Dump = True
    # acc.getService("ConditionStore").Dump = True
    # acc.printConfig(withDetails=True)
    # ConfigFlags.dump()

    # Execute and finish
    sc = acc.run(maxEvents=-1)

    # Success should be 0
    sys.exit(not sc.isSuccess())    
