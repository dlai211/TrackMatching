#include "AthenaBaseComps/AthHistogramAlgorithm.h"
#include "GeneratorObjects/McEventCollection.h"
#include "TrackerSimEvent/FaserSiHitCollection.h"
#include "TrackerRawData/FaserSCT_RDO_Container.h"
#include "TrackerSimData/TrackerSimDataCollection.h"
#include <TH1.h>
#include <math.h>
#include <TProfile.h>

/* RDORead reading example - Ryan Rice-Smith, UC Irvine */

class RDOReadAlg : public AthHistogramAlgorithm
{
    public:
    RDOReadAlg(const std::string& name, ISvcLocator* pSvcLocator);

    virtual ~RDOReadAlg();

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

    private:
    TH1* m_hist;  // Example histogram
    TH1* m_incAnglHist;
    TProfile* m_hprof;

    // Read handle keys for data containers
    // Any other event data can be accessed identically
    // Note the key names ("BeamTruthEvent" or "SCT_Hits") are Gaudi properties and can be configured at run-time
    SG::ReadHandleKey<McEventCollection> m_mcEventKey       { this, "McEventCollection", "TruthEvent" };
    SG::ReadHandleKey<FaserSiHitCollection> m_faserSiHitKey { this, "FaserSiHitCollection", "SCT_Hits" };
    SG::ReadHandleKey<FaserSCT_RDO_Container> m_faserRdoKey { this, "FaserSCT_RDO_Container", "SCT_RDOs"};
    SG::ReadHandleKey<TrackerSimDataCollection> m_sctMap {this, "TrackerSimDataCollection", "SCT_SDO_Map"};
};