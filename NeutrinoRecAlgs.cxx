#include "NeutrinoRecAlgs.h"
#include "NeutrinoSimEvent/NeutrinoHit.h"
#include "NeutrinoReadoutGeometry/NeutrinoDetectorElement.h"
#include "NeutrinoReadoutGeometry/EmulsionDetectorManager.h"
#include "StoreGate/StoreGateSvc.h"
#include "NeutrinoIdentifier/EmulsionID.h"
#include "GeoPrimitives/CLHEPtoEigenConverter.h"

NeutrinoRecAlgs::NeutrinoRecAlgs(const std::string& name, ISvcLocator* pSvcLocator)
: AthHistogramAlgorithm(name, pSvcLocator) { m_hist = nullptr; }

NeutrinoRecAlgs::~NeutrinoRecAlgs() { }

StatusCode NeutrinoRecAlgs::initialize()
{


    ATH_CHECK(detStore()->retrieve(m_emulsiom, "Emulsion"));
    ATH_CHECK(detStore()->retrieve(m_sID, "EmulsionID"));

    //Prepare tree for emulsion hit data
    m_NuHit_tree = new TTree("m_NuHit_tree","Tree containing emulsion hits data");
    m_NuHit_tree->Branch("m_runnumber",&m_runnumber,"m_runnumber/I");
    m_NuHit_tree->Branch("m_event_id",&m_event_id,"m_event_id/I");
    m_NuHit_tree->Branch("m_x_start",&m_x_start,"m_x_start/F");
    m_NuHit_tree->Branch("m_y_start",&m_y_start,"m_y_start/F");
    m_NuHit_tree->Branch("m_z_start",&m_z_start,"m_z_start/F");
    m_NuHit_tree->Branch("m_x_end",&m_x_end,"m_x_end/F");
    m_NuHit_tree->Branch("m_y_end",&m_y_end,"m_y_end/F");
    m_NuHit_tree->Branch("m_z_end",&m_z_end,"m_z_end/F");
    m_NuHit_tree->Branch("m_slope_x",&m_slope_x,"m_slope_x/F");
    m_NuHit_tree->Branch("m_slope_y",&m_slope_y,"m_slope_y/F");
    m_NuHit_tree->Branch("m_slope_z",&m_slope_z,"m_slope_z/F");
    m_NuHit_tree->Branch("m_track_id_hit",&m_track_id_hit,"m_track_id_hit/I");
    m_NuHit_tree->Branch("m_track_id_MC",&m_track_id_MC,"m_track_id_MC/I");
    m_NuHit_tree->Branch("m_pdg_MC",&m_pdg_MC,"m_pdg_MC/I");
    m_NuHit_tree->Branch("m_px_MC",&m_px_MC,"m_px_MC/D");
    m_NuHit_tree->Branch("m_py_MC",&m_py_MC,"m_py_MC/D");
    m_NuHit_tree->Branch("m_pz_MC",&m_pz_MC,"m_pz_MC/D");
    m_NuHit_tree->Branch("m_energy_MC",&m_energy_MC,"m_energy_MC/D");
    m_NuHit_tree->Branch("m_kinetic_energy_MC",&m_kinetic_energy_MC,"m_kinetic_energy_MC/D");
    m_NuHit_tree->Branch("m_mass_MC",&m_mass_MC,"m_mass_MC/D");

    m_NuHit_tree->Branch("m_vx_prod_hit",&m_vx_prod_hit,"m_vx_prod_hit/F");
    m_NuHit_tree->Branch("m_vy_prod_hit",&m_vy_prod_hit,"m_vy_prod_hit/F");
    m_NuHit_tree->Branch("m_vz_prod_hit",&m_vz_prod_hit,"m_vz_prod_hit/F");
    m_NuHit_tree->Branch("m_vx_decay_hit",&m_vx_decay_hit,"m_vx_decay_hit/F");
    m_NuHit_tree->Branch("m_vy_decay_hit",&m_vy_decay_hit,"m_vy_decay_hit/F");
    m_NuHit_tree->Branch("m_vz_decay_hit",&m_vz_decay_hit,"m_vz_decay_hit/F");
    m_NuHit_tree->Branch("m_x_start_global",&m_x_start_global,"m_x_start_global/F");
    m_NuHit_tree->Branch("m_y_start_global",&m_y_start_global,"m_y_start_global/F");
    m_NuHit_tree->Branch("m_z_start_global",&m_z_start_global,"m_z_start_global/F");
    m_NuHit_tree->Branch("m_x_end_global",&m_x_end_global,"m_x_end_global/F");
    m_NuHit_tree->Branch("m_y_end_global",&m_y_end_global,"m_y_end_global/F");
    m_NuHit_tree->Branch("m_z_end_global",&m_z_end_global,"m_z_end_global/F");
    m_NuHit_tree->Branch("m_slope_x_global",&m_slope_x_global,"m_slope_x_global/F");
    m_NuHit_tree->Branch("m_slope_y_global",&m_slope_y_global,"m_slope_y_global/F");
    m_NuHit_tree->Branch("m_slope_z_global",&m_slope_z_global,"m_slope_z_global/F");


    // m_NuHit_tree->Branch("m_slope_x_global",&m_slope_x_global,"m_slope_x_global/F");
    // m_NuHit_tree->Branch("m_slope_y_global",&m_slope_y_global,"m_slope_y_global/F");
    // m_NuHit_tree->Branch("m_slope_z_global",&m_slope_z_global,"m_slope_z_global/F");
    m_NuHit_tree->Branch("m_emeloss",&m_emeloss,"m_emeloss/F");
    m_NuHit_tree->Branch("m_plate",&m_plate,"m_plate/I");
    m_NuHit_tree->Branch("m_base",&m_base,"m_base/I");
    m_NuHit_tree->Branch("m_emmodule",&m_emmodule,"m_emmodule/I");
    m_NuHit_tree->Branch("m_film",&m_film,"m_film/I");


    m_NuMCTruth_tree = new TTree("m_NuMCTruth_tree","Tree containing simulated truth data");
    m_NuMCTruth_tree->Branch("m_runnumber",&m_runnumber,"m_runnumber/I");
    m_NuMCTruth_tree->Branch("m_event_id_MC",&m_event_id_MC,"m_event_id_MC/I");
    m_NuMCTruth_tree->Branch("m_track_id",&m_track_id,"m_track_id/I");
    m_NuMCTruth_tree->Branch("m_pdg_id",&m_pdg_id,"m_pdg_id/I");
    m_NuMCTruth_tree->Branch("m_vx_prod",&m_vx_prod,"m_vx_prod/F");
    m_NuMCTruth_tree->Branch("m_vy_prod",&m_vy_prod,"m_vy_prod/F");
    m_NuMCTruth_tree->Branch("m_vz_prod",&m_vz_prod,"m_vz_prod/F");
    m_NuMCTruth_tree->Branch("m_vx_decay",&m_vx_decay,"m_vx_decay/F");
    m_NuMCTruth_tree->Branch("m_vy_decay",&m_vy_decay,"m_vy_decay/F");
    m_NuMCTruth_tree->Branch("m_vz_decay",&m_vz_decay,"m_vz_decay/F");
    m_NuMCTruth_tree->Branch("m_px",&m_px,"m_px/D");
    m_NuMCTruth_tree->Branch("m_py",&m_py,"m_py/D");
    m_NuMCTruth_tree->Branch("m_pz",&m_pz,"m_pz/D");
    m_NuMCTruth_tree->Branch("m_energy",&m_energy,"m_energy/D");
    m_NuMCTruth_tree->Branch("m_kinetic_energy",&m_kinetic_energy,"m_kinetic_energy/D");
    m_NuMCTruth_tree->Branch("m_mass",&m_mass,"m_mass/D");
    
    m_NuMCTruth_tree->Branch("m_num_in_particle",&m_num_in_particle,"m_num_in_particle/I");
    m_NuMCTruth_tree->Branch("m_num_out_particle",&m_num_out_particle,"m_num_out_particle/I");

    m_NuMCTruth_tree->Branch("m_trackid_begin_in_particle",&m_trackid_begin_in_particle,"m_trackid_begin_in_particle/I");
    m_NuMCTruth_tree->Branch("m_trackid_begin_out_particle",&m_trackid_begin_out_particle,"m_trackid_begin_out_particle/I");
    m_NuMCTruth_tree->Branch("m_trackid_end_in_particle",&m_trackid_end_in_particle,"m_trackid_end_in_particle/I");
    m_NuMCTruth_tree->Branch("m_trackid_end_out_particle",&m_trackid_end_out_particle,"m_trackid_end_out_particle/I");

    m_NuMCTruth_tree->Branch("m_pdg_in_particle",&m_pdg_in_particle);
    m_NuMCTruth_tree->Branch("m_pdg_out_particle",&m_pdg_out_particle);

    m_NuMCTruth_tree->Branch("m_trackid_in_particle",&m_trackid_in_particle);
    m_NuMCTruth_tree->Branch("m_trackid_out_particle",&m_trackid_out_particle);

    m_NuMCTruth_tree->Branch("m_status",&m_status,"m_status/I");



    // initialize a histogram 
    // letter at end of TH1 indicated variable type (D double, F float etc)
    m_hist = new TH1D("eLoss", "SCT Hit Energy Loss", 100, 0, 1); //first string is root object name, second is histogram title
    m_module = new TH2D("module", "SCT Hit Module", 3, -1.5, 1.5, 4, -0.5, 3.5 );
    m_moduleSide1 = new TH2D("moduleSide1", "SCT Hit Module", 3, -1.5, 1.5, 4, -0.5, 3.5 );
    m_moduleSide2 = new TH2D("moduleSide2", "SCT Hit Module", 3, -1.5, 1.5, 4, -0.5, 3.5 );
    ATH_CHECK(histSvc()->regHist("/HIST/eloss", m_hist));
    ATH_CHECK(histSvc()->regHist("/HIST/modules", m_module));
    ATH_CHECK(histSvc()->regHist("/HIST/modulesSide1", m_moduleSide1));
    ATH_CHECK(histSvc()->regHist("/HIST/modulesSide2", m_moduleSide2));

    m_plate_preshower = new TH2D("plate", "Scint Hit Plate", 3, -1, 1, 4, 0, 1 );
    m_plate_trigger = new TH2D("plate", "Scint Hit Plate", 3, -1, 1, 4, 0, 1 );
    m_plate_veto = new TH2D("plate", "Scint Hit Plate", 3, -1, 1, 4, 0, 1 );
    ATH_CHECK(histSvc()->regHist("/HIST/plate_preshower", m_plate_preshower));
    ATH_CHECK(histSvc()->regHist("/HIST/plate_trigger", m_plate_trigger));
    ATH_CHECK(histSvc()->regHist("/HIST/plate_veto", m_plate_veto));

    m_ecalEnergy = new TH1D("ecalEnergy", "Ecal Energy Fraction", 100, 0.0, 0.20);
    ATH_CHECK(histSvc()->regHist("/HIST/ecal_energy", m_ecalEnergy));

    ATH_CHECK(histSvc()->regTree("/HIST/m_NuHit_tree", m_NuHit_tree));
    ATH_CHECK(histSvc()->regTree("/HIST/m_NuMCTruth_tree", m_NuMCTruth_tree));


    // initialize data handle keys
    ATH_CHECK( m_mcEventKey.initialize() );
    // ATH_CHECK( m_mcG4EventKey.initialize() );

    ATH_CHECK( m_faserSiHitKey.initialize() );

    ATH_CHECK( m_emulsionHitKey.initialize() );

    ATH_CHECK( m_preshowerHitKey.initialize() );
    ATH_CHECK( m_triggerHitKey.initialize() );
    ATH_CHECK( m_vetoHitKey.initialize() );

    ATH_CHECK( m_ecalHitKey.initialize() );

    ATH_MSG_INFO( "Using GenEvent collection with key " << m_mcEventKey.key());
    // ATH_MSG_INFO( "Using G4 collection with key " << m_mcG4EventKey.key());
    ATH_MSG_INFO( "Using Faser SiHit collection with key " << m_faserSiHitKey.key());

    ATH_MSG_INFO( "Using Emulsion NeutrinoHit collection with key " << m_emulsionHitKey.key() );

    ATH_MSG_INFO( "Using Preshower ScintHit collection with key " << m_preshowerHitKey.key());
    ATH_MSG_INFO( "Using Trigger ScintHit collection with key " << m_triggerHitKey.key());
    ATH_MSG_INFO( "Using Veto ScintHit collection with key " << m_vetoHitKey.key());

    ATH_MSG_INFO( "Using CaloHit collection with key " << m_ecalHitKey.key());
    return StatusCode::SUCCESS;


}

StatusCode NeutrinoRecAlgs::execute()
{
    // Handles created from handle keys behave like pointers to the corresponding container
    SG::ReadHandle<McEventCollection> h_mcEvents(m_mcEventKey);
    ATH_MSG_INFO("Read McEventContainer with " << h_mcEvents->size() << " events");
    if (h_mcEvents->size() == 0) return StatusCode::FAILURE;


    SG::ReadHandle<NeutrinoHitCollection> h_emulsionHits(m_emulsionHitKey);
    if (h_emulsionHits.isValid())
    {
        ATH_MSG_INFO("Read NeutrinoHitCollection with " << h_emulsionHits->size() << " hits");
    }
    SG::ReadHandle<FaserSiHitCollection> h_siHits(m_faserSiHitKey);
    ATH_MSG_INFO("Read FaserSiHitCollection with " << h_siHits->size() << " hits");
    SG::ReadHandle<ScintHitCollection> h_preshowerHits(m_preshowerHitKey);
    ATH_MSG_INFO("Read ScintHitCollection/Preshower with " << h_preshowerHits->size() << " hits");
    SG::ReadHandle<ScintHitCollection> h_triggerHits(m_triggerHitKey);
    ATH_MSG_INFO("Read ScintHitCollection/Trigger with " << h_triggerHits->size() << " hits");
    SG::ReadHandle<ScintHitCollection> h_vetoHits(m_vetoHitKey);
    ATH_MSG_INFO("Read ScintHitCollectionVeto with " << h_vetoHits->size() << " hits");    
    SG::ReadHandle<CaloHitCollection> h_ecalHits(m_ecalHitKey);

    // Since we have no pile-up, there should always be a single GenEvent in the container
    const HepMC::GenEvent* ev = (*h_mcEvents)[0];
    if (ev == nullptr) 
    {
        ATH_MSG_FATAL("GenEvent pointer is null");
        return StatusCode::FAILURE;
    }
    ATH_MSG_INFO("Event contains " << ev->particles_size() << " truth particles" );
    

    double ePrimary = 0;
    if (ev->particles_size() > 0)
    {
        HepMC::GenEvent::particle_const_iterator p = ev->particles_begin();
        ePrimary = (*p)->momentum().e();

    }

    m_event_id = Gaudi::Hive::currentContext().eventID().event_number();
    m_runnumber = Gaudi::Hive::currentContext().eventID().run_number();

    //Fill tree with simulated particles (Genie)
    for (const HepMC::GenEvent* event : *h_mcEvents) {

    event->print();

    int count_prod = 0;
    int count_decay = 0;
    m_num_in_particle = 0;
    m_num_out_particle = 0;

    for (const HepMC::GenParticle* particle : event->particle_range()) {

    const HepMC::FourVector& momentum = particle->momentum();
    m_track_id = particle->barcode();
    m_px = momentum.x();
    m_py = momentum.y();
    m_pz = momentum.z();
    m_energy = momentum.t();
    m_mass = particle->generated_mass();
    m_pdg_id = particle->pdg_id();
    m_kinetic_energy = m_energy - m_mass;      

    m_event_id_MC = m_event_id;

    m_trackid_begin_in_particle = 0;
    m_trackid_begin_out_particle = 0;
    m_trackid_end_in_particle = 0;
    m_trackid_end_out_particle = 0;
    m_num_in_particle = 0;
    m_num_out_particle = 0;

    m_pdg_in_particle.clear();
    m_pdg_out_particle.clear();
    m_trackid_in_particle.clear();
    m_trackid_out_particle.clear();


    m_status = particle->status();
    int raw_num=0;
    /* Loop over particles at end vertex of particle */

    if (particle->production_vertex()){
    const auto vertex = particle->production_vertex()->position();
    m_vx_prod = vertex.x();
    m_vy_prod = vertex.y();
    m_vz_prod = vertex.z();
    m_num_in_particle = particle->production_vertex()->particles_in_size();
    int in_count=0;
    std::cout<<"number of vertex going in : "<< m_num_in_particle<<std::endl;    
    // for(HepMC::GenVertex::particles_in_const_iterator partin =  particle->production_vertex()->particles_in_const_begin();partin != particle->production_vertex()->particles_in_const_end();partin++){
    for ( HepMC::GenVertex::particle_iterator partin= particle->production_vertex()->particles_begin(HepMC::parents);
    partin != particle->production_vertex()->particles_end(HepMC::parents);
    ++partin ){
        if ( (*partin) ) {
            m_pdg_in_particle.emplace_back((*partin)->pdg_id());
            m_trackid_in_particle.emplace_back((*partin)->barcode());
            (*partin)->print();
            in_count++;
            std::cout<<"Particles : "<<particle->barcode()<<std::endl;
            std::cout<<"Incoming particles : "<<in_count<<std::endl;
            std::cout<<"barcode : "<<(*partin)->barcode()<<std::endl;
            raw_num++;

        }//if ( !(*partin)->production_vertex() )
    }//particle iterator
    }//if (particle->production_vertex())

    std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
    std::cout<<"                              "<<std::endl;
    std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;


    if (particle->end_vertex()){
    const auto vertex_decay = particle->end_vertex()->position();
    m_vx_decay = vertex_decay.x();
    m_vy_decay = vertex_decay.y();
    m_vz_decay = vertex_decay.z();
    m_num_out_particle = particle->end_vertex()->particles_out_size();
    int out_count=0;
    std::cout<<"number of vertex going out : "<<m_num_out_particle<<std::endl;
    // for(HepMC::GenVertex::particles_out_const_iterator partout = particle->end_vertex()->particles_out_const_begin();partout != particle->end_vertex()->particles_out_const_end();partout++){
    for ( HepMC::GenVertex::particle_iterator partout= particle->end_vertex()->particles_begin(HepMC::children);
    partout != particle->end_vertex()->particles_end(HepMC::children);
    ++partout ){

        if ( (*partout) ) {
            m_pdg_out_particle.emplace_back((*partout)->pdg_id());
            m_trackid_out_particle.emplace_back((*partout)->barcode()); 
            (*partout)->print();
            out_count++;
            std::cout<<"Particles : "<<particle->barcode()<<std::endl;
            std::cout<<"Outgoing particles : "<<out_count<<std::endl;
            std::cout<<"barcode : "<<(*partout)->barcode()<<std::endl;
            raw_num++;

    }//if ( !(*partout)->end_vertex() )
    }//GenVetex iterator
    }//if (particle->end_vertex())

    m_NuMCTruth_tree->Fill();    

    std::cout<<"##############################"<<std::endl;
    std::cout<<"raw_num = "<<raw_num<<std::endl;
    std::cout<<"##############################"<<std::endl;
    std::cout<<"##############################"<<std::endl;
    std::cout<<"##############################"<<std::endl;
    }//HepMC::GenParticle* particles
    }//HepMC::GenEvent* event : *h_mcEvents


    if (h_emulsionHits.isValid() && h_emulsionHits->size()!=0)
    {
        //Loop over all hits; print
        for (const NeutrinoHit& hit : *h_emulsionHits)
        {
            // hit.print();

        Double_t thickness_base = 0.21;//mm
        Double_t thickness_plate = 1.;//mm
        Double_t thickness_film = 0.07;//mm

        //module_num * base_num = 770 plates
        //Int_t module_num = 35;
        Int_t base_num = 22;


        //x, y : local coordinate
	    //z : converted to the global coorfinate.	    
        m_x_start = hit.localStartPosition()[0];
        m_y_start = hit.localStartPosition()[1];	      
        m_z_start = (thickness_plate+thickness_base+(2*thickness_film))*(base_num*hit.getModule() + hit.getBase())+(thickness_base+thickness_film)*hit.getFilm() + hit.localStartPosition()[2];

        m_x_end = hit.localEndPosition()[0];
        m_y_end = hit.localEndPosition()[1];	      
        m_z_end = (thickness_plate+thickness_base+(2*thickness_film))*(base_num*hit.getModule() + hit.getBase())+(thickness_base+thickness_film)*hit.getFilm() + hit.localEndPosition()[2];

        Double_t len_mirco_track = pow((hit.localEndPosition()[0] - hit.localStartPosition()[0])*(hit.localEndPosition()[0] - hit.localStartPosition()[0])+(hit.localEndPosition()[1] - hit.localStartPosition()[1])*(hit.localEndPosition()[1] - hit.localStartPosition()[1])+(hit.localEndPosition()[2] - hit.localStartPosition()[2])*(hit.localEndPosition()[2] - hit.localStartPosition()[2]),0.5);

	    //3D vector in local coordinate
        m_slope_x = (hit.localEndPosition()[0] - hit.localStartPosition()[0])/len_mirco_track;
        m_slope_y = (hit.localEndPosition()[1] - hit.localStartPosition()[1])/len_mirco_track;
        m_slope_z = (hit.localEndPosition()[2] - hit.localStartPosition()[2])/len_mirco_track;


        //PlateID
        m_plate = base_num*hit.getModule() + hit.getBase();
        m_film = hit.getFilm();
        m_base = hit.getBase();
        m_emmodule = hit.getModule();        
	
        m_track_id_hit = hit.trackNumber();
        m_emeloss = hit.energyLoss();
        
        if (!hit.particleLink().isValid())
            continue;
        HepMcParticleLink partLink = hit.particleLink();
	    const HepMC::GenParticle * truthParticle = partLink.cptr();
        m_pdg_MC = truthParticle->pdg_id();
        m_track_id_MC = truthParticle->barcode();

        m_px_MC = truthParticle->momentum().x();
        m_py_MC = truthParticle->momentum().y();
        m_pz_MC = truthParticle->momentum().z();
        m_energy_MC = truthParticle->momentum().t();
        m_mass_MC = truthParticle->generated_mass();
        m_kinetic_energy_MC = m_energy_MC - m_mass_MC;      

        if(truthParticle->production_vertex()){
        const auto truth_vertex_hit = truthParticle->production_vertex()->position();
        m_vx_prod_hit = truth_vertex_hit.x();
        m_vy_prod_hit = truth_vertex_hit.y();
        m_vz_prod_hit = truth_vertex_hit.z();
        }

        if(truthParticle->end_vertex()){
        const auto truth_vertex_hit_decay = truthParticle->end_vertex()->position();
        m_vx_decay_hit = truth_vertex_hit_decay.x();
        m_vy_decay_hit = truth_vertex_hit_decay.y();
        m_vz_decay_hit = truth_vertex_hit_decay.z();
        }

        Identifier id;
        const NeutrinoDD::NeutrinoDetectorElement *geoelement=NULL;
        id = m_sID->film_id(hit.getModule(), hit.getBase(), hit.getFilm());
        geoelement = m_emulsiom->getDetectorElement(id);      

        if (!geoelement)
            continue;

        const HepGeom::Point3D<double> globalStartPos = Amg::EigenTransformToCLHEP(geoelement->transformHit()) * HepGeom::Point3D<double>(hit.localStartPosition());
        const HepGeom::Point3D<double> globalEndPos   = Amg::EigenTransformToCLHEP(geoelement->transformHit()) * HepGeom::Point3D<double>(hit.localEndPosition());

        m_x_start_global = globalStartPos.x();
        m_y_start_global = globalStartPos.y();
        m_z_start_global = globalStartPos.z();

        m_x_end_global = globalEndPos.x();
        m_y_end_global = globalEndPos.y();
        m_z_end_global = globalEndPos.z();


        Double_t len_mirco_track_global = pow((m_x_end_global - m_x_start_global)*(m_x_end_global - m_x_start_global)+(m_y_end_global - m_y_start_global)*(m_y_end_global - m_y_start_global)+(m_z_end_global - m_z_start_global)*(m_z_end_global - m_z_start_global),0.5);

	    //3D vector in local coordinate
        m_slope_x_global = (hit.localEndPosition()[0] - hit.localStartPosition()[0])/len_mirco_track_global;
        m_slope_y_global = (hit.localEndPosition()[1] - hit.localStartPosition()[1])/len_mirco_track_global;
        m_slope_z_global = (hit.localEndPosition()[2] - hit.localStartPosition()[2])/len_mirco_track_global;

        m_NuHit_tree->Fill();

        }        
    }

    m_NuHit_tree->Write();    

    // The hit container might be empty because particles missed the wafers
    //if (h_siHits->size() == 0) return StatusCode::SUCCESS;

    if (h_siHits->size()!=0){
        // Loop over all hits; print and fill histogram
        for (const FaserSiHit& hit : *h_siHits)
        {
            // hit.print();
            m_hist->Fill( hit.energyLoss() );
            m_module->Fill( hit.getModule(), hit.getRow() );
            if (hit.getSensor() == 0)
            {
                m_moduleSide1->Fill( hit.getModule(), hit.getRow());
            }
            else
            {
                m_moduleSide2->Fill( hit.getModule(), hit.getRow());
            }

        }
    }

    if (h_preshowerHits->size()!=0){
        for (const ScintHit& hit : *h_preshowerHits)
        {
            // hit.print();
            m_hist->Fill(hit.energyLoss());
            m_plate_preshower->Fill(hit.getStation(),hit.getPlate(),hit.energyLoss());

        }
    }

    if (h_triggerHits->size()!=0){
        for (const ScintHit& hit : *h_triggerHits)
        {
            // hit.print();
            m_hist->Fill(hit.energyLoss());
            m_plate_trigger->Fill(hit.getStation(),hit.getPlate(),hit.energyLoss());

        }
    }

    if (h_vetoHits->size()!=0){
        for (const ScintHit& hit : *h_vetoHits)
        {
            // hit.print();
            m_hist->Fill(hit.energyLoss());
            m_plate_veto->Fill(hit.getStation(),hit.getPlate(),hit.energyLoss());

        }
    }

    if (h_ecalHits->size() != 0)
    {
        double ecalTotal = 0.0;
        for (const CaloHit& hit : *h_ecalHits)
        {
            ecalTotal += hit.energyLoss();
        }
        if (ePrimary > 0) m_ecalEnergy->Fill(ecalTotal/ePrimary);
    }
    return StatusCode::SUCCESS;
}


StatusCode NeutrinoRecAlgs::finalize()
{
    return StatusCode::SUCCESS;
}
