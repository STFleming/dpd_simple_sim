#include "SimSystem.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SimSystem::SimSystem(float N, unsigned D, unsigned verbosity) {
    _N = N;
    _D = D;
    _verbosity = verbosity;
    _t = 0.0; // simulation starts at t=0.0
    _ts = 0;
    _dt = 0.05; // guess
    
    // create the global list of particles
    _particles = new std::vector<Particle *>();

    if(_N == 0) { // cannot have a problem with no size
       std::runtime_error("Problem must have a size >0\n");
    }

    // calculate the size of each spatial unit
    _unit_size = _N / (float)D;

    // create the list of _cubes
    _cubes = new std::vector<SpatialUnit *>(); 
    for(int x=0; x<_D; x++) {
        for(int y=0; y<_D; y++) {
            for(int z=0; z<_D; z++) {
                float fpx = (float)x * _unit_size;
                float fpy = (float)y * _unit_size;
                float fpz = (float)z * _unit_size;

                // create a SpatialUnit and add it to the cubes
                if(_verbosity >= 3)
                    printf("Constructing a cube (size=%f) at x:%f, y:%f, z:%f\n", _unit_size, fpx, fpy, fpz);
                    _cubes->push_back(new SpatialUnit(_unit_size,fpx,fpy,fpz,_verbosity)); // Add a new cube of size _N/_D at x,y,z
            }
        }
    }
}

// populates the universe with particles from a JSON file
// no longer needed?
//void SimSystem::populateFromJSON(std::string jsonfile) {
//    std::ifstream particle_file(jsonfile);
//    Json::Value root;
//    Json::Reader reader;    
//
//    // read in the JSON data
//    bool isJsonOK = (reader.parse(particle_file, root));
//    
//    if(isJsonOK) {
//        // read the particles from the JSON file
//        Json::Value particles = root["particles"];
//
//        // loop through the particles
//        for(int i=0; i < particles.size(); ++i) {
//           uint32_t cid = particles[i].get("id", 0xFFFFFFFF).asUInt();
//           position_t cpos;
//           cpos.x = particles[i].get("x", 0.0).asFloat(); 
//           cpos.y = particles[i].get("y", 0.0).asFloat(); 
//           cpos.z = particles[i].get("z", 0.0).asFloat(); 
//
//           // create the new particle
//           Particle *p = new Particle(cpos); // initialise with the position
//           p->setID(cid); // set it's ID 
//   
//           // add it to the universe
//           _particles->push_back(p);
//           allocateParticleToSpatialUnit(p);
//        }
//    } else {
//        std::runtime_error("Error: the input JSON file could not be parsed\n");
//    }
//}

// emits the state of the system as a JSON file
void SimSystem::emitJSON(std::string jsonfile) {
    std::ofstream out;
    out.open(jsonfile);
    out <<"{\n";
    out << "  \"particles\":[\n";

    // iterate through the particles and write the JSON 
    for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
        Particle *p = *i;
        out << "\t{\"id\":"<<p->getID()<<", \"x\":"<<p->getPos().x<<", \"y\":"<<p->getPos().y<<", \"z\":"<<p->getPos().z<<"}";
        if(p->getID() != (_particles->size()-1) )
            out << ",\n";
        else
            out << "\n";
    }

    out << "]}\n";
    out.close();
    return;
}

// iterators
SimSystem::iterator SimSystem::begin() { return _cubes->begin(); }
SimSystem::iterator SimSystem::end() { return _cubes->end(); }
SimSystem::p_iterator SimSystem::p_begin() { return _particles->begin(); }
SimSystem::p_iterator SimSystem::p_end() { return _particles->end(); }

/**! Destructor: cleans up the _cubes */
SimSystem::~SimSystem() {
    for(iterator i=begin(), e=end(); i!=e; ++i){
        SpatialUnit *cur = *i;
        if(_verbosity>=4)
            printf("removing cube (size=%f) at x:%f, y:%f, z:%f\n", cur->getSize(), cur->getPos().x, cur->getPos().y, cur->getPos().z);
        delete cur;
    }  
    delete _cubes;
    delete _particles; // this should probably iterate through and delete them all not the spatial units...
}

// allocates a particles to spatial processing units
void SimSystem::allocateParticleToSpatialUnit(Particle *p) {
    for(iterator i=begin(), ie=end(); i!=ie; ++i){
       SpatialUnit *cur = *i;
       if(cur->checkPos(p->getPos())) {
          cur->addParticle(p);
          return;
       } 
    }    
    std::runtime_error("Could not find a SpatialUnit that could support the particle.\n");
}

// runs the simulation for a given number of timesteps
// emits the state of each particle given by the emitrate
void SimSystem::run(uint32_t period, float emitrate) {

    // TODO this is where the real computation needs to be done
    // currently it just moves the particles in a random direction (to test the interface)
    unsigned start_ts = _ts;
    clock_t last_emit = clock();
    while(_ts <= start_ts + period) { // run for period timesteps

        for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
            Particle *p = *i; 
            position_t cur_p = p->getPos();
            // add a bit to the position
            cur_p.x = cur_p.x + 0.001;
            if(cur_p.x >= _N)
                cur_p.x = cur_p.x - _N;

            cur_p.y = cur_p.y + 0.001;
            if(cur_p.y >= _N)
                cur_p.y = cur_p.y - _N;

            cur_p.z = cur_p.z + 0.001;
            if(cur_p.z >= _N)
                cur_p.z = cur_p.z - _N;

            //update the particle with the new position
            p->setPos(cur_p);

            if ((float(clock() - last_emit) / CLOCKS_PER_SEC) > emitrate) {
                emitJSON("state_frame.json");
                printf("update\n"); // update command set via stdout to nodejs server
                //printf("{\"id\":%d, \"x\":%.2f, \"y\":%.2f, \"z\":%.2f}\n",p->getID(), p->getPos().x, p->getPos().y, p->getPos().z); 
                fflush(stdout);
                last_emit = clock();
            }

        } 
 

        // update time
        _ts = _ts + 1;
        _t = _t + _dt;
    }

}

// adds a particle to the system
void SimSystem::addParticle(Particle *p){
    // add the particle to the global particles list
    p->setID(_particles->size()); // just use the current number of particles as a unique ID
    _particles->push_back(p);     
    allocateParticleToSpatialUnit(p); // allocate particle to a processor
}
