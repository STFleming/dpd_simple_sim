#include "SimSystem.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SimSystem::SimSystem(float N, float dt, float r_c, unsigned D, unsigned verbosity) {
    _N = N;
    _D = D;
    _verbosity = verbosity;
    _t = 0.0; // simulation starts at t=0.0
    _ts = 0;
    _dt = dt; 
    _r_c = r_c; // interaction cutoff 
    
    // create the global list of particles
    _particles = new std::vector<Particle *>();

    // create a global list of particle pairs used in the seq_run naive solver
    _seq_pairs = new PartPair();

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

// emits the state of the system as a JSON file
void SimSystem::emitJSON(std::string jsonfile) {
    std::ofstream out;
    out.open(jsonfile);
    out <<"{\n";
    out << "  \"particles\":[\n";

    // iterate through the particles and write the JSON 
    for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
        Particle *p = *i;
        // check to make sure the particle position makes sense
        if (!isnan(p->getPos().x()) && !isnan(p->getPos().y()) && !isnan(p->getPos().z())) {
            out << "\t{\"id\":"<<p->getID()<<", \"x\":"<<p->getPos().x()<<", \"y\":"<<p->getPos().y()<<", \"z\":"<<p->getPos().z()<<", \"type\":"<<p->getType()<<"}";
            if(p->getID() != (_particles->size()-1) )
                out << ",\n";
            else
                out << "\n";
        } else {
           printf("Error: NaN encountered when trying to export state of particle: %u\n", p->getID());
           exit(1);
        }
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
       if(cur->checkPos(vec2pos(p->getPos()))) {
          cur->addParticle(p);
          return;
       } 
    }    
    std::runtime_error("Could not find a SpatialUnit that could support the particle.\n");
}

// runs the simulation for a given number of timesteps sequentially there is no parallelism here
// emits the state of each particle given by the emitrate
// this is the naive implementation of this algorithm (horrifically slow)
void SimSystem::seq_run(uint32_t period, float emitrate) {

    unsigned start_ts = _ts;
    clock_t last_emit = clock();
    while(_ts <= start_ts + period) { // run for period timesteps

        for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
            for(p_iterator j=p_begin(), je=p_end(); j!=je; ++j) {
                Particle *p1 = *i;
                Particle *p2 = *j;
                if(p1 != p2) { // particles do not apply forces to themselves
                    if ( p1->getPos().toroidal_dist(p2->getPos(), _N) <= _r_c ) {
                        // this particle is in range
                        // make sure that we have not done this pairwise interaction already
                        bool already_pair = false;
                        for(PartPair::iterator pi=_seq_pairs->begin(), pie=_seq_pairs->end(); pi!=pie; ++pi) {
                            std::tuple<Particle *, Particle*> cur_pair = *pi;
                            bool pair_one = (p1->getID() == std::get<0>(cur_pair)->getID()) && (p2->getID() == std::get<1>(cur_pair)->getID());
                            bool pair_two = (p2->getID() == std::get<0>(cur_pair)->getID()) && (p1->getID() == std::get<1>(cur_pair)->getID());
                            if (pair_one || pair_two) {
                               already_pair = true;
                               break;
                            }
                        }
                        if(!already_pair) {
                            // do the force update
                            p1->callConservative(p2);
                            p1->callDrag(p2);
                            p1->callRandom(p2);
                            std::tuple<Particle *, Particle*> part_pair = std::make_tuple(p1, p2);
                            _seq_pairs->push_back(part_pair);
                        } 
                    } 
                }
            }
        } 


        // update v_i and r_i for each particle
        for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
             Particle *p = *i;
             float mass = p->getMass();
             Vector3D acceleration = p->getForce()/mass;
             Vector3D delta_v = (p->getForce()/mass) * _dt;
             // update velocity
             p->setVelo(p->getVelo() + delta_v); 

             //if((p->getID() == 0) && (p->getForce().mag() != 0.0)) {
             //}

             //const float speedlimit = 1.0*0.95;
             //if (p->getVelo().mag() >= speedlimit) {
             //      printf("[particle:%u] Speed limit has been exceeded velocity: %.8f \n",  p->getID(), p->getVelo().mag());
             //      Vector3D delta_r = p->getVelo()*_dt + acceleration*0.5*_dt*_dt;
             //      printf("------------------\n"); // dump the stats for this particle            
             //      printf("[particle:%u] force = %.8f\n", p->getID(), p->getForce().mag());
             //      printf("[particle:%u] acceleration = %.8f\n", p->getID(), acceleration.mag());
             //      printf("[particle:%u] delta_v = %.8f\n", p->getID(), delta_v.mag());
             //      printf("[particle:%u] velocity = %.8f\n", p->getID(), p->getVelo().mag());
             //      printf("[particle:%u] distance travlled = %.8f\n", p->getID(), delta_r.mag());
             //      float scale = sqrt((speedlimit*speedlimit)/((delta_v.x()*delta_v.x())+(delta_v.y()*delta_v.y())+(delta_v.z()*delta_v.z())));
             //      //delta_v.set(scale*delta_v.x(), scale*delta_v.y(), scale*delta_v.z());
             //      //printf("  adjusted to delta_v: %.2f  using scale factor:%.2f\n", delta_v.mag(), scale);
             //      exit(0);
      
             //}

             // euler update position & include wraparound
             //Vector3D point = p->getPos() +p->getVelo()*_dt;

             // velocity verlet
             Vector3D point = p->getPos() + p->getVelo()*_dt + acceleration*0.5*_dt*_dt; 


             // wraparound x direction
             if(point.x() < 0.0)
                  point.x(point.x() + _N); 
             if(point.x() >= _N)
                  point.x(point.x() - _N); 

             // wrapointaround y direction
             if(point.y() < 0.0)
                  point.y(point.y() + _N); 
             if(point.y() >= _N)
                  point.y(point.y() - _N); 

             // wrapointaround z direction
             if(point.z() < 0.0)
                  point.z(point.z() + _N); 
             if(point.z() >= _N)
                  point.z(point.z() - _N); 

             // update the position
             p->setPos(point); 

             // clear force 
             p->setForce(Vector3D(0.0, 0.0, 0.0));
        } 

        // cleanup: 
        _seq_pairs->clear(); // we no longer need to track the pairs

        // do we want to emit the state of the simulation
        if ((float(clock() - last_emit) / CLOCKS_PER_SEC) > emitrate) {
            emitJSON("state.json");
            printf("update\n"); // update command set via stdout to nodejs server
            fflush(stdout);
            last_emit = clock();
        }


        // update time
        _ts = _ts + 1;
        _t = _t + _dt;
    }
}

//void SimSystem::seq_run(uint32_t period, float emitrate) {
//
//    unsigned start_ts = _ts;
//    clock_t last_emit = clock();
//    while(_ts <= start_ts + period) { // run for period timesteps
//
//        for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
//            const float p_speed = 0.01;
//            Particle *p = *i; 
//            position_t cur_p = p->getPos();
//            // add a bit to the position
//            cur_p.x = cur_p.x + p_speed;
//            if(cur_p.x >= _N)
//                cur_p.x = cur_p.x - _N;
//
//            cur_p.y = cur_p.y + p_speed;
//            if(cur_p.y >= _N)
//                cur_p.y = cur_p.y - _N;
//
//            cur_p.z = cur_p.z + p_speed;
//            if(cur_p.z >= _N)
//                cur_p.z = cur_p.z - _N;
//
//            //update the particle with the new position
//            p->setPos(cur_p);
//
//            if ((float(clock() - last_emit) / CLOCKS_PER_SEC) > emitrate) {
//                emitJSON("state.json");
//                printf("update\n"); // update command set via stdout to nodejs server
//                //printf("{\"id\":%d, \"x\":%.2f, \"y\":%.2f, \"z\":%.2f}\n",p->getID(), p->getPos().x, p->getPos().y, p->getPos().z); 
//                fflush(stdout);
//                last_emit = clock();
//            }
//
//        } 
// 
//
//        // update time
//        _ts = _ts + 1;
//        _t = _t + _dt;
//    }
//
//}

// adds a particle to the system
void SimSystem::addParticle(Particle *p){
    // add the particle to the global particles list
    p->setID(_particles->size()); // just use the current number of particles as a unique ID
    _particles->push_back(p);     
    allocateParticleToSpatialUnit(p); // allocate particle to a processor
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

