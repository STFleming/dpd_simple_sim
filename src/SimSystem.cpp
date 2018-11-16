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

                spatial_unit_address_t c_addr = {x,y,z};
                _cubes->push_back(new SpatialUnit(_unit_size,fpx,fpy,fpz,c_addr, _verbosity)); // Add a new cube of size _N/_D at x,y,z
            }
        }
    }

    // construct the neighbours list for each SpatialUnit
    // not worrying too much about efficiency for the time being for this
    for(int x=0; x<_D; x++) {
        for(int y=0; y<_D; y++) {
            for(int z=0; z<_D; z++) {
                  
                // calculate offsets
                int x_neg, y_neg, z_neg; 
                int x_pos, y_pos, z_pos;
                
                // assign the x offsets
                if(x==0) {
                  x_neg = _D-1; 
                  x_pos = x+1;
                } else if (x == (_D-1)) {
                  x_neg = x-1;
                  x_pos = 0; 
                } else {
                  x_neg = x-1;
                  x_pos = x+1;
                }

                // assign the y offsets
                if(y==0) {
                  y_neg = _D-1; 
                  y_pos = y+1;
                } else if (y == (_D-1)) {
                  y_neg = y-1;
                  y_pos = 0; 
                } else {
                  y_neg = y-1;
                  y_pos = y+1;
                }

                // assign the z offsets
                if(z==0) {
                  z_neg = _D-1; 
                  z_pos = z+1;
                } else if (z == (_D-1)) {
                  z_neg = z-1;
                  z_pos = 0; 
                } else {
                  z_neg = z-1;
                  z_pos = z+1;
                }


                 // x ,y ,z
                 spatial_unit_address_t curr = {x, y, z};
                 SpatialUnit *curr_spu = getSpatialUnit(curr); 
                 // add all 27 neighbours to the neighbours list
                 // z = -1
                   // { -1,-1,-1 },  { -1,0,-1 },  { -1, +1,-1 }
                      curr.x = x_neg; curr.y = y_neg; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_neg; curr.y = y; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_neg; curr.y = y_pos; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                   // { 0,-1, -1 },  { 0, 0,-1 },  { 0, +1, -1 }
                      curr.x = x; curr.y = y_neg; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x; curr.y = y; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x; curr.y = y_pos; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                   // { +1,-1,-1 },  { +1,0,-1 },  { +1, +1,-1 }
                      curr.x = x_pos; curr.y = y_neg; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_pos; curr.y = y; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_pos; curr.y = y_pos; curr.z = z_neg;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                 // z = 0
                   // { -1,-1,0 },  { -1,0,0 },  { -1, +1,0 }
                      curr.x = x_neg; curr.y = y_neg; curr.z = z;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_neg; curr.y = y; curr.z = z;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_neg; curr.y = y_pos; curr.z = z;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                   // { 0,-1, 0 },  { 0, 0, 0 },  { 0, +1, 0 }
                      curr.x = x; curr.y = y_neg; curr.z = z;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      // skipping! one is not a neighbour of oneself
                      //curr.x = x; curr.y = y; curr.z = z;
                      //curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x; curr.y = y_pos; curr.z = z;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                   // { +1,-1, 0 },  { +1,0, 0 },  { +1, +1, 0 }
                      curr.x = x_pos; curr.y = y_neg; curr.z = z;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_pos; curr.y = y; curr.z = z;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_pos; curr.y = y_pos; curr.z = z;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 


                 // z = +1
                   // { -1,-1,+1 },  { -1,0,+1},  { -1, +1,+1 }
                      curr.x = x_neg; curr.y = y_neg; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_neg; curr.y = y; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_neg; curr.y = y_pos; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                   // { 0,-1, +1 },  { 0, 0, +1 },  { 0, +1, +1 }
                      curr.x = x; curr.y = y_neg; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x; curr.y = y; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x; curr.y = y_pos; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                   // { +1,-1, +1 },  { +1,0, +1 },  { +1, +1, +1 }
                      curr.x = x_pos; curr.y = y_neg; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_pos; curr.y = y; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

                      curr.x = x_pos; curr.y = y_pos; curr.z = z_pos;
                      curr_spu->addNeighbour(getSpatialUnit(curr)); 

            }
        }
    }

}

// returns a spatial unit when provided with a spatial unti address
SpatialUnit * SimSystem::getSpatialUnit(spatial_unit_address_t x) {
    for(SimSystem::iterator i=begin(); i!=end(); ++i){
        SpatialUnit *c = *i;
        spatial_unit_address_t c_addr = c->getAddr();
        if( (x.x == c_addr.x) && (x.y == c_addr.y) && (x.z == c_addr.z)) {
           return c;
        } 
    }
    return NULL;
}

// bonds particle i to particle j in the universe 
void SimSystem::bond(Particle *i, Particle *j, std::function<void(Particle *me, Particle *other)> bondf){
    assert(!i->isBonded());
    i->setBond(j, bondf);
    if(i->getPos().toroidal_dist(j->getPos(), _N) >= _r_c){
         printf("Error: particle %d (%s) and %d (%s) are too far away (dist=%.2f) to be bonded\n", i->getID(), i->getPos().str().c_str(), j->getID(), j->getPos().str().c_str(), i->getPos().toroidal_dist(j->getPos(), _N));
         exit(EXIT_FAILURE);
    } 
    return;
}

// emits the state of the system from the spatial units (relative addressing)
void SimSystem::emitJSONFromSU(std::string jsonfile) {
    std::ofstream out;
    out.open(jsonfile);
    out <<"{\n";
    out << "  \"beads\":[\n";

    // iterate through the spatial units and grab the offsets
    for(iterator i=begin(); i!=end(); ++i){
        SpatialUnit *s = *i;
        spatial_unit_address_t pos = s->getAddr(); 
    
        float x_off = pos.x * s->getSize(); 
        float y_off = pos.y * s->getSize(); 
        float z_off = pos.z * s->getSize(); 

        // iterate through all particles of this spatial unit
        for(SpatialUnit::iterator ip=s->begin(); ip!=s->end(); ++ip){
           Particle *cp = *ip;
           float cp_x = cp->getPos().x() + x_off;
           float cp_y = cp->getPos().y() + y_off;
           float cp_z = cp->getPos().z() + z_off;

           // print the position of the particles to see if they are actually stuck just on the 
           // display or in the actual simulation
           //printf("bead:%d pos:<%.4f, %.4f, %.4f>\n",cp->getID(), cp_x, cp_y, cp_z);

           // check to make sure that the particle position makes sense
           if (!isnan(cp_x) && !isnan(cp_y) && !isnan(cp_z)) {
                out << "\t{\"id\":"<<cp->getID()<<", \"x\":"<<cp_x<<", \"y\":"<<cp_y<<", \"z\":"<<cp_z<<", \"vx\":"<<cp->getVelo().x()<<", \"vy\":"<<cp->getVelo().y()<<", \"vz\":"<<cp->getVelo().z()<<", \"type\":"<<cp->getType()<<"}";
                out << ",\n";
            } else {
               printf("Error: NaN encountered when trying to export state of particle: %u\n", cp->getID());
               exit(EXIT_FAILURE);
            }
        }
    }
    
    // remove the last comma
    //out << "\b\b";
    //out << "  \n";
    out.seekp(-2,std::ios::end);
    out << "  \n";
    
    out << "]}\n";
    out.close();
    return;
}

// emits the state of the system as a JSON file
void SimSystem::emitJSON(std::string jsonfile) {
    std::ofstream out;
    out.open(jsonfile);
    out <<"{\n";
    out << "  \"beads\":[\n";

    // iterate through the particles and write the JSON 
    for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
        Particle *p = *i;
        // check to make sure the particle position makes sense
        if (!isnan(p->getPos().x()) && !isnan(p->getPos().y()) && !isnan(p->getPos().z())) {
            out << "\t{\"id\":"<<p->getID()<<", \"x\":"<<p->getPos().x()<<", \"y\":"<<p->getPos().y()<<", \"z\":"<<p->getPos().z()<<", \"vx\":"<<p->getVelo().x()<<", \"vy\":"<<p->getVelo().y()<<", \"vz\":"<<p->getVelo().z()<<", \"type\":"<<p->getType()<<"}";
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
   
    float cube_size = _N/_D;
    spatial_unit_address_t su;
    su.x = ceil(p->getPos().x()/cube_size) - 1;
    su.y = ceil(p->getPos().y()/cube_size) - 1;
    su.z = ceil(p->getPos().z()/cube_size) - 1;
 
    Vector3D relative_pos;
    relative_pos.x(p->getPos().x() - (su.x*cube_size));
    relative_pos.y(p->getPos().y() - (su.y*cube_size));
    relative_pos.z(p->getPos().z() - (su.z*cube_size));

    // make the particle have a position relative to it's spatial unit
    p->setPos(relative_pos);

    SpatialUnit *sup = getSpatialUnit(su);

    sup->addLocalParticle(p);

}

//! prints the number of particles allocated per spatial unit
void SimSystem::printSpatialAllocation() {
   for(SimSystem::iterator i=begin(), e=end(); i!=e; ++i){
       SpatialUnit *cur = *i;
       spatial_unit_address_t pos =  cur->getAddr();
       std::cout << "Number of particles for this spatial unit ("<<pos.x << "," <<pos.y<<","<<pos.z<<"): " << cur->numBeads() << "\n"; 
   }
}

// runs the simulation for a given period uses the SpatialUnits and passing beads between them to simulate the MPI
// based approach
void SimSystem::run(uint32_t period, float emitrate) {

 clock_t last_emit = clock();
 unsigned emit_cnt=0;

 for(uint32_t t=0; t<period; t++) {
     // for each spatial unit create their local view of the world based on their neighbours states
     for(iterator i=begin(); i!=end(); ++i) {
         SpatialUnit *s = *i; // the current spatial unit

         // iterate over itself and apply the forces to it's internal particles
         for(SpatialUnit::iterator csu_p1=s->begin(); csu_p1 !=s->end(); ++csu_p1){
             Particle *p1 = *csu_p1;
             for(SpatialUnit::iterator csu_p2=s->begin(); csu_p2!=s->end(); ++csu_p2){
               Particle *p2 = *csu_p2;
               if(p1->getID() != p2->getID()) {
                  if(p1->getPos().dist(p2->getPos()) <= _r_c){
                      // do the force update
                      p1->callConservative(p2);
                      p1->callDrag(p2);
                      p1->callRandom(p2);

                      Particle *bond1 = p1->getBondedParticle();
                      Particle *bond2 = p2->getBondedParticle();
                      if( (bond1 == p2) || (bond2 == p1) ) { // these particles are bonded 
                         p1->callBond();
                      }  
                   }
               } 
             }
          }


         // iterate over all neighbours and solve the forces for all the particles local to this spatial unit          
         for(SpatialUnit::n_iterator j=s->n_begin(); j!=s->n_end(); ++j) {
            SpatialUnit *neighbour = *j;

            // calculate the relative position offsets for this neighbour
            spatial_unit_address_t n_addr = neighbour->getAddr();
            spatial_unit_address_t this_addr = s->getAddr(); // we need wraparound for this

            // get relative positions from the abs spatial unit addresses
            int x_rel = n_addr.x - this_addr.x;
            int y_rel = n_addr.y - this_addr.y;
            int z_rel = n_addr.z - this_addr.z; 
            
            // periodic boundary adjusting
            if(x_rel > 1)
                  x_rel = -1;
            else if (x_rel < -1)
                  x_rel = 1;
            
            if(y_rel > 1)
                  y_rel = -1;
            else if (y_rel < -1)
                  y_rel = 1;

            if(z_rel > 1)
                  z_rel = -1;
            else if (z_rel < -1)
                  z_rel = 1;
                   

            float x_off = (float)(x_rel)*neighbour->getSize();
            float y_off = (float)(y_rel)*neighbour->getSize();
            float z_off = (float)(z_rel)*neighbour->getSize();

            // loop over all the particles in this neighbour and apply the forces to our particles
            for(SpatialUnit::iterator np=neighbour->begin(); np!=neighbour->end(); ++np){
              // move the position of this particle to be relative to our own (we will need to restore it once we are done) 
               Particle *fp = *np;
               Vector3D fp_pos = fp->getPos();
               Vector3D n_fp_pos(fp_pos.x() + x_off, fp_pos.y() + y_off, fp_pos.z() + z_off);
               fp->setPos(n_fp_pos); 
                
               // loop over all our particles and apply the forces and update
               for(SpatialUnit::iterator this_pi = s->begin(); this_pi != s->end(); ++this_pi){
                   Particle *this_p = *this_pi; // apply the forces to our particle
                   if(this_p->getPos().dist(fp->getPos()) <= _r_c) {
 
                       // these particles are in range of each other
                       // do the force update
                       this_p->callConservative(fp);
                       this_p->callDrag(fp);
                       this_p->callRandom(fp);

                       Particle *bond1 = fp->getBondedParticle();
                       Particle *bond2 = this_p->getBondedParticle();
                       if( (bond1 == this_p) || (bond2 == fp) ) { // these particles are bonded 
                          this_p->callBond();
                       }  

                  } 

              }
              fp->setPos(fp_pos); // restore the position
            }
         }
     }

     // we are done updating all the forces, now we can start updating the positions and doing particle migration
     // update v_i and r_i for each particle
     // iterate over all spatial units
     std::vector<Particle *> to_migrate; // list of particles we want to migrate
     std::vector<SpatialUnit *> src_migrate; // the source spatial unit for each migration
     std::vector<SpatialUnit *> dest_migrate; // the destination spatial unit for each migration
     for(iterator sui=begin(); sui!=end(); ++sui) {
         SpatialUnit *su = *sui;
         std::vector<Particle *> tmp_su_particles = su->copyOfParticles();
         for(SpatialUnit::iterator i=tmp_su_particles.begin(), ie=tmp_su_particles.end(); i!=ie; ++i) {
              Particle *p = *i;
              float mass = p->getMass();
              Vector3D acceleration = p->getForce()/mass;
              Vector3D delta_v = (p->getForce()/mass) * _dt;
              // update velocity
              p->setVelo(p->getVelo() + delta_v); 
              //p->setVelo(Vector3D(0.1,0.1,0.1)); 

              // euler update position & include wraparound
              //Vector3D point = p->getPos() +p->getVelo()*_dt;

              // velocity verlet
              Vector3D point = p->getPos() + p->getVelo()*_dt + acceleration*0.5*_dt*_dt; 

              // here we need to check for migration
              spatial_unit_address_t dest; // the spatial unit where we might be potentially sending this particle
              spatial_unit_address_t curr_addr = su->getAddr();
              bool migrating = false;

              // migration in the x direction
              if(point.x() >= su->getSize()) {
                 migrating = true;
                 if(curr_addr.x == (_D - 1))
                     dest.x = 0;
                 else
                     dest.x = curr_addr.x + 1; 
                 point.x(point.x() - su->getSize());
              } else if (point.x() < 0.0) {
                 migrating = true;
                 if(curr_addr.x == 0)
                     dest.x = _D - 1;
                 else
                     dest.x = curr_addr.x - 1; 
                 point.x(point.x() + su->getSize());
              } else {
                 dest.x = curr_addr.x; 
              } 

              // migration in the y direction
              if(point.y() >= su->getSize()) {
                 migrating = true;
                 if(curr_addr.y == (_D - 1))
                     dest.y = 0;
                 else
                     dest.y = curr_addr.y + 1; 
                 point.y(point.y() - su->getSize());
              } else if (point.y() < 0.0) {
                 migrating = true;
                 if(curr_addr.y == 0)
                     dest.y = _D - 1;
                 else
                     dest.y = curr_addr.y - 1; 
                 point.y(point.y() + su->getSize());
              } else {
                 dest.y = curr_addr.y; 
              } 

              // migration in the z direction
              if(point.z() >= su->getSize()) {
                 migrating = true;
                 if(curr_addr.z == (_D - 1))
                     dest.z = 0;
                 else
                     dest.z = curr_addr.z + 1; 
                 point.z(point.z() - su->getSize());
              } else if (point.z() < 0.0) {
                 migrating = true;
                 if(curr_addr.z == 0)
                     dest.z = _D - 1;
                 else
                     dest.z = curr_addr.z - 1; 
                 point.z(point.z() + su->getSize());
              } else {
                 dest.z = curr_addr.z; 
              } 

              // do we actually need to migrate?
              if(migrating) {
                 SpatialUnit *dest_su = getSpatialUnit(dest);
                 to_migrate.push_back(p);            
                 src_migrate.push_back(su);
                 dest_migrate.push_back(dest_su);
              }

              // update the position
              p->setPos(point); 

              // clear force 
              p->setForce(Vector3D(0.0, 0.0, 0.0));
         }
     }

     // perform all the migrations
     for(unsigned i=0; i<to_migrate.size(); i++){
          Particle *p = to_migrate[i];
          SpatialUnit *src = src_migrate[i];
          SpatialUnit *dst = dest_migrate[i];

          dst->addLocalParticle(p); 
          if(!src->removeParticle(p)) {
            printf("\n\nERROR: unable to remove a particle from a spatial unit\n"); 
            spatial_unit_address_t su_addr = src->getAddr();
            printf("particle being removed p:%d pos:%s velo:%s\n", p->getID(), p->getPos().str().c_str(), p->getVelo().str().c_str()); 
            exit(EXIT_FAILURE);
          }
     }

     // clear the migration lists 
     to_migrate.clear();
     src_migrate.clear();
     dest_migrate.clear();

     // do we want to emit the state of the simulation
     if ((float(clock() - last_emit) / CLOCKS_PER_SEC) > emitrate) {
         emit_cnt++;

         emitJSONFromSU("state.json"); // also replace the latest frame
         printf("u\n"); // update command set via stdout to nodejs server
         fflush(stdout);
         last_emit = clock();
     }

     // update time
     _ts = _ts + 1;
     _t = _t + _dt;

   } // period has elapsed

}


// runs the simulation for a given number of timesteps sequentially there is no parallelism here
// emits the state of each particle given by the emitrate
// this is the naive implementation of this algorithm (horrifically slow)
void SimSystem::seq_run(uint32_t period, float emitrate) {

    unsigned start_ts = _ts;
    clock_t last_emit = clock();
    unsigned emit_cnt=0;
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

                            Particle *bond1 = p1->getBondedParticle();
                            Particle *bond2 = p2->getBondedParticle();
                            if( (bond1 == p2) || (bond2 == p1) ) { // these particles are bonded 
                               p1->callBond();
                            }  

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


             // sanity assertion check for bonded particles
             //if(p->isBonded()) {
             //   if(p->getPos().toroidal_dist(p->getBondedParticle()->getPos(), _N) > _r_c){
             //        Particle *i = p;
             //        Particle *j = p->getBondedParticle();
             //        std::cerr << "Error! bonded beads: "<< i->getID()<<" <-> "<< j->getID()<<"  have a broken bond\n";
             //        std::cerr << "Pos of Bead: "<<i->getID() << " is "<< i->getPos().str() <<"\n"; 
             //        std::cerr << "Pos of Bead: "<<j->getID() << " is "<< j->getPos().str() <<"\n"; 
             //        std::cerr << "Previous pos of Bead: "<<i->getID() << " is "<< i->getPrevPos().str() <<"\n"; 
             //        std::cerr << "Previous pos of Bead: "<<j->getID() << " is "<< j->getPrevPos().str() <<"\n"; 
             //        std::cerr << "Velocity of Bead: "<<i->getID() << " is "<< i->getVelo().str() <<"\n"; 
             //        std::cerr << "Velocity of Bead: "<<j->getID() << " is "<< j->getVelo().str() <<"\n"; 
             //        std::cerr << "Distance between Beads: "<<i->getPos().toroidal_dist(j->getPos(), _N) <<"\n"; 
             //        std::cerr << "Force of beads: Bead = "<<i->getID()<< " (Force = " << i->getForce().str() << ")   Bead = "<< j->getID() <<" (Force = " << j->getForce().str() <<")\n";
             //        exit(EXIT_FAILURE);
             //   }
             //}             

            // clear force 
            p->setForce(Vector3D(0.0, 0.0, 0.0));

        } 


        // cleanup: 
        _seq_pairs->clear(); // we no longer need to track the pairs

        // do we want to emit the state of the simulation
        if ((float(clock() - last_emit) / CLOCKS_PER_SEC) > emitrate) {
            //std::stringstream curframe;
            //curframe << "frames/state_" << emit_cnt << ".json"; 
            //emitJSON(curframe.str()); //emit the frame into the frames folder for later retrieval
            emit_cnt++;

            emitJSON("state.json"); // also replace the latest frame
            printf("u\n"); // update command set via stdout to nodejs server
            fflush(stdout);
            last_emit = clock();
        }


        // update time
        _ts = _ts + 1;
        _t = _t + _dt;
    }
}

// adds a particle to the system with global coordinates
void SimSystem::addParticle(Particle *p){
    // add the particle to the global particles list
    p->setID(_particles->size()); // just use the current number of particles as a unique ID
    _particles->push_back(p);     
    allocateParticleToSpatialUnit(p); // allocate particle to a processor
}

// adds a particle to the system with local coordinates
void SimSystem::addLocalParticle(Particle *p){
    // add the particle to the global particles list
    p->setID(_particles->size()); // just use the current number of particles as a unique ID
    _particles->push_back(p);     
}


// populates the universe with particles from a JSON file
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

