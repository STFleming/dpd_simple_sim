#include "SimSystem.hpp"

#ifndef _SIMSYSTEM_IMPL
#define _SIMSYSTEM_IMPL

/**! Constructor: creates a problem of NxNxN cubes */
template <class S>
SimSystem<S>::SimSystem(S N, S dt, S r_c, unsigned D, unsigned verbosity) {
    _N = N;
    _D = D;
    _verbosity = verbosity;
    _t = 0.0; // simulation starts at t=0.0
    _ts = 0;
    _dt = dt; 
    _r_c = r_c; // interaction cutoff 
    
    // create the global list of particles
    _particles = new std::vector<Particle<S> *>();

    if(_N == 0) { // cannot have a problem with no size
       std::runtime_error("Problem must have a size >0\n");
    }

    // calculate the size of each spatial unit
    _unit_size = _N / (float)D;

    // create the list of _cubes
    _cubes = new std::vector<SpatialUnit<S> *>(); 
    for(int x=0; x<_D; x++) {
        for(int y=0; y<_D; y++) {
            for(int z=0; z<_D; z++) {
                S fpx = (float)x * _unit_size;
                S fpy = (float)y * _unit_size;
                S fpz = (float)z * _unit_size;

                // create a SpatialUnit and add it to the cubes
                if(_verbosity >= 3)
                    printf("Constructing a cube (size=%f) at x:%f, y:%f, z:%f\n", _unit_size, fpx, fpy, fpz);

                spatial_unit_address_t c_addr = {x,y,z};
                _cubes->push_back(new SpatialUnit<S>(_unit_size,fpx,fpy,fpz,c_addr, _verbosity)); // Add a new cube of size _N/_D at x,y,z
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
                 SpatialUnit<S> *curr_spu = getSpatialUnit(curr); 
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
template <class S>
SpatialUnit<S> * SimSystem<S>::getSpatialUnit(spatial_unit_address_t x) {
    for(SimSystem<S>::iterator i=begin(); i!=end(); ++i){
        SpatialUnit<S> *c = *i;
        spatial_unit_address_t c_addr = c->getAddr();
        if( (x.x == c_addr.x) && (x.y == c_addr.y) && (x.z == c_addr.z)) {
           return c;
        } 
    }
    return NULL;
}

// bonds particle i to particle j in the universe 
template<class S>
void SimSystem<S>::bond(Particle<S> *i, Particle<S> *j, std::function<void(Particle<S> *me, Particle<S> *other)> bondf){
    i->setOutBond(j, bondf);
    j->setInBond(i, bondf);
    return;
}

// emits the state of the system from the spatial units (relative addressing)
template<class S>
void SimSystem<S>::emitJSONFromSU(std::string jsonfile) {
    std::ofstream out;
    out.open(jsonfile);
    out <<"{\n";
    out << "  \"beads\":[\n";

    // iterate through the spatial units and grab the offsets
    for(iterator i=begin(); i!=end(); ++i){
        SpatialUnit<S> *s = *i;
        spatial_unit_address_t pos = s->getAddr(); 
    
        float x_off = pos.x * s->getSize(); 
        float y_off = pos.y * s->getSize(); 
        float z_off = pos.z * s->getSize(); 

        // iterate through all particles of this spatial unit
        for(auto ip=s->begin(); ip!=s->end(); ++ip){
           Particle<S> *cp = *ip;
           float cp_x = cp->getPos().x() + x_off;
           float cp_y = cp->getPos().y() + y_off;
           float cp_z = cp->getPos().z() + z_off;

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
template<class S>
void SimSystem<S>::emitJSON(std::string jsonfile) {
    std::ofstream out;
    out.open(jsonfile);
    out <<"{\n";
    out << "  \"beads\":[\n";

    // iterate through the particles and write the JSON 
    for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
        Particle<S> *p = *i;
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
template<class S>
typename SimSystem<S>::iterator SimSystem<S>::begin() { return _cubes->begin(); }
template<class S>
typename SimSystem<S>::iterator SimSystem<S>::end() { return _cubes->end(); }
template<class S>
typename SimSystem<S>::p_iterator SimSystem<S>::p_begin() { return _particles->begin(); }
template<class S>
typename SimSystem<S>::p_iterator SimSystem<S>::p_end() { return _particles->end(); }

/**! Destructor: cleans up the _cubes */
template<class S>
SimSystem<S>::~SimSystem() {
    for(iterator i=begin(), e=end(); i!=e; ++i){
        SpatialUnit<S> *cur = *i;
        if(_verbosity>=4)
            printf("removing cube (size=%f) at x:%f, y:%f, z:%f\n", cur->getSize(), cur->getPos().x, cur->getPos().y, cur->getPos().z);
        delete cur;
    }  
    delete _cubes;
    delete _particles; // this should probably iterate through and delete them all not the spatial units...
}

// allocates a particles to spatial processing units
template<class S>
void SimSystem<S>::allocateParticleToSpatialUnit(Particle<S> *p) {

    S cube_size = _N/_D;
    spatial_unit_address_t su;
    su.x = floor(p->getPos().x()/cube_size);
    su.y = floor(p->getPos().y()/cube_size);
    su.z = floor(p->getPos().z()/cube_size);
 
    Vector3D relative_pos;
    relative_pos.x(p->getPos().x() - (su.x*cube_size));
    relative_pos.y(p->getPos().y() - (su.y*cube_size));
    relative_pos.z(p->getPos().z() - (su.z*cube_size));

    // make the particle have a position relative to it's spatial unit
    p->setPos(relative_pos);
    p->setPrevPos(relative_pos);

    SpatialUnit<S> *sup = getSpatialUnit(su);

    sup->addLocalParticle(p);

}

//! prints the number of particles allocated per spatial unit
template<class S>
void SimSystem<S>::printSpatialAllocation() {
   for(SimSystem<S>::iterator i=begin(), e=end(); i!=e; ++i){
       SpatialUnit<S> *cur = *i;
       spatial_unit_address_t pos =  cur->getAddr();
       std::cout << "Number of particles for this spatial unit ("<<pos.x << "," <<pos.y<<","<<pos.z<<"): " << cur->numBeads() << "\n"; 
   }
}

// runs the simulation for a given period uses the SpatialUnits and passing beads between them to simulate the MPI
// based approach
template<class S>
void SimSystem<S>::run(uint32_t period, float emitrate) {

 clock_t last_emit = clock();
 unsigned emit_cnt=0;
 _grand = rand(); // get a new global random number 

 for(uint32_t t=0; t<period; t++) {

     // for each spatial unit create their local view of the world based on their neighbours states
     for(iterator i=begin(); i!=end(); ++i) {
         SpatialUnit<S> *s = *i; // the current spatial unit

         // iterate over itself and apply the forces to it's internal particles
         for(auto csu_p1=s->begin(); csu_p1 !=s->end(); ++csu_p1){
             Particle<S> *p1 = *csu_p1;
             for(auto csu_p2=s->begin(); csu_p2!=s->end(); ++csu_p2){
               Particle<S> *p2 = *csu_p2;
               if(p1->getID() != p2->getID()) {
                  if(p1->getPos().dist(p2->getPos()) <= _r_c){
                      // do the force update
                      p1->callConservative(p2);
                      p1->callDrag(p2);
                      p1->callRandom(_grand, p2);

                      // check to see if this bead is the Out Bond
                      Particle<S> *outBond = p1->getOutBondBead();
                      if(outBond == p2) {
                        p1->callBond();
                      }                       
 
                      // check to see we are the In Bond bead for the foreign particle fp
                      Particle<S> *inBond = p1->getInBondBead();
                      if(inBond == p2){
                         p2->callInverseBond();
                      }

                   } 
               } 
             }
          }


         // iterate over all neighbours and solve the forces for all the particles local to this spatial unit          
         for(auto j=s->n_begin(); j!=s->n_end(); ++j) {
            SpatialUnit<S> *neighbour = *j;

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
                   

            S x_off = (float)(x_rel)*neighbour->getSize();
            S y_off = (float)(y_rel)*neighbour->getSize();
            S z_off = (float)(z_rel)*neighbour->getSize();

            // loop over all the particles in this neighbour and apply the forces to our particles
            for(auto np=neighbour->begin(); np!=neighbour->end(); ++np){
              // move the position of this particle to be relative to our own (we will need to restore it once we are done) 
               Particle<S> *fp = *np;
               Vector3D fp_pos = fp->getPos();
               Vector3D n_fp_pos(fp_pos.x() + x_off, fp_pos.y() + y_off, fp_pos.z() + z_off);
               fp->setPos(n_fp_pos); 
 
               // also do the same thing for it's previous position
               Vector3D fp_prev_pos = fp->getPrevPos();
               Vector3D n_pfp_pos(fp_prev_pos.x() + x_off, fp_prev_pos.y() + y_off, fp_prev_pos.z() + z_off);
               fp->setPrevPos(n_pfp_pos); 
               
                
               // loop over all our particles and apply the forces and update
               for(auto this_pi = s->begin(); this_pi != s->end(); ++this_pi){
                   Particle<S> *this_p = *this_pi; // apply the forces to our particle
                   if(this_p->getPos().dist(fp->getPos()) <= _r_c) {
 
                       // these particles are in range of each other
                       // do the force update
                       this_p->callConservative(fp);
                       this_p->callDrag(fp);
                       this_p->callRandom(_grand, fp);
                        
                       // check to see if this bead is the Out Bond
                       Particle<S> *outBond = this_p->getOutBondBead();
                       if(outBond == fp) {
                         this_p->callBond();
                       }                       
 
                       // check to see we are the In Bond bead for the foreign particle fp
                       Particle<S> *inBond = this_p->getInBondBead();
                       if(inBond == fp){
                          fp->callInverseBond();
                       }

                  } else { 
                  }
                  
              }
              fp->setPos(fp_pos); // restore the position
              fp->setPrevPos(fp_prev_pos); // restore the position
            }
         }
     }

     // we are done updating all the forces, now we can start updating the positions and doing particle migration
     // update v_i and r_i for each particle
     // iterate over all spatial units
     std::vector<Particle<S> *> to_migrate; // list of particles we want to migrate
     std::vector<SpatialUnit<S> *> src_migrate; // the source spatial unit for each migration
     std::vector<SpatialUnit<S> *> dest_migrate; // the destination spatial unit for each migration
     for(iterator sui=begin(); sui!=end(); ++sui) {
         SpatialUnit<S> *su = *sui;
         std::vector<Particle<S> *> tmp_su_particles = su->copyOfParticles();
         for(auto i=tmp_su_particles.begin(), ie=tmp_su_particles.end(); i!=ie; ++i) {
              Particle<S> *p = *i;
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
                 SpatialUnit<S> *dest_su = getSpatialUnit(dest);
                 to_migrate.push_back(p);            
                 src_migrate.push_back(su);
                 dest_migrate.push_back(dest_su);
              }

              // update the position
              p->setPrevPos(p->getPos());
              p->setPos(point); 

              // clear force 
              p->setForce(Vector3D(0.0, 0.0, 0.0));
         }
     }

     // perform all the migrations
     for(unsigned i=0; i<to_migrate.size(); i++){
          Particle<S> *p = to_migrate[i];
          SpatialUnit<S> *src = src_migrate[i];
          SpatialUnit<S> *dst = dest_migrate[i];

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
     _grand = rand();

   } // period has elapsed

}


// runs the simulation for a given number of timesteps sequentially there is no parallelism here
// emits the state of each particle given by the emitrate
// this is the naive implementation of this algorithm (horrifically slow)
template<class S>
void SimSystem<S>::seq_run(uint32_t period, float emitrate) {

}

// adds a particle to the system with global coordinates
template<class S>
void SimSystem<S>::addParticle(Particle<S> *p){
    // add the particle to the global particles list
    p->setID(_particles->size()); // just use the current number of particles as a unique ID
    _particles->push_back(p);     
    allocateParticleToSpatialUnit(p); // allocate particle to a processor
}

// adds a particle to the system with local coordinates
template<class S>
void SimSystem<S>::addLocalParticle(Particle<S> *p){
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

#endif /* _SIMSYSTEM_IMPL */
