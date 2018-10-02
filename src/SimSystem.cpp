#include "SimSystem.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SimSystem::SimSystem(float N, unsigned D, unsigned verbosity) {
    _N = N;
    _D = D;
    _verbosity = verbosity;
    
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

// emits the state of the system as a JSON file
void SimSystem::emitJSON() {
    std::ofstream out;
    out.open("state.json");
    out <<"{\n";
    out << "  \"particles\":[\n";

    // iterate through the particles and write the JSON 
    for(p_iterator i=p_begin(), ie=p_end(); i!=ie; ++i) {
        Particle *p = *i;
        out << "\t{\"id\":\"p_"<<p->getID()<<"\", \"x\":"<<p->getPos().x<<", \"y\":"<<p->getPos().y<<", \"z\":"<<p->getPos().z<<"}";
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

// adds a particle to the system
void SimSystem::addParticle(Particle *p){
    // add the particle to the global particles list
    p->setID(_particles->size()); // just use the current number of particles as a unique ID
    _particles->push_back(p);     

    // determine which SpatialUnit should host this particle
    for(iterator i=begin(), ie=end(); i!=ie; ++i){
       SpatialUnit *cur = *i;
       if(cur->checkPos(p->getPos())) {
          cur->addParticle(p);
          return;
       } 
    }    
    std::runtime_error("Could not find a SpatialUnit that could support the particle.\n");
}
