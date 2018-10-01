#include "SimSystem.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SimSystem::SimSystem(unsigned N, unsigned D, unsigned verbosity) {
    _N = N;
    _D = D;
    _verbosity = verbosity;

    // _N must be an integer multiple of _D
    if(_N % _D != 0) {
        std::runtime_error("Problem must be discretised into an integer number of chunks i.e. N \% D == 0\n"); 
    }

    if(_N == 0) { // cannot have a problem with no size
       std::runtime_error("Problem must have a size >0\n");
    }

    // create the list of _cubes
    _cubes = new std::vector<SpatialUnit *>(); 
    for(int x=0; x<_D; x++) {
        for(int y=0; y<_D; y++) {
            for(int z=0; z<_D; z++) {
                // create a SpatialUnit and add it to the cubes
                if(_verbosity >= 3)
                    printf("Constructing a cube (size=%d) at x:%d, y:%d, z:%d\n", (_N/_D), x, y, z);
                _cubes->push_back(new SpatialUnit((_N/_D), x,y,z, _verbosity)); // Add a new cube of size _N/_D at x,y,z
            }
        }
    }
}

// iterators
SimSystem::iterator SimSystem::begin() { return _cubes->begin(); }
SimSystem::iterator SimSystem::end() { return _cubes->end(); }

/**! Destructor: cleans up the _cubes */
SimSystem::~SimSystem() {
    for(iterator i=begin(), e=end(); i!=e; ++i){
        SpatialUnit *cur = *i;
        if(_verbosity>=4)
            printf("removing cube (size=%d) at x:%d, y:%d, z:%d\n", cur->getSize(), cur->getPos().x, cur->getPos().y, cur->getPos().z);
        delete cur;
    }  
    delete _cubes;
}

// adds a particle to the system
void SimSystem::addParticle(Particle *p){
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
