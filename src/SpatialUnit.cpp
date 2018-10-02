#include "SpatialUnit.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SpatialUnit::SpatialUnit(float size, float x, float y, float z, unsigned verbosity) {
    _size = size;
    _pos.x = x;
    _pos.y = y;
    _pos.z = z;
    _verbosity = verbosity;
    
    // create the heap vector of particles
    _particles = new std::vector<Particle *>();
}

/**! Destructor: cleans up the particles */
SpatialUnit::~SpatialUnit() {
   // TODO iterate through all the particles and delete them
   // TODO delete the vector of particles

}

/**! returns true if this position is within this spatial unit */
bool SpatialUnit::checkPos(position_t p) {
    if( (p.x < _pos.x) || (p.x >= _pos.x + _size) )
        return false;
    if( (p.y < _pos.y) || (p.y >= _pos.y + _size) )
        return false;
    if( (p.z < _pos.z) || (p.z >= _pos.z + _size) )
        return false;
    return true;
} 

/**! add a particle to this spatial unit at position p */
void SpatialUnit::addParticle(Particle* p) {
   // first check to make sure that we can actually add a particle
   if(checkPos(p->getPos())) {
       _particles->push_back(p); // add the particle 
       if(_verbosity >= 2)
           printf("    spatial unit center=(x:%f,y:%f,z:%f) is hosting a particle at pos=(x:%f, y:%f, z:%f)\n", _pos.x, _pos.y, _pos.z, p->getPos().x, p->getPos().y, p->getPos().z);
   } else {
     std::runtime_error("Error: this particle cannot fit in within this cube!\n"); 
   } 
}

/**! returns the position of this cube */
position_t SpatialUnit::getPos() {
    return _pos;
}

/**! returns the size of this cube */
float SpatialUnit::getSize() {
    return _size;
}

/**! iterators for the start and end of the particle vector */
SpatialUnit::iterator SpatialUnit::begin(){ return _particles->begin(); }
SpatialUnit::iterator SpatialUnit::end(){ return _particles->end(); }
