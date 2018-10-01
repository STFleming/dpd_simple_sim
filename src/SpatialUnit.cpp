#include "SpatialUnit.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SpatialUnit::SpatialUnit(unsigned size, uint16_t x, uint16_t y, uint16_t z, unsigned verbosity) {
    _size = size;
    _center.x = x;
    _center.y = y;
    _center.z = z;
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
    if( (p.x >= _center.x + _size/2) || (p.x <= _center.x - _size/2) )
        return false;
    if( (p.y >= _center.y + _size/2) || (p.y <= _center.y - _size/2) )
        return false;
    if( (p.z >= _center.z + _size/2) || (p.z <= _center.z - _size/2) )
        return false;
    return true;
} 

/**! add a particle to this spatial unit at position p */
void SpatialUnit::addParticle(Particle* p) {
   // first check to make sure that we can actually add a particle
   if(checkPos(p->getPos())) {
       _particles->push_back(p); // add the particle 
       if(_verbosity >= 2)
           printf("    spatial unit center=(x:%d,y:%d,z:%d) is hosting a particle at pos=(x:%d, y:%d, z:%d)\n", _center.x, _center.y, _center.z, p->getPos().x, p->getPos().y, p->getPos().z);
   } else {
     std::runtime_error("Error: this particle cannot fit in within this cube!\n"); 
   } 
}

/**! returns the position of this cube */
position_t SpatialUnit::getPos() {
    return _center;
}

/**! returns the size of this cube */
unsigned SpatialUnit::getSize() {
    return _size;
}
