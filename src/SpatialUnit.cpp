#include "SpatialUnit.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SpatialUnit::SpatialUnit(unsigned size, uint16_t x, uint16_t y, uint16_t z) {
    _size = size;
    _center.x = x;
    _center.y = y;
    _center.z = z;
}

/**! Destructor: cleans up the particles */
SpatialUnit::~SpatialUnit() {

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
