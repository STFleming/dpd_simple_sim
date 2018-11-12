#include "SpatialUnit.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SpatialUnit::SpatialUnit(float size, float x, float y, float z, spatial_unit_address_t addr, unsigned verbosity) {
    _size = size;
    _pos.x = x;
    _pos.y = y;
    _pos.z = z;
    _verbosity = verbosity;
    _addr = addr;
    
    // create the heap vector of particles
    _particles = new std::vector<Particle *>();
    // create the heap vector of neighbours
    _neighbours = new std::vector<SpatialUnit *>();
}

/**! Destructor: cleans up the particles */
SpatialUnit::~SpatialUnit() {
   for(iterator i=begin(), ie=end(); i!=ie; ++i) {
       Particle *p = *i;
       delete p;
   }
   delete _particles;
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

/**! returns the address for this spatial unit */
spatial_unit_address_t SpatialUnit::getAddr() { return _addr; }

/**! adds a neighbour to the neighbour list*/
void SpatialUnit::addNeighbour(SpatialUnit *s){
    _neighbours->push_back(s); // it is the responsibility of the universe to make sure the neighbour makes sense
}

/**! add a particle to this spatial unit at position p */
void SpatialUnit::addParticle(Particle* p) {
   // first check to make sure that we can actually add a particle
   if(checkPos(vec2pos(p->getPos()))) {
       _particles->push_back(p); // add the particle 
       if(_verbosity >= 2)
           printf("    spatial unit center=(x:%f,y:%f,z:%f) is hosting a particle at pos=(x:%f, y:%f, z:%f)\n", _pos.x, _pos.y, _pos.z, p->getPos().x(), p->getPos().y(), p->getPos().z());
   } else {
     std::runtime_error("Error: this particle cannot fit in within this cube!\n"); 
   } 
}

/**! removes a particle from this spatial unit 
     returns true if the removal was successful
*/
bool SpatialUnit::removeParticle(Particle *p){
  for(SpatialUnit::iterator i=begin(); i!=end(); ++i){
      Particle *curr = *i;
      if(p->getID() == curr->getID()){
          // we've found our particle
          _particles->erase(i);
          return true;
      } 
  } 
  return false;
}

/**! returns the number of beads within this spatial unit */
unsigned SpatialUnit::numBeads() {
    return _particles->size();
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

/**! iterators for the start and end of the neighbours vector */
SpatialUnit::n_iterator SpatialUnit::n_begin(){ return _neighbours->begin(); }
SpatialUnit::n_iterator SpatialUnit::n_end(){ return _neighbours->end(); }
