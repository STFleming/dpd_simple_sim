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
   //for(iterator i=begin(), ie=end(); i!=ie; ++i) {
   //    Particle *p = *i;
   //    delete p;
   //}
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
    if(s==NULL) {
         printf("ERROR: trying to assign a NULL as a neighbour to a spatial unit!\n");
         exit(EXIT_FAILURE);
    }
    _neighbours->push_back(s); // it is the responsibility of the universe to make sure the neighbour makes sense
}

/**! returns the number of neighbours this spatial unit has */
unsigned SpatialUnit::numNeighbours(){ return _neighbours->size(); }

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

/**! adds a local particle to a spatial unit */
void SpatialUnit::addLocalParticle(Particle *p){
     // check to make sure the position is local  
     Vector3D p_pos = p->getPos();

     if( ( p_pos.x() < 0 ) || (p_pos.x() > _size) ){
         printf("Error adding particle p(%d): the x dimension for this particle is not local:%s\n", p->getID(), p_pos.str().c_str());
         fflush(stdout);
         exit(EXIT_FAILURE);
     }  
    
     if( ( p_pos.y() < 0 ) || (p_pos.y() > _size) ){
         printf("Error adding particle p(%d): the y dimension for this particle is not local pos:%s\n", p->getID(), p_pos.str().c_str());
         fflush(stdout);
         exit(EXIT_FAILURE);
     }  
     
     if( ( p_pos.z() < 0 ) || (p_pos.z() > _size) ){
         printf("Error adding particle p(%d): the z dimension for this particle is not local pos:%s\n", p->getID(), p_pos.str().c_str());
         fflush(stdout);
         exit(EXIT_FAILURE);
     }  
     _particles->push_back(p);
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

// returns a copy of the particles list (so things can be iterated over and removed)
std::vector<Particle *> SpatialUnit::copyOfParticles(){
    std::vector<Particle *> t;
    for(SpatialUnit::iterator i=begin(); i!=end(); ++i){
         Particle *p = *i;
         t.push_back(p);
    }
    return t;
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
