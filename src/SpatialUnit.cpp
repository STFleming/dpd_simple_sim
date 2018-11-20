#include "SpatialUnit.hpp"

#ifndef _SPATIAL_UNIT_IMPL
#define _SPATIAL_UNIT_IMPL

/**! Constructor: creates a problem of NxNxN cubes */
template<class S>
SpatialUnit<S>::SpatialUnit(S size, S x, S y, S z, spatial_unit_address_t addr, unsigned verbosity) {
    _size = size;
    _pos.x = x;
    _pos.y = y;
    _pos.z = z;
    _verbosity = verbosity;
    _addr = addr;
    
    // create the heap vector of particles
    _particles = new std::vector<Particle<S> *>();
    // create the heap vector of neighbours
    _neighbours = new std::vector<SpatialUnit<S> *>();
}

/**! Destructor: cleans up the particles */
template<class S>
SpatialUnit<S>::~SpatialUnit() {
   //for(iterator i=begin(), ie=end(); i!=ie; ++i) {
   //    Particle<S> *p = *i;
   //    delete p;
   //}
   delete _particles;
}

/**! returns true if this position is within this spatial unit */
template<class S>
bool SpatialUnit<S>::checkPos(position_t<S> p) {
    if( (p.x < _pos.x) || (p.x >= _pos.x + _size) )
        return false;
    if( (p.y < _pos.y) || (p.y >= _pos.y + _size) )
        return false;
    if( (p.z < _pos.z) || (p.z >= _pos.z + _size) )
        return false;
    return true;
} 

/**! returns the address for this spatial unit */
template<class S>
spatial_unit_address_t SpatialUnit<S>::getAddr() { return _addr; }

/**! adds a neighbour to the neighbour list*/
template<class S>
void SpatialUnit<S>::addNeighbour(SpatialUnit<S> *s){
    if(s==NULL) {
         printf("ERROR: trying to assign a NULL as a neighbour to a spatial unit!\n");
         exit(EXIT_FAILURE);
    }
    _neighbours->push_back(s); // it is the responsibility of the universe to make sure the neighbour makes sense
}

/**! returns the number of neighbours this spatial unit has */
template<class S>
unsigned SpatialUnit<S>::numNeighbours(){ return _neighbours->size(); }

/**! add a particle to this spatial unit at position p */
template<class S>
void SpatialUnit<S>::addParticle(Particle<S>* p) {
   // first check to make sure that we can actually add a particle
   if(checkPos(vec2pos(p->getPos()))) {
       _particles->push_back(p); // add the particle 
   } else {
     std::runtime_error("Error: this particle cannot fit in within this cube!\n"); 
   } 
}

/**! adds a local particle to a spatial unit */
template<class S>
void SpatialUnit<S>::addLocalParticle(Particle<S> *p){
     // check to make sure the position is local  
     Vector3D<S> p_pos = p->getPos();

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
template <class S>
bool SpatialUnit<S>::removeParticle(Particle<S> *p){
  for(SpatialUnit<S>::iterator i=begin(); i!=end(); ++i){
      Particle<S> *curr = *i;
      if(p->getID() == curr->getID()){
          // we've found our particle
          _particles->erase(i);
          return true;
      } 
  } 
  return false;
}

// returns a copy of the particles list (so things can be iterated over and removed)
template<class S>
std::vector<Particle<S> *> SpatialUnit<S>::copyOfParticles(){
    std::vector<Particle<S> *> t;
    for(SpatialUnit<S>::iterator i=begin(); i!=end(); ++i){
         Particle<S> *p = *i;
         t.push_back(p);
    }
    return t;
}

/**! returns the number of beads within this spatial unit */
template<class S>
unsigned SpatialUnit<S>::numBeads() {
    return _particles->size();
}

/**! returns the position of this cube */
template<class S>
position_t<S> SpatialUnit<S>::getPos() {
    return _pos;
}

/**! returns the size of this cube */
template<class S>
S SpatialUnit<S>::getSize() {
    return _size;
}

/**! iterators for the start and end of the particle vector */
template<class S>
typename SpatialUnit<S>::iterator SpatialUnit<S>::begin(){ return _particles->begin(); }
template<class S>
typename SpatialUnit<S>::iterator SpatialUnit<S>::end(){ return _particles->end(); }

/**! iterators for the start and end of the neighbours vector */
template<class S>
typename SpatialUnit<S>::n_iterator SpatialUnit<S>::n_begin(){ return _neighbours->begin(); }
template<class S>
typename SpatialUnit<S>::n_iterator SpatialUnit<S>::n_end(){ return _neighbours->end(); }

#endif /* _SPATIAL_UNIT_IMPL */
