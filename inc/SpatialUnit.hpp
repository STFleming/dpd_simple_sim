#ifndef __SPATIALUNIT_H
#define __SPATIALUNIT_H

#include <vector> 
#include <stdexcept>
#include <assert.h>
#include "Particle.hpp"
#include "utils.hpp"

// lets make the spatial units addressable
typedef struct spatial_unit_address_ {
    int x;
    int y;
    int z;
} spatial_unit_address_t;

//! SpatialUnit
//! A class used to define the system that is being simulated
template<class S>
class SpatialUnit {
    public:
        SpatialUnit(S size, S x, S y, S z, spatial_unit_address_t addr, unsigned verbosity=0); /**< constructor takes the size of this unit of space and it's center position */       
        ~SpatialUnit(); /**< Destructor */
        position_t<S> getPos(); /**< returns the position for this cube */
        S getSize(); /**< returns the size of this cube */
        void addParticle(Particle<S>* p); /**< Adds a particle to this spatial unit at (global) position p */ 
        void addLocalParticle(Particle<S> *p); /**< Adds a particle to this spatial unit at (local) position p */
        bool removeParticle(Particle<S> *p); /**< Removes a particle from this SpatialUnit */
        unsigned numBeads(); /**< returns the number of beads for this spatial unit */
        bool checkPos(position_t<S> p); /**< returns True if this position is within this spatial unit */
        void addNeighbour(SpatialUnit<S> *s); /**< Adds a neighbour to this spatial unit */
        unsigned numNeighbours(); /**, returns the number of neighbours this spatial unit has */
        spatial_unit_address_t getAddr(); /**< returns the address for this spatial unit */

        // make a copy of the spatial_unit list
        std::vector<Particle<S> *> copyOfParticles(); /**< returns a copy of the _particles list */

        // particle iterators
        typedef typename std::vector<Particle<S> *>::iterator iterator; /**< iterator type for the particles */
        iterator begin(); /**< returns the start of the particle vector */
        iterator end(); /**< returns the end of the particle vector */

        // neighbour iterators
        typedef typename std::vector<SpatialUnit<S> *>::iterator n_iterator; /**< iterator type for the neighbours */
        n_iterator n_begin(); /**< returns the start of the neighbours vector */
        n_iterator n_end(); /**< returns the end of the neighbours vector */
        
    private:
        S _size; /**< The size of one of the dimensions of this cube _size*_size*_size */
        position_t<S> _pos; /**< The start position of this cube */ 
        std::vector<Particle<S>*> *_particles; /**< the collection of particles being managed by this cube */ 
        unsigned _verbosity; /**< the verbosity of the output */
        std::vector<SpatialUnit<S> *> *_neighbours; /**< the list of neighbouring spatial units */
        spatial_unit_address_t _addr; /**< A unique address for this spatial unit (location based) */
};


#include "../src/SpatialUnit.cpp"

#endif /* __SPATIALUNIT_H */ 
