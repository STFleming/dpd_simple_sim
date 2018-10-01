#ifndef __SPATIALUNIT_H
#define __SPATIALUNIT_H

#include <vector> 
#include <stdexcept>
#include <assert.h>
#include "Particle.hpp"
#include "utils.hpp"

//! SpatialUnit
//! A class used to define the system that is being simulated
class SpatialUnit {
    public:
        SpatialUnit(unsigned size, uint16_t x, uint16_t y, uint16_t z, unsigned verbosity=0); /**< constructor takes the size of this unit of space and it's center position */       
        ~SpatialUnit(); /**< Destructor */
        position_t getPos(); /**< returns the position for this cube */
        unsigned getSize(); /**< returns the size of this cube */
        void addParticle(Particle* p); /**< Adds a particle to this spatial unit at (global) position p */ 
        bool checkPos(position_t p); /**< returns True if this position is within this spatial unit */

        // TODO add iterators for all the particles in a spatial unit
    private:
        unsigned _size; /**< The size of one of the dimensions of this cube _size*_size*_size */
        position_t _center; /**< The center of this cube */ 
        std::vector<Particle*> *_particles; /**< the collection of particles being managed by this cube */ 
        unsigned _verbosity; /**< the verbosity of the output */
};


#endif /* __SPATIALUNIT_H */ 
