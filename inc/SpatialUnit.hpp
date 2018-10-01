#ifndef __SPATIALUNIT_H
#define __SPATIALUNIT_H

#include <vector> 
#include <stdexcept>
#include <assert.h>
#include "Particle.hpp"
#include "utils.h"

//! SpatialUnit
//! A class used to define the system that is being simulated
class SpatialUnit {
    public:
        SpatialUnit(unsigned size, uint16_t x, uint16_t y, uint16_t z); /**< constructor takes the size of this unit of space and it's center position */       
        ~SpatialUnit(); /**< Destructor */
        position_t getPos(); /**< returns the position for this cube */
        unsigned getSize(); /**< returns the size of this cube */
    private:
        unsigned _size; /**< The size of one of the dimensions of this cube _size*_size*_size */
        position_t _center; /**< The center of this cube */ 
        std::vector<Particle*> *_particles; /**< the collection of particles being managed by this cube */ 
};


#endif /* __SPATIALUNIT_H */ 
