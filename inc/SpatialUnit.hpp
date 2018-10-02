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
        SpatialUnit(float size, float x, float y, float z, unsigned verbosity=0); /**< constructor takes the size of this unit of space and it's center position */       
        ~SpatialUnit(); /**< Destructor */
        position_t getPos(); /**< returns the position for this cube */
        float getSize(); /**< returns the size of this cube */
        void addParticle(Particle* p); /**< Adds a particle to this spatial unit at (global) position p */ 
        bool checkPos(position_t p); /**< returns True if this position is within this spatial unit */

        // particle iterators
        typedef std::vector<Particle *>::iterator iterator; /**< iterator type for the particles */
        iterator begin(); /**< returns the start of the particle vector */
        iterator end(); /**< returns the end of the particle vector */
        
    private:
        float _size; /**< The size of one of the dimensions of this cube _size*_size*_size */
        position_t _pos; /**< The start position of this cube */ 
        std::vector<Particle*> *_particles; /**< the collection of particles being managed by this cube */ 
        unsigned _verbosity; /**< the verbosity of the output */
};


#endif /* __SPATIALUNIT_H */ 
