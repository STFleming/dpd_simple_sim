#ifndef __SIM_SYSTEM_H
#define __SIM_SYSTEM_H

#include "SpatialUnit.hpp" // contains the particle class
#include <vector> 
#include <stdexcept>
#include <assert.h>
#include <cstdio>

#include "Particle.hpp"
#include "utils.hpp"

//! SimSystem
//! A class used to define the system that is being simulated
class SimSystem {
    public:
        SimSystem(float N, unsigned D, unsigned verbosity=0); /**< constructor takes the size of the problem as an argument */       
        ~SimSystem(); /**< Destructor */

        // Iterator type
        typedef std::vector<SpatialUnit *>::iterator iterator;
        iterator begin(); /**< returns the begining of the _cubes vector */
        iterator end(); /**< returns the end of the _cubes vector */

        // Particle management
        void addParticle(Particle *p); /**< Add a particle to the system */

    private:
        float _N; /**< the size of the problem is NxNxN */
        unsigned _D; /**< the number of discrete cubes along an axis */
        unsigned _verbosity = 0; /**< the verbosity of the simulation*/
        float _unit_size; /**< each spatial unit has a size _unit_size x _unit_size x _unit_size */
        std::vector<SpatialUnit *>* _cubes; /**< contains the cubes of SpatialUnits (chunks) of the problem space */ 
};


#endif /* __SIM_SYSTEM_H */ 
