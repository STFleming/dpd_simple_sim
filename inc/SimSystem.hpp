#ifndef __SIM_SYSTEM_H
#define __SIM_SYSTEM_H

#include "SpatialUnit.hpp" // contains the particle class
#include <vector> 
#include <stdexcept>
#include <assert.h>
#include <cstdio>
#include <string>

#include "Particle.hpp"
#include "utils.hpp"

#include <iostream> // for emitting JSON files
#include <fstream>
#include <jsoncpp/json/value.h> // for parsing JSON files
#include <jsoncpp/json/reader.h> // for parsing JSON files

//! SimSystem
//! A class used to define the system that is being simulated
class SimSystem {
    public:
        SimSystem(float N, unsigned D, unsigned verbosity=0); /**< constructor takes the size of the problem as an argument */       
        ~SimSystem(); /**< Destructor */

        // Iterator type
        // default iterators is for the spatial units
        typedef std::vector<SpatialUnit *>::iterator iterator;
        iterator begin(); /**< returns the begining of the _cubes vector */
        iterator end(); /**< returns the end of the _cubes vector */
        // but we also have particle iterators for the global particles list
        typedef std::vector<Particle *>::iterator p_iterator;
        p_iterator p_begin(); /**< returns the start of the global particle list */
        p_iterator p_end(); /**< returns the start of the global particle list */

        // Particle management
        void addParticle(Particle *p); /**< Add a particle to the system */
        void populateFromJSON(std::string jsonfile); /**< populates the system with particles contained within a JSON file*/
        void allocateParticleToSpatialUnit(Particle *p); /**< Allocates all the particles to a spatial processing unit */

        // exporting
        void emitJSON(); /**< emits global particle IDs and position as one JSON file (used for initial values) */

    private:
        float _N; /**< the size of the problem is NxNxN */
        unsigned _D; /**< the number of discrete cubes along an axis */
        unsigned _verbosity = 0; /**< the verbosity of the simulation*/
        float _unit_size; /**< each spatial unit has a size _unit_size x _unit_size x _unit_size */
        std::vector<SpatialUnit *>* _cubes; /**< contains the cubes of SpatialUnits (chunks) of the problem space */ 
        std::vector<Particle *>* _particles; /**< the global list of particles */
};


#endif /* __SIM_SYSTEM_H */ 
