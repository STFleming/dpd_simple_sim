#ifndef __SIM_SYSTEM_H
#define __SIM_SYSTEM_H

#include "SpatialUnit.hpp" // contains the particle class
#include <vector> 
#include <tuple> 
#include <stdexcept>
#include <assert.h>
#include <cstdio>
#include <string>
#include <sstream>

#include "Particle.hpp"
#include "utils.hpp"

#include <iostream> // for emitting JSON files
#include <fstream>
#include <functional>
#include <jsoncpp/json/value.h> // for parsing JSON files
#include <jsoncpp/json/reader.h> // for parsing JSON files

#include <ctime> // for timing the emit output

//! SimSystem
//! A class used to define the system that is being simulated
template <class S> // S is the type for this simulation for example a fixedap type or double/float
class SimSystem {
    public:
        SimSystem(S N, S dt, S r_c,  unsigned D, unsigned verbosity=0); /**< constructor takes the size of the problem as an argument */       
        ~SimSystem(); /**< Destructor */

        // Iterator type
        // default iterators is for the spatial units
        typedef typename std::vector<SpatialUnit<S> *>::iterator iterator;
        iterator begin(); /**< returns the begining of the _cubes vector */
        iterator end(); /**< returns the end of the _cubes vector */
        // but we also have particle iterators for the global particles list
        typedef typename std::vector<Particle<S> *>::iterator p_iterator;
        p_iterator p_begin(); /**< returns the start of the global particle list */
        p_iterator p_end(); /**< returns the start of the global particle list */

        // simulation management
        void run(uint32_t period, float emitrate); /**< runs the simulation for period timesteps */
        void seq_run(uint32_t period, float emitrate); /**< runs the simulation sequentially (no parallelism) for period timesteps emitting it's state every emitrate timesteps */

        // Particle management
        void addParticle(Particle<S> *p); /**< Adds a particle to the system with a global position */
        void addLocalParticle(Particle<S> *p); /**< Add a particle to the system with position relative to its spatial unti */
        //void populateFromJSON(std::string jsonfile); /**< populates the system with particles contained within a JSON file*/
        void allocateParticleToSpatialUnit(Particle<S> *p); /**< Allocates all the particles to a spatial processing unit */
        void bond(Particle<S> *i, Particle<S> *j, std::function<void(Particle<S>* me, Particle<S> *other)> bondf); /**< Bonds particles i and j with bond function bondf */
        SpatialUnit<S> * getSpatialUnit(spatial_unit_address_t x); /**< returns a spatial unit given a spatial unit address */

        // exporting
        void emitJSON(std::string jsonfile); /**< emits global particle IDs and position as one JSON file (used for initial values) */
        void emitJSONFromSU(std::string jsonfile); /**< emits the global particle IDs and the position from each spatial unit as a JSON file */
        // debug & testing
        // print the number of particles per spatial unit
        void printSpatialAllocation();

    private:
        S _t; /**< the current simulation time */
        S _dt; /**< simulation timestep size */
        S _r_c; /**< the cutoff distance */
        unsigned _ts; /**< the current integer timestep */
        S _N; /**< the size of the problem is NxNxN */
        unsigned _D; /**< the number of discrete cubes along an axis */
        unsigned _verbosity = 0; /**< the verbosity of the simulation*/
        S _unit_size; /**< each spatial unit has a size _unit_size x _unit_size x _unit_size */
        std::vector<SpatialUnit<S> *>* _cubes; /**< contains the cubes of SpatialUnits (chunks) of the problem space */ 
        std::vector<Particle<S> *>* _particles; /**< the global list of particles */
        uint32_t _grand; /**< a global random number used for dt10's has based random number tech , updated every simulation step and passed to the rand force function*/
        
};

#include "../src/SimSystem.cpp"

#endif /* __SIM_SYSTEM_H */ 
