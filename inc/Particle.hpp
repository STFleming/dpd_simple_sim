#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "utils.hpp"
#include <cstdint>
#include <functional>
#include "Vector3D.hpp"

// Particle class
template<class S>
class Particle {
    public:
        Particle(Vector3D<S> pos, uint32_t type, S mass, std::function<void(Particle<S> *me, Particle<S> *other)> conservative, std::function<void(Particle<S> *me, Particle<S> *other)> drag, std::function<void(uint32_t grand, Particle<S> *me, Particle<S> *other)> random); /**< Constructor that sets the initial position of the particle */
        ~Particle(); /**< destructor that destroys this particle */ 
        Vector3D<S> getPos(); /**< get the position of this particle */
        void setPos(Vector3D<S> npos); /**< set a new position for this particle */
        void setPrevPos(Vector3D<S> npos); /**< set a new previous position for this particle and rewrite history */
        Vector3D<S> getPrevPos(); /**< get the previous position of this particle */
        Vector3D<S> getVelo(); /**< gets the velocity for this particle */
        void setVelo(Vector3D<S> velo); /**< sets the velocity for this particle */
        S getMass(); /**< return the mass of this particle */
        uint32_t getID(); /**< returns the ID for this particle */
        void setID(uint32_t id); /**< sets teh ID of the particle */
        uint32_t getType(); /**< gets the type of this particle */
  
        // for updating the force of this particle
        void setForce(Vector3D<S> force); /**< sets a new force for this particle */
        Vector3D<S> getForce(void); /**< returns the current force of this particle */

        // for calling the separate force functions
        void callConservative(Particle<S> *other); /**< calls the conservative force and applied it to this */
        void callDrag(Particle<S> *other); /**< calls the drag force and applied to this */
        void callRandom(uint32_t grand, Particle<S> *other); /**< calls the random force function applied to this */
        void callBond(); /**< calls the force on the bonded particle from this particles perspective */
        void callInverseBond(); /**< calls the force on the bond from the other particles perspective */

        // sets the bond force
        void setInBond(Particle<S> *p, std::function<void(Particle<S> *me, Particle<S> *other)> bondf); /**< bonds this particle to another particle */
        void setOutBond(Particle<S> *p, std::function<void(Particle<S> *me, Particle<S> *other)> bondf); /**< bonds this particle to another particle */
        bool isBonded(); 
        Particle<S> * getInBondBead(); /**< returns a pointer to the bonded particle */
        Particle<S> * getOutBondBead(); /**< returns a pointer to the bonded particle */

    private:
        Vector3D<S> _velocity; /**< the current velocity of the particle */
        Vector3D<S> _pos; /**< the current position of this particle*/ 
        Vector3D<S> _prev_pos; /**< the current position of this particle*/ 
        S _mass; /**< the mass of this particle */
        uint32_t _id; /**< the unique ID for this particle */
        Vector3D<S> _force; /**< the forces accumulated on this particle for this timestep */
        uint32_t _type; /**< identifier for the type of this particle (user configurable) */
        bool _isBonded; /**< True if this particle is bonded to another particle */
        Particle<S>* _inBond; /**< the particle that this particle may or may not be bonded to */
        Particle<S>* _outBond; /**< the particle that this particle may or may not be bonded to */

        // Force functions
        std::function<void(Particle<S> * me, Particle<S> * other)> _conservative; /**< the pairwise conservative force function */ 
        std::function<void(Particle<S> * me, Particle<S> * other)> _drag; /**< the pairwise drag force function */
        std::function<void(uint32_t grand, Particle<S> * me, Particle<S> * other)> _random; /**< the pairwise random force function */
        std::function<void(Particle<S> *me, Particle<S> *other)> _bond; /**< the pairwise bonded force (may or may not exist) */
};

#include "../src/Particle.cpp"

#endif /* __PARTICLE_H */
