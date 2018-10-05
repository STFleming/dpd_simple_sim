#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "utils.hpp"
#include <cstdint>
#include <functional>
#include "Vector3D.hpp"

// Particle class
class Particle {
    public:
        Particle(Vector3D pos, float mass, std::function<void(Particle *me, Particle *other)> conservative, std::function<void(Particle *me, Particle *other)> drag, std::function<void(Particle *me, Particle *other)> random); /**< Constructor that sets the initial position of the particle */
        ~Particle(); /**< destructor that destroys this particle */ 
        Vector3D getPos(); /**< get the position of this particle */
        void setPos(Vector3D npos); /**< set a new position for this particle */
        Vector3D getVelo(); /**< gets the velocity for this particle */
        void setVelo(Vector3D velo); /**< sets the velocity for this particle */
        float getMass(); /**< return the mass of this particle */
        uint32_t getID(); /**< returns the ID for this particle */
        void setID(uint32_t id); /**< sets teh ID of the particle */
  
        // for updating the force of this particle
        void setForce(Vector3D force); /**< sets a new force for this particle */
        Vector3D getForce(void); /**< returns the current force of this particle */

        // for calling the separate force functions
        void callConservative(Particle *other); /**< calls the conservative force and applied it to this */
        void callDrag(Particle *other); /**< calls the drag force and applied to this */
        void callRandom(Particle *other); /**< calls the random force function applied to this */

    private:
        Vector3D _velocity; /**< the current velocity of the particle */
        Vector3D _pos; /**< the current position of this particle*/ 
        float _mass; /**< the mass of this particle */
        uint32_t _id; /**< the unique ID for this particle */
        Vector3D _force; /**< the forces accumulated on this particle for this timestep */

        // Force functions
        std::function<void(Particle * me, Particle * other)> _conservative; /**< the pairwise conservative force function */ 
        std::function<void(Particle * me, Particle * other)> _drag; /**< the pairwise drag force function */
        std::function<void(Particle * me, Particle * other)> _random; /**< the pairwise random force function */
};


#endif /* __PARTICLE_H */
