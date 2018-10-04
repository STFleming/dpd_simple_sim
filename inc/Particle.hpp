#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "utils.hpp"
#include <cstdint>
#include <functional>

// Particle class
class Particle {
    public:
        Particle(position_t pos, std::function<void(Particle *me, Particle *other)> conservative, std::function<void(Particle *me, Particle *other)> drag, std::function<void(Particle *me, Particle *other)> random); /**< Constructor that sets the initial position of the particle */
        ~Particle(); /**< destructor that destroys this particle */ 
        position_t getPos(); /**< get the position of this particle */
        void setPos(position_t npos); /**< set a new position for this particle */
        uint32_t getID(); /**< returns the ID for this particle */
        void setID(uint32_t id); /**< sets teh ID of the particle */

    private:
        float _velocity; /**< the current velocity of the particle */
        position_t _pos; /**< the current position of this particle*/ 
        float _mass; /**< the mass of this particle */
        uint32_t _id; /**< the unique ID for this particle */

        // Force functions
        std::function<void(Particle * me, Particle * other)> _conservative; /**< the pairwise conservative force function */ 
        std::function<void(Particle * me, Particle * other)> _drag; /**< the pairwise drag force function */
        std::function<void(Particle * me, Particle * other)> _random; /**< the pairwise random force function */
};


#endif /* __PARTICLE_H */
