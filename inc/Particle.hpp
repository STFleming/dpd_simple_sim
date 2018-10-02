#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "utils.hpp"
#include <cstdint>

// Particle class
class Particle {
    public:
        Particle(position_t pos); /**< Constructor that sets the initial position of the particle */
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
};


#endif /* __PARTICLE_H */
