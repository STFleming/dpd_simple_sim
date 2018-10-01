#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "utils.h"

// Particle class
class Particle {
    public:
        Particle(position_t pos); /**< Constructor that sets the initial position of the particle */
        ~Particle(); /**< destructor that destroys this particle */ 
        position_t getPos(); /**< get the position of this particle */
        void setPos(position_t npos); /**< set a new position for this particle */

    private:
        uint16_t _velocity; /**< the current velocity of the particle */
        position_t _pos; /**< the current position of this particle*/ 
        uint16_t _mass; /**< the mass of this particle */
};


#endif /* __PARTICLE_H */
