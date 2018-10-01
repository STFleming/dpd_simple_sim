#include "Particle.hpp"

// constructor (sets the initial position of the particle)
Particle::Particle(position_t pos) {
   _pos = pos;
}

// destructor (destroys this particle)
Particle::~Particle() {
}

// sets a new position for this particle
void Particle::setPos(position_t pos) {
    _pos = pos;
}

// gets the current position of the particle
position_t Particle::getPos(){
   return _pos; 
}
