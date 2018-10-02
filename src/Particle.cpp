#include "Particle.hpp"

// constructor (sets the initial position of the particle)
Particle::Particle(position_t pos) {
   _pos = pos;
   _velocity = 0.0;
   _mass = 0.0;
   _id = 0xFFFFFFFF;
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

// gets the ID of this particle
uint32_t Particle::getID() {
    return _id;
}

// sets the ID of this particle
void Particle::setID(uint32_t id) {
    _id = id;
    return;
} 
