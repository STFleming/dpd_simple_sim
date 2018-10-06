#include "Particle.hpp"

// constructor (sets the initial position of the particle)
Particle::Particle(Vector3D pos, uint32_t type, float mass, std::function<void(Particle *me, Particle *other)> conservative, std::function<void(Particle *me, Particle *other)> drag, std::function<void(Particle *me, Particle *other)> random) {
   _pos = pos;
   _velocity = Vector3D(0.0, 0.0, 0.0);
   _mass = mass;
   _id = 0xFFFFFFFF;
   _force = Vector3D(0.0, 0.0, 0.0);
   _type = type;

   // assign the force functions
   _conservative = conservative;
   _drag = drag;
   _random = random; 
}

// destructor (destroys this particle)
Particle::~Particle() {
}

// sets a new position for this particle
void Particle::setPos(Vector3D npos) {
    _pos = npos;
}

// gets the current position of the particle
Vector3D Particle::getPos(){
   return _pos; 
}

// gets the velocity for the particle
Vector3D Particle::getVelo() {
    return _velocity;
}

// sets the velocity of this particle
void Particle::setVelo(Vector3D velo) {
    _velocity = velo;
}

// get the mass of this particle
float Particle::getMass() {
    return _mass;
}

// gets the ID of this particle
uint32_t Particle::getID() {
    return _id;
}

// gets the type of this particle
uint32_t Particle::getType() {
    return _type;
}

// sets the ID of this particle
void Particle::setID(uint32_t id) {
    _id = id;
    return;
} 

// getters and setters for the force of this particle
void Particle::setForce(Vector3D force) {
    _force = force;
}

Vector3D Particle::getForce(void) {
    return _force;
}

// calls the conservative force function
void Particle::callConservative(Particle *other) {
    _conservative(this, other);
}

// calls the drag force function 
void Particle::callDrag(Particle *other) {
    _drag(this, other);
}

// calls the random force function
void Particle::callRandom(Particle *other) {
    _random(this, other);
}
