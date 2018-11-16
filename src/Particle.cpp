#include "Particle.hpp"

// constructor (sets the initial position of the particle)
Particle::Particle(Vector3D pos, uint32_t type, float mass, std::function<void(Particle *me, Particle *other)> conservative, std::function<void(Particle *me, Particle *other)> drag, std::function<void(uint32_t grand, Particle *me, Particle *other)> random) {
   _pos = pos;
   _prev_pos = pos+0.0001; // small addition to avoid initial condition problems
   _velocity = Vector3D(0.0, 0.0, 0.0);
   _mass = mass;
   _id = 0xFFFFFFFF;
   _force = Vector3D(0.0, 0.0, 0.0);
   _type = type;

   // by default particles are not bonded to any other particle
   _isBonded = false;
   _bondParticle = NULL; // the particle that this could be bonded to
   _bond = NULL; // the force function that is called between the bonds

   // assign the force functions
   _conservative = conservative;
   _drag = drag;
   _random = random; 
}

// destructor (destroys this particle)
Particle::~Particle() {
}

// used to set particle bonds
bool Particle::isBonded() { return _isBonded;}

// used to set the particle this should be bonded to
void Particle::setBond(Particle *p, std::function<void(Particle *me, Particle *other)> bondf){
    _isBonded = true;
    _bondParticle = p;
    _bond = bondf;  
}

// returns the pointer to the particle that this is bonded to
Particle * Particle::getBondedParticle(){ return _bondParticle; }

// sets a new position for this particle
void Particle::setPos(Vector3D npos) {
    _prev_pos = _pos;
    _pos = npos;
}

// gets the current position of the particle
Vector3D Particle::getPos(){
   return _pos; 
}

// gets the prev position of the particle
Vector3D Particle::getPrevPos(){
   return _prev_pos; 
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
void Particle::callRandom(uint32_t grand, Particle *other) {
    _random(grand, this, other);
}

// calls the bond force function
void Particle::callBond() {
    if(_isBonded){
      _bond(this, getBondedParticle());
    }
}
