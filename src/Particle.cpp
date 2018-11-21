#include "Particle.hpp"

#ifndef _PARTICLE_IMPL
#define _PARTICLE_IMPL

// constructor (sets the initial position of the particle)
template<class S>
Particle<S>::Particle(Vector3D<S> pos, uint32_t type, S mass, std::function<void(Particle<S> *me, Particle<S> *other)> conservative, std::function<void(Particle<S> *me, Particle<S> *other)> drag, std::function<void(uint32_t grand, Particle<S> *me, Particle<S> *other)> random) {
   _pos = pos;
   _prev_pos = pos+S(0.01); // small addition to avoid initial condition problems
   _velocity = Vector3D<S>(S(0.0), S(0.0), S(0.0));
   _mass = mass;
   _id = 0xFFFFFFFF;
   _force = Vector3D<S>(S(0.0), S(0.0), S(0.0));
   _type = type;

   // by default particles are not bonded to any other particle
   _isBonded = false;
   _inBond = NULL; // the particle that this could be bonded to
   _outBond = NULL; // the particle that this could be bonded to
   _bond = NULL; // the force function that is called between the bonds

   // assign the force functions
   _conservative = conservative;
   _drag = drag;
   _random = random; 
}

// destructor (destroys this particle)
template<class S>
Particle<S>::~Particle() {
}

// used to set particle bonds
template<class S>
bool Particle<S>::isBonded() { return _isBonded;}

// used to set the particle this should be bonded to
template<class S>
void Particle<S>::setInBond(Particle<S> *p, std::function<void(Particle<S> *me, Particle<S> *other)> bondf){
    _isBonded = true;
    _inBond = p;
    _bond = bondf;  
}

// used to set the particle this should be bonded to
template<class S>
void Particle<S>::setOutBond(Particle<S> *p, std::function<void(Particle<S> *me, Particle<S> *other)> bondf){
    _isBonded = true;
    _outBond = p;
    _bond = bondf;  
}

// returns the pointer to the particle that this is bonded to
template<class S>
Particle<S> * Particle<S>::getInBondBead(){ return _inBond; }
template<class S>
Particle<S> * Particle<S>::getOutBondBead(){ return _outBond; }

// sets a new position for this particle
template<class S>
void Particle<S>::setPos(Vector3D<S> npos) {
    _pos = npos;
}

// sets a new previous position for this particle
template<class S>
void Particle<S>::setPrevPos(Vector3D<S> npos) {
    _prev_pos = npos;
}

// gets the current position of the particle
template<class S>
Vector3D<S> Particle<S>::getPos(){
   return _pos; 
}

// gets the prev position of the particle
template<class S>
Vector3D<S> Particle<S>::getPrevPos(){
   return _prev_pos; 
}


// gets the velocity for the particle
template<class S>
Vector3D<S> Particle<S>::getVelo() {
    return _velocity;
}

// sets the velocity of this particle
template<class S>
void Particle<S>::setVelo(Vector3D<S> velo) {
    _velocity = velo;
}

// get the mass of this particle
template<class S>
S Particle<S>::getMass() {
    return _mass;
}

// gets the ID of this particle
template<class S>
uint32_t Particle<S>::getID() {
    return _id;
}

// gets the type of this particle
template<class S>
uint32_t Particle<S>::getType() {
    return _type;
}

// sets the ID of this particle
template<class S>
void Particle<S>::setID(uint32_t id) {
    _id = id;
    return;
} 

// getters and setters for the force of this particle
template<class S>
void Particle<S>::setForce(Vector3D<S> force) {
    _force = force;
}

template<class S>
Vector3D<S> Particle<S>::getForce(void) {
    return _force;
}

// calls the conservative force function
template<class S>
void Particle<S>::callConservative(Particle<S> *other) {
    _conservative(this, other);
}

// calls the drag force function 
template<class S>
void Particle<S>::callDrag(Particle<S> *other) {
    _drag(this, other);
}

// calls the random force function
template<class S>
void Particle<S>::callRandom(uint32_t grand, Particle<S> *other) {
    _random(grand, this, other);
}

// calls the bond force function
template<class S>
void Particle<S>::callBond() {
    if(_isBonded){
      _bond(this, getOutBondBead());
    }
}

// calls the bond force function from the other way around
template<class S>
void Particle<S>::callInverseBond() {
   if(_isBonded) {
      _bond(getOutBondBead(), this);
   }
}

#endif /* _PARTICLE_IMPL */
