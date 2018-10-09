#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>

#define DELTA_T 0.00001 
#define UNISIZE_D 10.0 // the size of a single dimension of the universe
#define R_C 1.0 

//const float A[2][2] = { {0.0, 0.0}, {0.0, 0.0}}; // interaction matrix
const float A[2][2] = { {2.0, 2.0}, {2.0, 2.0}}; // interaction matrix

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){

    const float a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const float r_c = R_C; // the interaction cutoff

    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.toroidal_dist(r_j, UNISIZE_D); // get the distance
    Vector3D r_ij = r_j - r_i; // vector between the points
 
    // Equation 8.5 in the dl_meso manual
    Vector3D force = (r_ij/r_ij_dist) * -1.0*a_ij * (1.0 - (r_ij_dist/r_c));

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    other->setForce( other->getForce() + force*-1.0); 
    return;
}

// drag pairwie force declaration
void dragF(Particle *me, Particle *other) {
    return;
}

// random pairwise force declaration
void randF(Particle *me, Particle *other) {
   return;
}

// Test program
int main() {
   // size of the universe
   const float unisize = UNISIZE_D;

   // mass of the particles
   const float mass_p0 = 1.0;
   const float mass_p1 = 1.0;

   SimSystem universe(unisize, DELTA_T, R_C, 10, 0);
   
   // Add some particles to the system
   Particle *p1 = new Particle(Vector3D(unisize/2.0,0.0,0.0), 0, mass_p0, conF, dragF, randF);
   universe.addParticle(p1);
   Particle *p2 = new Particle(Vector3D(unisize/2.0 + 0.8,0.0,0.0), 0, mass_p0, conF, dragF, randF);
   universe.addParticle(p2);
 
   // emit the initial state (read by the web renderer interface)
   universe.emitJSON("state.json");

   // run the universe on a single thread for many timesteps, emitting it's value every 0.07 seconds 
   universe.seq_run(1000000000, 0.07);

   std::cout << "done\n";
}
