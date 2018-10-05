#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){

    const float a_ij = 0.05; // the interaction strength
    const float r_c = 1.5; // the interaction cutoff

    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D r_ij = r_j - r_i; // vector between the points
 
    // Equation 8.5 in the dl_meso manual
    Vector3D force = (r_ij/r_ij_dist) * a_ij * (1.0 - (r_ij_dist/r_c));

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    other->setForce( other->getForce() + force); 
    return;
}

// drag pairwie force declaration
void dragF(Particle *me, Particle *other) {
   // do the pairwise drag force calculation
}

// random pairwise force declaration
void randF(Particle *me, Particle *other) {
   // do the pairwise random force calculation
}

// Test program
int main() {
   // size of the universe
   const unsigned unisize = 3000;

   // number of particles (beads) in the universe
   const unsigned n = 500;

   SimSystem universe(unisize, 10, 0);
   
   // Add some particles to the system
   for(unsigned i=0; i<n; i++) {
       Particle *p = new Particle(randPos(unisize), conF, dragF, randF);
       universe.addParticle(p);
   }
 
   // emit the initial state (read by the web renderer interface)
   universe.emitJSON("state.json");

   // run the universe on a single thread for many timesteps, emitting it's value every 0.07 seconds 
   universe.seq_run(1000000000, 0.07);

   std::cout << "done\n";
}
