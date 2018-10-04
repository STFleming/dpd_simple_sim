#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"

#define A_IJ 0.05 // the interaction strength
#define R_C 50.0 // the interaction cutoff

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){
    // do the pairwise conservative force
    //float f = 0.5 * A_IJ * (1 - (dist(me->getPos(), other->getPos()) / R_C ));
    float r_ij = dist(me->getPos(), other->getPos());
    float tmp_f = A_IJ*(1 - r_ij/R_C)/r_ij; 
    
    // update the forces acting on the two particles
    me->setForce( me->getForce() + f); 
    other->setForce( other->getForce() + f); 
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
