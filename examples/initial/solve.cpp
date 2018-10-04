#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){
    // do the pairwise conservative force
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
   const unsigned unisize = 200;

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

   // run the universe for 10000 timesteps, emitting it's value every 0.5 seconds 
   universe.run(1000000000, 0.1);

   std::cout << "done\n";
}
