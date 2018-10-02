#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// generates some particles program
int main() {
   SimSystem universe(500, 10, 0);

   // Add some particles to the system
   const unsigned n = 1000;
   for(unsigned i=0; i<n; i++){
       Particle *p = new Particle(randPos(100));
       universe.addParticle(p);
   }
   
   universe.emitJSON("state.json");
   std::cout << "done\n";
}
