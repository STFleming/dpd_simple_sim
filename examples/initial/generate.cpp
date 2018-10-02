#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// generates some particles program
int main() {
   SimSystem universe(100, 10, 0);

   // Add some particles to the system
   const unsigned n = 100;
   for(unsigned i=0; i<n; i++){
       Particle *p = new Particle(randPos(100));
       universe.addParticle(p);
   }
   
   universe.emitJSON();
   std::cout << "done\n";
}
