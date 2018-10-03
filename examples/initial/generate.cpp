#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// generates some particles program
int main() {
   const unsigned unisize = 200;

   SimSystem universe(unisize, 10, 0);

   // Add some particles to the system
   const unsigned n = 500;
   for(unsigned i=0; i<n; i++){
       Particle *p = new Particle(randPos(unisize));
       universe.addParticle(p);
   }
   
   universe.emitJSON("state.json");
   std::cout << "done\n";
}
