#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// Test program
int main() {
   SimSystem sys(100, 10, 3);

   // Add some particles to the system
   const unsigned n = 100;
   for(unsigned i=0; i<n; i++){
       Particle *p = new Particle(randPos(100));
       sys.addParticle(p);
       printf("Adding a new particle to the system at x:%d, y:%d, z:%d\n", p->getPos().x, p->getPos().y, p->getPos().z);
   }

   std::cout << "done\n";
}
