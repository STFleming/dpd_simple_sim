#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// Test program
int main() {
   SimSystem sys(100, 10, 1);

   // Add some particles to the system
   Particle *p = new Particle(randPos(100));
   printf("Adding a new particle to the system at x:%d, y:%d, z:%d\n", p->getPos().x, p->getPos().y, p->getPos().z);

   std::cout << "done\n";
}
