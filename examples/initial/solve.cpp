#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// Test program
int main() {
   SimSystem universe(500, 10, 0);
   
   // read in the particles and initial state
   universe.populateFromJSON("state.json");
 
   // run the universe for 10000 timesteps, emitting it's value every 0.5 seconds 
   universe.run(1000000000, 0.1);

   std::cout << "done\n";
}
