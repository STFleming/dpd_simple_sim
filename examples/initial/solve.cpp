#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// Test program
int main() {
   SimSystem universe(100, 10, 4);
   
   // read in the particles and initial state
   universe.populateFromJSON("state.json");
 
   // run the universe for 10000 timesteps, emitting it's value every 100 timesteps
   universe.run(10000, 100);

   std::cout << "done\n";
}
