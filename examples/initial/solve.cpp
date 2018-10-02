#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"

// Test program
int main() {
   SimSystem universe(100, 10, 4);
   
   // read in the particles and initial state
   universe.populateFromJSON("state.json");
 
   std::cout << "done\n";
}
