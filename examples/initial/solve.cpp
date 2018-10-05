#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>

#define DELTA_T 0.001
#define R_C 10.5 

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){

    const float a_ij = 0.005; // the interaction strength
    const float r_c = R_C; // the interaction cutoff

    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D r_ij = r_j - r_i; // vector between the points
 
    // Equation 8.5 in the dl_meso manual
    Vector3D force = (r_ij/r_ij_dist) * a_ij * (1.0 - (r_ij_dist/r_c));

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    other->setForce( other->getForce() + force); 
    return;
}

// drag pairwie force declaration
void dragF(Particle *me, Particle *other) {
    
    const float a_ij = 0.005; // the interaction strength
    const float r_c = R_C; // the interaction cutoff
    const float drag_coef = 0.005; // the drag coefficient (no idea what to set this at)

    // position and distance
    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D r_ij = r_j - r_i; // vector between the points
 
    // switching function
    float w_d = (1 - r_ij_dist / r_c)*(1 - r_ij_dist / r_c);

    // velocities
    Vector3D v_i = me->getVelo(); 
    Vector3D v_j = other->getVelo();
    Vector3D v_ij = v_i - v_j; // relative velocity
     
    Vector3D force = (r_ij / (r_ij_dist * r_ij_dist)) * w_d * r_ij.dot(v_ij) * (-1.0 * drag_coef); 

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    other->setForce( other->getForce() + force); 
    return;
}

// random pairwise force declaration
void randF(Particle *me, Particle *other) {

   std::default_random_engine generator;
   std::normal_distribution<float> rand_var(0.0,1.0);

   const float K_B = 0.001; // no idea
   const float T = 0.0003; // also no idea
   const float drag_coef = 0.005; // the drag coefficient (no idea what to set this at)
   const float sigma_ij = 2*drag_coef*K_B*T; // the temperature coef
   const float dt = DELTA_T; // another guess
   const float r_c = R_C; // the interaction cutoff

   // position and distance
   Vector3D r_i = me->getPos();
   Vector3D r_j = other->getPos();
   float r_ij_dist = r_i.dist(r_j); // get the distance
   Vector3D r_ij = r_j - r_i; // vector between the points

   // switching function
   float w_r = (1 - r_ij_dist/r_c)*(1 - r_ij_dist/r_c)*(1 - r_ij_dist/r_c)*(1 - r_ij_dist/r_c);
      
   // force calculation
   Vector3D force = (r_ij / r_ij_dist)*sqrt(dt)*rand_var(generator)*w_r*sigma_ij;  

   // update the forces acting on the two particles
   me->setForce( me->getForce() + force); 
   other->setForce( other->getForce() + force); 
   return;
}

// Test program
int main() {
   // size of the universe
   const unsigned unisize = 3000;

   // number of particles (beads) in the universe
   const unsigned n = 500;

   // mass of the particles
   const float mass = 0.00001;

   SimSystem universe(unisize, DELTA_T, R_C, 10, 0);
   
   // Add some particles to the system
   for(unsigned i=0; i<n; i++) {
       Particle *p = new Particle(randPos(unisize), mass, conF, dragF, randF);
       universe.addParticle(p);
   }
 
   // emit the initial state (read by the web renderer interface)
   universe.emitJSON("state.json");

   // run the universe on a single thread for many timesteps, emitting it's value every 0.07 seconds 
   universe.seq_run(1000000000, 0.07);

   std::cout << "done\n";
}
