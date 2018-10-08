#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>

#define DELTA_T 1.0/128 
#define UNISIZE_D 10.0 // the size of a single dimension of the universe
#define R_C 1.0 

//const float A[2][2] = { {0.25, 16.0}, {16.0, 0.25}}; // interaction matrix
const float A[2][2] = { {0.0, 0.0}, {0.0, 0.0}}; // interaction matrix

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){

    const float a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const float r_c = R_C; // the interaction cutoff

    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.toroidal_dist(r_j, UNISIZE_D); // get the distance
    Vector3D r_ij = r_j - r_i; // vector between the points
 
    // Equation 8.5 in the dl_meso manual
    Vector3D force = (r_ij/r_ij_dist) * a_ij * (1.0 - (r_ij_dist/r_c));

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    other->setForce( other->getForce() + force*-1.0); 
    return;
}

// drag pairwie force declaration
void dragF(Particle *me, Particle *other) {
    
    const float r_c = R_C; // the interaction cutoff
    const float drag_coef = 0.998; // the drag coefficient (no idea what to set this at)

    // position and distance
    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.toroidal_dist(r_j, UNISIZE_D); // get the distance
    //float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D r_ij = r_j - r_i; // vector between the points

    // switching function
    float w_d = (1.0 - r_ij_dist / r_c)*(1.0 - r_ij_dist / r_c)*(1.0 - r_ij_dist / r_c)*(1.0 - r_ij_dist / r_c);

    // velocities
    Vector3D v_i = me->getVelo(); 
    Vector3D v_j = other->getVelo();
    Vector3D v_ij = v_i - v_j; // relative velocity
     
    // Equation 8.7 in dl_meso manual
    Vector3D force = (r_ij / (r_ij_dist * r_ij_dist)) * w_d * r_ij.dot(v_ij) * (-1.0 * drag_coef); 

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    other->setForce( other->getForce() + force*-1.0); 
    return;
}

// random pairwise force declaration
void randF(Particle *me, Particle *other) {

   float rand_var = (rand() / (float)RAND_MAX * 1.0);

   const float K_BT = 20.0;
   const float drag_coef = 0.998; // the drag coefficient (no idea what to set this at)
   const float sigma_ij = sqrt(2*drag_coef*K_BT); // the temperature coef
   const float dt = DELTA_T; 
   const float r_c = R_C; // the interaction cutoff

   // position and distance
   Vector3D r_i = me->getPos();
   Vector3D r_j = other->getPos();
   float r_ij_dist = r_i.toroidal_dist(r_j, UNISIZE_D); // get the distance
   Vector3D r_ij = r_j - r_i; // vector between the points

   // switching function
   float w_r = (1.0 - r_ij_dist/r_c)*(1.0 - r_ij_dist/r_c);
      
   // Equation 8.6 in dl_meso for random force 
   Vector3D force = (r_ij / r_ij_dist)*sqrt(dt)*rand_var*w_r*sigma_ij;  

   // update the forces acting on the two particles
   me->setForce( me->getForce() + force); 
   other->setForce( other->getForce() + force*-1.0); 
   return;
}

// Test program
int main() {
   // size of the universe
   const float unisize = UNISIZE_D;

   // number of particles (beads) in the universe
   const unsigned n = 1000;

   // mass of the particles
   const float mass_p0 = 1.0;
   const float mass_p1 = 1.0;

   SimSystem universe(unisize, DELTA_T, R_C, 10, 0);
   
   // Add some particles to the system
   for(unsigned i=0; i<n/2; i++) {
       // add a particle with type 0
       Particle *p1 = new Particle(randPos(unisize), 0, mass_p0, conF, dragF, randF);
       universe.addParticle(p1);
       // add a particle with type 1 
       Particle *p2 = new Particle(randPos(unisize), 1, mass_p1, conF, dragF, randF);
       universe.addParticle(p2);
   }
 
   // emit the initial state (read by the web renderer interface)
   universe.emitJSON("state.json");

   // run the universe on a single thread for many timesteps, emitting it's value every 0.07 seconds 
   universe.seq_run(1000000000, 0.07);

   std::cout << "done\n";
}
