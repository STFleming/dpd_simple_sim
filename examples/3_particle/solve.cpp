#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>

#define DELTA_T 0.01 
#define UNISIZE_D 12.0 // the size of a single dimension of the universe
#define R_C 1.0 

#define W 0.75
#define WW 0.5

//const float A[2][2] = { {0.0, 0.0}, {0.0, 0.0}}; // interaction matrix
const float A[3][3] = { {WW, W, W}, {W, 1.0, 0.05},  {W, 0.05, 0.01}}; // interaction matrix

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){

    const float a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const float r_c = R_C; // the interaction cutoff

    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.toroidal_dist(r_j, UNISIZE_D); // get the distance
    //Vector3D r_ij = r_i - r_j; // vector between the points
    Vector3D r_ij = r_i.toroidal_subtraction(r_j, UNISIZE_D, R_C);
 
    // Equation 8.5 in the dl_meso manual
    Vector3D force = (r_ij/r_ij_dist) * (a_ij * (1.0 - (r_ij_dist/r_c)));

    assert(r_ij_dist <= R_C);

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    other->setForce( other->getForce() + force*-1.0); 
    return;
}

// drag pairwie force declaration
void dragF(Particle *me, Particle *other) {
    
    const float r_c = R_C; // the interaction cutoff
    const float drag_coef = 0.998; // the drag coefficient (no idea what to set this at)

   // // position and distance
    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.toroidal_dist(r_j, UNISIZE_D); // get the distance
    //Vector3D r_ij = r_j - r_i; // vector between the points
    Vector3D r_ij = r_i.toroidal_subtraction(r_j, UNISIZE_D, R_C);

    // switching function
    float w_d = (1.0 - r_ij_dist / r_c)*(1.0 - r_ij_dist / r_c);

    // velocities
    Vector3D v_i = me->getVelo(); 
    Vector3D v_j = other->getVelo();
    Vector3D v_ij = v_i - v_j; // relative velocity
     
    Vector3D force = (r_ij / (r_ij_dist * r_ij_dist)) * w_d * r_ij.dot(v_ij) * (-1.0 * drag_coef); 

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    other->setForce( other->getForce() + force*-1.0); 
    return;
}

// random pairwise force declaration
void randF(Particle *me, Particle *other) {

   const float K_BT = 0.05;
   const float drag_coef = 0.998; // the drag coefficient (no idea what to set this at)
   const float sigma_ij = sqrt(2*drag_coef*K_BT); // the temperature coef
   const float dt = DELTA_T; 
   const float r_c = R_C; // the interaction cutoff

   // position and distance
   Vector3D r_i = me->getPos();
   Vector3D r_j = other->getPos();
   float r_ij_dist = r_i.toroidal_dist(r_j, UNISIZE_D); // get the distance
   //Vector3D r_ij = r_j - r_i; // vector between the points
   Vector3D r_ij = r_i.toroidal_subtraction(r_j, UNISIZE_D, R_C);

   // switching function
   float w_r = (1.0 - r_ij_dist/r_c);
      
   // force calculation
   float r = (rand() / (float)RAND_MAX * 1.0);
   Vector3D force = (r_ij / r_ij_dist)*sqrt(dt)*r*w_r*sigma_ij;  

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
   const unsigned n = 1728;

   // mass of the particles
   const float mass_p0 = 1.0;
   const float mass_p1 = 1.0;
   const float mass_p2 = 1.0;

   SimSystem universe(unisize, DELTA_T, R_C, 10, 0);
   
   // Add lots of water
   for(unsigned i=0; i<(n-50); i++) {
       Particle *p1 = new Particle(randPos(unisize), 0, mass_p0, conF, dragF, randF);
       universe.addParticle(p1);
   }

  // Add the first particle
  for(unsigned i=0; i<20; i++) {
       Particle *p2 = new Particle(randPos(unisize/2) + unisize/4, 1, mass_p1, conF, dragF, randF);
       universe.addParticle(p2);
  }
  for(unsigned i=0; i<60; i++) {
       Particle *p3 = new Particle(randPos(unisize/2) + unisize/4, 2, mass_p2, conF, dragF, randF);
       universe.addParticle(p3);
  }

   // emit the initial state (read by the web renderer interface)
   universe.emitJSON("state.json");

   // run the universe on a single thread for many timesteps, emitting it's value every 0.07 seconds 
   universe.seq_run(1000000000, 0.07);

   std::cout << "done\n";
}
