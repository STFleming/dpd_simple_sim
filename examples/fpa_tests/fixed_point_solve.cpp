// Oil and water DPD simulation
// 3 bead types: water, red oil, and yellow oil

#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>
#include "fixed_ap.h"

#define DELTA_T 0.02 
#define UNISIZE_D 10.0 // the size of a single dimension of the universe
#define R_C 1.0 

const fixap<int32_t,28> A[3][3] = {  {fixap<int32_t,28>(25.0), fixap<int32_t,28>(75.0), fixap<int32_t,28>(35.0)}, 
                                     {fixap<int32_t,28>(75.0), fixap<int32_t,28>(25.0), fixap<int32_t,28>(50.0)},  
                                     {fixap<int32_t,28>(35.0), fixap<int32_t,28>(50.0), fixap<int32_t,28>(25.0)}};

Vector3D<fixap<int32_t,28>> randPos(unsigned N){
    Vector3D<fixap<int32_t,28>> t_pos;
    fixap<int32_t,28> x(rand() / (float)RAND_MAX * N);
    fixap<int32_t,28> y(rand() / (float)RAND_MAX * N);
    fixap<int32_t,28> z(rand() / (float)RAND_MAX * N);
    t_pos.set(x,y,z);
    return t_pos;
}


//! generates a random position within a given space (NxN)
Vector3D<fixap<int32_t,28>> rand2DPos(unsigned N){
    Vector3D<fixap<int32_t,28>> t_pos;
    fixap<int32_t,28> x(rand() / (float)RAND_MAX * N);
    fixap<int32_t,28> y(rand() / (float)RAND_MAX * N);
    fixap<int32_t,28> z(0.0);
    t_pos.set(x,y,z);
    return t_pos;
}


// conservative pairwise force declaration
void conF(Particle<fixap<int32_t,28>> *me, Particle<fixap<int32_t,28>> *other){

    fixap<int32_t,28> a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const fixap<int32_t,28> r_c(R_C); // the interaction cutoff

    Vector3D<fixap<int32_t,28>> r_i = me->getPos();
    Vector3D<fixap<int32_t,28>> r_j = other->getPos();
    fixap<int32_t,28> r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<fixap<int32_t,28>> r_ij = r_i - r_j;
 
    // Equation 8.5 in the dl_meso manual
    Vector3D<fixap<int32_t, 28>> force = (r_ij/r_ij_dist) * (a_ij * (fixap<int32_t,28>(1.0) - (r_ij_dist/r_c)));

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    return;
}

// drag pairwie force declaration
void dragF(Particle<fixap<int32_t,28>> *me, Particle<fixap<int32_t,28>> *other) {
    
    const fixap<int32_t,28> r_c(R_C); // the interaction cutoff
    const fixap<int32_t,28> drag_coef(4.5); // the drag coefficient (no idea what to set this at)

    // position and distance
    Vector3D<fixap<int32_t,28>> r_i = me->getPos();
    Vector3D<fixap<int32_t,28>> r_j = other->getPos();
    fixap<int32_t,28> r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<fixap<int32_t,28>> r_ij = r_j - r_i; // vector between the points

    // switching function
    fixap<int32_t,28> w_d = (fixap<int32_t,28>(1.0) - r_ij_dist / r_c)*(fixap<int32_t,28>(1.0) - r_ij_dist / r_c);

    // velocities
    Vector3D<fixap<int32_t,28>> v_i = me->getVelo(); 
    Vector3D<fixap<int32_t,28>> v_j = other->getVelo();
    Vector3D<fixap<int32_t,28>> v_ij = v_i - v_j; // relative velocity
     
    Vector3D<fixap<int32_t,28>> force = (r_ij / (r_ij_dist * r_ij_dist)) * w_d * r_ij.dot(v_ij) * (fixap<int32_t,28>(-1.0) * drag_coef); 

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    //printf("p:%d dragF new force=%s\n", me->getID(), me->getForce().str().c_str());
    return;
}

// dt10's hash based random num gen
uint32_t pairwise_rand(uint32_t grand, uint32_t pid1, uint32_t pid2){
    uint32_t la=std::min(pid1, pid2);
    uint32_t lb=std::max(pid1, pid2);
    uint32_t s0 = (pid1 ^ grand)*pid2;
    uint32_t s1 = (pid2 ^ grand)*pid1;
    return s0 + s1;
}

// random pairwise force declaration
void randF(uint32_t grand, Particle<fixap<int32_t,28>> *me, Particle<fixap<int32_t,28>> *other) {

   const fixap<int32_t,28> K_BT(1.0);
   const fixap<int32_t,28> drag_coef(4.5); // the drag coefficient (no idea what to set this at)
   const fixap<int32_t,28> sigma_ij(sqrt(2*drag_coef*K_BT)); // the temperature coef
   const fixap<int32_t,28> dt(DELTA_T); 
   const fixap<int32_t,28> r_c(R_C); // the interaction cutoff

   // position and distance
   Vector3D<fixap<int32_t,28>> r_i = me->getPos();
   Vector3D<fixap<int32_t,28>> r_j = other->getPos();
   fixap<int32_t,28> r_ij_dist = r_i.dist(r_j); // get the distance
   Vector3D<fixap<int32_t,28>> r_ij = r_j - r_i; // vector between the points

   // switching function
   fixap<int32_t,28> w_r = (fixap<int32_t,28>(1.0) - r_ij_dist/r_c);
      
   // random number generation
   fixap<int32_t,28> r((pairwise_rand(grand, me->getID(), other->getID()) / (float)(RAND_MAX)) * 0.5);

   // force calculation
   Vector3D<fixap<int32_t,28>> force = (r_ij / r_ij_dist)*sqrt(dt)*r*w_r*sigma_ij;  

   // update the forces acting on the two particles
   me->setForce( me->getForce() + force*fixap<int32_t,28>(-1.0)); 
   return;
}

// Test program
int main() {
   // size of the universe
   const float unisize = UNISIZE_D;

   // mass of the particles
   const fixap<int32_t,28> mass_p0(1.0);

   const unsigned num_cubes = 10;
   const fixap<int32_t,28> cube_size(unisize/num_cubes);

   SimSystem<fixap<int32_t,28>> universe(fixap<int32_t,28>(unisize), fixap<int32_t,28>(DELTA_T), fixap<int32_t,28>(R_C), num_cubes, 0);

   // add water
   for(int w=0; w<500; w++){
       Particle<fixap<int32_t,28>> *p = new Particle<fixap<int32_t,28>>(rand2DPos(unisize), 0, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // add orange oil 
   for(int o_o=0; o_o<300; o_o++){
       Particle<fixap<int32_t, 28>> *p = new Particle<fixap<int32_t,28>>(rand2DPos(unisize), 1, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // add green oil 
   for(int g_o=0; g_o<200; g_o++){
       Particle<fixap<int32_t,28>> *p = new Particle<fixap<int32_t,28>>(rand2DPos(unisize), 2, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");
   
   // run the simulation universe for some timesteps 
   universe.run(-1, 0.1);

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");

}
