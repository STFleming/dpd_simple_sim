// Oil and water DPD simulation
// 3 bead types: water, red oil, and yellow oil

#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include "fixed_ap.h"
#include <random>

#define DELTA_T 0.02 
#define UNISIZE_D 10.0 // the size of a single dimension of the universe
#define R_C 1.0 
typedef float p_type;
//typedef fixap<int32_t, 27> p_type;

const p_type A[3][3] = { {p_type(25.0), p_type(75.0), p_type(35.0)}, 
                        {p_type(75.0), p_type(25.0), p_type(50.0)},  
                        {p_type(35.0), p_type(50.0), p_type(25.0)}}; // interaction matrix

Vector3D<p_type> randPos(unsigned N){
    Vector3D<p_type> t_pos;
    p_type x(rand() / (float)RAND_MAX * N);
    p_type y(rand() / (float)RAND_MAX * N);
    p_type z(rand() / (float)RAND_MAX * N);
    t_pos.set(x,y,z);
    return t_pos;
}


//! generates a random position within a given space (NxN)
Vector3D<p_type> rand2DPos(unsigned N){
    Vector3D<p_type> t_pos;
    p_type x(rand() / (float)RAND_MAX * N);
    p_type y(rand() / (float)RAND_MAX * N);
    p_type z(0.0);
    t_pos.set(x,y,z);
    return t_pos;
}


// conservative pairwise force declaration
void conF(Particle<p_type> *me, Particle<p_type> *other){

    p_type a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const p_type r_c(R_C); // the interaction cutoff

    Vector3D<p_type> r_i = me->getPos();
    Vector3D<p_type> r_j = other->getPos();
    p_type r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<p_type> r_ij = r_i - r_j;
 
    // Equation 8.5 in the dl_meso manual
    Vector3D<p_type> force = (r_ij/r_ij_dist) * (a_ij * (p_type(1.0) - (r_ij_dist/r_c)));

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    return;
}

// drag pairwie force declaration
void dragF(Particle<p_type> *me, Particle<p_type> *other) {
    
    const p_type r_c(R_C); // the interaction cutoff
    const p_type drag_coef(4.5); // the drag coefficient

    // position and distance
    Vector3D<p_type> r_i = me->getPos();
    Vector3D<p_type> r_j = other->getPos();
    p_type r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<p_type> r_ij = r_j - r_i; // vector between the points

    // switching function
    p_type w_d = (p_type(1.0) - r_ij_dist / r_c)*(p_type(1.0) - r_ij_dist / r_c);

    // velocities
    Vector3D<p_type> v_i = me->getVelo(); 
    Vector3D<p_type> v_j = other->getVelo();
    Vector3D<p_type> v_ij = v_i - v_j; // relative velocity
     
    Vector3D<p_type> force = (r_ij / (r_ij_dist * r_ij_dist)) * w_d * r_ij.dot(v_ij) * (p_type(-1.0) * drag_coef); 

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
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
void randF(uint32_t grand, Particle<p_type> *me, Particle<p_type> *other) {

   const p_type K_BT(1.0);
   const p_type drag_coef(4.5); // the drag coefficient (no idea what to set this at)
   const p_type sigma_ij(sqrt(2*drag_coef*K_BT)); // the temperature coef
   const p_type dt(DELTA_T); 
   const p_type r_c(R_C); // the interaction cutoff

   // position and distance
   Vector3D<p_type> r_i = me->getPos();
   Vector3D<p_type> r_j = other->getPos();
   p_type r_ij_dist = r_i.dist(r_j); // get the distance
   Vector3D<p_type> r_ij = r_j - r_i; // vector between the points

   // switching function
   p_type w_r = (p_type(1.0) - r_ij_dist/r_c);
      
   // random number generation
   p_type r((pairwise_rand(grand, me->getID(), other->getID()) / (float)(RAND_MAX)) * 0.5);

   // force calculation
   Vector3D<p_type> force = (r_ij / r_ij_dist)*sqrt(dt)*r*w_r*sigma_ij;  

   // update the forces acting on the two particles
   me->setForce( me->getForce() + force*p_type(-1.0)); 
   return;
}

// Test program
int main() {
   // size of the universe
   const float unisize = UNISIZE_D;

   // mass of the particles
   const p_type mass_p0(1.0);

   const unsigned num_cubes = 2;
   const p_type cube_size(unisize/num_cubes);

   SimSystem<p_type> universe(p_type(unisize), p_type(DELTA_T), p_type(R_C), num_cubes, 0);

   // add water
   for(int w=0; w<500; w++){
       Particle<p_type> *p = new Particle<p_type>(rand2DPos(unisize), 0, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // add orange oil 
   //for(int o_o=0; o_o<300; o_o++){
   //    Particle<p_type> *p = new Particle<p_type>(rand2DPos(unisize), 1, mass_p0, conF, dragF, randF);
   //    universe.addParticle(p);
   //}

   //// add green oil 
   //for(int g_o=0; g_o<200; g_o++){
   //    Particle<p_type> *p = new Particle<p_type>(rand2DPos(unisize), 2, mass_p0, conF, dragF, randF);
   //    universe.addParticle(p);
   //}

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");
   
   // run the simulation universe for some timesteps 
   universe.run(-1, 0.1);

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");

}
