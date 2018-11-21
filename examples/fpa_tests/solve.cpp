// Oil and water DPD simulation
// 3 bead types: water, red oil, and yellow oil

#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>
#include "fixed_ap.h"

#define DELTA_T 0.0002 
//#define DELTA_T 0.02 
#define UNISIZE_D 10.0 // the size of a single dimension of the universe
#define R_C 1.0 

typedef fixap<int32_t, 27> fix_p;
//typedef float fix_p;

const fix_p A[3][3] = {  {fix_p(25.0), fix_p(75.0), fix_p(35.0)}, 
                                     {fix_p(75.0), fix_p(25.0), fix_p(50.0)},  
                                     {fix_p(35.0), fix_p(50.0), fix_p(25.0)}};

Vector3D<fix_p> randPos(unsigned N){
    Vector3D<fix_p> t_pos;
    fix_p x(rand() / (float)RAND_MAX * N);
    fix_p y(rand() / (float)RAND_MAX * N);
    fix_p z(rand() / (float)RAND_MAX * N);
    t_pos.set(x,y,z);
    return t_pos;
}


//! generates a random position within a given space (NxN)
Vector3D<fix_p> rand2DPos(unsigned N){
    Vector3D<fix_p> t_pos;
    fix_p x(rand() / (float)RAND_MAX * N);
    fix_p y(rand() / (float)RAND_MAX * N);
    fix_p z(0.0);
    t_pos.set(x,y,z);
    return t_pos;
}


// conservative pairwise force declaration
void conF(Particle<fix_p> *me, Particle<fix_p> *other){

    
    fix_p a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const fix_p r_c(R_C); // the interaction cutoff

    Vector3D<fix_p> r_i = me->getPos();
    Vector3D<fix_p> r_j = other->getPos();
    fix_p r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<fix_p> r_ij = r_i - r_j;
 
    // Equation 8.5 in the dl_meso manual
    Vector3D<fix_p> force = (r_ij/r_ij_dist) * (a_ij * (fix_p(1.0) - (r_ij_dist/r_c)));

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 

    return;
}

// drag pairwie force declaration
void dragF(Particle<fix_p> *me, Particle<fix_p> *other) {

    
    const fix_p r_c(R_C); // the interaction cutoff
    const fix_p drag_coef(4.5); // the drag coefficient (no idea what to set this at)

    // position and distance
    Vector3D<fix_p> r_i = me->getPos();
    Vector3D<fix_p> r_j = other->getPos();
    fix_p r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<fix_p> r_ij = r_j - r_i; // vector between the points

    // switching function
    fix_p w_d = (fix_p(1.0) - r_ij_dist / r_c)*(fix_p(1.0) - r_ij_dist / r_c);


    // velocities
    Vector3D<fix_p> v_i = me->getVelo(); 
    Vector3D<fix_p> v_j = other->getVelo();
    Vector3D<fix_p> v_ij = v_i - v_j; // relative velocity

    Vector3D<fix_p> force = (r_ij / (r_ij_dist * r_ij_dist)) * w_d * r_ij.dot(v_ij) * (fix_p(-1.0) * drag_coef); 

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
void randF(uint32_t grand, Particle<fix_p> *me, Particle<fix_p> *other) {

   const fix_p K_BT(1.0);
   const fix_p drag_coef(4.5); // the drag coefficient (no idea what to set this at)
   const fix_p sigma_ij(sqrt(2*drag_coef*K_BT)); // the temperature coef
   const fix_p dt(DELTA_T); 
   const fix_p r_c(R_C); // the interaction cutoff

   // position and distance
   Vector3D<fix_p> r_i = me->getPos();
   Vector3D<fix_p> r_j = other->getPos();
   fix_p r_ij_dist = r_i.dist(r_j); // get the distance
   Vector3D<fix_p> r_ij = r_j - r_i; // vector between the points

   // switching function
   fix_p w_r = (fix_p(1.0) - r_ij_dist/r_c);
      
   // random number generation
   fix_p r((pairwise_rand(grand, me->getID(), other->getID()) / (float)(RAND_MAX)) * 0.5);

   // force calculation
   Vector3D<fix_p> force = (r_ij / r_ij_dist)*sqrt(dt)*r*w_r*sigma_ij;  

   // update the forces acting on the two particles
   me->setForce( me->getForce() + force*fix_p(-1.0)); 

   return;
}

// Test program
int main() {
   // size of the universe
   const float unisize = UNISIZE_D;

   // mass of the particles
   const fix_p mass_p0(1.0);

   const unsigned num_cubes = 10;
   const fix_p cube_size(unisize/num_cubes);

   SimSystem<fix_p> universe(fix_p(unisize), fix_p(DELTA_T), fix_p(R_C), num_cubes, 0);

   Particle<fix_p> *p0 = new Particle<fix_p>(rand2DPos(unisize), 0, mass_p0, conF, dragF, randF);
   Vector3D<fix_p> offset(fix_p(0.25),fix_p(0.0),fix_p(0.0));
   Vector3D<fix_p> p1_pos = p0->getPos() + offset;
   Particle<fix_p> *p1 = new Particle<fix_p>(p1_pos, 0, mass_p0, conF, dragF, randF);
   Particle<fix_p> *p2 = new Particle<fix_p>(p1_pos + Vector3D<fix_p>(fix_p(0.0), fix_p(0.75), fix_p(0.0)), 0, mass_p0, conF, dragF, randF);
   universe.addParticle(p0);
   universe.addParticle(p1);
   //universe.addParticle(p2);

   //// add water
   //for(int w=0; w<500; w++){
   //    Particle<fix_p> *p = new Particle<fix_p>(rand2DPos(unisize), 0, mass_p0, conF, dragF, randF);
   //    universe.addParticle(p);
   //}

   //// add orange oil 
   //for(int o_o=0; o_o<300; o_o++){
   //    Particle<fix_p> *p = new Particle<fix_p>(rand2DPos(unisize), 1, mass_p0, conF, dragF, randF);
   //    universe.addParticle(p);
   //}

   //// add green oil 
   //for(int g_o=0; g_o<200; g_o++){
   //    Particle<fix_p> *p = new Particle<fix_p>(rand2DPos(unisize), 2, mass_p0, conF, dragF, randF);
   //    universe.addParticle(p);
   //}

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");
   
   // run the simulation universe for some timesteps 
   universe.run(-1, 0.1);

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");

}
