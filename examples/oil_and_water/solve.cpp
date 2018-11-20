// Oil and water DPD simulation
// 3 bead types: water, red oil, and yellow oil

#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>

#define DELTA_T 0.02 
//#define DELTA_T 0.001 
#define UNISIZE_D 10.0 // the size of a single dimension of the universe
#define R_C 1.0 

const float A[3][3] = { {25.0, 75.0, 35.0}, 
                        {75.0, 25.0, 50.0},  
                        {35.0, 50.0, 25.0}}; // interaction matrix

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){

    const float a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const float r_c = R_C; // the interaction cutoff

    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D r_ij = r_i - r_j;
 
    // Equation 8.5 in the dl_meso manual
    Vector3D force = (r_ij/r_ij_dist) * (a_ij * (1.0 - (r_ij_dist/r_c)));

    assert(r_ij_dist <= R_C);

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    return;
}

// drag pairwie force declaration
void dragF(Particle *me, Particle *other) {
    
    const float r_c = R_C; // the interaction cutoff
    const float drag_coef = 4.5; // the drag coefficient (no idea what to set this at)

    // position and distance
    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D r_ij = r_j - r_i; // vector between the points

    // switching function
    float w_d = (1.0 - r_ij_dist / r_c)*(1.0 - r_ij_dist / r_c);

    // velocities
    Vector3D v_i = me->getVelo(); 
    Vector3D v_j = other->getVelo();
    Vector3D v_ij = v_i - v_j; // relative velocity
     
    Vector3D force = (r_ij / (r_ij_dist * r_ij_dist)) * w_d * r_ij.dot(v_ij) * (-1.0 * drag_coef); 

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
void randF(uint32_t grand, Particle *me, Particle *other) {

   const float K_BT = 1.0;
   const float drag_coef = 4.5; // the drag coefficient (no idea what to set this at)
   const float sigma_ij = sqrt(2*drag_coef*K_BT); // the temperature coef
   const float dt = DELTA_T; 
   const float r_c = R_C; // the interaction cutoff

   // position and distance
   Vector3D r_i = me->getPos();
   Vector3D r_j = other->getPos();
   float r_ij_dist = r_i.dist(r_j); // get the distance
   Vector3D r_ij = r_j - r_i; // vector between the points

   // switching function
   float w_r = (1.0 - r_ij_dist/r_c);
      
   // random number generation
   float r = ((pairwise_rand(grand, me->getID(), other->getID()) / (float)(RAND_MAX)) * 0.5);

   // force calculation
   Vector3D force = (r_ij / r_ij_dist)*sqrt(dt)*r*w_r*sigma_ij;  

   //printf("p:%d <-> p:%d  r=%u   f=%s\n", me->getID(), other->getID(), pairwise_rand(grand, me->getID(), other->getID()), force.str().c_str()); 

   // update the forces acting on the two particles
   me->setForce( me->getForce() + force*-1.0); 
   //printf("p:%d randF new force=%s\n", me->getID(), me->getForce().str().c_str());
   return;
}

// Test program
int main() {
   // size of the universe
   const float unisize = UNISIZE_D;

   // mass of the particles
   const float mass_p0 = 1.0;

   const unsigned num_cubes = 2;
   const float cube_size = unisize/num_cubes;

   SimSystem universe(unisize, DELTA_T, R_C, num_cubes, 0);

   // add water
   for(int w=0; w<600; w++){
       Particle *p = new Particle(randPos(unisize), 0, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // add orange oil 
   for(int o_o=0; o_o<300; o_o++){
       Particle *p = new Particle(randPos(unisize), 1, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // add green oil 
   for(int g_o=0; g_o<100; g_o++){
       Particle *p = new Particle(randPos(unisize), 2, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");
   
   // run the simulation universe for some timesteps 
   universe.run(-1, 0.25);

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");

}
