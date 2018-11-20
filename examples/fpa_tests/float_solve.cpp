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

//! generates a random position within a given space (NxNxN)
Vector3D<float> randPos(float N){
    Vector3D<float> t_pos;
    float x = (rand() / (float)RAND_MAX * N);
    float y = (rand() / (float)RAND_MAX * N);
    float z = (rand() / (float)RAND_MAX * N);
    t_pos.set(x,y,z);
    return t_pos;
}

//template<class C, unsigned F>
//Vector3D<fixap<C,F>> randPos(fixap<C,F> N){
//    Vector3D<fixap<C,F>> t_pos;
//    fixap<C,F> x(rand() / (float)RAND_MAX * N);
//    fixap<C,F> y(rand() / (float)RAND_MAX * N);
//    fixap<C,F> z(rand() / (float)RAND_MAX * N);
//    t_pos.set(x,y,z);
//    return t_pos;
//}


//! generates a random position within a given space (NxN)
Vector3D<float> rand2DPos(float N){
    Vector3D<float> t_pos;
    float x = (rand() / (float)RAND_MAX * N);
    float y = (rand() / (float)RAND_MAX * N);
    float z = 0.0;
    t_pos.set(x,y,z);
    return t_pos;
}

//template<class C, unsigned F>
//Vector3D<fixap<C,F>> rand2DPos(fixap<C,F> N){
//    Vector3D<fixap<C,F>> t_pos;
//    fixap<C,F> x(rand() / (float)RAND_MAX * N);
//    fixap<C,F> y(rand() / (float)RAND_MAX * N);
//    fixap<C,F> z(0.0);
//    t_pos.set(x,y,z);
//    return t_pos;
//}


// conservative pairwise force declaration
void conF(Particle<float> *me, Particle<float> *other){

    const float a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const float r_c = R_C; // the interaction cutoff

    Vector3D<float> r_i = me->getPos();
    Vector3D<float> r_j = other->getPos();
    float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<float> r_ij = r_i - r_j;
 
    // Equation 8.5 in the dl_meso manual
    Vector3D<float> force = (r_ij/r_ij_dist) * (a_ij * (1.0 - (r_ij_dist/r_c)));

    assert(r_ij_dist <= R_C);

    // update the forces acting on the two particles
    me->setForce( me->getForce() + force); 
    return;
}

// drag pairwie force declaration
void dragF(Particle<float> *me, Particle<float> *other) {
    
    const float r_c = R_C; // the interaction cutoff
    const float drag_coef = 4.5; // the drag coefficient (no idea what to set this at)

    // position and distance
    Vector3D<float> r_i = me->getPos();
    Vector3D<float> r_j = other->getPos();
    float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<float> r_ij = r_j - r_i; // vector between the points

    // switching function
    float w_d = (1.0 - r_ij_dist / r_c)*(1.0 - r_ij_dist / r_c);

    // velocities
    Vector3D<float> v_i = me->getVelo(); 
    Vector3D<float> v_j = other->getVelo();
    Vector3D<float> v_ij = v_i - v_j; // relative velocity
     
    Vector3D<float> force = (r_ij / (r_ij_dist * r_ij_dist)) * w_d * r_ij.dot(v_ij) * (-1.0 * drag_coef); 

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
void randF(uint32_t grand, Particle<float> *me, Particle<float> *other) {

   const float K_BT = 1.0;
   const float drag_coef = 4.5; // the drag coefficient (no idea what to set this at)
   const float sigma_ij = sqrt(2*drag_coef*K_BT); // the temperature coef
   const float dt = DELTA_T; 
   const float r_c = R_C; // the interaction cutoff

   // position and distance
   Vector3D<float> r_i = me->getPos();
   Vector3D<float> r_j = other->getPos();
   float r_ij_dist = r_i.dist(r_j); // get the distance
   Vector3D<float> r_ij = r_j - r_i; // vector between the points

   // switching function
   float w_r = (1.0 - r_ij_dist/r_c);
      
   // random number generation
   float r = ((pairwise_rand(grand, me->getID(), other->getID()) / (float)(RAND_MAX)) * 0.5);

   // force calculation
   Vector3D<float> force = (r_ij / r_ij_dist)*sqrt(dt)*r*w_r*sigma_ij;  

   // update the forces acting on the two particles
   me->setForce( me->getForce() + force*-1.0); 
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

   SimSystem<float> universe(unisize, DELTA_T, R_C, num_cubes, 0);

   // add water
   for(int w=0; w<500; w++){
       Particle<float> *p = new Particle<float>(rand2DPos(unisize), 0, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // add orange oil 
   for(int o_o=0; o_o<300; o_o++){
       Particle<float> *p = new Particle<float>(rand2DPos(unisize), 1, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // add green oil 
   for(int g_o=0; g_o<200; g_o++){
       Particle<float> *p = new Particle<float>(rand2DPos(unisize), 2, mass_p0, conF, dragF, randF);
       universe.addParticle(p);
   }

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");
   
   // run the simulation universe for some timesteps 
   universe.run(-1, 0.1);

   // emit the state of the simulation
   universe.emitJSONFromSU("state.json");

}
