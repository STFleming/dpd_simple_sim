#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>
#include "fixed_ap.h"

//#define DELTA_T 0.02 
#define DELTA_T 0.01 
#define UNISIZE_D 10.0 // the size of a single dimension of the universe
#define R_C 1.0 

typedef float p_type;
//typedef fixap<int32_t, 27> p_type;

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

const p_type A[2][2] = { {p_type(25.0), p_type(25.0)}, {p_type(25.0), p_type(25.0)}}; // interaction matrix

// conservative pairwise force declaration
void conF(Particle<p_type> *me, Particle<p_type> *other){

    p_type a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const p_type r_c(R_C); // the interaction cutoff

    Vector3D<p_type> r_i = me->getPos();
    Vector3D<p_type> r_j = other->getPos();
    p_type r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D<p_type> r_ij = r_i - r_j; // vector between the points
 
    // Equation 8.5 in the dl_meso manual
    Vector3D<p_type> force = (r_ij/r_ij_dist) * (a_ij * (p_type(1.0) - (r_ij_dist/r_c)));

    // update the forces acting on the two particles
    //me->setForce( me->getForce() + force); 

    return;
}

// drag pairwie force declaration
void dragF(Particle<p_type> *me, Particle<p_type> *other) {
    
    const p_type r_c(R_C); // the interaction cutoff
    const p_type drag_coef(4.5); // the drag coefficient (no idea what to set this at)

   // // position and distance
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
}


// this is a hookean harmonic bond 
void bondF(Particle<p_type> *me, Particle<p_type> *other) {

   const p_type K_div_2(128.0/2); // spring constant guess at a value 
   const p_type r_o(0.5); // the equilibrium bond length 

   // variables for the current position
   Vector3D<p_type> r_i = me->getPos();
   Vector3D<p_type> r_j = other->getPos();
   p_type r_ij_dist = r_i.dist(r_j); // get the distance
   Vector3D<p_type> r_ij = r_i - r_j;

   // variables for the previous position
   Vector3D<p_type> p_r_i = me->getPrevPos();
   Vector3D<p_type> p_r_j = other->getPrevPos();
   p_type p_r_ij_dist = p_r_i.dist(p_r_j);
   Vector3D<p_type> p_r_ij = p_r_i - p_r_j;

   p_type U_r_ij = K_div_2 * (r_ij_dist - r_o)*(r_ij_dist - r_o); 
   p_type p_U_r_ij = K_div_2 * (p_r_ij_dist - r_o)*(p_r_ij_dist - r_o); 

   Vector3D<p_type> b_force;
   if(r_ij_dist == p_r_ij_dist) {
      b_force.set(p_type(0.0),p_type(0.0),p_type(0.0));
   } else {
      b_force = r_ij * (p_type(-1.0)/r_ij_dist) * ((U_r_ij - p_U_r_ij)/(r_ij_dist - p_r_ij_dist));  
   }

   // update the force values
   me->setForce( me->getForce() + b_force);
}

// returns a safe random position that will not break a polymer bond
Vector3D<p_type> safe_water_pos(std::vector<Particle<p_type> *> polymer_beads, float unisize){
    bool safe=false;
    Vector3D<p_type> w_pos;
    while(!safe) {
      safe = true;
      w_pos = rand2DPos(unisize);
      for(auto i=polymer_beads.begin(); i!=polymer_beads.end(); ++i){
          Particle<p_type>* p = *i;
          if(p->getPos().dist(w_pos) <= 0.75) {
             safe = false;
          }  
      }
    }
    return w_pos; 
}

// Test program
int main() {
   // size of the universe
   const float unisize = UNISIZE_D;

   // mass of the particles
   const p_type mass_w = 1.0;
   const p_type mass_p = 1.0;

   const unsigned num_cubes = 10;

   SimSystem<p_type> universe(p_type(unisize), p_type(DELTA_T), p_type(R_C), num_cubes, 0);
   
   std::vector<Particle<p_type> *> polymer_beads;

   // add a polymer
   // create the polymer at a random position
   for(int num_p=0; num_p <5; num_p++) {
       Vector3D<p_type> offset(p_type(0.0),p_type(0.5),p_type(0.0)); //how far apart each bond particle is initially placed
       Vector3D<p_type> start = rand2DPos(unisize);
       Particle<p_type>* p0 = new Particle<p_type>(start, 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p1 = new Particle<p_type>(p0->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p2 = new Particle<p_type>(p1->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p3 = new Particle<p_type>(p2->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p4 = new Particle<p_type>(p3->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p5 = new Particle<p_type>(p4->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p6 = new Particle<p_type>(p5->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p7 = new Particle<p_type>(p6->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p8 = new Particle<p_type>(p7->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
       Particle<p_type>* p9 = new Particle<p_type>(p8->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 

       // keep a list of polymer beads to make sure we don't place a water bead between bonds
       polymer_beads.push_back(p0);
       polymer_beads.push_back(p1);
       polymer_beads.push_back(p2);
       polymer_beads.push_back(p3);
       polymer_beads.push_back(p4);
       polymer_beads.push_back(p5);
       polymer_beads.push_back(p6);
       polymer_beads.push_back(p7);
       polymer_beads.push_back(p8);
       polymer_beads.push_back(p9);

       // add the polymer to the simulation universe
       universe.addParticle(p0);
       universe.addParticle(p1);
       universe.addParticle(p2);
       universe.addParticle(p3);
       universe.addParticle(p4);
       universe.addParticle(p5);
       universe.addParticle(p6);
       universe.addParticle(p7);
       universe.addParticle(p8);
       universe.addParticle(p9);

       // add the bonds between the polymer links
       universe.bond(p0,p1, bondF);
       universe.bond(p1,p2, bondF);
       universe.bond(p2,p3, bondF);
       universe.bond(p3,p4, bondF);
       universe.bond(p4,p5, bondF);
       universe.bond(p5,p6, bondF);
       universe.bond(p6,p7, bondF);
       universe.bond(p7,p8, bondF);
       universe.bond(p8,p9, bondF);
   } 

   // add water
   for(int w=0; w<125; w++){
       // we want to add water, but not in between a polymer bond
       Particle<p_type> *p = new Particle<p_type>(safe_water_pos(polymer_beads, unisize), 0, mass_w, conF, dragF, randF);
       universe.addParticle(p);
   }

   // run the universe on a single thread for many timesteps, emitting it's value every 0.07 seconds 
   universe.run(-1, 0.1);

   // emit the initial state (read by the web renderer interface)
   universe.emitJSONFromSU("state.json");

   std::cout << "done\n";
}
