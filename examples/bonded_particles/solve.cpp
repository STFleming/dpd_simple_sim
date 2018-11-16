#include <iostream>
#include <cstdio>
#include "SimSystem.hpp"
#include "Particle.hpp"
#include "utils.hpp"
#include <random>

#define DELTA_T 0.02 
#define UNISIZE_D 10.0 // the size of a single dimension of the universe
#define R_C 1.0 

const float A[2][2] = { {25.0, 25.0}, {25.0, 25.0}}; // interaction matrix

// conservative pairwise force declaration
void conF(Particle *me, Particle *other){

    const float a_ij = A[me->getType()][other->getType()]; // the interaction strength
    const float r_c = R_C; // the interaction cutoff

    Vector3D r_i = me->getPos();
    Vector3D r_j = other->getPos();
    float r_ij_dist = r_i.dist(r_j); // get the distance
    Vector3D r_ij = r_i - r_j; // vector between the points
 
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

   // // position and distance
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

    return;
}

// random pairwise force declaration
void randF(Particle *me, Particle *other) {

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
      
   // force calculation
   float r = (rand() / (float)RAND_MAX * 1.0);
   Vector3D force = (r_ij / r_ij_dist)*sqrt(dt)*r*w_r*sigma_ij;  

   // update the forces acting on the two particles
   me->setForce( me->getForce() + force); 

   return;
}

// this is a hookean harmonic bond 
void bondF(Particle *me, Particle *other) {

   const float K_div_2 = 128.0/2; // spring constant guess at a value 
   const float r_o = 0.5; // the equilibrium bond length 

   // variables for the current position
   Vector3D r_i = me->getPos();
   Vector3D r_j = other->getPos();
   float r_ij_dist = r_i.dist(r_j); // get the distance
   Vector3D r_ij = r_i - r_j;

   // variables for the previous position
   Vector3D p_r_i = me->getPrevPos();
   Vector3D p_r_j = other->getPrevPos();
   float p_r_ij_dist = p_r_i.dist(p_r_j);
   Vector3D p_r_ij = p_r_i - p_r_j;

   float U_r_ij = K_div_2 * (r_ij_dist - r_o)*(r_ij_dist - r_o); 
   float p_U_r_ij = K_div_2 * (p_r_ij_dist - r_o)*(p_r_ij_dist - r_o); 

   Vector3D b_force;
   if(r_ij_dist == p_r_ij_dist) {
      b_force.set(0.0,0.0,0.0);
   } else {
      b_force = r_ij * (-1.0/r_ij_dist) * ((U_r_ij - p_U_r_ij)/(r_ij_dist - p_r_ij_dist));  
   }

   // update the force values
   me->setForce( me->getForce() + b_force);
}

// Test program
int main() {
   // size of the universe
   const float unisize = UNISIZE_D;

   // number of particles (beads) in the universe
   const unsigned n = 1000;

   // the number of polymers 
   const unsigned nP = 1; 
   // pLen the length of the polymer 
   const unsigned pLen = 10;

   // mass of the particles
   const float mass_w = 1.0;
   const float mass_p = 1.0;

   const unsigned num_cubes = 10;

   SimSystem universe(unisize, DELTA_T, R_C, num_cubes, 0);
   
   // add water
   for(int w=0; w<1000; w++){
       Particle *p = new Particle(randPos(unisize), 0, mass_w, conF, dragF, randF);
       universe.addParticle(p);
   }

   // add a polymer
   // create the polymer at a random position
   Vector3D offset(0.0,0.5,0.0); //how far apart each bond particle is initially placed
   //Vector3D start = randPos(unisize);
   Vector3D start(5.1,5.1,5.1);
   Particle* p0 = new Particle(start, 1, mass_p, conF, dragF, randF); 
   Particle* p1 = new Particle(p0->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
   //Particle* p2 = new Particle(p1->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
   //Particle* p3 = new Particle(p2->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
   //Particle* p4 = new Particle(p3->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
   //Particle* p5 = new Particle(p4->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
   //Particle* p6 = new Particle(p5->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
   //Particle* p7 = new Particle(p6->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
   //Particle* p8 = new Particle(p7->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 
   //Particle* p9 = new Particle(p8->getPos().modulo_add(offset, unisize), 1, mass_p, conF, dragF, randF); 

   // add the polymer to the simulation universe
   universe.addParticle(p0);
   universe.addParticle(p1);
   //universe.addParticle(p2);
   //universe.addParticle(p3);
   //universe.addParticle(p4);
   //universe.addParticle(p5);
   //universe.addParticle(p6);
   //universe.addParticle(p7);
   //universe.addParticle(p8);
   //universe.addParticle(p9);

   // add the bonds between the polymer links
   universe.bond(p0,p1, bondF);
   //universe.bond(p1,p2, bondF);
   //universe.bond(p2,p3, bondF);
   //universe.bond(p3,p4, bondF);
   //universe.bond(p4,p5, bondF);
   //universe.bond(p5,p6, bondF);
   //universe.bond(p6,p7, bondF);
   //universe.bond(p7,p8, bondF);
   //universe.bond(p8,p9, bondF);
   

   // run the universe on a single thread for many timesteps, emitting it's value every 0.07 seconds 
   universe.run(-1, 0.1);

   // emit the initial state (read by the web renderer interface)
   universe.emitJSONFromSU("state.json");

   std::cout << "done\n";
}
