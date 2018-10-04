#ifndef __UTILS_H
#define __UTILS_H

#include <stdlib.h> // rand()
//#include <cstdint> // uint16_t
#include <math.h> // sqrt
//#include "Particle.hpp"

//! A typedef used for the position of a cube center or a particle
typedef struct _position_t { float x; float y; float z; } position_t;

//! generates a random position within a given space (NxNxN)
position_t randPos(float N);

//! computes the euclidean distance between two particles
float dist(position_t a, position_t b);

#endif /* __UTILS_H */
