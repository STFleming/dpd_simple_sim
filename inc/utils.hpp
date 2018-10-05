#ifndef __UTILS_H
#define __UTILS_H

#include <stdlib.h> // rand()
//#include <cstdint> // uint16_t
#include <math.h> // sqrt
//#include "Particle.hpp"
#include "Vector3D.hpp"

typedef struct _position_t { float x; float y; float z; } position_t;

//! generates a random position within a given space (NxNxN)
Vector3D randPos(float N);

//! convert between a Vector3D and a point (position_t)
position_t vec2pos(Vector3D a);

#endif /* __UTILS_H */
