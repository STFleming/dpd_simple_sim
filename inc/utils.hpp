#ifndef __UTILS_H
#define __UTILS_H

#include <stdlib.h> // rand()
//#include <cstdint> // uint16_t
#include <math.h> // sqrt
//#include "Particle.hpp"
#include "Vector3D.hpp"

template<class S>
struct position_t { S x; S y; S z; };

//! generates a random position within a given space (NxNxN)
Vector3D randPos(float N);

//! generates a random position within a given space (NxN)
Vector3D rand2DPos(float N);

//! convert between a Vector3D and a point (position_t)
template<class S>
position_t<S> vec2pos(Vector3D a){
    position_t<S> t = {(S)a.x(), (S)a.y(), (S)a.z()};
    return t;
}

#endif /* __UTILS_H */
