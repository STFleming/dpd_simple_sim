#ifndef __UTILS_H
#define __UTILS_H

#include <stdlib.h> // rand()
//#include <cstdint> // uint16_t
#include <math.h> // sqrt
//#include "Particle.hpp"
#include "Vector3D.hpp"
#include "fixed_ap.h"

template<class S>
struct position_t { S x; S y; S z; };

////! generates a random position within a given space (NxNxN)
//Vector3D<float> randPos(float N){ 
//    Vector3D<float> t_pos;
//    float x = (rand() / (float)RAND_MAX * N); 
//    float y = (rand() / (float)RAND_MAX * N); 
//    float z = (rand() / (float)RAND_MAX * N); 
//    t_pos.set(x,y,z);
//    return t_pos; 
//}
//
//template<class C, unsigned F>
//Vector3D<fixap<C,F>> randPos(fixap<C,F> N){ 
//    Vector3D<fixap<C,F>> t_pos;
//    fixap<C,F> x(rand() / (float)RAND_MAX * N); 
//    fixap<C,F> y(rand() / (float)RAND_MAX * N); 
//    fixap<C,F> z(rand() / (float)RAND_MAX * N); 
//    t_pos.set(x,y,z);
//    return t_pos; 
//}
//
//
////! generates a random position within a given space (NxN)
//Vector3D<float> rand2DPos(float N){ 
//    Vector3D<float> t_pos;
//    float x = (rand() / (float)RAND_MAX * N); 
//    float y = (rand() / (float)RAND_MAX * N); 
//    float z = 0.0;
//    t_pos.set(x,y,z);
//    return t_pos; 
//}
//
//template<class C, unsigned F>
//Vector3D<fixap<C,F>> rand2DPos(fixap<C,F> N){ 
//    Vector3D<fixap<C,F>> t_pos;
//    fixap<C,F> x(rand() / (float)RAND_MAX * N); 
//    fixap<C,F> y(rand() / (float)RAND_MAX * N); 
//    fixap<C,F> z(0.0);
//    t_pos.set(x,y,z);
//    return t_pos; 
//}

//! convert between a Vector3D and a point (position_t)
template<class S>
position_t<S> vec2pos(Vector3D<S> a){
    position_t<S> t = {(S)a.x(), (S)a.y(), (S)a.z()};
    return t;
}

#include "../src/utils.cpp"

#endif /* __UTILS_H */
