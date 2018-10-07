#include "utils.hpp"

Vector3D randPos(float N){
    Vector3D t_pos;
    float x = (rand() / (float)RAND_MAX * N);
    float y = (rand() / (float)RAND_MAX * N);
    float z = (rand() / (float)RAND_MAX * N);
    t_pos.set(x,y,z);
    return t_pos; 
}

Vector3D rand2DPos(float N){
    Vector3D t_pos;
    float x = (rand() / (float)RAND_MAX * N);
    float y = (rand() / (float)RAND_MAX * N);
    float z = 0.0;
    t_pos.set(x,y,z);
    return t_pos; 
}

// conversion function
position_t vec2pos(Vector3D a) {
    position_t t = {a.x(), a.y(), a.z()};
    return t;
}
