#include "utils.hpp"

position_t randPos(float N){
    position_t t_pos;
    t_pos.x = (rand() / (float)RAND_MAX * N);
    t_pos.y = (rand() / (float)RAND_MAX * N);
    t_pos.z = (rand() / (float)RAND_MAX * N);
    return t_pos; 
}


// distance function for computing the euclidean distance between two particles
float dist(position_t a, position_t b) {
    // calculate the position
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z));
}
