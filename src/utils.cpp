#include "utils.hpp"

position_t randPos(float N){
    position_t t_pos;
    t_pos.x = (rand() / (float)RAND_MAX * N);
    t_pos.y = (rand() / (float)RAND_MAX * N);
    t_pos.z = (rand() / (float)RAND_MAX * N);
    return t_pos; 
}

