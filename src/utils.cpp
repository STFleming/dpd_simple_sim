#include "utils.hpp"

position_t randPos(float N){
    position_t t_pos;
    t_pos.x = (rand() / (float)RAND_MAX * N) + 1;
    t_pos.y = (rand() / (float)RAND_MAX * N) + 1;
    t_pos.z = (rand() / (float)RAND_MAX * N) + 1;
    return t_pos; 
}

