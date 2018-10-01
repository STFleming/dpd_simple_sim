#include "utils.hpp"

position_t randPos(uint16_t N){
    position_t t_pos;
    t_pos.x = rand() % N; 
    t_pos.y = rand() % N; 
    t_pos.z = rand() % N; 
    return t_pos; 
}

