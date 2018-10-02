#ifndef __UTILS_H
#define __UTILS_H

#include <stdlib.h> // rand()
//#include <cstdint> // uint16_t

//! A typedef used for the position of a cube center or a particle
typedef struct _position_t { float x; float y; float z; } position_t;

//! generates a random position within a given space (NxNxN)
position_t randPos(float N);

#endif /* __UTILS_H */
