#include "SpatialUnit.hpp"

/**! Constructor: creates a problem of NxNxN cubes */
SpatialUnit::SpatialUnit(unsigned size, uint16_t x, uint16_t y, uint16_t z) {
    _size = size;
    _center.x = x;
    _center.y = y;
    _center.z = z;
}

/**! Destructor: cleans up the particles */
SpatialUnit::~SpatialUnit() {

}

/**! returns the position of this cube */
position_t SpatialUnit::getPos() {
    return _center;
}

/**! returns the size of this cube */
unsigned SpatialUnit::getSize() {
    return _size;
}
