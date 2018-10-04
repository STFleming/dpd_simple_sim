// implementation file for the 3DVector class
#include "3DVector.hpp"

// constructor
3DVector::3DVector(float x, float y, float z){
   _x = x; 
   _y = y; 
   _z = z;
}

// destructor
3DVector::~3DVector() {

}

// copy constructor
3DVector::3DVector(const 3DVector &in) {
   _x = in.x;
   _y = in.y;
   _z = in.z;
}

// clears the vector
void 3DVector::clear() {
   _x = 0.0;
   _y = 0.0;
   _z = 0.0;
}

// getters
float 3DVector::x() { return _x; }
float 3DVector::y() { return _y; }
float 3DVector::z() { return _z; }

// setter
void 3DVector::set(float x, float y, float z) {
    _x = x;
    _y = y;
    _z = z;
}
