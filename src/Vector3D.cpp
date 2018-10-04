// implementation file for the Vector3D class
#include "Vector3D.hpp"

// constructor
Vector3D::Vector3D(float x, float y, float z){
   _x = x; 
   _y = y; 
   _z = z;
}

// default constructor
Vector3D::Vector3D() {
  _x = 0.0;
  _y = 0.0;
  _z = 0.0;
}

// destructor
Vector3D::~Vector3D() {

}

// copy constructor
Vector3D::Vector3D(const Vector3D &in) {
   _x = in._x;
   _y = in._y;
   _z = in._z;
}

// clears the vector
void Vector3D::clear() {
   _x = 0.0;
   _y = 0.0;
   _z = 0.0;
}

// getters
float Vector3D::x() { return _x; }
float Vector3D::y() { return _y; }
float Vector3D::z() { return _z; }

// setter
void Vector3D::set(float x, float y, float z) {
    _x = x;
    _y = y;
    _z = z;
}

// operators
//----------------------
// vector * vector
Vector3D Vector3D::operator*(Vector3D a) {
   Vector3D c;
   float x = a.x() * _x;
   float y = a.y() * _y;
   float z = a.z() * _z;
   c.set(x,y,z);
   return c;
}

// vector * scalar
Vector3D Vector3D::operator*(float a) {
   Vector3D c;
   float x = a * _x;
   float y = a * _y;
   float z = a * _z;
   c.set(x,y,z);
   return c;
}
