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
   float x =  _x * a.x();
   float y =  _y * a.y();
   float z =  _z * a.z();
   c.set(x,y,z);
   return c;
}

// vector * scalar
Vector3D Vector3D::operator*(float a) {
   Vector3D c;
   float x = _x * a;
   float y = _y * a;
   float z = _z * a;
   c.set(x,y,z);
   return c;
}

// vector + vector
Vector3D Vector3D::operator+(Vector3D a) {
   Vector3D c;
   float x =  _x  + a.x();
   float y =  _y  + a.y();
   float z =  _z  + a.z();
   c.set(x,y,z);
   return c;
}

// vector + scalar
Vector3D Vector3D::operator+(float a) {
   Vector3D c;
   float x = _x + a;
   float y = _y + a;
   float z = _z + a;
   c.set(x,y,z);
   return c;
}

// vector - vector
Vector3D Vector3D::operator-(Vector3D a) {
   Vector3D c;
   float x = _x - a.x();
   float y = _y - a.y();
   float z = _z - a.z();
   c.set(x,y,z);
   return c;
}

// vector - scalar
Vector3D Vector3D::operator-(float a) {
   Vector3D c;
   float x = _x - a;
   float y = _y - a;
   float z = _z - a;
   c.set(x,y,z);
   return c;
}

// vector / scalar
Vector3D Vector3D::operator/(float a) {
   Vector3D c;
   float x = _x / a;
   float y = _y / a;
   float z = _z / a;
   c.set(x,y,z);
   return c;
}

// dot product
float Vector3D::dot(Vector3D a) {
    return (_x * a.x()) + (_y * a.y()) + (_z * a.z());
}

// cross product
Vector3D Vector3D::cross(Vector3D a){
   Vector3D c;
   float x = (_y*a.z()) - (_z*a.y()); 
   float y = (_z*a.x()) - (_x*a.z()); 
   float z = (_x*a.y()) - (_y*a.x()); 
   c.set(x,y,z);
   return c;
}

// mag
float Vector3D::mag(){
   return sqrt(_x*_x + _y*_y + _z*_z);
} 

// dist
float Vector3D::dist(Vector3D a) {
   Vector3D c = *this - a; 
   return c.mag();
}
