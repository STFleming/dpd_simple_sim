//! A 3D vector class for computing forces etc.. in the simple DPD particle simulator
#ifndef __VECTOR_3D_H
#define __VECTOR_3D_H

#include <math.h>

class Vector3D {
   public:
       Vector3D(float x, float y, float z); /**< constructor */
       Vector3D(); /**< default constructor */
       ~Vector3D(); /**< destructor */ 
       Vector3D(const Vector3D &in); /**< copy constructor */

       // various getters and setters
       void clear(); /**< clears the current vector values (sets them all to zero) */
       float x(); /**< returns the x value of this vector */
       float y(); /**< returns the x value of this vector */
       float z(); /**< returns the x value of this vector */
       void set(float x, float y, float z); /**< sets the vector value */       

       //         operations 
       // ----------------------------
       // multiplication
       Vector3D operator*(Vector3D a); /**< multiple a vector to this */
       Vector3D operator*(float a); /**< multiple a scalar to this */
       // addition
       Vector3D operator+(Vector3D a); /**< add a vector to this */
       Vector3D operator+(float a); /**< add a scalar to this */
       // subtraction 
       Vector3D operator-(Vector3D a); /**< subtract a vector from this */
       Vector3D operator-(float a); /**< subtract a scalar to this */
       // scalar division
       Vector3D operator/(float a); /**< divide the vector by a scalar value */

       // dot product
       float dot(Vector3D a); /**< computes the dot product between this and vector a */

       // cross product
       Vector3D cross(Vector3D a); /**< computes the cross product between this and vector a */
       
       // magnitude
       float mag(); /**< calculates the magnitude of this vector */
   
       // distance
       float dist(Vector3D a); /**< calculates the euclidean distance */

   private:
       float _x;
       float _y;
       float _z; 
};

#endif /* __VECTOR_3D_H */
