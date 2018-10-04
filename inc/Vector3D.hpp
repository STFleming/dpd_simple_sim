//! A 3D vector class for computing forces etc.. in the simple DPD particle simulator

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
       Vector3D operator*(Vector3D a); /**< multiple a vector to this */
       Vector3D operator*(float a); /**< multiple a scalar to this */

   private:
       float _x;
       float _y;
       float _z; 
};
