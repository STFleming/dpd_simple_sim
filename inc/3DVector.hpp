//! A 3D vector class for computing forces etc.. in the simple DPD particle simulator

class 3DVector {
   public:
       3DVector(float x, float y, float z); /**< constructor */
       ~3DVector(); /**< destructor */ 
       3DVector(const 3DVector &in); /**< copy constructor */

       // various getters and setters
       void clear(); /**< clears the current vector values (sets them all to zero) */
       float x(); /**< returns the x value of this vector */
       float y(); /**< returns the x value of this vector */
       float z(); /**< returns the x value of this vector */
       void set(float x, float y, float z); /**< sets the vector value */       

       //         operations 
       // ----------------------------

       // vec and scalar

       // scalar and vec

       // vec and vec

   private:
       float _x;
       float _y;
       float _z; 
};
