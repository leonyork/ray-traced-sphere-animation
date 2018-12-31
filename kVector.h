//
// File:   kVector.h
// Author: leon
//
// Created on 08 July 2007, 23:48
//

#ifndef _KVECTOR_H
#define	_KVECTOR_H

#include <iostream>

#define real float
#define EPSILON 1.0e-05

namespace kBase {
    /**
     * 3D Vector
     * Vector with x,y and z components
     */
    class v3D {
    public:
        real x,     //!< x component of vector
                y,  //!< y component of vector
                z;  //!< z component of vector
        
        v3D(real x, real y, real z);    //!< Creates a vector (x,y,z)
        v3D();                          //!< Creates a zero vector
        
        void set(real x, real y, real z);   //!< Sets the vector to (x,y,z)
        void set(const v3D &v);                   //!< Sets the vector to v
        void set();                         //!< Sets the vector to the zero vector
        
        real getSqrLength() const;                //!< Gets the square of the length
        real getLength() const;                   //!< Gets the length of the vector
        bool isZero() const;                      //!< Is this vector within Epsilon of the zero vector?
        
        v3D operator+ (const v3D &v) const;       //!< Adds the vector v to this vector
        void operator+= (const v3D &v);           //!< @see v3D operator+ (v3D &v)
        v3D operator- (const v3D &v) const;       //!< Subtracts the vector v from this vector
        void operator-= (const v3D &v);           //!< @see v3D operator- (v3D &v)        
        void operator*= (real a);                 //!< Multiplies by a scalar
        bool operator== (v3D &v) const;                 //!< Are the two vectors with Epsilon of each other. Does not negate direction.
        
        friend v3D add(const v3D &u, const v3D &v);                //!< Adds the two vectors
        friend v3D scalarMultiply(real a, const v3D &v);     //!< Multiplies the vector by a scalar
        friend real dotProduct(const v3D &u, const v3D &v);        //!< Returns the dot product of two vectors
        friend v3D crossProduct(const v3D &u, const v3D &v);       //!< Returns the cross product of two vectors
        friend v3D componentProduct(const v3D &u, const v3D &v);   //!< Multiplies corresponding components to create a new vector
        
        friend v3D operator* (real a, const v3D &v);         //!< Multiplication by scalars
        friend v3D operator* (const v3D &v, real a);         //!< Multiplication by scalars
        
        friend std::ostream &operator<< (std::ostream &os, const v3D &v);    //!< Output vector in the form (x,y,z)
        friend std::istream &operator>> (std::istream &is, v3D &v);    //!< Input vector in the form x >> y >> z
    };
    
    v3D add(const v3D &u, const v3D &v);                //!< Adds the two vectors
    v3D scalarMultiply(real a, const v3D &v);     //!< Multiplies the vector by a scalar
    real dotProduct(const v3D &u, const v3D &v);        //!< Returns the dot product of two vectors
    v3D crossProduct(const v3D &u, const v3D &v);       //!< Returns the cross product of two vectors
    v3D componentProduct(const v3D &u, const v3D &v);   //!< Multiplies corresponding components to create a new vector
    
    
    v3D operator* (real a, const v3D &v);         //!< Multiplication by scalars
    v3D operator* (const v3D &v, real a);         //!< Multiplication by scalars
    
    std::ostream &operator<< (std::ostream &os, const v3D &v);    //!< Output vector in the form (x,y,z)
    std::istream &operator>> (std::istream &is, v3D &v);    //!< Input vector in the form x >> y >> z
    
    /**
     * 2D Vector
     * Vector with x and y components
     */
    class v2D {
    public:
        real x,     //!< x component of vector
                y;  //!< y component of vector
        
        v2D(real x, real y);            //!< Creates a vector (x,y)
        v2D();                          //!< Creates a zero vector
        
        void set(real x, real y);           //!< Sets the vector to (x,y)
        void set(const v2D &v);                   //!< Sets the vector to v
        void set();                         //!< Sets the vector to the zero vector
        
        real getSqrLength() const;          //!< Gets the square of the length
        real getLength() const;             //!< Gets the length of the vector
        bool isZero() const;                //!< Is the vector within Epsilon of the zero vector
        
        v2D operator+ (const v2D &v) const;             //!< Adds the vector v to this vector
        void operator+= (const v2D &v);           //!< @see v2D operator+ (v2D &v)
        v2D operator- (const v2D &v) const;             //!< Subtracts the vector v from this vector
        void operator-= (const v2D &v);           //!< @see v2D operator- (v2D &v)        
        void operator*= (real a);           //!< Multiplies by a scalar
        bool operator== (v2D &v) const;
        
        friend v2D add(const v2D &u, const v2D &v);                //!< Adds the two vectors
        friend v2D scalarMultiply(real a, const v2D &v);     //!< Multiplies the vector by a scalar
        friend real dotProduct(const v2D &u, const v2D &v);        //!< Returns the dot product of two vectors
        friend real determinant(const v2D &u, const v2D &v);       //!< Returns the determinant of the matrix (u v) (u,v are column vectors).
        friend v2D componentProduct(const v2D &u, const v2D &v);   //!< Multiplies corresponding components to create a new vector
        
        friend v2D operator* (real a, const v2D &v);         //!< Multiplication by scalars
        friend v2D operator* (const v2D &v, real a);         //!< Multiplication by scalars
        
        friend std::ostream &operator<< (std::ostream &os, const v2D &v);    //!< Output vector in the form (x,y)
        friend std::istream &operator>> (std::istream &is, v2D &v);    //!< Input vector in the form x >> y
    };
    
    v2D add(const v2D &u, const v2D &v);                //!< Adds the two vectors
    v2D scalarMultiply(real a, const v2D &v);     //!< Multiplies the vector by a scalar
    real dotProduct(const v2D &u, const v2D &v);        //!< Returns the dot product of two vectors
    real determinant(const v2D &u, const v2D &v);       //!< Returns the determinant of the matrix (u v) (u,v are column vectors).
    v2D componentProduct(const v2D &u, const v2D &v);   //!< Multiplies corresponding components to create a new vector
    
    v2D operator* (real a, const v2D &v);         //!< Multiplication by scalars
    v2D operator* (const v2D &v, real a);         //!< Multiplication by scalars
    
    std::ostream &operator<< (std::ostream &os, const v2D &v);    //!< Output vector in the form (x,y)
    std::istream &operator>> (std::istream &is, const v2D &v);    //!< Input vector in the form x >> y
    
    /**
     * 3D Ray
     *  Structure with displacement and velocity (e.g a Ray of light)
     */
    class r3D {
    public:
        v3D x,      //!< Displacement of Ray
                v;  //!< Velocity of Ray
        
        r3D();      //!< Creates ray with zero displacement and velocity
        r3D(v3D x, v3D v);  //!< Creates ray with displacement x and velocity y
        
        void set();                     //!< Sets the displacement and velocity of the ray to zero
        void set(const v3D &x, const v3D &v);         //!< Sets the displacement and velocity of the ray
        void set(const r3D &r);         //!< Sets the ray to be the same as another ray
        
        r3D operator+ (real t) const;         //!< Advances the ray's displacement by applying the velocity over a time t
        void operator+= (real t);       //!< Advances the ray's displacement by applying the velocity over a time t
        
        friend std::ostream &operator<< (std::ostream &os, const r3D &r);    //!< Output ray in the form [displacement,velocity]
        friend std::istream &operator>> (std::istream &is, r3D &r);    //!< Input ray in the form x.x >> x.y >> x.z >> v.x >> v.y >> v.z
        
    };
    
    std::ostream &operator<< (std::ostream &os, const r3D &r);    //!< Output ray in the form [displacement -> velocity]
    std::istream &operator>> (std::istream &is, r3D &r);    //!< Input ray in the form x.x >> x.y >> x.z >> v.x >> v.y >> v.z
    
    /**
     * PointStack
     * A stack of at most maxLength (specified in constructor) 3D points in order of addition.
     * If more than maxLength points are added, then older points are deleted to make room
     * for the new points.
     */
    class kPointStack {
    protected:
        v3D* points;
        unsigned char length;
        unsigned char maxLength;
        unsigned char curPos;
        
    public:
        kPointStack(unsigned char maxLength);                              //!< Set up stack with no points, but which will have maxLength number of points
        kPointStack(unsigned char maxLength, const v3D &firstPoint);       //!< Set up stack with the first point, with a maximum of maxLength number of points.
        kPointStack(unsigned char maxLength, const v3D *points, unsigned char numOfPoints);           //!< Set up stack with the points specified and a maximum of maxLength number of points.
        
        void newPoint(const v3D &point);    //!< Add a point
        
        unsigned char getLength() const;    //!< Get current number of points in stack
        unsigned char getMaxLength() const; //!< Get maximum number of points stack can have
        v3D getPoint(unsigned char n) const;       //!< Get the n'th point (with 0 being the newest, MaxLength being the oldest)
    };
    
    /**
     * PointStack
     * A stack of at most maxLength (specified in constructor) 2D points in order of addition.
     * If more than maxLength points are added, then older points are deleted to make room
     * for the new points.
     */
    class kPointStack2D {
    protected:
        v2D* points;
        unsigned char length;
        unsigned char maxLength;
        unsigned char curPos;
        
    public:
        kPointStack2D(unsigned char maxLength);                              //!< Set up stack with no points, but which will have maxLength number of points
        kPointStack2D(unsigned char maxLength, const v2D &firstPoint);       //!< Set up stack with the first point, with a maximum of maxLength number of points.
        kPointStack2D(unsigned char maxLength, const v2D *points, unsigned char numOfPoints);           //!< Set up stack with the points specified and a maximum of maxLength number of points.
        
        void newPoint(const v2D &point);    //!< Add a point
        
        unsigned char getLength() const;    //!< Get current number of points in stack
        unsigned char getMaxLength() const; //!< Get maximum number of points stack can have
        v2D getPoint(unsigned char n) const;       //!< Get the n'th point (with 0 being the newest, MaxLength being the oldest)
    }; 
    
    /**
     * Particle
     * A particle with constant acceleration
     */    
    class kParticle
    {
    public:
        v3D x;  //!< Displacement
        v3D v;  //!< Velocity
        v3D a;  //!< Acceleration        
    
        kParticle();    //!< Constructs Particle with initial displacement, velocity and acceleration zero
        kParticle(const v3D &x);   //!< Constructs Particle with initial displacement x, velocity and acceleration zero
        kParticle(const v3D &x, const v3D &v); //!< Constructs Particle with initial displacement x, velocity v and acceleration zero
        kParticle(const v3D &x, const v3D &v, const v3D &a);  //!< Constructs Particle with initial displacement x, velocity v and acceleration a
        
        kParticle operator+(real t) const;        //!< Move the particle through the timestep t
        void operator+=(real t);            //!< Move the particle through the timestep t
    };
    
}

#endif	/* _KVECTOR_H */

