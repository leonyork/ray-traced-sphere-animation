// 
// File:   spherelit.h
// Author: leon
//
// Created on 12 July 2007, 12:18
//

#ifndef _SPHERELIT_H
#define	_SPHERELIT_H

#include "kVector.h"

using kBase::v3D;

/**
 * Camera
 * Camera with fixed output width, height and distance from viewer to screen.
 * The camera row by row (i.e. runs horizontal scanlines). Calculates the 
 * directions of light vectors from the viewer through the pixel points (x,y)
 * on construction. So as long as the direction of viewing is not changed,
 * these only need to be computed once.
 * TODO: Add rotations (nb. directly modify lightDirections array by rotating
 * each vector - this will be faster than recalculating them all). Also add
 * errors for 0 direction.
 */
class FixedCamera
{
protected:
    v3D x;          //!< The position of the camera
    v3D d;          //!< The direction it points (cannot be zero)
    unsigned int dvs;        //!< The distance from the viewer to the screen (in pixels)
    unsigned int width;      //!< Width of the resulting image in pixels
    unsigned int height;     //!< Height of the resulting image in pixels
    
    real rotation;          //!< Rotation from horizontal (i.e. the angle the horizontal of the camera makes with x-z plane)
    
    v3D *lightDirections;     //!< An array with size width*height to store the vectors giving the direction light travels through a pixel x,y.
public:
    FixedCamera(const v3D &x, const v3D &d, unsigned int width, unsigned int height, unsigned int dvs, real rotation=0); //!< Set the camera position and direction (should not be zero), the width and the height (in pixels) of the final image, and the distance of the viewer from the screen (in pixels!!!). (Optional: set the rotation of the horizontal axis of the camera from the x-z plane. Default =0)
    
    void moveCameraTo(const v3D &x);    //!< Moves the camera to x (note: v.quick as lightDirections do not need to be computed
    const v3D &getPosition() const;     //!< Get the position of the light
    
    const v3D &getVectorThrough(unsigned int x, unsigned int y) const; //!< Gives the direction of the light passing through the (x,y) pixel of the camera.
    
    unsigned int getWidth() const;          //!< Get the width of the camera output
    unsigned int getHeight() const;         //!< Get the height of the camera output    
    
    friend class FixedSphere;
    friend class FixedSphereTransparent;
};

/**
 * Sphere
 * A Sphere with fixed radius and position. This is done 
 * so cast rays from certain pixels will always hit the same place
 */
class FixedSphere
{
protected:
    v3D x;          //!< Center of the sphere
    real r;         //!< Radius of the sphere
    real *intersectTimes; //!< Array of size width*height to hold the lambda which gives us the length of the vector passing through pixel (x,y) and "landing" on the sphere
    v3D *reflection;      //!< Array of size width*height to hold the reflection vector of light coming from the pixel (x,y) of the camera
    const FixedCamera *const camera;    //!< Pointer to the camera to use    
    FixedSphere(const FixedCamera &camera);                      //!< Here to allow for inheritance    
public:    

    FixedSphere(const v3D &center, real radius, const FixedCamera &camera); //!< Sphere with given radius and centre. Camera at position given in camera and pointing in direction r3D.v   
    
    const v3D& getCentre() const;                                     //!< Return the centre of the sphere
    const real getRadius() const;                                     //!< Return the radius of the sphere
    
    v3D getIntersect(unsigned int x, unsigned int y) const;                  //!< Get the vector that intersects the sphere from the pixel x,y of the camera   
    real getIntersectTime(unsigned int x, unsigned int y) const;            //!< Get the "time" that the light through point (x,y) of the camera screen would have to travel before it hit the sphere,

    const v3D &getReflection(unsigned int x, unsigned int y) const;         //!< Get the reflection vector for light hitting the sphere from the pixel (x,y) of the camera
};

class FixedSphereTransparent: public FixedSphere
{
protected:
    v3D *refraction; //!< Array of size width*height to store the second point on the sphere the light hits if the first ray of light was refracted.
    v3D *refractionReflection;   //!< Array of size width*height to store the reflection (across the sphere norm) of the refraction.
    v3D *refractionRefraction ;      //!< Array of size width*height to store the direction of the light after the second refraction (i.e. coming out the other side of the sphere at the point given in the refraction array)
    v3D *refractionReflectionPoint;     //!< Array of size width*height holding the point the light hits in the sphere after being refracted, then reflected off the inside of the sphere
    v3D *refractionReflectionRefraction;    //!< Array of size width*height holding the direction of the light after refraction, then reflection and refraction to leave the sphere (at point refractionReflectionPoint)

    real *reflectionIntensity;
    real *refractionRefractionIntensity;
    real *refractionReflectionRefractionIntensity;
    
public:
    FixedSphereTransparent(const v3D &center, real radius, const FixedCamera &camera); //!< Sphere with given radius and centre. Camera at position given in camera and pointing in direction r3D.v   
    
    const v3D &getRefractionPoint(unsigned int x, unsigned int y) const;      //!< Get the point on the sphere, the refracted light coming from the (x,y) pixel of the camera will hit. 
    const v3D &getRefractionReflection(unsigned int x, unsigned int y) const; //!< Get the reflection of the light ray after refraction
    const v3D &getRefractionRefraction(unsigned int x, unsigned int y) const; //!< Get the direction of the light ray as it leaves the other end of the sphere (at the point getRefractionPoint()) @see getRefractionPoint()
    const v3D &getRefractionReflectionPoint(unsigned int x, unsigned int y) const; //!< Get the position of the light ray as it leaves after entering the sphere and being reflected.
    const v3D &getRefractionReflectionRefraction(unsigned int x, unsigned int y) const; //!< Get the direction of the light ray as it leaves the sphere after being reflected inside the sphere (at the point getRefractionPoint()) @see getRefractionReflectionPoint()
    
    real getReflectionIntensity(unsigned int x, unsigned int y) const;              //!< Get the intensity of the first reflection
    real getRefractionRefractionIntensity(unsigned int x, unsigned int y) const;  //!< Get the intensity of the ray (between 0.0 and 1.0) after two refractions as it leaves the sphere
    real getRefractionReflectionRefractionIntensity(unsigned int x, unsigned int y) const; //!< Get the intensity of the ray (between 0.0 and 1.0) after a reflection inside the sphere and then a refraction.

    real pathBlocked(const v3D& start, const v3D& finish) const;        //!< Does the line from start to finish intersect the sphere? (False for start or finish within epsilon of the sphere unless the line intersects elsewhere). Note that this used to be boolean, however to get a smoother image, we return a real number which tells us how closely the ray from start to finish was to being blocked
};

#endif	/* _SPHERELIT_H */

