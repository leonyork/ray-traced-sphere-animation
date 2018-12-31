#include "spherelit.h"
#include "kVector.h"
#include <math.h>

using kBase::v3D;
using kBase::dotProduct;
using kBase::crossProduct;

/*************************************************
 * Fixed Camera
 *************************************************/

FixedCamera::FixedCamera(const v3D &x, const v3D &d, unsigned int width, unsigned int height, unsigned int dvs, real rotation):
    x(x), d(d), width(width), height(height), dvs(dvs), rotation(rotation) {
        // TODO: produce error
        if (d.isZero()) {
            this->d = v3D(1, 1, 1);
        }
        
        //ensure d is a unit vector
        this->d*=1/(this->d).getLength();
        
        // Create the array for the light directions
        lightDirections = new v3D[width*height];
        
        // We first create two axis that we use as our x and y axis for the image
        // These are the vectors before a rotation of [rotation].
        
        v3D xVector =  crossProduct(v3D(0, 1, 0), d);
        v3D yVector =  crossProduct(d, xVector);
        
        real sintheta = sin(rotation);
        real costheta = cos(rotation);
        
        // We apply the rotation about the d vector of xVector and yVector
        // Note that both X and Y will be unit vectors, so we let 1.0 float
        // correspond to 1 pixel.
        v3D X = sintheta*yVector + costheta*xVector;
        v3D Y = costheta*yVector - sintheta*xVector;
        
        
        // Go through each pixel and create a vector from -d*dvs to xi*X+yi*Y
        // (note: we forget about x, as it doesn't matter where the camera is,
        // the directions will be identical)
        
        //std::cout << X << Y << d;
        
        // The X ,Y and Z components of the vector wrt X, Y and d
        v3D Ycomp;
        v3D Xcomp;
        v3D Zcomp = (real)dvs*d;
        
        // Centre of the camera screen
        real cX =  ((real)(width+1))/2.0;
        real cY =  ((real)(height+1))/2.0;
        
        // Number which gives us the index number for the lightDirections array for the row we are on
        unsigned int yTimesWidth;
        
        for (unsigned int yi=0; yi<height; yi++) {
            Ycomp = ((cY-(real)yi))*Y;
            yTimesWidth = width*yi;
            
            for (unsigned int xi=0; xi<width; xi++) {
                Xcomp = (((real)xi-cX))*X;
                lightDirections[yTimesWidth+xi]=Xcomp+Ycomp+Zcomp;
                
                // Make the vector a unit vector (possibly not necessary)
                lightDirections[yTimesWidth+xi]*=1/lightDirections[yTimesWidth+xi].getLength();
            }
        }
        //std::cout << lightDirections[20*150+20];
    }
    
// Don't need to change the direction vectors for this one.
    void FixedCamera::moveCameraTo(const v3D &x) {
        this->x=x;
    }
    
    const v3D &FixedCamera::getPosition() const {
        return x;
    }
    
// TODO: Add checking for x,y within range
    const v3D &FixedCamera::getVectorThrough(unsigned int x, unsigned int y) const {
        return lightDirections[y*width+x];
    }
    
    unsigned int FixedCamera::getWidth() const {
        return width;
    }
    
    unsigned int FixedCamera::getHeight() const {
        return height;
    }
    
    /**************************************************************
     * Fixed Sphere. A sphere which is meant to stay fixed
     * TODO: Make reflection a smaller array by only making it
     * the size it needs to be - only add vectors if we actually
     * hit the sphere. Currently reflection is an array of size
     * width*height but it could be smaller.
     **************************************************************/
    
    FixedSphere::FixedSphere(const FixedCamera &camera)
    : camera(&camera) {
        
    }
    
    FixedSphere::FixedSphere(const v3D &center, real radius, const FixedCamera &camera)
    :x(center), r(radius), camera(&camera) {
        
        intersectTimes = new real[camera.width*camera.height];
        reflection = new v3D[camera.width*camera.height];
        
        real vdotxsqr;
        real vdotx;
        real vsqr;
        real sqroot;
        unsigned int yTimesWidth;
        v3D v;
        v3D xp = x-camera.getPosition();
        
        // Vector for the norm vector from the sphere
        v3D norm;
        
        for (unsigned int yi = 0; yi<camera.height; yi++) {
            yTimesWidth=yi*camera.width;
            for (unsigned int xi = 0; xi<camera.width; xi++) {
                v=camera.getVectorThrough(xi, yi);
                
                vdotx=kBase::dotProduct(v, xp);
                vdotxsqr=vdotx*vdotx;
                
                vsqr=v.getSqrLength();
                sqroot = vsqr*(r*r-xp.getSqrLength())+(vdotxsqr);
                if (sqroot >= EPSILON) {
                    if (vdotxsqr-sqroot >= EPSILON) {
                        intersectTimes[yTimesWidth+xi]=(vdotx-sqrt(sqroot))/vsqr;
                    }
                    else {
                        intersectTimes[yTimesWidth+xi]=(vdotx+sqrt(sqroot))/vsqr;
                    }
                    
                    // This is the calculation of the reflection vector
                    
                    // Find the normal and make it a unit vector
                    norm = (v*intersectTimes[yTimesWidth+xi])+camera.getPosition()-x;
                    norm *= 1/norm.getLength();
                    
                    // Calculate the vector and normalise it
                    reflection[yTimesWidth+xi]= -2*(dotProduct(v, norm))*norm+v;
                    reflection[yTimesWidth+xi]*=1/reflection[yTimesWidth+xi].getLength();
                }
                else {
                    intersectTimes[yTimesWidth+xi]=0.0;
                    
                }
            }
        }
        std::cout << reflection[75*150+75] << "\n";
    }
    
    const v3D& FixedSphere::getCentre() const {
        return x;
    }
    
    const real FixedSphere::getRadius() const {
        return r;
    }
    
    v3D FixedSphere::getIntersect(unsigned int x, unsigned int y) const {
        return intersectTimes[y*camera->width+x]*camera->getVectorThrough(x, y);
    }
    
    real FixedSphere::getIntersectTime(unsigned int x, unsigned int y) const {
        return intersectTimes[y*camera->width+x];
    }
    
    const v3D &FixedSphere::getReflection(unsigned int x, unsigned int y) const {
        return reflection[y*camera->width+x];
    }
    
    
    /************************************************
     * Transparent Sphere
     * A fixed sphere with calculations for the
     * refraction of light through the sphere
     ************************************************/
    
    FixedSphereTransparent::FixedSphereTransparent(const v3D &center, real radius, const FixedCamera &camera)
    :FixedSphere(camera) {
        x=center;
        r=radius;
        intersectTimes = new real[camera.width*camera.height];
        reflection = new v3D[camera.width*camera.height];
        refraction = new v3D[camera.width*camera.height];
        refractionReflection = new v3D[camera.width*camera.height];
        refractionRefraction = new v3D[camera.width*camera.height];
        refractionReflectionPoint = new v3D[camera.width*camera.height];
        refractionReflectionRefraction = new v3D[camera.width*camera.height];
        
        reflectionIntensity = new real[camera.width*camera.height];
        refractionRefractionIntensity = new real[camera.width*camera.height];
        refractionReflectionRefractionIntensity = new real[camera.width*camera.height];
        
        real vdotxsqr;
        real vdotx;
        real vsqr;
        real sqroot;
        unsigned int yTimesWidth;
        v3D v;
        
        // dotProduct((xp+lambda*v),(xp+lambda*v))=r*r
        v3D xp = camera.getPosition()-x;
        
        // Vector for the norm vector from the sphere
        v3D norm;
        
        // dotProduct((x0+lambda*v),(x0+lambda*v))=r*r
        v3D x0;
        
        // Useful for calculating refraction
        real cost1, cost2;
        
        const real n1 = 1.0;
        const real n2 = 1.5;
        const real n1_n2=  n1/n2;     // refractive index ratio (air approx 1, glass approx 1.5 => ratio = 2/3);
        const real n2_n1= n2/n1;
        const real n1minusn2=n1-n2;
        const real n2minusn1=n2-n1;
        
        real intensity;
        
        for (unsigned int yi = 0; yi<camera.height; yi++) {
            yTimesWidth=yi*camera.width;
            for (unsigned int xi = 0; xi<camera.width; xi++) {
                v=camera.getVectorThrough(xi, yi);
                
                vdotx=dotProduct(v, xp);
                vdotxsqr=vdotx*vdotx;
                
                vsqr=v.getSqrLength();
                sqroot = vsqr*(r*r-xp.getSqrLength())+(vdotxsqr);
                if (sqroot >= 0) {
                    sqroot=sqrt(sqroot);
                    if (vdotx+sqroot < EPSILON) {
                        intersectTimes[yTimesWidth+xi]=-1*(vdotx+sqroot)/vsqr;
                    }
                    else {
                        intersectTimes[yTimesWidth+xi]=(sqroot-vdotx)/vsqr;
                    }
                    
                    refractionRefractionIntensity[yTimesWidth+xi]=1.0; //Start with intensity equal to 1
                    refractionReflectionRefractionIntensity[yTimesWidth+xi]=1.0;
                    
                    // This is the calculation of the reflection vector
                    v*=1/v.getLength();
                    
                    // Find the normal and make it a unit vector
                    norm = (v*intersectTimes[yTimesWidth+xi])+camera.x-this->x;
                    norm *= 1/this->r;
                    
                    // Calculate the vector and normalise it
                    reflection[yTimesWidth+xi]= -2*(dotProduct(v, norm))*norm+v;
                    
                    // This is the calculation of the refraction vectors
                    
                    //x0 relative to the camera
                    x0 = (v*intersectTimes[yTimesWidth+xi])+xp;
                    
                    // Calculations for refraction
                    
                    cost1 = 1*dotProduct(v, norm);
                    cost2 = sqrt(1-n1_n2*n1_n2*(1-cost1*cost1));
                    
                    // Reuse v as the new v, i.e. the refracted light ray
                    v=n1_n2*v+(cost2-n1_n2*cost1)*norm;
                    
                    // Work out what intensity of light is refracted (this is updated again later)
                    intensity=0.5*(n1minusn2*cost1+n1minusn2*cost2)/(n1*cost1+n2*cost2);
                    if (intensity<0) {
                        intensity*=-1;
                    }
                    reflectionIntensity[yTimesWidth+xi]=intensity;
                    refractionRefractionIntensity[yTimesWidth+xi]*=1.0-intensity;
                    
                    // Calculate where it hits the sphere again
                    
                    vdotx=dotProduct(v, x0);
                    vdotxsqr=vdotx*vdotx;
                    
                    vsqr=v.getSqrLength();
                    sqroot = vsqr*(r*r-x0.getSqrLength())+(vdotxsqr);
                    
                    
                    // We know the light is definitly going to hit the sphere again
                    // So there's no need to check the square root is postive
                    // We also know it's going to be the larger of the two solutions
                    // as the other one is the first point the light hit the sphere
                    // (i.e. lambda = 0)
                    
                    refraction[yTimesWidth+xi]=(sqrt(sqroot)-vdotx)/vsqr*v+x0;
                    
                    norm = refraction[yTimesWidth+xi] - x;
                    norm *= 1/norm.getLength();
                    
                    refractionReflection[yTimesWidth+xi]= -2*(dotProduct(v, norm))*norm+v;
                    refractionReflection[yTimesWidth+xi]*= 1/refractionReflection[yTimesWidth+xi].getLength();
                    
                    
                    // Calculations for refraction
                    v*=1/v.getLength();
                    
                    cost1 = dotProduct(v, norm);
                    cost2 = sqrt(1-n2_n1*n2_n1*(1-cost1*cost1));
                    
                    // Reuse v as the new v, i.e. the refracted light ray
                    v=n2_n1*v-(cost2+n2_n1*cost1)*norm;
                    
                    refractionRefraction[yTimesWidth+xi] = v;
                    refractionRefraction[yTimesWidth+xi]*= 1/refractionRefraction[yTimesWidth+xi].getLength();
                    
                    // Work out what intensity of light is refracted (this is updated again later)
                    intensity=0.5*(n2minusn1*cost1+n2minusn1*cost2)/(n2*cost1+n1*cost2);
                    if (intensity<0) {
                        intensity*=-1;
                    }
                    refractionRefractionIntensity[yTimesWidth+xi]*=1.0-intensity;
                    refractionReflectionRefractionIntensity[yTimesWidth+xi]*=intensity;
                    
                    // Now we work out where the reflected light goes
                    x0=refraction[yTimesWidth+xi]-x;
                    v=refractionReflection[yTimesWidth+xi];
                    
                    vdotx=dotProduct(v, x0);
                    vdotxsqr=vdotx*vdotx;
                    
                    vsqr=v.getSqrLength();
                    sqroot = vsqr*(r*r-x0.getSqrLength())+(vdotxsqr);
                    
                    refractionReflectionPoint[yTimesWidth+xi]=(vdotx+sqrt(sqroot))/vsqr*v+x0;
                    
                    // Calculate the refraction
                    
                    norm = refractionReflectionPoint[yTimesWidth+xi]-x;
                    norm *= 1/norm.getLength();
                    
                    cost1 = dotProduct(v, norm);
                    cost2 = sqrt(1-n2_n1*n2_n1*(1-cost1*cost1));
                    
                    v=n2_n1*v-(cost2+n2_n1*cost1)*norm;
                    
                    
                    refractionReflectionRefraction[yTimesWidth+xi] = v;
                    refractionReflectionRefraction[yTimesWidth+xi]*=1/refractionReflectionRefraction[yTimesWidth+xi].getLength();
                    
                    
                    // The intensity of this light:
                    intensity=0.5*(n2minusn1*cost1+n2minusn1*cost2)/(n2*cost1+n1*cost2);
                    if (intensity<0) {
                        intensity*=-1;
                    }
                    
                    
                    refractionReflectionRefractionIntensity[yTimesWidth+xi]*=1.0-intensity;
                }
                else {
                    intersectTimes[yTimesWidth+xi]=0.0;
                    /*                reflection[yTimesWidth+xi]=v3D();
                refractionReflection[yTimesWidth+xi]=v3D();
                refractionRefraction[yTimesWidth+xi]=v3D();    */
                }
                
                
            }
        }
        int a=75;
        int b=30;
        std::cout << intersectTimes[b*150+a]*camera.getVectorThrough(a,b)+camera.getPosition() << "\n";
        std::cout << reflection[b*150+a] << "\n";
        std::cout << refraction[b*150+a] << "\n";
        std::cout << refractionRefraction[b*150+a] << "\n";
        std::cout << refractionReflection[b*150+a] << "\n";
        std::cout << refractionReflectionPoint[b*150+a] << "\n";
        std::cout << refractionReflectionRefraction[b*150+a] << "\n";
        std::cout << refractionRefractionIntensity[b*150+a] << "\n";
    }
    
    const v3D &FixedSphereTransparent::getRefractionPoint(unsigned int x, unsigned int y) const {
        return refraction[y*camera->width+x];
    }
    
    const v3D &FixedSphereTransparent::getRefractionReflection(unsigned int x, unsigned int y) const {
        return refractionReflection[y*camera->width+x];
    }
    
    const v3D &FixedSphereTransparent::getRefractionRefraction(unsigned int x, unsigned int y) const {
        return refractionRefraction[y*camera->width+x];
    }
    
    const v3D &FixedSphereTransparent::getRefractionReflectionPoint(unsigned int x, unsigned int y) const {
        return refractionReflectionPoint[y*camera->width+x];
    }
    
    const v3D &FixedSphereTransparent::getRefractionReflectionRefraction(unsigned int x, unsigned int y) const {
        return refractionReflectionRefraction[y*camera->width+x];
    }
    
    real FixedSphereTransparent::getReflectionIntensity(unsigned int x, unsigned int y) const
    {
        return reflectionIntensity[y*camera->width+x];
    }
    
    real FixedSphereTransparent::getRefractionRefractionIntensity(unsigned int x, unsigned int y) const {
        return refractionRefractionIntensity[y*camera->width+x];
    }
    
    real FixedSphereTransparent::getRefractionReflectionRefractionIntensity(unsigned int x, unsigned int y) const {
        return refractionReflectionRefractionIntensity[y*camera->width+x];
    }
    
    real FixedSphereTransparent::pathBlocked(const v3D& start, const v3D& finish) const {
        v3D v = finish-start;
        v3D x0 = start-this->x;
        real vdotx=dotProduct(v, x0);
        real vdotxsqr=vdotx*vdotx;
        real vsqr=v.getSqrLength();
        real sqroot = vsqr*(r*r-x0.getSqrLength())+(vdotxsqr);
        
        real lambda;
        real MYEPSILON = EPSILON*10;
        
        if (sqroot >= MYEPSILON) {
            sqroot=sqrt(sqroot);
            if (vdotx+sqroot <= MYEPSILON) {
                lambda=-1*(vdotx+sqroot)/vsqr;
            }
            else {
                lambda=(sqroot-vdotx)/vsqr;
            }
            
            // We now make a scale between -MYEPSILON and MYEPSILON at each of lambda = 0 and lambda = 1
            // And return a value depending on how close to zero the lambda is.
            
            if ( ( lambda >= MYEPSILON ) && ( lambda <= 1.0-MYEPSILON) )
            {
                return 0.0;
            }
            else
            {
                if (lambda<0.5)
                {
                    // lambda is close to zero
                    lambda = (lambda-MYEPSILON);
                    if (lambda<=0)
                        return 1.0;                        
                    else
                        return lambda/MYEPSILON;
                }
                else
                {
                    // lambda is close to one
                    lambda = (lambda-1.0+MYEPSILON);
                    if (lambda>=0)
                        return 1.0;                        
                    else
                        return -1*lambda/MYEPSILON; 
                }
            }
        }
        else
            return false;
    }
