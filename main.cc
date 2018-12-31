//
// File:   main.cc
// Author: leon
//
// Created on 11 July 2007, 18:43
//

#include <stdlib.h>
#include <iostream>
#include "animation.h"
#include <math.h>
#include "kVector.h"
#include "spherelit.h"

using kBase::v3D;
using kBase::kParticle;
using kBase::dotProduct;

// The array of lights
kParticle *p;

// The number of lights
const int num = 3;

// Width and height of animation
const int width = 600;
const int height = 400;

// Pointers to the sphere and to the camera
FixedSphereTransparent *spr;
FixedCamera *cam;

// Variables used in DrawScene. Saves a little time as we don't create them each time DrawScene is run
float r, g, b;
float prod1, prod2, prod3, prod;
int x, y, k;
v3D v, w;


real t=0;
real sintheta;
real costheta;
real plu=3.14/2;
real fact=7;
real dist=80;

void DrawScene(Surface &drawSurface) {

    // TODO: Make faster with sintheta, costheta pre-computed for each of the angles
    for (k=0; k < num; k++) {
        // Angles to determine which plane the lights rotate in

        sintheta=sin(k*2.0/((real)num)*3.14+plu);
        costheta=cos(k*2.0/((real)num)*3.14+plu);



        real add = 4*k*3.14/num;

        p[k].x=v3D(dist*costheta*cos(t/fact + add), dist*sin(t/fact + add), dist*sintheta*cos(t/fact + add))-spr->getCentre();

    }
    t++;

    for(y=0;y<height;y++) {
        for(x=0;x<width;x++) {

            if (spr->getIntersectTime(x, y)<EPSILON) {
                drawSurface.drawPixel(x, y, 0, 0, 0);
            }
            else {
                v = spr->getIntersectTime(x, y)*cam->getVectorThrough(x, y);
                r=0.0;
                g=0.0;
                b=0.0;
                for (int k=0; k < num; k++) {

                    // For each light source we calculate what it contributes to the pixels colour.
                    // prod1 gives

                    w=(p[k].x-(v+cam->getPosition()));
                    prod1 = dotProduct(spr->getReflection(x, y), w)*1/(w.getLength());   // Basic light calculation
                    if (prod1<0)
                        prod1=0;
                    prod1*=prod1;
                    //prod1*=prod1;
                    //prod1*=spr->getReflectionIntensity(x, y)+0.5;
                    //prod1*=prod1;
                   // prod1 *= spr->pathBlocked(p[k].x, v+cam->getPosition());                // Times by whether the sphere is blocked or not




                    w=(p[k].x-spr->getRefractionPoint(x, y));
                    prod2 = dotProduct(spr->getRefractionRefraction(x, y), w)*1/(w.getLength());
                    if (prod2<0)
                        prod2=0;

                    prod2*=prod2;
                    prod2*=prod2;
                    prod2*=spr->getRefractionRefractionIntensity(x, y);
                    //prod2*=prod2;



                    prod3=prod2*spr->pathBlocked(p[k].x, spr->getRefractionPoint(x, y));


                    prod = (0.99*prod1)+(0.0*prod2)+(0.0*prod3);

                    if (k%3==2)
                        r+=prod;
                    if (k%3==1)
                        g+=prod;
                    if (k%3==0)
                        b+=prod;
                }
                if (r>1) r=1;
                if (g>1) g=1;
                if (b>1) b=1;
                drawSurface.drawPixel(x, y, (Uint8)(r*255), (Uint8)(g*255), (Uint8)(b*255));
            }
        }
    }
}


//
//
//
int main(int argc, char** argv) {

    p=new kBase::kParticle[num];

    FixedCamera c = FixedCamera(v3D(-7, 0, -7), v3D(1/sqrt(2), 0, 1/sqrt(2)), width, height, 800);

    cam = &c;

    FixedSphereTransparent s = FixedSphereTransparent(v3D(0, 0, 0), 2, c);

    spr=&s;

    Animation a=Animation(width, height, DrawScene);
    return (EXIT_SUCCESS);

}
