// 
// File:   animation.h
// Author: leon
//
// Created on 11 July 2007, 20:41
//

#ifndef _ANIMATION_H
#define	_ANIMATION_H

#include <SDL/SDL.h>

class Surface
{
protected:
    int w, h;
    SDL_Surface *screen; 
    
    void sLock();
    void sUnlock();  
    
public:
    
    Surface(int width, int height);
    void drawPixel(int x, int y, Uint8 R, Uint8 G, Uint8 B);
    int width();
    int height();
    
    friend class Animation;   
};

class Animation
{
protected:
    //void (*DrawFunc)(Surface&);    
    
public: 
    Surface drawSurface;
    Animation(int width, int height, void (*DrawFunc)(Surface&));       
};

#endif	/* _ANIMATION_H */

