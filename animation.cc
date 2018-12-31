#include "animation.h"

Surface::Surface(int width, int height)
{
    this->w=width;
    this->h=height;
    
    if ( SDL_Init(SDL_INIT_VIDEO) < 0 ) {
        //printf("Unable to init SDL: %s\n", SDL_GetError());
        return;
    }
    atexit(SDL_Quit);
    screen=SDL_SetVideoMode(width, height, 32, SDL_HWSURFACE|SDL_DOUBLEBUF);
    if ( screen == NULL ) {
        //printf("Unable to set "+width+"x"+height+" video: %s\n", SDL_GetError());
        return;
    }   
}

void Surface::drawPixel(int x, int y,
        Uint8 R, Uint8 G, Uint8 B) {
    Uint32 color = SDL_MapRGB(screen->format, R, G, B);
    switch (screen->format->BytesPerPixel) {
        case 1: // Assuming 8-bpp
        {
            Uint8 *bufp;
            bufp = (Uint8 *)screen->pixels + y*screen->pitch + x;
            *bufp = color;
        }
        break;
        case 2: // Probably 15-bpp or 16-bpp
        {
            Uint16 *bufp;
            bufp = (Uint16 *)screen->pixels + y*screen->pitch/2 + x;
            *bufp = color;
        }
        break;
        case 3: // Slow 24-bpp mode, usually not used
        {
            Uint8 *bufp;
            bufp = (Uint8 *)screen->pixels + y*screen->pitch + x * 3;
            if(SDL_BYTEORDER == SDL_LIL_ENDIAN) {
                bufp[0] = color;
                bufp[1] = color >> 8;
                bufp[2] = color >> 16;
            } else {
                bufp[2] = color;
                bufp[1] = color >> 8;
                bufp[0] = color >> 16;
            }
        }
        break;
        case 4: // Probably 32-bpp
        {
            Uint32 *bufp;
            bufp = (Uint32 *)screen->pixels + y*screen->pitch/4 + x;
            *bufp = color;
        }
        break;
    }
}

void Surface::sLock() {
    if ( SDL_MUSTLOCK(screen) ) {
        if ( SDL_LockSurface(screen) < 0 ) {
            return;
        }
    }
}

void Surface::sUnlock() {
    if ( SDL_MUSTLOCK(screen) ) {
        SDL_UnlockSurface(screen);
    }
}

int Surface::width()
{
    return w;
}

int Surface::height()
{
    return h;
}

Animation::Animation(int width, int height, void (*DrawFunc)(Surface&)):drawSurface(width,height)
{
   // this->DrawFunc=DrawFunc;
    
   
    int done=0;
    
    while(done == 0) {
        SDL_Event event;
        
        while ( SDL_PollEvent(&event) ) {
            if ( event.type == SDL_QUIT )  {  done = 1;  }
            
            if ( event.type == SDL_KEYDOWN ) {
                if ( event.key.keysym.sym == SDLK_ESCAPE ) { done = 1; }
            }
        }
        drawSurface.sLock();        
        DrawFunc(this->drawSurface);
        drawSurface.sUnlock();
        SDL_Flip(drawSurface.screen);
    }    
}
