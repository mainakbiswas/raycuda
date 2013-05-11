//
//  Framework for a raytracer
//  File: light.h
//
//  Created for the Computer Science course "Introduction Computer Graphics"
//  taught at the University of Groningen by Tobias Isenberg.
//
//  Authors:
//    Maarten Everts
//    Jasper van de Gronde
//
//  This framework is inspired by and uses code of the raytracer framework of 
//  Bert Freudenberg that can be found at
//  http://isgwww.cs.uni-magdeburg.de/graphik/lehre/cg2/projekt/rtprojekt.html 
//

#ifndef D_LIGHT_H_PG2BAJRA
#define D_LIGHT_H_PG2BAJRA

#include "d_triple.h"

class Light
{
public:
    __device__ Light(Point pos,Color c) : position(pos), color(c)
    { }

    Point position;
    Color color;
};

#endif /* end of include guard: LIGHT_H_PG2BAJRA */
