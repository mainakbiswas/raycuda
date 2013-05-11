//
//  Framework for a raytracer
//  File: hit.h
//
//  Created for the Computer Science course "Introduction Computer Graphics"
//  taught at the University of Groningen by Tobias Isenberg.
//
//  Authors:
//    Maarten Everts
//    Jasper van de Gronde
//    Matthew van der Zwan
//
//  This framework is inspired by and uses code of the raytracer framework of 
//  Bert Freudenberg that can be found at
//  http://isgwww.cs.uni-magdeburg.de/graphik/lehre/cg2/projekt/rtprojekt.html 
//

#ifndef D_HIT_H_
#define D_HIT_H_

#include <limits>
#include "d_triple.h"

#include <math_constants.h>    //for cudart_nan_f
#define NO_HIT Hit(CUDART_NAN_F,Vector(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F))

class Hit
{
public:
    double t;
    Vector N;

    __device__ Hit(const double t, const Vector &normal)
        : t(t), N(normal)
    { }

/*
    static const Hit NO_HIT() { static Hit no_hit(std::numeric_limits<double>::quiet_NaN(),Vector(std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN())); return no_hit; }
*/
};

#endif /* end of include guard: HIT_H_ */
