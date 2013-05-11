//
//  Framework for a raytracer
//  File: triangle.h
//
//  Created for the Computer Science course "Introduction Computer Graphics"
//  taught at the University of Groningen by Tobias Isenberg.
//
//  Authors:
//    Mainak Biswas
//
//  This framework is inspired by and uses code of the raytracer framework of 
//  Bert Freudenberg that can be found at
//  http://isgwww.cs.uni-magdeburg.de/graphik/lehre/cg2/projekt/rtprojekt.html 
//

#ifndef D_TRIANGLE_H
#define D_TRIANGLE_H

#include "d_hit.h"
#include "d_boundbox.h"
#include "d_triple.h"
#include "d_face.h"
#include "d_ray.h"
#include "d_material.h"

#include <math.h>
#include <stdlib.h>


class Triangle
{
public:
    __device__ Triangle(Point p1, Point p2, Point p3) : 
    vertex1(p1), vertex2(p2), vertex3(p3){ }

    __device__ Hit intersect(const Ray &ray);
    
    __device__ BoundBox getBoundBox();

    Point vertex1;
    Point vertex2;
    Point vertex3;
    Material material;
};

__device__ Hit Triangle::intersect(const Ray &ray)
{
    Vector e1, e2, h, s, q;
    double a,f,u,v;
    e1 = vertex2 - vertex1;
    e2 = vertex3 - vertex1;
    
    h = ray.D.cross(e2);
    a = e1.dot(h);
    
    if(a > -0.000001 && a < 0.000001)
        return NO_HIT;
    
    f = 1 / a;
    s = ray.O - vertex1;
    u = f * s.dot(h);
    
    if(u < 0.0 || u > 1.0)
        return NO_HIT;
    
    q = s.cross(e1);
    v = f * ray.D.dot(q);
    
    if(v < 0.0 || u+v > 1.0)
        return NO_HIT;
    
    double t;
    t = f * e2.dot(q);
    
    if(t <= 0.000001)
        return NO_HIT;
    
    Vector N;
    N = (vertex2 - vertex1).cross(vertex3 - vertex1);
    N.normalize();
    
    return Hit(t, N);
}

__device__ BoundBox Triangle::getBoundBox(){
    double x1 = vertex1.x;
    double y1 = vertex1.y;
    double z1 = vertex1.z;
    double x2 = vertex2.x;
    double y2 = vertex2.y;
    double z2 = vertex2.z;
    double x3 = vertex3.x;
    double y3 = vertex3.y;
    double z3 = vertex3.z;
    double maxx = x1 > x2? (x1 > x3? x1 : x3) : (x2 > x3 ? x2: x3);
    double minx = x1 < x2? (x1 < x3? x1 : x3) : (x2 < x3 ? x2: x3);
    double maxy = y1 > y2? (y1 > y3? y1 : y3) : (y2 > y3 ? y2: y3);
    double miny = y1 < y2? (y1 < y3? y1 : y3) : (y2 < y3 ? y2: y3);
    double maxz = z1 > z2? (z1 > z3? z1 : z3) : (z2 > z3 ? z2: z3);
    double minz = z1 < z2? (z1 < z3? z1 : z3) : (z2 < z3 ? z2: z3);
    
    return BoundBox(
            Face('x', minx),
            Face('x', maxx),
            Face('y', miny),
            Face('y', maxy),
            Face('z', minz),
            Face('z', maxz)
              );
}

#endif /* end of include guard: SPHERE_H_115209AE */
