#ifndef D_BOX_H
#define D_BOX_H

#include "d_face.h"
#include "d_triangle.h"
#include "d_triple.h"
#include "d_hit.h"
#include <stdlib.h>
#include <math_constants.h>

class Box
{
public:
    Face f1, f2, f3, f4, f5, f6;
    Triangle* objects[300];
    unsigned int nObj;
    
    __device__ Box(Face f1, Face f2, Face f3, Face f4, Face f5, Face f6) : 
    f1(f1), f2(f2), f3(f3), f4(f4), f5(f5), f6(f6) { 
        nObj = 0;
    }
    
    __device__ Box() { nObj = 0;}
    
    __device__ void addObject(Triangle *obj){
        unsigned int place = atomicInc(&nObj,299);
        objects[place] = obj;
    }
    
    //returns closest face where ray intersects this box
    __device__ double intersect(const Ray &ray){
        Face* f[]={&f1, &f2, &f3, &f4, &f5, &f6};
        double dist = CUDART_INF_F;
        for(int i=0;i<6;i++){
            double d = computeDist(ray, f[i]);
            if(!isnan(d) && d < dist){
                dist = d;
            }
        }
        return dist;
    }
    
    __device__ bool contains(Point p){
        if(p.x < f1.value - 0.00001)
            return false;
        if(p.x > f2.value + 0.00001)
            return false;
        if(p.y < f3.value - 0.00001)
            return false;
        if(p.y > f4.value + 0.00001)
            return false;
        if(p.z < f5.value - 0.00001)
            return false;
        if(p.z > f6.value + 0.00001)
            return false;
        return true;
    }
    
    __device__ double computeDist(const Ray &ray, Face *f){
        Hit h = f->intersect(ray);
        if(!isnan(h.t)){
            Point p = ray.at(h.t);
            if(contains(p))
                return h.t;
            else
                return CUDART_NAN_F;
        }
        return h.t;
    }
};

#endif
