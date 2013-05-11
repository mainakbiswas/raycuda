#ifndef D_FACE_H
#define D_FACE_H

#include "d_triple.h"
#include "d_ray.h"
#include "d_hit.h"

//face can only be along x, y or z plane
class Face
{
public:         //eqn is axis = value
    char axis;
    double value; 
    
    __device__ Face(char a='i', double v=0.0) : axis(a), value(v) { }
    __device__ Face(const Face &f) : axis(f.axis), value(f.value) { }
    
    __device__ Hit intersect(const Ray &ray){
        double v1, v2, v3;
        
        v1 = value;
        switch(axis){
            case 'x':
                v2 = ray.O.x;
                v3 = ray.D.x;
                break;
            case 'y':
                v2 = ray.O.y;
                v3 = ray.D.y;
                break;
            case 'z':
                v2 = ray.O.z;
                v3 = ray.D.z;
                break;
        }        
        
        if(v3 == 0.0) return NO_HIT; //parallel
        
        double dist = (v1 - v2) / v3;
        if(dist <= 0.0)
            return NO_HIT;
        else
            return Hit(dist, Vector(0,0,0));  //normal not used, and I want null objects in c++
    }
    
    __device__ bool operator==(const Face &f){
        return (axis == f.axis && value == f.value);
    }
};

#endif
