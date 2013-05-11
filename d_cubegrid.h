#ifndef D_CUBE_GRID_H
#define D_CUBE_GRID_H

#include "d_box.h"
#include "d_triangle.h"
#include "d_boundbox.h"
#include "d_triple.h"
#include "d_hit.h"
#define NCUBES 100

class CubeGrid
{
public:
    //Box *boxes; //[NCUBES*NCUBES*NCUBES];
    Box bound;
        
    //__device__ void storeObjects(Triangle objects[], int nObj);
    
    __device__ Triangle* getNearest(const Ray &ray, Box *boxes);
};
/*
__device__ void CubeGrid::storeObjects(Triangle *objects, int nObj){
    if(nObj == 0) return;
    double maxx, minx, maxy, miny, maxz, minz;
    BoundBox bb = objects[0].getBoundBox();
    minx = bb.f1.value;
    maxx = bb.f2.value;
    miny = bb.f3.value;
    maxy = bb.f4.value;
    minz = bb.f5.value;
    maxz = bb.f6.value;
    unsigned int i,j,k;
    for(i=1; i < nObj; i++){
        bb = objects[i].getBoundBox();
        double x1,x2,y1,y2,z1,z2;
        x1 = bb.f1.value;
        x2 = bb.f2.value;
        y1 = bb.f3.value;
        y2 = bb.f4.value;
        z1 = bb.f5.value;
        z2 = bb.f6.value;
        
        minx = x1 < minx? x1 : minx;
        maxx = x2 > maxx? x2 : maxx;
        miny = y1 < miny? y1 : miny;
        maxy = y2 > maxy? y2 : maxy;
        minz = z1 < minz? z1 : minz;
        maxz = z2 > maxz? z2 : maxz;
    }
    bound.f1.axis = 'x';
    bound.f1.value = minx;
    bound.f2.axis = 'x';
    bound.f2.value = maxx;
    bound.f3.axis = 'y';
    bound.f3.value = miny;
    bound.f4.axis = 'y';
    bound.f4.value = maxy;
    bound.f5.axis = 'z';
    bound.f5.value = minz;
    bound.f6.axis = 'z';
    bound.f6.value = maxz;
    
    double xstep = (maxx - minx) / NCUBES;
    double ystep = (maxy - miny) / NCUBES;
    double zstep = (maxz - minz) / NCUBES;
    
    for(i=0;i<NCUBES;i++)
        for(j=0;j<NCUBES;j++)
            for(k=0;k<NCUBES;k++){
                boxes[i][j][k].f1.axis = 'x';
                boxes[i][j][k].f1.value = minx + xstep * i;
                boxes[i][j][k].f2.axis = 'x';
                boxes[i][j][k].f2.value = minx + xstep * (i+1);
                boxes[i][j][k].f3.axis = 'y';
                boxes[i][j][k].f3.value = miny + ystep * j;
                boxes[i][j][k].f4.axis = 'y';
                boxes[i][j][k].f4.value = miny + ystep * (j+1);
                boxes[i][j][k].f5.axis = 'z';
                boxes[i][j][k].f5.value = minz + zstep * k;
                boxes[i][j][k].f6.axis = 'z';
                boxes[i][j][k].f6.value = minz + zstep * (k+1);
            }
    
    unsigned int xmin,xmax, ymin, ymax, zmin, zmax;
    for(i=0;i < nObj; i++){
        BoundBox bb = objects[i].getBoundBox();
        xmin = (int)((bb.f1.value - minx) / xstep);
        xmax = (int)((bb.f2.value - minx) / xstep);
        ymin = (int)((bb.f3.value - miny) / ystep);
        ymax = (int)((bb.f4.value - miny) / ystep);
        zmin = (int)((bb.f5.value - minz) / zstep);
        zmax = (int)((bb.f6.value - minz) / zstep);
        
        unsigned int l;
        for(j=xmin;j<=xmax && j<NCUBES;j++)
            for(k=ymin;k<=ymax && k<NCUBES;k++)
                for(l=zmin;l<=zmax && l<NCUBES;l++)
                    boxes[j][k][l].addObject(objects+i);
    }
}
*/

__device__ Triangle* CubeGrid::getNearest(const Ray& ray, Box *boxes){
    //if ray does not intersect boundbox return NULL
    int firstI, firstJ, firstK;
    double minx = bound.f1.value;
    double maxx = bound.f2.value;
    double miny = bound.f3.value;
    double maxy = bound.f4.value;
    double minz = bound.f5.value;
    double maxz = bound.f6.value;
    
    double xstep = (maxx - minx) / NCUBES;
    double ystep = (maxy - miny) / NCUBES;
    double zstep = (maxz - minz) / NCUBES;
    
    Point p;
    Ray currentRay(p, ray.D);
    if(bound.contains(ray.O)){
        firstI = (int)((ray.O.x - minx) / xstep);
        firstJ = (int)((ray.O.y - miny) / ystep);
        firstK = (int)((ray.O.z - minz) / zstep);
        currentRay.O = ray.O;
    }
    else{
        double t = bound.intersect(ray);
        if(t == CUDART_INF_F)
            return NULL;             //ray does not even hot the big bounding box
        else{
            p = ray.at(t);
            firstI = (int)((p.x - minx) / xstep);
            firstJ = (int)((p.y - miny) / ystep);
            firstK = (int)((p.z - minz) / zstep);
            currentRay.O = p;
        }    
    }
    
    int i = firstI==NCUBES?NCUBES-1:firstI, j = firstJ==NCUBES?NCUBES-1:firstJ, k = firstK==NCUBES?NCUBES-1:firstK;
    
    while( 
            i >= 0 && i <= NCUBES - 1 &&
            j >= 0 && j <= NCUBES - 1 &&
            k >= 0 && k <= NCUBES - 1
         ){
        //for all triangles in boxes[i][j][k] find the one with min distance and having hit point in this box
        Triangle *obj = NULL;
        Hit min_hit(CUDART_INF_F, Vector());
        for(unsigned int l=0; l < boxes[i*NCUBES*NCUBES + j*NCUBES + k].nObj; l++){
            Hit hit((boxes[i*NCUBES*NCUBES + j*NCUBES + k].objects[l])->intersect(ray));
            if (!isnan(hit.t) && hit.t<min_hit.t && boxes[i*NCUBES*NCUBES + j*NCUBES + k].contains(ray.at(hit.t))) {
                min_hit = hit;
                obj = boxes[i*NCUBES*NCUBES + j*NCUBES + k].objects[l];
            }
        }
        
        if(obj)
            return obj;
        
        double dist[3];
        if(ray.D.x < 0){
            dist[0] = (boxes[i*NCUBES*NCUBES + j*NCUBES + k].f1.intersect(ray)).t;
        }
        else{
            dist[0] = (boxes[i*NCUBES*NCUBES + j*NCUBES + k].f2.intersect(ray)).t;
        }
        if(ray.D.y < 0){
            dist[1] = (boxes[i*NCUBES*NCUBES + j*NCUBES + k].f3.intersect(ray)).t;
        }
        else{
            dist[1] = (boxes[i*NCUBES*NCUBES + j*NCUBES + k].f4.intersect(ray)).t;
        }
        if(ray.D.z < 0){
            dist[2] = (boxes[i*NCUBES*NCUBES + j*NCUBES + k].f5.intersect(ray)).t;
        }
        else{
            dist[2] = (boxes[i*NCUBES*NCUBES + j*NCUBES + k].f6.intersect(ray)).t;
        }
        
        int l=4; double min = CUDART_INF_F;
        for(int m = 0; m<3;m++){
            if(!isnan(dist[m]) && dist[m] < min){
                min = dist[m];
                l = m;
            }
        }
        switch(l){
            case 0:
                if(ray.D.x < 0) i--; else i++; break;
            case 1:
                if(ray.D.y < 0) j--; else j++; break;
            case 2:
                if(ray.D.z < 0) k--; else k++; break;
        }
    }
    
    return NULL;
}

#endif
