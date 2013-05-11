#include "d_triangle.h"
#include "d_light.h"
#include "d_ray.h"
#include "d_cubegrid.h"
#include "d_triple.h"
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

//double array serialization of scene
//eye, (nlight,) lpos, lcolor, ..., (ntri,) triP1, triP2, triP3, col, ka, kd, ks, n, ... ; () items reqd but not in array

// shadow,nlight,ntri will be passed to function for tracing

//double array is copied to mem and memory for device objects allocated
__global__ void createObjects(double *d_arr, Point *d_eye, int nLight, Light *d_lights, int nTri, Triangle *d_triangles, CubeGrid *d_cubeDS){
    if(blockIdx.x==0 && threadIdx.x==0){    //create eye and lights
        double eyex = d_arr[0];
        double eyey = d_arr[1];
        double eyez = d_arr[2];
        new (d_eye) Point(eyex, eyey, eyez);
        int cur_pos = 3;
        for(int i=0; i < nLight; i++){
            double lp1,lp2,lp3,lc1,lc2,lc3;
            lp1 = d_arr[cur_pos++];
            lp2 = d_arr[cur_pos++];
            lp3 = d_arr[cur_pos++];
            lc1 = d_arr[cur_pos++];
            lc2 = d_arr[cur_pos++];
            lc3 = d_arr[cur_pos++];
            new (d_lights + i) Light(Point(lp1,lp2,lp3), Color(lc1,lc2,lc3));
        }
        new (d_cubeDS) CubeGrid();
    }
    
    //create triangles
    double *tri = d_arr + 3 + (6*nLight);
    int myTri = nTri / (gridDim.x * blockDim.x) + 1;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    double *tristart = tri + id * myTri * 16;
    int curpos = 0;
    Triangle *myTriStart = d_triangles + id * myTri;
    for(int i = 0; i < myTri && (id*myTri+i < nTri); i++){
        double d[16];
        for(int j=0;j<16;j++)
            d[j] = tristart[curpos++];
        Point p1(d[0], d[1], d[2]);
        Point p2(d[3], d[4], d[5]);
        Point p3(d[6], d[7], d[8]);
        Color c(d[9], d[10], d[11]);
        double ka = d[12];
        double kd = d[13];
        double ks = d[14];
        double n = d[15];
        new (myTriStart + i) Triangle(p1,p2,p3);
        myTriStart[i].material.color = c;
        myTriStart[i].material.ka = ka;
        myTriStart[i].material.kd = kd;
        myTriStart[i].material.ks = ks;
        myTriStart[i].material.n = n;
    }
}

//createcubegrid must be called <<<1,1>>> fashion
__global__ void createCubeGridBoxes(Box *d_boxes){
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    new (d_boxes + id) Box();
}

/*
__global__ void createCubeGrid(CubeGrid *d_cubeDS){
        new (d_cubeDS) CubeGrid();
        //d_cubeDS->boxes = d_boxes;
}
*/

__global__ void createBoxesofCubeGrid(CubeGrid *d_cubeDS, double maxx, double maxy, double maxz,
                                        double minx, double miny, double minz, 
                                        double xstep, double ystep, double zstep, Box *boxes){
    if(blockIdx.x == 0 && threadIdx.x == 0){
        d_cubeDS->bound.f1.axis = 'x';
        d_cubeDS->bound.f1.value = minx;
        d_cubeDS->bound.f2.axis = 'x';
        d_cubeDS->bound.f2.value = maxx;
        d_cubeDS->bound.f3.axis = 'y';
        d_cubeDS->bound.f3.value = miny;
        d_cubeDS->bound.f4.axis = 'y';
        d_cubeDS->bound.f4.value = maxy;
        d_cubeDS->bound.f5.axis = 'z';
        d_cubeDS->bound.f5.value = minz;
        d_cubeDS->bound.f6.axis = 'z';
        d_cubeDS->bound.f6.value = maxz;
    }
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int k = id % NCUBES;
    int j = (id / NCUBES) % NCUBES;
    int i = (id / (NCUBES * NCUBES)) % NCUBES;
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f1.axis = 'x';
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f1.value = minx + xstep * i;
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f2.axis = 'x';
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f2.value = minx + xstep * (i+1);
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f3.axis = 'y';
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f3.value = miny + ystep * j;
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f4.axis = 'y';
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f4.value = miny + ystep * (j+1);
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f5.axis = 'z';
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f5.value = minz + zstep * k;
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f6.axis = 'z';
    /*d_cubeDS->*/boxes[i*NCUBES*NCUBES + j*NCUBES + k].f6.value = minz + zstep * (k+1); 
}

__global__ void insertIntoCubeGrid(CubeGrid *d_cubeDS, Triangle *d_triangles, int nTri, double minx, double miny, double minz,
                                                                                        double xstep, double ystep, double zstep, Box *boxes){
    int myTri = nTri / (gridDim.x * blockDim.x) + 1;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    Triangle *myTriStart = d_triangles + id * myTri;
    unsigned int xmin,xmax, ymin, ymax, zmin, zmax;
    for(int i = 0; i < myTri && (id*myTri+i < nTri); i++){
        BoundBox bb = myTriStart[i].getBoundBox();        
        xmin = (int)((bb.f1.value - minx) / xstep);
        xmax = (int)((bb.f2.value - minx) / xstep);
        ymin = (int)((bb.f3.value - miny) / ystep);
        ymax = (int)((bb.f4.value - miny) / ystep);
        zmin = (int)((bb.f5.value - minz) / zstep);
        zmax = (int)((bb.f6.value - minz) / zstep);
        for(unsigned int j=xmin;j<=xmax && j<NCUBES;j++)
            for(unsigned int k=ymin;k<=ymax && k<NCUBES;k++)
                for(unsigned int l=zmin;l<=zmax && l<NCUBES;l++)
                    /*d_cubeDS->*/boxes[j*NCUBES*NCUBES + k*NCUBES + l].addObject(myTriStart+i);
    }
}

#define MAXRECDEPTH 3
#define SHADOWS 1

__device__ Color computeShading(Light l, Point p, Vector n, Vector view, Material m, CubeGrid *cubeDS, int isReflectedLight, Box *box){
    Color color;
    Vector lightVector = -(l.position - p).normalized(); // vector from light source to point
    Ray shadowRay(p, -lightVector); //ray from point to light source
    double NdotL = n.dot(-lightVector); //dot product of surface normal and vector from point to light source
    Vector R = (n*2.0*(NdotL) + lightVector).normalized();
    double VdotR = view.dot(R);

    Triangle *obj = cubeDS->getNearest(shadowRay, box);
    if(!SHADOWS || isReflectedLight || obj==NULL){  // no obstruction on the way to light source
        //diffuse colour
        NdotL = NdotL>0.0?NdotL:0.0;//max(NdotL, 0.0);
        Color dmc = m.color * m.kd;
        if(!isReflectedLight)
            color += (l.color * dmc * NdotL);
            
        //specular colour
        VdotR = VdotR>0.0?VdotR:0.0;
        VdotR = pow(VdotR, m.n);
        color += l.color * m.ks * VdotR;
    }
    if(!isReflectedLight)
        color += l.color * m.color * m.ka; //ambient light
    return color;
}  

template<int numRefls> __device__ Color
trace(Ray &ray, CubeGrid *cubeDS, Light* d_lights, int nLight, Box *box){
    Triangle *obj = cubeDS->getNearest(ray, box);
    if (!obj) return Color(0.0, 0.0, 0.0);
    Triangle tri = *(obj);
    Hit min_hit = tri.intersect(ray);
    Material material = tri.material;              //the hit objects material
    Point hit = ray.at(min_hit.t);                 //the hit point
    Vector N = min_hit.N;                          //the normal at hit point
    Vector V = -ray.D;                             //the view vector
    
    Color color, reflColor;
    double NdotView = N.dot(ray.D);
    Vector reflVector = -(N*2.0*(NdotView) - ray.D); // reflected vector from hit point away from object
    Ray reflRay(hit, reflVector);
    if(material.ks > 0.0)
        reflColor = trace<numRefls -1>(reflRay, cubeDS, d_lights, nLight, box);
    
    for(int i = 0; i < nLight; i++){
        color += computeShading(d_lights[i], hit, N, V, material, cubeDS, 0, box);
    }
    Light reflLight(reflRay.at(1.0), reflColor);
    color += computeShading(reflLight, hit, N, V, material, cubeDS, 1, box);
    return color;
}

template<> __device__ Color 
trace<-1>(Ray &ray, CubeGrid *cubeDS, Light* d_lights, int nLight, Box *box){
    return Color(0.0, 0.0, 0.0);
} 

//OUR PICTURE IS 400*400

//should be called <<<25*25,16*16>>> fashion ,1 thread per ray, eye and lights shared try cubeDS too in shared
__global__ void computePixel(Point *d_eye, int nLight, Light *d_lights, CubeGrid *d_cubeDS, double *d_pixels, Box *d_boxes){
    extern __shared__ char mem[];
    Point* eye = (Point*)mem;
    Light* lights = (Light*) (mem + sizeof(Point));
    if(threadIdx.x == 0){
        new (eye) Point(*d_eye);
        for(int i=0;i<nLight;i++){
            new (lights+i) Light(*(d_lights+i));
        }
    }
    __syncthreads(); 
    //from thrd.x and thrd.y compute moton i,j
    //int i = threadIdx.x;
    //int j = threadIdx.y;
    //i = ((i&8)<<3) | ((i&4)<<2) | ((i&2)<<1) | (i&1);
    //j = ((j&8)<<3) | ((j&4)<<2) | ((j&2)<<1) | (j&1);
    //i = (i<<1) | j;
    //int a = (i&(15<<4))>>4;
    //int b = (i&15);
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    Point pixel(x, 400-1-y, 0);
    Ray ray(*eye, (pixel-(*eye)).normalized());
    Color col = trace<MAXRECDEPTH>(ray, d_cubeDS, lights, nLight, d_boxes);
    int pos = (y * 400 + x) * 3;
    d_pixels[pos] = col.r;
    d_pixels[pos+1] = col.g;
    d_pixels[pos+2] = col.b;
}

//parse file and create double array on host by mallocing
typedef struct{
    double *sceneArr;
    int nLight;
    int nTri;
}Scene;

Scene createScene(const char *fileName){
    Scene scene;
    double eye[3];
    FILE* fp = fopen(fileName, "r");
    fscanf(fp, "%lf %lf %lf", &eye[0],&eye[1],&eye[2]);
    fscanf(fp, "%d %d", &scene.nLight, &scene.nTri);
    size_t nDouble = scene.nLight * 6 + scene.nTri * 16;
    scene.sceneArr = (double*)malloc(sizeof(double) * (nDouble + 3));
    scene.sceneArr[0] = eye[0];
    scene.sceneArr[1] = eye[1];
    scene.sceneArr[2] = eye[2];
    size_t curpos = 3;
    for(size_t i=0; i<nDouble; i++){
        fscanf(fp,"%lf",&scene.sceneArr[curpos]);
        curpos++;
    }
    return scene;
}

//call gpu funcs, alloc mamory, copy things
double *getPixelsByRender(Scene scene){
    Point *d_eye;
    Light *d_lights;
    Triangle *d_triangles;
    CubeGrid *d_cubeDS;
    double *d_sceneArr;
    
    double *result, *d_result;
    Box *d_boxes;
    
    cudaError_t err;
    
    err = cudaMalloc((void**)&d_eye, sizeof(Point));
    if(err!=cudaSuccess){
        printf("\n1 %s\n", cudaGetErrorString(err));
    }
    err = cudaMalloc((void**)&d_lights, sizeof(Light) * scene.nLight);
    if(err!=cudaSuccess){
        printf("\n2 %s\n", cudaGetErrorString(err));
    }
    err = cudaMalloc((void**)&d_triangles, sizeof(Triangle) * scene.nTri);
    if(err!=cudaSuccess){
        printf("\n3 %s\n", cudaGetErrorString(err));
    }
    err = cudaMalloc((void**)&d_cubeDS, sizeof(CubeGrid));
    if(err!=cudaSuccess){
        printf("\n4 %s\n", cudaGetErrorString(err));
    }
    err = cudaMalloc((void**)&d_sceneArr, sizeof(double) * (scene.nLight * 6 + scene.nTri * 16 + 3));
    if(err!=cudaSuccess){
        printf("\n5 %s\n", cudaGetErrorString(err));
    }
    err = cudaMalloc((void**)&d_result, sizeof(double) * 400 * 400 * 3); 
    if(err!=cudaSuccess){
        printf("\n6 %s\n", cudaGetErrorString(err));
    }
    err = cudaMalloc((void**)&d_boxes, sizeof(Box) * NCUBES * NCUBES * NCUBES);
    if(err!=cudaSuccess){
        printf("\n7 %s\n", cudaGetErrorString(err));
    }

    
    cudaMemcpy(d_sceneArr, scene.sceneArr, sizeof(double) * (scene.nLight * 6 + scene.nTri * 16 + 3), cudaMemcpyHostToDevice);
    result = (double*)malloc(sizeof(double) * 400 * 400 * 3);
    struct timeval t_start,end,diff;
    gettimeofday(&t_start,NULL);
    createObjects<<<512,256>>>(d_sceneArr, d_eye, scene.nLight, d_lights, scene.nTri, d_triangles, d_cubeDS);
    //err = cudaGetLastError();
    //if(err!=cudaSuccess){
    //    printf("\n8 %s\n", cudaGetErrorString(err));
    //}
    //cudaStreamSynchronize(0);
    //gettimeofday(&end,NULL);
    //timersub(&end, &t_start, &diff);
    //printf("Created objects in %llu usecs\n", diff.tv_sec*1000000+diff.tv_usec );
    createCubeGridBoxes<<<(NCUBES*NCUBES*NCUBES/200), 200>>>(d_boxes); /*(NCUBES*NCUBES*NCUBES/200),200*/
    //err = cudaGetLastError();
    //if(err!=cudaSuccess){
    //    printf("\n9 %s\n", cudaGetErrorString(err));
    //}
    //createCubeGrid<<<1,1>>>(d_cubeDS);
    //err = cudaGetLastError();
    //if(err!=cudaSuccess){
    //    printf("\n10 %s\n", cudaGetErrorString(err));
    //}
    //cudaStreamSynchronize(0);
    //gettimeofday(&end,NULL);
    //timersub(&end, &t_start, &diff);
    //printf("Call to new cubegrid completed in %llu usecs\n", diff.tv_sec*1000000+diff.tv_usec );
    double maxx,minx,maxy,miny,maxz,minz;
    int start = 3 + 6 * scene.nLight;
    maxx = scene.sceneArr[start];
    minx = maxx;
    maxy = scene.sceneArr[start+1];
    miny = maxy;//scene.sceneArr[start+1];
    maxz = scene.sceneArr[start+2];
    minz = maxz; //scene.sceneArr[start+2];
    for(int i=0; i<scene.nTri; i++){
        double x = scene.sceneArr[start+i*16];
        double y = scene.sceneArr[start+i*16+1];
        double z = scene.sceneArr[start+i*16+2];
        maxx = max(x, maxx);
        minx = min(x, minx);
        maxy = max(y, maxy);
        miny = min(y, miny);
        maxz = max(z, maxz);
        minz = min(z, minz);
        x = scene.sceneArr[start+i*16+3];
        y = scene.sceneArr[start+i*16+4];
        z = scene.sceneArr[start+i*16+5];
        maxx = max(x, maxx);
        minx = min(x, minx);
        maxy = max(y, maxy);
        miny = min(y, miny);
        maxz = max(z, maxz);
        minz = min(z, minz); 
        x = scene.sceneArr[start+i*16+6];
        y = scene.sceneArr[start+i*16+7];
        z = scene.sceneArr[start+i*16+8];
        maxx = max(x, maxx);
        minx = min(x, minx);
        maxy = max(y, maxy);
        miny = min(y, miny);
        maxz = max(z, maxz);
        minz = min(z, minz);
    }
    
    double xstep = (maxx - minx) / NCUBES;
    double ystep = (maxy - miny) / NCUBES;
    double zstep = (maxz - minz) / NCUBES;
    //gettimeofday(&end,NULL);
    //timersub(&end, &t_start, &diff);
    //printf("Computed max, min in %llu usecs\n", diff.tv_sec*1000000+diff.tv_usec );    
 
    createBoxesofCubeGrid<<<(NCUBES*NCUBES*NCUBES/200), 200>>>(d_cubeDS,maxx,maxy,maxz,minx,miny,minz,xstep,ystep,zstep, d_boxes);
    //err = cudaGetLastError();
    //if(err!=cudaSuccess){
    //    printf("\n11 %s\n", cudaGetErrorString(err));
    //}
    //cudaStreamSynchronize(0);
    //gettimeofday(&end,NULL);
    //timersub(&end, &t_start, &diff);
    //printf("Created boxes in %llu usecs\n", diff.tv_sec*1000000+diff.tv_usec );
    insertIntoCubeGrid<<<512,256>>>(d_cubeDS, d_triangles, scene.nTri, minx, miny, minz, xstep, ystep, zstep, d_boxes);
    //err = cudaGetLastError();
    //if(err!=cudaSuccess){
    //    printf("\n12 %s\n", cudaGetErrorString(err));
    //}
    //cudaStreamSynchronize(0);
    //gettimeofday(&end,NULL);
    //timersub(&end, &t_start, &diff);
    //printf("Created cubeDS in %llu usecs\n", diff.tv_sec*1000000+diff.tv_usec );
    //computepixel 25*25,16*16 with (Point *d_eye, int nLight, Light *d_lights, CubeGrid *d_cubeDS, double *d_pixels) remember shared mem
    //can material be stored in shared?
    dim3 grdDim(25,25);
    dim3 blkDim(16,16);
    computePixel<<<grdDim,blkDim,(sizeof(Point) + scene.nLight * sizeof(Light))>>>(d_eye, scene.nLight, d_lights, d_cubeDS, d_result, d_boxes);
    err = cudaGetLastError();
    if(err!=cudaSuccess){
        printf("\n13 %s\n", cudaGetErrorString(err));
    }
    //cudaStreamSynchronize(0);
    cudaMemcpy(result, d_result, sizeof(double) * 400 * 400 * 3, cudaMemcpyDeviceToHost);
    gettimeofday(&end,NULL);
    timersub(&end, &t_start, &diff);
    printf("Computed pixels in %llu usecs\n", diff.tv_sec*1000000+diff.tv_usec );
    
    return result;
}

void createPPM(double *pixels, const char *ppmFile){
    FILE* fp = fopen(ppmFile,"w");
    fprintf(fp,"P3\n# CREATOR: RAYTRACER\n400 400\n255\n");
    for(int i=0;i<400*400*3;i++){
        double val = pixels[i]>1.0?1.0:pixels[i]; //clamping here
        fprintf(fp,"%.0lf\n", val*255);
    }
    fclose(fp);
}


int main(int argc, char *argv[]){
    if(argc<3) printf("usage: %s inputFile outputFile\n",argv[0]);
    //printf("%u\n", sizeof(CubeGrid));
    Scene scene = createScene(argv[1]);
    struct timeval start,end,diff;
    gettimeofday(&start,NULL);
    double *pixels = getPixelsByRender(scene);
    gettimeofday(&end,NULL);
    timersub(&end, &start, &diff);
    //printf("Rendered image in %llu usecs\n", diff.tv_sec*1000000+diff.tv_usec );
    createPPM(pixels, argv[2]);
    return 0;
}
