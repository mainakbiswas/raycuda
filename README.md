raycuda
=======

Ray tracer with uniform space division acceleration structure in CUDA (50x speedup compared to rayserial)

Compile - nvcc -arch=sm_20 -O5 raycuda.cu -o ray

run - ./ray input.txt output.ppm

sample input - dragon.txt

input format
eyex eyey eyez
nLight nTri
Lightx Lighty Lightz
Lightcolr Lightcolg Lightcolb
... nLight times
Trip1x Trip1y Trip1z
Trip2x Trip2y Trip2z
Trip3x Trip3y Trip3z
TriColr Tricolg TriColb
Ka Kd Ks n
... nTri times
