#ifndef CUDA_TEST_H
#define CUDA_TEST_H
int device_setup(int num_pts, int num_voxels,   float* flattenXYZ,
                 float* flattenRGB,int* voxel_offset);
#endif
