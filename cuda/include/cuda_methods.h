#ifndef CUDA_METHODS_H
#define CUDA_METHODS_H

void cuda_resampling(int num_pts, int num_voxels,   float *flattenXYZ,
		     float *flattenRGB, int *voxel_offset, int *neighbor_ids,
		     int x_idx, int y_idx, int z_idx,int num_samples,uint *sample_arr,float *pdens);

void cuda_segmentation(int num_pts, int num_voxels, float* pdens,  float *flattenXYZ,
		       float *flattenRGB, int *voxel_offset, int *neighbor_ids,
		       int x_idx, int y_idx, int z_idx,int* parents_ptr);
#endif
