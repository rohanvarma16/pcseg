#include <stdio.h>
#include <vector>

#include <algorithm>
#include <functional>
#include <math.h>       /* pow */
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
__global__
void saxpy(int n, float a, float *x, float *y)
{
  if( blockIdx.x==1 && threadIdx.x ==1){
  printf("hello in kernel! \n");
  }
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}


__global__
void sampling(float* device_xyz, float* device_rgb, float* device_offset){
  int block_i = blockIdx.x;
  int block_j = blockIdx.y;
  int block_k = blockIDx.z;

  int blockId = ;
  // do pre-computation of neighbor offsets
  //find all 6 neighbors and offset:
  // be careful for last voxel.
  int neighborhood_size;
  int pts_voxel = device_offset[i+1] - device_offset[i];
  
  // initialize importance weight in shared memory for that voxel.
  // shared memory for importance is of size voxel_stride
  for(int global_iter = 0 ; i < pts_voxel ; i+= voxel_stride)
    
  for(int i = 0 ; i < neighborhood_size; i+=256){
    //load 256 neighbors information into shared memory.    
    for(int j =0 ; j < voxel_stride; j+=256){
      int pt_idx = j + threadIdx.x;
      for (int k = 0 ; k < 256 ; k++){
	float dist = (your_x - nbr_x)^2 +(your_y - nbr_y)^2 + (your_z -nbr_z)^2;
	float ew = exp(-1* dist/sigma_sq);
	float your_feature[0] += ew * your_x
	  float your_feature[1] += ew * your_x
	  float your_feature[2] += ew * your_z
	  ..
      }
    }

    for(int j = 0 ; j < voxelstrid ; j+=256){
    imp_wt[j+threadIdx.x] = norm(x- partial_sum computed)
      }
     
    //finished computing importanve weight for voxel stride
  }


      imp_wt[pt_idx] += ();

    }


  }


//2nd: kernel
//maybe use thrust to find local and global sums

// weighted sampling kernel
void weighted_sampling_kenrel{

  int num_samples = local_sum/global_sum * total_samples;
  bin
  for (int i = 0 ; i < bins ; i+=256){
    //    local sum over 256 strides
    for ( int j = 0 ; j < num_samples* ; j++)
      //each thread will generate a random number
      //      and check if it falls in this particular bin[i]
      // and if it does, then break.
  }

}



  }


  
  

}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) 
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
}
int device_setup(int num_pts, int num_voxels,   float* flattenXYZ,
		 float* flattenRGB,int* voxel_offset){


  printf("HELLO IN DEVICE SETUP!\n");

  float *device_xyz, *device_rgb;
  int* device_offset;
  cudaMalloc(&device_xyz, num_pts*3*sizeof(float));
  cudaMalloc(&device_rgb, num_pts*3*sizeof(float));
  cudaMalloc(&device_offset, num_voxels*3*sizeof(int));
  
  gpuErrchk(cudaMemcpy(device_offset,voxel_offset, num_voxels*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_xyz, flattenRGB, num_pts*3*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_rgb, flattenRGB, num_pts*3*sizeof(float),cudaMemcpyHostToDevice));



  float x_grid = 0.5;
  float y_grid = 0.5;
  float z_grid = 0.5;
  float max_xyz = std::max(max_x-min_x,max_y-min_y);
  max_xyz = std::max(max_xyz, max_z-min_z);
  int xy_idx = (max_xyz/x_grid)*(max_xyz/y_grid);
  int x_idx = max_xyz/x_grid;
  int y_idx  = (max_xyz/y_grid);
  int z_idx = (max_xyz/z_grid);

  dim3 gridDim(x_idx,y_idx,z_ids);
  dim3 blockDim(256);

  sampling<<<gridDim,blockDim>>>(device_xyz,device_rgb,device_offset);


}
/*{
  int N = 1<<20;
  float *x, *y, *d_x, *d_y;
  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));

  gpuErrchk(cudaMalloc(&d_x, N*sizeof(float)));
  gpuErrchk(cudaMalloc(&d_y, N*sizeof(float)));

  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  gpuErrchk(cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice));

  // Perform SAXPY on 1M elements
  saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y);

  gpuErrchk(cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost));

  float maxError = 0.0f;
  for (int i =0 ; i < 4 ; i++){
    printf("y= %0.4f \n", y[i]);
  }

  for (int i = 0; i < N; i++)
    maxError = max(maxError, abs(y[i]-4.0f));
  printf("Max error: %f\n", maxError);

  cudaFree(d_x);
  cudaFree(d_y);
  free(x);
  free(y);
  return(1);
*/




