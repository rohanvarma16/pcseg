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
void sampling(float* device_xyz, float* device_rgb, int* device_offset, int* neighbor_id, int xy_idx, int y_idx,float* imp_wt, float* pdensity){
  int block_i = blockIdx.x;
  int block_j = blockIdx.y;
  int block_k = blockIdx.z;

  int blockId = xy_idx * block_i + y_idx* block_j + block_k ;

  if(threadIdx.x ==0 ){
    //    printf("hello!, in block: %d \n", blockId);
  }
  int num_threads = 256;

  float pdensity_sum = 0.0;
  float nbr_feature[6];
  float norm_sum = 0.0;
  int my_num_pts = device_offset[blockId + 1 ] - device_offset[blockId];
  float3 nbr_rgb;
  float3 nbr_xyz;
  float3 my_rgb;
  float3 my_xyz;
  float sigma_sq = 0.00005;
  if(threadIdx.x == 0){
    //    printf("my_num_pts: %d , blockId: %d \n",my_num_pts,blockId );
  }
  for( int i = 0 ; i < my_num_pts ; i+= num_threads){
    pdensity_sum = 0.0;
    norm_sum = 0.0;
    nbr_feature[0] = 0.0;
    nbr_feature[1] = 0.0;
    nbr_feature[2] = 0.0;
    nbr_feature[3] = 0.0;
    nbr_feature[4] = 0.0;
    nbr_feature[5] = 0.0;
    
    if(threadIdx.x + i > my_num_pts){
      break;
    }

    my_xyz = *(float3*) &device_xyz[device_offset[blockId] + i + threadIdx.x];
    my_rgb = *(float3*) &device_rgb[device_offset[blockId] + i + threadIdx.x];  
    
    for(int j = 0 ; j < 7 ; j++){
      if(neighbor_id[7*blockId+j] == -1){
	continue;
      }
      int nbr_num_pts = device_offset[neighbor_id[7*blockId +j] + 1 ] - device_offset[neighbor_id[7*blockId + j] ];

      for(int k = 0 ; k < nbr_num_pts; k++){

	nbr_rgb = *(float3*) &device_rgb[device_offset[neighbor_id[blockId*7 + j]+k]];
	nbr_xyz = *(float3*) &device_xyz[device_offset[neighbor_id[blockId*7 + j]+k]];
	 float xyz_dist = pow(my_xyz.x - nbr_xyz.x,2.0)+ 
						     pow((my_xyz.y - nbr_xyz.y),2.0)
						     +pow((my_xyz.z - nbr_xyz.z),2.0);
						     float rgb_dist = pow((my_rgb.x - nbr_rgb.x),2.0)+
                                                     pow((my_rgb.y - nbr_rgb.y),2.0)
						     +pow((my_rgb.z - nbr_rgb.z),2.0);
						     
						     float Aij_ew = exp(-1.0 * xyz_dist/sigma_sq);
						     float Aij_ew2 = exp(-1.0 * (xyz_dist + rgb_dist)/sigma_sq);						 
						     
						     pdensity_sum += Aij_ew2;
						     nbr_feature[0] += Aij_ew * nbr_xyz.x;
						     nbr_feature[1] += Aij_ew* nbr_xyz.y;
						     nbr_feature[2] += Aij_ew* nbr_xyz.z;
						     nbr_feature[3] += Aij_ew* nbr_rgb.x;
						     nbr_feature[4] += Aij_ew* nbr_rgb.y;
						     nbr_feature[5] += Aij_ew* nbr_rgb.z; 
      }
    }
    
    norm_sum += pow(my_xyz.x - nbr_feature[0],2.0);
    norm_sum += pow(my_xyz.y - nbr_feature[1],2.0);
    norm_sum += pow(my_xyz.z - nbr_feature[2],2.0);
    norm_sum += pow(my_xyz.x - nbr_feature[3],2.0);
    norm_sum += pow(my_xyz.y - nbr_feature[4],2.0);
    norm_sum += pow(my_xyz.z - nbr_feature[5],2.0);
    pdensity[device_offset[blockId] + i + threadIdx.x] = pdensity_sum;
    imp_wt[device_offset[blockId] + i + threadIdx.x ] = norm_sum;
    if(threadIdx.x==0){
      printf("imp_wt: %f\n",norm_sum);
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
		 float* flattenRGB,int* voxel_offset, int* neighbor_ids,int x_idx,int y_idx,int z_idx){


  printf("HELLO IN DEVICE SETUP!\n");

  float *device_xyz, *device_rgb;
  int* device_offset, *device_neighbor_ids;
  float *imp_wt, *pdensity;

    cudaMalloc(&device_xyz, num_pts*3*sizeof(float));
  
  cudaMalloc(&device_xyz, num_pts*3*sizeof(float));
  cudaMalloc(&device_rgb, num_pts*3*sizeof(float));
  cudaMalloc(&device_offset, (num_voxels+1)*sizeof(int));
  cudaMalloc(&device_neighbor_ids , num_voxels * 7 * sizeof(int));

  cudaMalloc(&imp_wt, num_pts*sizeof(float));
  cudaMalloc(&pdensity, num_pts*sizeof(float));

  printf("finished mallocing!\n");
  gpuErrchk(cudaMemcpy(device_offset,voxel_offset, (num_voxels+1)*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_xyz, flattenXYZ, num_pts*3*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_rgb, flattenRGB, num_pts*3*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_neighbor_ids, neighbor_ids, num_voxels*7*sizeof(int),cudaMemcpyHostToDevice));
  int xy_idx = x_idx * y_idx;

  dim3 gridDim(x_idx,y_idx,z_idx);
  dim3 blockDim(256,1,1);
  printf("about to call kernel\n");
  sampling<<<gridDim,blockDim>>>(device_xyz,device_rgb,device_offset,device_neighbor_ids,xy_idx,y_idx,imp_wt,pdensity);
  printf("finished sampling!\n");
  cudaThreadSynchronize();
  float* host_imp_wt;
  host_imp_wt = (float*) malloc(num_pts*sizeof(float));
  gpuErrchk(cudaMemcpy(host_imp_wt,imp_wt,num_pts*sizeof(float),cudaMemcpyDeviceToHost));

  for(int i = 0 ; i < 100; i++){
    printf("imp_wt[%d] = %f\n",i,host_imp_wt[i]);
  }

  return(1);
}




