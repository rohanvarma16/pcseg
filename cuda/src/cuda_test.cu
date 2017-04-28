#include <stdio.h>
#include <vector>

#include <algorithm>
#include <functional>
#include <cuda_runtime.h>
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include <math.h>
#include <cuda.h>
#include <float.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/sort.h> 
#include <thrust/gather.h>
#include <thrust/execution_policy.h>
#include <iostream> 
#include <thrust/binary_search.h>
#include <thrust/random.h>
#include <curand.h>
#include <time.h>
#include <stdlib.h>

#define THREADS_PER_BLOCK 256

__global__
void saxpy(int n, float a, float *x, float *y)
{
  if( blockIdx.x==1 && threadIdx.x ==1){
  printf("hello in kernel! \n");
  }
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}

 __device__ __inline__ float my_exp(float a) {
  float b;
  //  b = exp(a);
  b = 2*a;
  return b;
}


struct GenUnifRands
{
    __device__
    float operator () (int idx)
  {
    thrust::default_random_engine randEng;
    thrust::uniform_real_distribution<float> uniDist;
    randEng.discard(idx);
    return uniDist(randEng);
  }
};


__global__
void segmentation_kernel(float* device_xyz,float* device_rgb,int* device_offset,
			 int* neighbor_id,float* device_pdens,int yz_idx,
			 int z_idx,int* parents,float* distances){
  int block_i = blockIdx.x;
  int block_j = blockIdx.y;
  int block_k = blockIdx.z;

  int blockId = yz_idx * block_i + z_idx* block_j + block_k ;
  int threadId = threadIdx.x;

  if( blockId==0 && threadId==0){
    printf("hello in kernel! \n");
  }

  int num_threads = THREADS_PER_BLOCK;
  int my_num_pts = device_offset[blockId + 1] - device_offset[blockId];
  float3 nbr_xyz;
  float3 my_xyz;
  int my_idx;
  int nbr_idx;
  // Iterate over all points in voxel, num_threads points at a time.


  for(int i = 0 ; i < my_num_pts ; i+= num_threads){

    if(threadId + i > my_num_pts){
      return;
    }
    my_idx = device_offset[blockId] + i + threadId;
    //my_idx = 0;
    my_xyz = *(float3*) &device_xyz[3*my_idx];
   

    // Iterate over the neighbouring blocks including yourself.
    float min_distance = 10000000.0;
    float current_parent = my_idx;
    
    for(int j = 0 ; j < 7 ; j++){
      if(neighbor_id[7*blockId+j] == -1){
	continue;
      }

      int nbr_num_pts = device_offset[neighbor_id[7*blockId +j] + 1] - device_offset[neighbor_id[7*blockId + j]];
      
      
      // Iterate over all the points in the neighbouring block.
      for(int k = 0; k < nbr_num_pts; k++){

	nbr_idx = device_offset[neighbor_id[blockId*7 + j]] + k;
	//	nbr_idx = 1;
	nbr_xyz = *(float3*) &device_xyz[3*nbr_idx];

	if(nbr_idx == my_idx){
	  continue;
	}

	float xyz_dist =  pow(my_xyz.x - nbr_xyz.x,2.0) + pow((my_xyz.y - nbr_xyz.y),2.0) +
	  pow((my_xyz.z - nbr_xyz.z),2.0);
	
	if(device_pdens[nbr_idx] > device_pdens[my_idx] && xyz_dist < min_distance){
	  min_distance = xyz_dist;
	  current_parent = nbr_idx;
	  //	  printf("I am here!");
	}
	
      }

    }
      
        distances[my_idx] = min_distance;
    parents[my_idx] = current_parent;
  }

}


__global__
void sampling(float *device_xyz, float *device_rgb, int *device_offset, 
              int *neighbor_id, int yz_idx, int z_idx, float *imp_wt, 
              float* pdensity)
{
  int block_i = blockIdx.x;
  int block_j = blockIdx.y;
  int block_k = blockIdx.z;

  int blockId = yz_idx * block_i + z_idx* block_j + block_k ;
  int threadId = threadIdx.x;
  
  if( blockId==0 && threadId==0){
    printf("hello in kernel! \n");
  }

  int num_threads = THREADS_PER_BLOCK;

  float pdensity_sum = 0.0;
  float nbr_feature[6];
  float norm_sum = 0.0;
  int my_num_pts = device_offset[blockId + 1] - device_offset[blockId];
  float3 nbr_rgb;
  float3 nbr_xyz;
  float3 my_rgb;
  float3 my_xyz;
  float sigma_sq = 0.00005;
  
  float Aij_ew;
  float Aij_ew2;
  // Iterate over all points in voxel, num_threads points at a time. 
    int sum_nbr_num_pts = 0;
  for(int i = 0 ; i < my_num_pts ; i+= num_threads){
    pdensity_sum = 0.0;
    norm_sum = 0.0;
    nbr_feature[0] = 0.0;
    nbr_feature[1] = 0.0;
    nbr_feature[2] = 0.0;
    nbr_feature[3] = 0.0;
    nbr_feature[4] = 0.0;
    nbr_feature[5] = 0.0;

    if(threadId + i > my_num_pts){
      return;
    }

    my_xyz = *(float3*) &device_xyz[3*(device_offset[blockId] + i + threadId)];
    my_rgb = *(float3*) &device_rgb[3*(device_offset[blockId] + i + threadId)];
    
    /*if(blockId == 32 && threadId == 0)
    {
        printf("my x %f, my y %f, my z %f\n", my_xyz.x, my_xyz.y, my_xyz.z);
        printf("my r %f, my g %f, my b %f\n", my_rgb.x, my_rgb.y, my_rgb.z);
    }*/

    // Iterate over the neighbouring blocks including yourself. 
    if(threadId == 0){
    sum_nbr_num_pts = 0;
    }
    for(int j = 0 ; j < 7 ; j++){
      if(neighbor_id[7*blockId+j] == -1){
	    continue;
      }

      int nbr_num_pts = device_offset[neighbor_id[7*blockId +j] + 1] - device_offset[neighbor_id[7*blockId + j]];
      
      if(threadId ==0){
      sum_nbr_num_pts += nbr_num_pts;
      }
      // Iterate over all the points in the neighbouring block.
      for(int k = 0; k < nbr_num_pts; k++){
        
	    nbr_xyz = *(float3*) &device_xyz[3*(device_offset[neighbor_id[blockId*7 + j]] + k)];
	    nbr_rgb = *(float3*) &device_rgb[3*(device_offset[neighbor_id[blockId*7 + j]] + k)];
    
        /*if(blockId == 32 && threadId == 0)
        {
        
  printf("neighbour block = %d, offset = %d\n", neighbor_id[blockId*7 + j],
                                3*(device_offset[neighbor_id[blockId*7+j]] + k)); 
          printf("nbr x %f, nbr y %f, nbr z %f\n", nbr_xyz.x, nbr_xyz.y, nbr_xyz.z);
          printf("nbr r %f, nbr g %f, nbr b %f\n", nbr_rgb.x, nbr_rgb.y, nbr_rgb.z);
        }*/

	    float xyz_dist =  pow(my_xyz.x - nbr_xyz.x,2.0) + pow((my_xyz.y - nbr_xyz.y),2.0) +
                         pow((my_xyz.z - nbr_xyz.z),2.0);
 
        float rgb_dist =  pow((my_rgb.x - nbr_rgb.x),2.0) + pow((my_rgb.y - nbr_rgb.y),2.0) +
						 pow((my_rgb.z - nbr_rgb.z),2.0);
						     
	Aij_ew =  __expf(-1.0f * (xyz_dist/sigma_sq));
	Aij_ew2 = __expf(-1.0f * (xyz_dist + rgb_dist)/sigma_sq); 
	//Aij_ew = pow(xyz_dist,2.0);
	//Aij_ew2 = pow(xyz_dist + rgb_dist,2.0);
	    if(threadId == 1 && blockId==4){
	      printf("xyz_dst: %f , rgb_dst: %f, Aij_ew: %f, Aij_ew2: %f\n",xyz_dist,rgb_dist,Aij_ew, Aij_ew2);
	    }
	    pdensity_sum += Aij_ew2;
		nbr_feature[0] += Aij_ew * nbr_xyz.x;
	    nbr_feature[1] += Aij_ew * nbr_xyz.y;
		nbr_feature[2] += Aij_ew * nbr_xyz.z;
		nbr_feature[3] += Aij_ew * nbr_rgb.x;
	    nbr_feature[4] += Aij_ew * nbr_rgb.y;
	    nbr_feature[5] += Aij_ew * nbr_rgb.z; 
      }
    }
    
    norm_sum += pow(my_xyz.x - nbr_feature[0],2.0);
    norm_sum += pow(my_xyz.y - nbr_feature[1],2.0);
    norm_sum += pow(my_xyz.z - nbr_feature[2],2.0);
    norm_sum += pow(my_rgb.x - nbr_feature[3],2.0);
    norm_sum += pow(my_rgb.y - nbr_feature[4],2.0);
    norm_sum += pow(my_rgb.z - nbr_feature[5],2.0);
    
    pdensity[device_offset[blockId] + i + threadId] = pdensity_sum;

    imp_wt[device_offset[blockId] + i + threadId] = norm_sum;
    
  }
  /*
  if((threadId == 0)){
    printf("my_num_pts: %d ,my_sum_nbr_pts:%d, blockId: %d \n", my_num_pts,sum_nbr_num_pts, blockId );
    }
  */
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


int device_setup(int num_pts, int num_voxels,   float *flattenXYZ,
		 float *flattenRGB, int *voxel_offset, int *neighbor_ids,
		 int x_idx, int y_idx, int z_idx,int num_samples,uint *sample_arr,float *pdens)
{
  printf("HELLO IN DEVICE SETUP!\n");

  float *device_xyz, *device_rgb;
  int *device_offset, *device_neighbor_ids;

  //  cudaMalloc(&device_xyz, num_pts*3*sizeof(float));
  cudaMalloc(&device_xyz, num_pts*3*sizeof(float));
  cudaMalloc(&device_rgb, num_pts*3*sizeof(float));
  cudaMalloc(&device_offset, (num_voxels+1)*sizeof(int));
  cudaMalloc(&device_neighbor_ids , num_voxels * 7 * sizeof(int));
  //  cudaMalloc(&imp_wt, num_pts*sizeof(float));
  //cudaMalloc(&pdensity, num_pts*sizeof(float));
  //host_imp_wt = (float*) malloc(num_pts*sizeof(float));
  //host_pdensity = (float*) malloc(num_pts*sizeof(float));

  thrust::device_vector<float> dev_imp_wt(num_pts);
  thrust::device_vector<float> dev_pdensity(num_pts);
  thrust::host_vector<float> host_imp_wt(num_pts);
  thrust::host_vector<float> host_pdensity(num_pts);
  
  float* imp_wt = thrust::raw_pointer_cast(dev_imp_wt.data());
  float* pdensity = thrust::raw_pointer_cast(dev_pdensity.data());
  printf("finished mallocing!\n");

  


  gpuErrchk(cudaMemcpy(device_offset,voxel_offset, (num_voxels+1)*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_xyz, flattenXYZ, num_pts*3*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_rgb, flattenRGB, num_pts*3*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_neighbor_ids, neighbor_ids, num_voxels*7*sizeof(int),cudaMemcpyHostToDevice));
  
  int yz_idx = y_idx * z_idx;

  dim3 gridDim(x_idx,y_idx,z_idx);
  dim3 blockDim(THREADS_PER_BLOCK,1,1);
  
  printf("about to call kernel\n");
  sampling<<<gridDim,blockDim>>>(device_xyz,device_rgb,device_offset,device_neighbor_ids,yz_idx,z_idx,imp_wt,pdensity);
  printf("finished sampling!\n");
  
  cudaDeviceSynchronize();
  
  thrust::copy(dev_imp_wt.begin(),dev_imp_wt.end(),host_imp_wt.begin());
  thrust::copy(dev_pdensity.begin(),dev_pdensity.end(),host_pdensity.begin());


  //  gpuErrchk(cudaMemcpy(host_imp_wt,imp_wt,num_pts*sizeof(float),cudaMemcpyDeviceToHost));
  //gpuErrchk(cudaMemcpy(host_pdensity,pdensity,num_pts*sizeof(float),cudaMemcpyDeviceToHost));
  
  for(int i = 0 ; i < 10; i++){
    printf("imp_wt[%d] = %0.9f pdensity[%d] = %f\n",i,host_imp_wt[i],i,host_pdensity[i]);
  }
  

  /*
  thrust::host_vector<int> h_vec(20); thrust::generate(h_vec.begin(), h_vec.end(), rand);
  for(int i = 0; i < 20 ; i++){
    printf("i:%d, %d\n" , i, h_vec[i]);
  }
  // transfer data to the device
  thrust::device_vector<int> d_vec = h_vec;
  // sort data on the device (805 Mkeys/sec on GeForce GTX 480)
  thrust::sort(d_vec.begin(), d_vec.end()); // transfer data back to host
  thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
  for(int i = 0 ; i < 20 ; i++){
    printf("i:%d, %d\n" , i, h_vec[i]);
}
  */




  // WEIGHTED SAMPLING


  
  printf("====================\n");
  printf("going to do weighted sampling\n");
  // step 1: normalize: (TODO: use functor to make fast )
  float norm_sum = thrust::reduce(dev_imp_wt.begin(),dev_imp_wt.end());
  float norm_factor = 1.0/norm_sum;
  printf("norm_sum: %f\n",norm_sum);

  thrust::device_vector<float> temp(dev_imp_wt.size());
  thrust::fill(temp.begin(),temp.end(),norm_factor);
  thrust::transform(dev_imp_wt.begin(),dev_imp_wt.end(),temp.begin(),dev_imp_wt.begin(),thrust::multiplies<float>());

  thrust::copy(dev_imp_wt.begin(),dev_imp_wt.end(),host_imp_wt.begin());
  printf("\n normalized weights \n");
    for(int i = 0 ; i < 5 ; i++){
      printf("h_wt[%d] = %f , ",i,host_imp_wt[i]);
    }

    // step 2: compute prefix sum (clusive scan)
    thrust::device_vector<float> wts_rs(dev_imp_wt.size());
    thrust::inclusive_scan(dev_imp_wt.begin(),dev_imp_wt.end(),wts_rs.begin());
    printf("\n rolling_sum \n");
    thrust::copy(wts_rs.begin(),wts_rs.end(),host_imp_wt.begin());
    for(int i = 0 ; i < 5 ; i++){
      printf("h_wt[%d] = %0.9f , ",i,host_imp_wt[i]);
    }

    // step 3: generate uniform random numbers:
    srand(time(NULL));
    int seed = rand();
    printf("\n get random samples \n");

    thrust::device_vector<float> d_unifrands(num_samples);
    thrust::transform( thrust::make_counting_iterator(seed), thrust::make_counting_iterator(seed + num_samples),
                       d_unifrands.begin(),GenUnifRands());


    // step 4 : generate (weighted) random samples
    thrust::device_vector<unsigned int> samples(num_samples);
    thrust::lower_bound(wts_rs.begin(),wts_rs.end(),d_unifrands.begin(),d_unifrands.end(),samples.begin());
    thrust::host_vector<unsigned int> h_samples(num_samples);
    thrust::copy(samples.begin(),samples.end(),h_samples.begin());
    for(int i = 0 ; i < 10 ; i++){
      printf("h_samples[%d] = %d " , i , h_samples[i]);
    }
   
    printf("\n");
    uint* pc_samples = thrust::raw_pointer_cast(h_samples.data());
    //uint* samples_arr = (uint*) malloc(num_samples* sizeof(uint));
    memcpy(sample_arr, pc_samples,num_samples*sizeof(uint));
    
    
    thrust::device_vector<float> dev_pdensity_rs(num_samples);
    thrust::host_vector<float> host_pdensity_rs(num_samples);
    thrust::gather(thrust::device,samples.begin(),samples.end(),dev_pdensity.begin(),dev_pdensity_rs.begin());
    thrust::copy(dev_pdensity_rs.begin(),dev_pdensity_rs.end(),host_pdensity_rs.begin());
    float* host_pdensity_rs_ptr = thrust::raw_pointer_cast(host_pdensity_rs.data());
    memcpy(pdens, host_pdensity_rs_ptr,num_samples*sizeof(float));
    
    /*
    for(int i = 0 ; i < 10; i++){
      printf("samples_arr[%d] = %d , " , i ,samples_arr[i] );
    }
    */

    /*
    thrust::host_vector<float> h_unifrands(num_samples);
    thrust::copy(d_unifrands.begin(),d_unifrands.end(),h_unifrands.begin());
    for(int i = 0 ; i < 10 ; i++){
      printf("unifrands[%d] = %f " , i , h_unifrands[i]);
    }
    */

    // step 5: use a gather operator to "compress" xyz rgb into appropriate form:

  return(1);
}






int segmentation(int num_pts, int num_voxels, float* pdens,  float *flattenXYZ,
                 float *flattenRGB, int *voxel_offset, int *neighbor_ids,
                 int x_idx, int y_idx, int z_idx,int* parents_ptr)
{
  printf("HELLO IN SEGMENTATION KERNEL!\n");

  float *device_xyz, *device_rgb, *device_pdens;
  int *device_offset, *device_neighbor_ids;


  cudaMalloc(&device_xyz, num_pts*3*sizeof(float));
  cudaMalloc(&device_rgb, num_pts*3*sizeof(float));
  cudaMalloc(&device_offset, (num_voxels+1)*sizeof(int));
  cudaMalloc(&device_neighbor_ids , num_voxels * 7 * sizeof(int));
  cudaMalloc(&device_pdens, num_pts*sizeof(float));

  thrust::device_vector<int> dev_parents(num_pts);
  thrust::host_vector<int> host_parents(num_pts);
  thrust::device_vector<float> dev_distances(num_pts);
  thrust::host_vector<float> host_distances(num_pts);

  int* parents = thrust::raw_pointer_cast(dev_parents.data());
  float* distances = thrust::raw_pointer_cast(dev_distances.data());
  printf("finished mallocing!\n");

  gpuErrchk(cudaMemcpy(device_offset,voxel_offset, (num_voxels+1)*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_xyz, flattenXYZ, num_pts*3*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_rgb, flattenRGB, num_pts*3*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_neighbor_ids, neighbor_ids, num_voxels*7*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(device_pdens, pdens, num_pts * sizeof(float), cudaMemcpyHostToDevice));

  int yz_idx = y_idx * z_idx;
  dim3 gridDim(x_idx,y_idx,z_idx);
  dim3 blockDim(THREADS_PER_BLOCK,1,1);

  printf("about to call kernel\n");
  segmentation_kernel<<<gridDim,blockDim>>>(device_xyz,device_rgb,device_offset,device_neighbor_ids,device_pdens,yz_idx,z_idx,parents,distances);
  printf("finished segmentation!\n");
  
  int num_gather = 6;
  thrust::device_vector<int> temp_1(num_pts);
  thrust::device_vector<int> temp_2(num_pts);

  thrust::copy(dev_parents.begin(),dev_parents.end(),temp_1.begin());
  for(int i  = 0 ; i < num_gather ; i++){
    thrust::gather(thrust::device, dev_parents.begin(),dev_parents.end(),temp_1.begin(),temp_2.begin());
    thrust::copy(temp_2.begin(),temp_2.end(),temp_1.begin());
    thrust::copy(temp_2.begin(),temp_2.end(),dev_parents.begin());
  }
  

  printf("finished tree cutting!\n");
  thrust::copy(dev_parents.begin(),dev_parents.end(),host_parents.begin());
   int* host_parents_ptr = thrust::raw_pointer_cast(host_parents.data());
   memcpy(parents_ptr, host_parents_ptr,num_pts*sizeof(int));
   
  for(int i = 0 ; i < 5 ; i++){
    printf("parents[%d] = %d \n", i ,host_parents[i]);
  }
  
  
  printf("DONE!!!!");
  cudaDeviceSynchronize();
  return(1);
}
