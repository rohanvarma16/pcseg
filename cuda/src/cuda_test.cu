#include <stdio.h>

__global__
void saxpy(int n, float a, float *x, float *y)
{
  if( blockIdx.x==1 && threadIdx.x ==1){
  printf("hello in kernel! \n");
  }
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
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
int run_main(void)
{
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
}



