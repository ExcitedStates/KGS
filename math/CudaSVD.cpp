#include "CudaSVD.h"

#include <algorithm>
#include <iostream>

using namespace std;

#ifdef __CUDA_SVD

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cusolverDn.h>
#include <cuda_runtime_api.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_double.h>


void gpuAssert(cudaError_t code, char *file, int line, bool abort=true);
void gpuErrchk(cudaError_t ans);


void CudaSVD::UpdateFromMatrix()
{

  int O = std::min(m,n);

  // Setting the device matrix and moving the host matrix to the device
  double *d_A;
  gpuErrchk(cudaMalloc(&d_A,     m * n * sizeof(double)));
  gpuErrchk(cudaMemcpy(d_A, matrix->data, m * n * sizeof(double), cudaMemcpyHostToDevice));

  // Device side SVD workspace and matrices
  int work_size = 0;

  int *devInfo;  gpuErrchk(cudaMalloc(&devInfo,      sizeof(int)));
  double *d_U;   gpuErrchk(cudaMalloc(&d_U,  m * m * sizeof(double)));
  double *d_V;   gpuErrchk(cudaMalloc(&d_V,  n * n * sizeof(double)));
  double *d_S;   gpuErrchk(cudaMalloc(&d_S,      O * sizeof(double)));

  cusolverStatus_t stat;

  // CUDA solver initialization
  cusolverDnHandle_t solver_handle;
  cusolverDnCreate(&solver_handle);

  stat = cusolverDnDgesvd_bufferSize(solver_handle, m, n, &work_size);
  if(stat != CUSOLVER_STATUS_SUCCESS ) std::cout << "Initialization of cuSolver failed."<<endl;

  double *work;    gpuErrchk(cudaMalloc(&work, work_size * sizeof(double)));

  // CUDA SVD execution
  stat = cusolverDnDgesvd(solver_handle, 'A', 'A', m, n, d_A, m, d_S, d_U, m, d_V, n, work, work_size, NULL, devInfo);
  cudaDeviceSynchronize();

  int devInfo_h = 0;
  gpuErrchk(cudaMemcpy(&devInfo_h, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
  std::cout << "devInfo = " << devInfo_h << "\n";

  switch(stat){
    case CUSOLVER_STATUS_SUCCESS:           std::cout << "SVD computation success\n";                       break;
    case CUSOLVER_STATUS_NOT_INITIALIZED:   std::cout << "Library cuSolver not initialized correctly\n";    break;
    case CUSOLVER_STATUS_INVALID_VALUE:     std::cout << "Invalid parameters passed\n";                     break;
    case CUSOLVER_STATUS_INTERNAL_ERROR:    std::cout << "Internal operation failed\n";                     break;
    case CUSOLVER_STATUS_ALLOC_FAILED:      std::cout << "Internal operation failed\n";                     break;
    case CUSOLVER_STATUS_ARCH_MISMATCH:     std::cout << "Internal operation failed\n";                     break;
    case CUSOLVER_STATUS_MAPPING_ERROR:     std::cout << "Internal operation failed\n";                     break;
    case CUSOLVER_STATUS_EXECUTION_FAILED:  std::cout << "Internal operation failed\n";                     break;
    default: std::cout<<"Unknown CUDA error"<<endl;
  }

  if (devInfo_h == 0 && stat == CUSOLVER_STATUS_SUCCESS) std::cout    << "SVD successful\n\n";

  // Moving the results from device to host
  gpuErrchk(cudaMemcpy(S->data, d_S,     O * sizeof(double), cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(V->data, d_V, n * n * sizeof(double), cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(U->data, d_U, m * m * sizeof(double), cudaMemcpyDeviceToHost));

  cusolverDnDestroy(solver_handle);


  // Free memory
  // ...
}

void gpuAssert(cudaError_t code, char *file, int line, bool abort)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) { exit(code); }
   }
}
void gpuErrchk(cudaError_t ans) { gpuAssert((ans), __FILE__, __LINE__); }


#else


void CudaSVD::UpdateFromMatrix()
{
  cerr<< "CudaSVD::updateFromMatrix error! CUDA not supported. Install CUDA and compile with -D__CUDA_SVD"<<endl;
  throw "CudaSVD::updateFromMatrix error! CUDA not supported. Install CUDA and compile with -D__CUDA_SVD";
}

#endif

