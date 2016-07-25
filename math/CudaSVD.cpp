#include "CudaSVD.h"

#include <algorithm>
#include <iostream>

using namespace std;

#ifdef __CUDA_SVD

#include <cuda_runtime.h>
#include <cusolverDn.h>


void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true);
void gpuErrchk(cudaError_t ans);


void CudaSVD::UpdateFromMatrix()
{

  int O = std::min(m,n);

  gsl_matrix* matrixTranspose = gsl_matrix_alloc(n,m);
  gsl_matrix_transpose_memcpy(matrixTranspose, matrix);

  // Setting the device matrix and moving the host matrix to the device
  double *d_A;
  gpuErrchk(cudaMalloc(&d_A,     m * n * sizeof(double)));
  gpuErrchk(cudaMemcpy(d_A, matrixTranspose->data, m * n * sizeof(double), cudaMemcpyHostToDevice));

//  double *d_C;
//  gpuErrchk(cudaMalloc(&d_C,     m * n * sizeof(double)));
//  double zero = 0;
//  double one = 1;
//
//
//  //Transpose d_A as cusolver assumes col-major order
//  cublasStatus_t transStat;
//  cublasHandle_t trans_handle;
//  cublasCreate(&trans_handle);
//  double *d_alpha;   gpuErrchk(cudaMalloc(&d_alpha, sizeof(double)));
//  double *d_beta;    gpuErrchk(cudaMalloc(&d_beta , sizeof(double)));
//  gpuErrchk(cudaMemcpy(d_alpha, &one , sizeof(double), cudaMemcpyHostToDevice));
//  gpuErrchk(cudaMemcpy(d_beta,  &zero, sizeof(double), cudaMemcpyHostToDevice));
//  transStat = cublasDgeam(trans_handle, CUBLAS_OP_T, CUBLAS_OP_N, m, n, d_alpha, d_A, m, d_beta, NULL, m, d_C, m);

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

  stat = cusolverDnDgesvd(solver_handle, 'A', 'A', m, n, d_A, m, d_S, d_U, m, d_V, n, work, work_size, NULL, devInfo);
  cudaDeviceSynchronize();

  int devInfo_h = 0;
  gpuErrchk(cudaMemcpy(&devInfo_h, devInfo, sizeof(int), cudaMemcpyDeviceToHost));

//  switch(transStat){
//    case CUBLAS_STATUS_SUCCESS: std::cout<<"Transpose succesful\n"; break;
//    case CUBLAS_STATUS_INVALID_VALUE: std::cout<<"Transpose failed: Invalid value\n"; break;
//    default: std::cout<<"Transpose failed (code "<<transStat<<")\n";
//  }
  switch(stat){
    case CUSOLVER_STATUS_SUCCESS:           break; //std::cout << "SVD computation success\n";                       break;
    case CUSOLVER_STATUS_NOT_INITIALIZED:   std::cout << "Library cuSolver not initialized correctly\n";    break;
    case CUSOLVER_STATUS_INVALID_VALUE:     std::cout << "Invalid parameters passed\n";                     break;
    case CUSOLVER_STATUS_INTERNAL_ERROR:    std::cout << "Internal operation failed\n";                     break;
    case CUSOLVER_STATUS_ALLOC_FAILED:      std::cout << "Internal operation failed\n";                     break;
    case CUSOLVER_STATUS_ARCH_MISMATCH:     std::cout << "Internal operation failed\n";                     break;
    case CUSOLVER_STATUS_MAPPING_ERROR:     std::cout << "Internal operation failed\n";                     break;
    case CUSOLVER_STATUS_EXECUTION_FAILED:  std::cout << "Internal operation failed\n";                     break;
    default: std::cout<<"Unknown CUDA error"<<endl;
  }

//  if (devInfo_h == 0 && stat == CUSOLVER_STATUS_SUCCESS) std::cout    << "SVD successful\n\n";

  // Moving the results from device to host
  gpuErrchk(cudaMemcpy(S->data, d_S,     O * sizeof(double), cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(V->data, d_V, n * n * sizeof(double), cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(U->data, d_U, m * m * sizeof(double), cudaMemcpyDeviceToHost));

  cusolverDnDestroy(solver_handle);

  gsl_matrix_transpose(U); //In-place transpose as U comes out in col-major format

  // Free memory
  gsl_matrix_free(matrixTranspose);
}

void gpuAssert(cudaError_t code, const char *file, int line, bool abort)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) { exit(code); }
   }
}
void gpuErrchk(cudaError_t ans) {
  const string& f = __FILE__;
  gpuAssert((ans), f.c_str(), __LINE__);
}


#else


void CudaSVD::UpdateFromMatrix()
{
  cerr<< "CudaSVD::updateFromMatrix error! CUDA not supported. Install CUDA and compile with -D__CUDA_SVD"<<endl;
  throw "CudaSVD::updateFromMatrix error! CUDA not supported. Install CUDA and compile with -D__CUDA_SVD";
}

#endif

