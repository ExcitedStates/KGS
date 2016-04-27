
#include "gsl_helpers.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_double.h>

#include "Logger.h"
#include "MathUtility.h"

using namespace std;

void gsl_matrix_cout (const gsl_matrix *m) {
  for (int i=0; i<m->size1; ++i) {
    for (int j=0; j<m->size2; ++j){
      log() << gsl_matrix_get(m,i,j);
      if(j == (m->size2-1))
        log() << ";" << endl;
      else
        log() << ",\t";
    }
  }
}

void gsl_matrix_outtofile (const gsl_matrix *m, const string& filename) {
  ofstream output( filename.c_str() );
  if(!output.is_open()) {
    cerr<<"Cannot write to "<<filename<<". You might need to create output directory first"<<endl;
    exit(-1);
  }
  for (int i=0; i<m->size1; ++i) {
    for (int j=0; j<m->size2; ++j) {
      output << setprecision(14) << gsl_matrix_get(m,i,j); //<< setw(20)
      //output << gsl_matrix_get(m,i,j) << " ";
      if(j == (m->size2-1))
        output << endl;
      else
        output << " ";
    }
  }
  output.close();
}

void gsl_vector_cout (const gsl_vector *v) {
  for (int i=0; i<v->size; ++i)
    cout << gsl_vector_get(v,i) << "\t";
  cout << endl;
  //    log() << gsl_vector_get(v,i) << "\t";
  //log() << endl;
}

void gsl_vector_outtofile (const gsl_vector *v, const string& filename) {
    ofstream output(filename.c_str());
    for (int i=0; i<v->size; ++i)
        output << setprecision(14) << gsl_vector_get(v,i) << endl;
    //output << gsl_vector_get(v,i) << endl;
    output.close();
}


double gsl_vector_length(const gsl_vector *v) {
  return gsl_blas_dnrm2(v);
//  double len = 0;
//
//  for (int i=0; i<v->size; ++i)
//    len += gsl_vector_get(v,i)*gsl_vector_get(v,i);
//
//  len = sqrt(len);
//
//  return len;
}

void gsl_vector_normalize(gsl_vector *v) {
  double norm = gsl_vector_length(v);
  gsl_vector_scale(v,1.0/norm);
}

void gsl_vector_scale_to_length(gsl_vector* ret, double length)
{
  gsl_vector_scale(ret, length/gsl_blas_dnrm2(ret));
}

void gsl_vector_scale_max_component(gsl_vector* v, double maxComponent)
{
  assert(maxComponent>=0);
  double curMax = gsl_vector_max(v);
  double curMin = gsl_vector_min(v);
  curMax = std::max(curMax,-curMin);
  if(curMax>maxComponent)
    gsl_vector_scale(v, maxComponent/curMax);
}

gsl_matrix* gsl_matrix_trans(gsl_matrix* A){
  int M = A->size1;
  int N = A->size2;
  gsl_matrix* ret = gsl_matrix_alloc(A->size2, A->size1);
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      gsl_matrix_set(ret,j,i, gsl_matrix_get(A, i,j));
    }
  }
  return ret;
}

gsl_vector* gsl_vector_copy(gsl_vector* v){
  gsl_vector* ret = gsl_vector_alloc(v->size);
  gsl_vector_memcpy(ret, v);
  return ret;
}

gsl_vector* gsl_matrix_vector_mul(gsl_matrix* A, gsl_vector* v){
  int MA = A->size1;
  int NA = A->size2;
  int Mv = v->size;
  if(NA!=Mv){
    cerr<<"Mismatching sizes: A: "<<MA<<"x"<<NA<<" , v: "<<Mv<<"x1"<<endl;
    exit(-1);
  }
  gsl_vector* ret = gsl_vector_alloc(MA);
  for(int i=0;i<MA;i++){
    double val = 0;
    for(int j=0;j<NA;j++){
      val+=gsl_matrix_get(A,i,j)*gsl_vector_get(v,j);
    }
    gsl_vector_set(ret, i, val);
  }

  return ret;
}
gsl_matrix* gsl_matrix_mul(gsl_matrix* A, gsl_matrix* B){
  int MA = A->size1;
  int NA = A->size2;
  int MB = B->size1;
  int NB = B->size2;
  if(NA!=MB){
    cerr<<"Mismatching sizes: A: "<<MA<<"x"<<NA<<" , B: "<<MB<<"x"<<NB<<endl;
    exit(-1);
  }
  gsl_matrix* ret = gsl_matrix_alloc(MA,NB);
  for(int i=0; i<MA; i++){
    for(int j=0; j<NB; j++){
      double val = 0;
      for(int k=0; k<NA; k++)
        val += gsl_matrix_get(A,i,k)*gsl_matrix_get(B,k,j);
      //c[i][j] = c[i][j] + a[i][k]*b[k][j];

      gsl_matrix_set(ret,i,j, val);
    }
  }
  return ret;
}


double gsl_vector_mean (const gsl_vector *v) {
  double sum = 0;
  for (int i=0; i<v->size; ++i)
    sum += gsl_vector_get(v,i);
  return sum/v->size;
}
//---------------------------------------------------------
gsl_matrix* pca (gsl_matrix* sample_matrix) {
  int m = sample_matrix->size1;
  int n = sample_matrix->size2;
  // subtract the column mean from each column
  for (int j=0; j<n; ++j) {
//		gsl_vector_view column_j = gsl_matrix_column(sample_matrix,j); // get run-time error possibily due to DLL setting
    gsl_vector* column_j = gsl_vector_calloc(m);
    for (int i=0; i<m; ++i)
      gsl_vector_set(column_j,i,gsl_matrix_get(sample_matrix,i,j));
    double column_mean = gsl_vector_mean(column_j);
    gsl_vector_add_constant(column_j,-column_mean);
    for (int i=0; i<m; ++i)
      gsl_matrix_set(sample_matrix,i,j,gsl_vector_get(column_j,i));
  }
  // compute the covariance matrix
  gsl_matrix* covariance_matrix = gsl_matrix_calloc(n,n);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,sample_matrix,sample_matrix,0.0,covariance_matrix);
  double ratio = 1.0/(m-1);
  gsl_matrix_scale(covariance_matrix,ratio);
  // do SVD on the covariance matrix
  gsl_matrix *V = gsl_matrix_calloc(n,n);
  gsl_vector *S = gsl_vector_calloc(n);
  gsl_vector *work = gsl_vector_calloc(n);
  gsl_linalg_SV_decomp(covariance_matrix,V,S,work);
  gsl_vector_free(S);
  gsl_vector_free(work);
  return V;
}

//gsl_vector* RandomUnitVectorInNullSpace (NullSpaceRet* nullspace) {
//	//int dimension = nullspace->ns->size;
//	int dimension = nullspace->nullspaceSize;
//	gsl_vector* random_vector = gsl_vector_calloc(dimension);
//	double *random_scaler = new double[dimension];
//	for (int i=0; i<dimension; ++i) {
//		random_scaler[i] = Random01();
//	}
//	for (int i=0; i<dimension; ++i) {
//		double value = 0;
//		for (int j=0; j<dimension; ++j) {
//			BUG alert: This should probably be:
//			BUG alert: value += gsl_matrix_get(nullspace->V,i,nullspace->n-dimension+j)*random_scaler[j];
//			value += gsl_matrix_get(nullspace->V,i,j)*random_scaler[j];
//			//value += gsl_matrix_get(nullspace->Svec,i,j)*random_scaler[j];
//		}
//		gsl_vector_set(random_vector,i,value);
//	}
//
//	// scale to length 1
//	double norm = gsl_vector_length(random_vector);
//	gsl_vector_scale(random_vector,1.0/norm);
//
//	return random_vector;
//}
//---------------------------------------------------------

gsl_vector* RandomUnitVector (int size) {
  gsl_vector* random_vector = gsl_vector_calloc(size);
  for (int i=0; i<size; ++i) {
    gsl_vector_set(random_vector,i,Random01());
  }
  gsl_vector_normalize(random_vector);
  return random_vector;
}

double frobenius_norm (const gsl_matrix *m) {
  double square_sum = 0;
  for (int i=0; i<m->size1; ++i) {
    for (int j=0; j<m->size2; ++j) {
      square_sum += pow(gsl_matrix_get(m,i,j),2);
    }
  }
  return pow(square_sum,0.5);
}

