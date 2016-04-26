/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

        Permission is hereby granted, free of charge, to any person obtaining a copy of
        this software and associated documentation files (the "Software"), to deal in
        the Software without restriction, including without limitation the rights to
        use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
        of the Software, and to permit persons to whom the Software is furnished to do
        so, subject to the following conditions:

        This entire text, including the above copyright notice and this permission notice
        shall be included in all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
        OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
        FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
        IN THE SOFTWARE.


*/

#include "JacobianRelated.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_blas.h>

#include "Logger.h"
#include "CTKTimer.h"
#include "Util.h"
#include "math/MathUtility.h"
// Intel MKL
#ifdef __INTEL_MKL
#include <mkl_lapack.h>
#include <math/gsl_helpers.h>
#include "freeMKL.h"
#endif
// CULA
#ifdef __GPU_CULA
#include <cula.hpp>
#endif

using namespace std;

gsl_matrix* NullSpaceRet::V = nullptr;
gsl_vector* NullSpaceRet::singularValues = nullptr;
//gsl_matrix* NullSpaceRet::Ut= nullptr;

NullSpaceRet::NullSpaceRet (int input_m, int input_n) : 
							m (input_m),
							n (input_n),
							m_rigidAngles( nullptr ),
							m_rigidHBonds( nullptr ),
							m_numCoordinated( 0 ),
							m_numRigid( 0 ),
							m_numRigidHBonds( 0 ),
							numOwners_( 0 ) 
{
	m_nullspaceBasis = nullptr;
	if( singularValues == nullptr){
		singularValues = gsl_vector_alloc(std::min(m,n));
	}else if( singularValues->size == std::min(m,n) ){
		//correct size already
	}else{
		gsl_vector_free( singularValues );
		singularValues = gsl_vector_alloc(std::min(m,n));
	}
	//V = gsl_matrix_alloc(n,n);//Works with old nr_svd
	//V = gsl_matrix_alloc(n,std::min(n,m));
	if(V==nullptr){
		V = gsl_matrix_alloc(n,n);//FULL
	}else if( V->size1==n && V->size2==n ){//FULL
	}else{
		gsl_matrix_free(V);
		V = gsl_matrix_alloc(n,n);//FULL
	}

	if(m_rigidAngles==nullptr){
		m_rigidAngles = gsl_vector_alloc(n);
	}
	if(m_rigidHBonds==nullptr){
		m_rigidHBonds = gsl_vector_alloc(n);
	}
}

NullSpaceRet::~NullSpaceRet () {
	gsl_matrix_free(m_nullspaceBasis);
	
	gsl_vector_free(m_rigidAngles);
	gsl_vector_free(m_rigidHBonds);
	//gsl_matrix_free(V);//Optimized away
}


// Accessor
unsigned int NullSpaceRet::numOwners () const {
	return numOwners_;
}

// Mutator
void NullSpaceRet::numOwners (unsigned int numOwners) {
	numOwners_ = numOwners;
}

void NullSpaceRet::NullSpacePostCompute()
{
	if(nullspaceSize == 0){
		cerr<<"No nullspace available for motion, returning."<<endl;
		m_nullspaceBasis = nullptr;
		//exit(-1);
	}
	else{
		gsl_matrix_view nullspaceBasis = gsl_matrix_submatrix(V,0,V->size2-nullspaceSize,n,nullspaceSize);
		m_nullspaceBasis = gsl_matrix_calloc(n,nullspaceSize);
		gsl_matrix_memcpy (m_nullspaceBasis, &nullspaceBasis.matrix);
	}
}

void NullSpaceRet::print () {
//	log()<< "n_ns = " << n_ns << endl;
//	log()<< "ns (" << ns->size << "): ";
//	    gsl_vector_cout(ns);
//	log() << "Sval (" << Sval->size << "): ";
//	    gsl_vector_cout(Sval);
//	log() << "Svec (" << Svec->size1 << "x" << Svec->size2 << "):" << endl;
//	gsl_matrix_cout(Svec);
//	log()<< "U (" << U->size1 << "x" << U->size2 << "):" << endl;
//	    gsl_matrix_cout(U);
//	log() << "V (" << V->size1 << "x" << V->size2 << "):" << endl;
//	    gsl_matrix_cout(V);
}

//---This file analyzes which dihedrals and hydrogen bonds are rigidified by constraints
void NullSpaceRet::RigidityAnalysis(gsl_matrix* HBondJacobian)
{
  /// First, check the dihedral angles for rigidity
  ///Nullspace dimension
  gsl_vector* currentRow = gsl_vector_alloc(n);
  double val=0.0;

  gsl_vector_set_zero(m_rigidAngles);
  gsl_vector_set_zero(m_rigidHBonds);

  bool moving=false;
  m_numCoordinated = 0;
  m_numRigid = 0;

  for(int i=0; i<n; i++){
    gsl_matrix_get_row(currentRow,V,i);

    for(int j=n-nullspaceSize; j<n; j++){
      val = fabs( gsl_vector_get(currentRow,j) );
      if( val > RIGID_TOL ){
        moving = true;
      }
    }

    /// at least one entry was greater than threshold --> dihedral is coordinated
    if( moving ) {
      m_numCoordinated++;
    }
    else {
      gsl_vector_set(m_rigidAngles,i,1); /// binary list of rigid dihedrals, complement to coordinated version
      m_numRigid++;
    }
    /// reset for the next dihedral / row
    moving = false;
  }

  log("constraints")<<"There are "<<m_numRigid << " rigidified and " << m_numCoordinated << " coordinated dihedrals" << endl;

	///Free memory
	gsl_vector_free(currentRow);

	/// Now, check the hydrogen Bonds for rigidity
	int numHBonds = HBondJacobian->size1;
	gsl_matrix* hBondNullspace = gsl_matrix_alloc(numHBonds,nullspaceSize);
	gsl_vector* currentHBondRow = gsl_vector_alloc(nullspaceSize);

	///Calculate the "In-Nullspace" Rotation of the hBonds
	gsl_matrix_view nullspaceBasis = gsl_matrix_submatrix(V,0,V->size2-nullspaceSize,n,nullspaceSize); //Matrix, top-left element, num rows, num cols
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, HBondJacobian, &nullspaceBasis.matrix, 0.0, hBondNullspace);

	//Only for test analysis: write hBondNullspace=J_h*N to file
//	string outNullH="KINARI_Comparison/output/hBondnullSpace.txt";
//	gsl_matrix_outtofile(hBondNullspace, outNullH);

	double maxVal;
	double minVal;
	m_numRigidHBonds = 0;

	for(int i=0; i<numHBonds; i++){
			gsl_matrix_get_row(currentHBondRow,hBondNullspace,i);
			gsl_vector_minmax(currentHBondRow, &minVal, &maxVal);

			if (minVal > -RIGID_TOL && maxVal < RIGID_TOL){
				gsl_vector_set(m_rigidHBonds,i,1);
				m_numRigidHBonds++;
			}
	}

//	log("constraints")<<"There are "<<m_numRigidHBonds<<" rigid out of "<<numHBonds<<" hydrogen bonds!"<<endl;
	cout<<"JacobianRelated.cpp: There are "<<m_numRigidHBonds<<" rigid out of "<<numHBonds<<" hydrogen bonds!"<<endl;

	gsl_vector_free(currentHBondRow);
	gsl_matrix_free(hBondNullspace);
}
//---------------------------------------------------------
void NullSpaceRet::ProjectOnNullSpace (gsl_vector *to_project, gsl_vector *after_project) {

	//The projection equation is:   after_project = N * N^T * to_project
	//Dimensions are: [n x 1] = [n x (n-r) ] [(n-r) x n] [n x 1]
	//Computations: N^T * to_project:	(n-r) * n  --> firstResult [(n-r) x 1]
	//Computations: N * firstResult:	(n) * (n-r)  --> finalResult [n x 1]
	//Computations overall: 2*n*(n-r)

	gsl_vector* firstResult = gsl_vector_alloc(nullspaceSize);
	gsl_blas_dgemv (CblasTrans, 1.0, m_nullspaceBasis, to_project, 0.0, firstResult);
	gsl_blas_dgemv (CblasNoTrans, 1.0, m_nullspaceBasis, firstResult, 0.0, after_project);
	gsl_vector_free(firstResult);


	//Old version
//	gsl_blas_dgemv (CblasNoTrans, 1.0, nullspace->P, to_project, 0.0, after_project);
}

#ifdef __INTEL_MKL
//---------------------------------------------------------
// INTEL MKL SVD implementation
//---------------------------------------------------------
void nr_svd (gsl_matrix* jac, gsl_vector* s_out, gsl_matrix* v_out) {
	const int m = jac->size1;
	const int n = jac->size2;

	const char* jobu  = "A";    // 'N'o columns of U are computed and returned in array U
	const char* jobvt = "A";    // 'A'll n rows of VT are returned in array VT

	double *a = new double [m*n];
	double *s= new double [std::min(m,n)];
	double *u= new double [m*m];
	double *vt =new double [n*n];

	// Start timer
	double start_time, end_time, svd_time; 
	CTKTimer timer;
	start_time = timer.getTimeNow();
	
	// Copy Jac from row-major to column-major format in a[]
	for ( int i=0; i<m; ++i )
		for ( int j=0; j<n; ++j )
			a[i+j*m] = gsl_matrix_get (jac, i, j);

        // Initialize MKL workspace
	int info = 0;
	double wkopt = 0;
	int lwrk = -1;
	dgesvd( jobu, jobvt, &m, &n, a, &m, s, u, &m, vt, &n, &wkopt, &lwrk, &info );
	lwrk = (int)wkopt;
	double *wrk = new double [lwrk];
	
	// Do SVD
	dgesvd( jobu, jobvt, &m, &n, a, &m, s, u, &m, vt, &n, wrk, &lwrk, &info);
	
	// End timer
	end_time = timer.getTimeNow();
	svd_time = end_time - start_time;
	log("dimitar") << "InfoSVD:\t svd_time   = " << svd_time   << " seconds."<< endl;
	
	// Populate s_out with the s singular values
    log("dimitar") << "InfoSVD)\t s_out Vector (m=" << s_out->size << ")" << endl;
	//for (int i=0;i<m;++i)
	for (int i=0;i<s_out->size;++i)
		gsl_vector_set ( s_out, i, s[i] ); 

	if (s_out->size<std::min(m,n)) {
        log("dimitar") << "InfoSVD)\t row_dim (m) > col_dim (n) of Jacobian." << endl;
        log("dimitar") << "InfoSVD)\t Print remaining singular values: " << endl;
		for (int i=s_out->size;i<std::min(m,n);++i)
            log("dimitar") << "\t s[" << i << "]=" << s[i] << " ";
        log("dimitar") << endl;
	}

//	// Populate ut_out with the v matrix values (vt column-major=v row-major)
//    log("dimitar") << "InfoSVD)\t ut_out Matrix (rows=" << ut_out->size1 << " columns=" << ut_out->size2 << ")" << endl;
//	for ( int i=0; i<m; ++i )
//		for ( int j=0; j<m; ++j ){//FULL
//			gsl_matrix_set (ut_out, i, j, u[i*m+j] );
//        }

	// Populate v_out with the v matrix values (vt column-major=v row-major)
    log("dimitar") << "InfoSVD)\t v_out Matrix (rows=" << v_out->size1 << " columns=" << v_out->size2 << ")" << endl;
	for ( int i=0; i<n; ++i )
		//for ( int j=0; j<std::min(m,n); ++j )//THIN
		for ( int j=0; j<n; ++j )//FULL
			gsl_matrix_set (v_out, i, j, vt[i*n+j] );
			
	// Free memory
	delete[] a;
	delete[] s;
	delete[] u;
	delete[] vt;
	delete[] wrk;

	my_mkl_free_buffers();

}

#elif __GPU_CULA
//---------------------------------------------------------
// GPU CULA SVD implementation
//---------------------------------------------------------
void printCulaError(culaStatus status) {
    char buf[256];

    if(!status) return;

    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    log("dimitar") << endl << buf << endl;

    exit(1);
}

int meetsMinimumCulaRequirements() {
    int cudaMinimumVersion = culaGetCudaMinimumVersion();
    int cudaRuntimeVersion = culaGetCudaRuntimeVersion();
    int cudaDriverVersion = culaGetCudaDriverVersion();
    int cublasMinimumVersion = culaGetCublasMinimumVersion();
    int cublasRuntimeVersion = culaGetCublasRuntimeVersion();

    if(cudaRuntimeVersion < cudaMinimumVersion) {
        log("dimitar")<<"InfoCula)\t CUDA runtime version is insufficient; version "<<cudaMinimumVersion<<" or greater is required."<<endl;
        return 0;
    }

    if(cudaDriverVersion < cudaMinimumVersion) {
        log("dimitar")<<"InfoCula)\t CUDA driver version is insufficient; version "<<cudaMinimumVersion<<" or greater is required."<<endl;
        return 0;
    }

    if(cublasRuntimeVersion < cublasMinimumVersion) {
        log("dimitar")<<"InfoCula)\t CUDA runtime version is insufficient; version "<<cublasMinimumVersion<<" or greater is required."<<endl;
        return 0;
    }

    return 1;
}

void nr_svd (gsl_matrix* jac, gsl_vector* s_out, gsl_matrix* v_out) {
	const int m = jac->size1; 
	const int n = jac->size2; 

        const char jobu  = 'N';    // 'N'o columns of U are computed and returned in array U
        const char jobvt = 'A';    // 'A'll n rows of VT are returned in array VT

        const int lda = m;         // "ld"=leading/physical dimension of array in mem
        const int ldu = m;        
        const int ldvt = n;       

	double *s = new double [std::min(m,n)];   // singular values of A, S(i) >= S(i+1), dimension=min(m,n)
	double *u = new double [m*m];
	double *vt = new double [n*n];

	// Start timer
        double start_time, end_time, svd_time;
	CTKTimer timer;
	start_time = timer.getTimeNow();

	// Jac[row-major] = Jac_t[column-major] format
	gsl_matrix* jac_t = gsl_matrix_alloc(n,m);
	gsl_matrix_transpose_memcpy(jac_t,jac);
	// Assume lda=jac_t->tda=m (which should be)!
        double *a = jac_t->data; 

        // Initialize CULA
        culaStatus status;
	if(meetsMinimumCulaRequirements()) {
        	status = culaInitialize();
        log("dimitar") << "InfoCula)\t Initializing CULA with status... " << culaGetStatusString(status) << endl;
		if(status != culaNoError) {
            log("dimitar") << "InfoCula)\t CULA Initialization error. Exit. " << endl;
			exit(1);
		}
	} else {
        log("dimitar") << "InfoCula)\t CULA Minimum requirements not met. Exit. " << endl;
		exit(1);
	}

        if(!a || !s || !u || !vt) {
        log("dimitar") << "InfoCula)\t Host side allocation error." << endl;
            	status = culaInsufficientMemory;
		exit(1);
        }
	
        // Do SVD
        status = culaGesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt);

	// End timer
        end_time = timer.getTimeNow();
        svd_time = end_time - start_time;
        log("dimitar") << "InfoSVD:\t svd_time   = " << svd_time   << " seconds."<< endl;

        // Check CULA errors
        if(status != culaNoError) printCulaError(status);

	// Populate s_out with the S singular values
    log("dimitar") << "InfoSVD)\t s_out Vector (m=" << s_out->size << ")" << endl;
	//for (int i=0;i<m;++i)
	for (int i=0;i<s_out->size;++i)
		gsl_vector_set ( s_out, i, s[i] ); 
	if (s_out->size<std::min(m,n)) {
        log("dimitar") << "InfoSVD)\t row_dim (m) > col_dim (n) of Jacobian." << endl;
        log("dimitar") << "InfoSVD)\t Print remaining singular values: " << endl;
		for (int i=s_out->size;i<std::min(m,n);++i)
            log("dimitar") << "\t s[" << i << "]=" << s[i] << " ";
        log("dimitar") << endl;
	}
	// Populate v_out with the V matrix values (vt column-major=v row-major)
    log("dimitar") << "InfoSVD)\t v_out Matrix (rows=" << v_out->size1 << " columns=" << v_out->size2 << ")" << endl;
	for ( int i=0; i<n; ++i )
         	for ( int j=0; j<n; ++j )
        		gsl_matrix_set (v_out, i, j, vt[i*n+j] );

	// Free memory
	gsl_matrix_free(jac_t);
	delete[] s;
	delete[] u;
	delete[] vt;

        // Shutdown
        culaShutdown();
    log("dimitar") << "InfoCula)\t Shutting down CULA with status... " << culaGetStatusString(status) << endl;

}
#else
//---------------------------------------------------------
// GSL SVD implementation
//---------------------------------------------------------


void nr_svd(gsl_matrix* jac, gsl_vector* s_out, gsl_matrix* v_out){ // we assume m<=n always //No we don't @RFonseca
	const int m = jac->size1;
	const int n = jac->size2;

	double start_time, end_time, svd_time; 
	CTKTimer timer;
	start_time = timer.getTimeNow();

	//if (m<n) {
	//	// make a transpose of Jac, because gsl_linalg_SV_decomp requires row_dim>=col_dim
	//	gsl_matrix* jac_t = gsl_matrix_alloc(n,m);
	//	gsl_matrix_transpose_memcpy(jac_t,jac);
	//	// do svd
	//	gsl_vector *S = gsl_vector_alloc(m);
	//	gsl_matrix *V = gsl_matrix_alloc(m,m);
	//	gsl_vector *work = gsl_vector_alloc(m);
	//	gsl_linalg_SV_decomp(jac_t,V,S,work);
	//	// copy the singular values into s_out
	//	for (int i=0; i<m; ++i)
	//		gsl_vector_set(s_out,i,gsl_vector_get(S,i));

	//	// do QR decomposition
	//	gsl_vector *tau = gsl_vector_alloc(m);
	//	// gsl_matrix_transpose_memcpy(jac_t,jac); // set jac_t again since it is modified during SVD computation.
	//	// It turns out it is unnecessary to restore jac_t. That is tricky.
	//	gsl_linalg_QR_decomp(jac_t,tau);
	//	// do QR unpack
	//	gsl_matrix *R = gsl_matrix_alloc(n,m);
	//	gsl_linalg_QR_unpack(jac_t,tau,v_out,R);

	//	// free memory
	//	gsl_matrix_free(jac_t);
	//	gsl_matrix_free(V);
	//	gsl_vector_free(S);
	//	gsl_vector_free(work);
	//	gsl_matrix_free(R);
	//	gsl_vector_free(tau);
	//}
	//else {  // This seems to modify jac, which should not happen?! @Dim
	//	gsl_vector *work = gsl_vector_alloc(n);
	//	gsl_linalg_SV_decomp(jac,v_out,s_out,work);
	//	gsl_vector_free(work);
	//}

	const int o = std::min(m,n);
	gsl_matrix* U = gsl_matrix_alloc(m,o);
	gsl_vector* work = gsl_vector_alloc(o);

	if (m<n) {
		gsl_matrix_transpose_memcpy(v_out,jac);
		gsl_linalg_SV_decomp(v_out,U,s_out,work);

	} else {
		gsl_matrix_memcpy(U, jac);
		gsl_linalg_SV_decomp(U,v_out,s_out,work);
	}

	gsl_vector_free(work);
	gsl_matrix_free(U);


	end_time = timer.getTimeNow();
	svd_time = end_time - start_time;
	log("dimitar") << "InfoSVD:\t svd_time   = " << svd_time   << " seconds."<< endl;

    log("dimitar") << "InfoSVD)\t s_out Vector (m=" << s_out->size << ")" << endl;
    log("dimitar") << "InfoSVD)\t v_out Matrix (rows=" << v_out->size1 << " columns=" << v_out->size2 << ")" << endl;

}
#endif

/*
 * Given a Jacobian matrix, computes its nullspace.
 */
void ComputeNullSpace(gsl_matrix* Jac, NullSpaceRet* Ret){
  nr_svd(Jac,NullSpaceRet::singularValues,NullSpaceRet::V);
  double maxSingularValue = gsl_vector_get(NullSpaceRet::singularValues,0);

  //Case with m < n and all singular values non-zero
  Ret->nullspaceSize = max((int)(NullSpaceRet::V->size2 - Ret->m),0);

  for (int i=0; i<NullSpaceRet::singularValues->size ; ++i){
    if ( gsl_vector_get(NullSpaceRet::singularValues,i) / maxSingularValue < SINGVAL_TOL ) {
      Ret->nullspaceSize = NullSpaceRet::V->size2-i;
      break;
    }
  }

  Ret->NullSpacePostCompute();
}

//---------------------------------------------------------
/**
 * Find the update to the dTheta vector that moves the end-effectors along the e-vector. It is assumed that 
 * $ \Delta\theta = \sum_{i=1}^r \frac{\sigma_i}{\sigma_i^2+\lambda^2} v_iu_i^T e $
 * where v_i and u_i are the columns of U and V, \sigma_i is the diagonal of D and \lambda is a damping constant.
 * For details see the least squared method in http://math.ucsd.edu/~sbuss/ResearchWeb/ikmethods/iksurvey.pdf
 * (This method might not work because nullspace->V is not calculated properly)
 */
//void TorsionUpdate(NullSpaceRet* nullspace, gsl_vector *e, double lambda, gsl_vector *dTheta) {
//	gsl_matrix* D = gsl_matrix_calloc(nullspace->n, nullspace->n);
//	for(int i=0;i<nullspace->n;i++)
//		gsl_matrix_set(D, i,i, gsl_vector_get(nullspace->Sval, i));
//
//
//	gsl_matrix* prod1 = gsl_matrix_calloc(nullspace->n,nullspace->m);
//	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, D, nullspace->Svec, 0.0, prod1);
//	gsl_matrix* prod2 = gsl_matrix_calloc(nullspace->n,nullspace->n);
//	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, nullspace->V, prod1, 0.0, prod2);
//	gsl_blas_dgemv(CblasNoTrans, 1.0, prod2, e, 1.0, dTheta);
//
//}



