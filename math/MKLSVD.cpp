#include "MKLSVD.h"

#include <algorithm>

#ifdef __INTEL_MKL
#include <mkl_lapack.h>

void MKLSVD::UpdateFromMatrix()
{
    int lwrk;
    const int mn = m*n;
    const int mm = m*m;
    const int nn = n*n;
    double *a = new double [mn];for(int i=0;i<mn;i++) a[i]=0;
    double *s= new double [m];  for(int i=0;i<m ;i++) s[i]=0;
    double *u= new double [mm]; for(int i=0;i<mm;i++) u[i]=0;
    double *vt =new double [nn];for(int i=0;i<nn;i++) vt[i]=0;
    const char* Ustr = "A"; //All info: Full SVD
    const char* Vstr = "A"; //All info: Full SVD

    // Copy Jac from row-major to column-major format in a[]
    for ( int i=0; i<m; ++i )
        for ( int j=0; j<n; ++j )
            a[i+j*m] = gsl_matrix_get (matrix, i, j );

    int info = 0;
    double wkopt = 0;
    lwrk = -1;
    dgesvd( Ustr, Vstr, &m, &n, a, &m, s, u, &m, vt, &n, &wkopt, &lwrk, &info );
    lwrk = (int)wkopt;
    double *wrk = new double [lwrk]; for(int i=0;i<lwrk;i++) wrk[i]=0;

    dgesvd(Ustr, Vstr, &m, &n, a, &m, s, u, &m, vt, &n, wrk, &lwrk, &info);


    // Convert from MKL to GSL matrices
    //TODO: Do this with gsl views on the arrays .. way faster

    // Populate the S singular values
    for (int i=0;i<std::min(m,n);++i)
        gsl_vector_set ( S, i, s[i] );

    // Populate the V matrix values (vt column-major=v row-major)
    for ( int i=0; i<n; ++i )
        for ( int j=0; j<n; ++j )
            gsl_matrix_set (V, i, j, vt[i*n+j] );

    // Populate the U matrix values (vt column-major=v row-major)
    for ( int i=0; i<m; ++i )
        for ( int j=0; j<m; ++j )
            gsl_matrix_set (U, j, i, u[i*m+j] );

    // Free memory
    delete[] a;
    delete[] s;
    delete[] u;
    delete[] vt;
    delete[] wrk;

}

#else

void MKLSVD::UpdateFromMatrix()
{
    throw "MKLSVD::update_from_matrix error! MKL not supported. Install MKL and compile with -D__INTEL_MKL";
}

#endif

