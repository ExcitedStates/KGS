#ifndef KGS_GSL_HELPERS_H
#define KGS_GSL_HELPERS_H

#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/** Print matrix to cout */
void gsl_matrix_cout (      const gsl_matrix *m);

/** Print matrix to file specified by `filename` */
void gsl_matrix_outtofile ( const gsl_matrix *m, const std::string& filename    );

/** Print vector to file specified by `filename` */
void gsl_vector_outtofile ( const gsl_vector *v, const std::string& filename    );

/** Print vector to cout */
void gsl_vector_cout (      const gsl_vector *v );

/** Return length of vector */
double gsl_vector_length(   const gsl_vector *v );

/** Scale vector to length 1.0 */
void gsl_vector_normalize ( gsl_vector *v );

/** Scale vector to specified length */
void gsl_vector_scale_to_length(gsl_vector* ret, double length);

/** Scale vector so largest absolute value of any component is at most `maxComponent`.
 * If the values of all components are less than `maxComponent` no scaling is performed. */
void gsl_vector_scale_max_component(gsl_vector* v, double maxComponent);


gsl_matrix* gsl_matrix_trans(gsl_matrix* A);
gsl_matrix* gsl_matrix_mul(gsl_matrix* A, gsl_matrix* B);
gsl_vector* gsl_matrix_vector_mul(gsl_matrix* A, gsl_vector* v);
gsl_matrix* pca (gsl_matrix* sample_matrix);
gsl_vector* RandomUnitVector (int size);
double frobenius_norm (const gsl_matrix *m);


#endif //KGS_GSL_HELPERS_H
