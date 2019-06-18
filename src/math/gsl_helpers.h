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

/** Print vector to output-stream */
void gsl_vector_out (const gsl_vector *v, std::ostream& os);

/** Print vector to cout */
void gsl_vector_cout (const gsl_vector *v);

/** Return length of vector */
double gsl_vector_length(   const gsl_vector *v );

/** Scale vector to length 1.0 */
void gsl_vector_normalize ( gsl_vector *v );

/** Scale vector to specified length */
void gsl_vector_scale_to_length(gsl_vector* ret, double length);

/** Scale vector so largest absolute value of any component is at most `maxComponent`.
 * If the values of all components are less than `maxComponent` no scaling is performed. */
void gsl_vector_scale_max_component(gsl_vector* v, double maxComponent);

/** Make a copy of the vector */
gsl_vector* gsl_vector_copy(gsl_vector*);

gsl_matrix* gsl_matrix_trans(gsl_matrix* A);

gsl_matrix* gsl_matrix_mul(gsl_matrix* A, gsl_matrix* B);

gsl_vector* gsl_matrix_vector_mul(gsl_matrix* A, gsl_vector* v);

gsl_matrix* pca (gsl_matrix* sample_matrix);

gsl_vector* RandomUnitVector (int size);

double frobenius_norm (const gsl_matrix *m);

/** Compute Shannon entropy (as in JCIM paper: fraction of significant contributors) for given input vector**/
double fractionOfSignificantContributors(const gsl_vector *v);

/** Compute Shannon entropy unnormalized: number of significant contributors of given input vector**/
double significantContributors(const gsl_vector *v);

/** Compute Shannon entropy from information theory in bits for given input vector**/
double shannonEntropyInBits(const gsl_vector* v);

/** Compute Shannon entropy from information theory in bits for given input vector**/
double shannonEntropyUnnormalizedInBits(const gsl_vector* v);


#endif //KGS_GSL_HELPERS_H
