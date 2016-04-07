
#include "gsl_helpers.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_blas.h>
#include <assert.h>

#include "Logger.h"

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
