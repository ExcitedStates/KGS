#ifndef RMSD_H_
#define RMSD_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <list>

#include "Metric.h"
#include "core/Configuration.h"
#include "core/Molecule.h"

namespace metrics{


class RMSD: public Metric{
 public:
  RMSD();
  RMSD(std::vector<Atom*>* atomsRMSD);

  double distance(Configuration*, Configuration*);
  static double distance_noOptimization(Configuration *c1, Configuration *c2);

  static double align(Molecule * other, Molecule * base);

 private:
  std::vector<Atom*>* m_atomsRMSD;

};

typedef struct
{
  float m[4][4];
} MATRIX;

#define vdiff2(a,b) ( ((a)[0]-(b)[0]) * ((a)[0]-(b)[0]) +	\
		((a)[1]-(b)[1]) * ((a)[1]-(b)[1]) + \
		((a)[2]-(b)[2]) * ((a)[2]-(b)[2]) )

double alignedrmsd(float *v1, float *v2, int N);
void centroid(float *ret, float *v, int N);
int getalignmtx(float *v1, float *v2, int N, MATRIX *mtx);
void crossproduct(float *ans, float *pt1, float *pt2);
void mtx_root(MATRIX *mtx);
int almostequal(MATRIX *a, MATRIX *b);
void mulpt(MATRIX *mtx, float *pt);
void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y);
void mtx_identity(MATRIX *mtx);
void mtx_trans(MATRIX *mtx, float x, float y, float z);
int mtx_invert(float *mtx, int N);
float absmaxv(float *v, int N);

/*
   calculate rmsd between two structures
Params: v1 - first set of points
v2 - second set of points
N - number of points
mtx - return for transfrom matrix used to align structures
Returns: rmsd score
Notes: mtx can be null. Transform will be rigid. Inputs must
be previously aligned for sequence alignment
*/
double rmsd(float *v1, float *v2, int N, float *mtx);


}
#endif
