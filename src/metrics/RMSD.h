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
  RMSD(Selection& selection);

  double distance(Configuration*, Configuration*);
  double distance_noOptimization(Configuration *c1, Configuration *c2);

  double align(Molecule * other, Molecule * base);

};

typedef struct
{
  double m[4][4];
} MATRIX;

#define vdiff2(a,b) ( ((a)[0]-(b)[0]) * ((a)[0]-(b)[0]) +	\
		((a)[1]-(b)[1]) * ((a)[1]-(b)[1]) + \
		((a)[2]-(b)[2]) * ((a)[2]-(b)[2]) )

double alignedrmsd(double *v1, double *v2, int N);
void centroid(double *ret, double *v, int N);
int getalignmtx(double *v1, double *v2, int N, MATRIX *mtx);
void crossproduct(double *ans, double *pt1, double *pt2);
void mtx_root(MATRIX *mtx);
int almostequal(MATRIX *a, MATRIX *b);
void mulpt(MATRIX *mtx, double *pt);
void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y);
void mtx_identity(MATRIX *mtx);
bool isRightHanded(MATRIX *mtx);
void mtx_trans(MATRIX *mtx, double x, double y, double z);
int mtx_invert(double *mtx, int N);
double absmaxv(double *v, int N);

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
double rmsd(double *v1, double *v2, int N, double *mtx);


}
#endif
