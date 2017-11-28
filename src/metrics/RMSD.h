/*

Excited States software: KGS
Contributors: See CONTRIBUTORS.txt
Contact: kgs-contact@simtk.org

Copyright (C) 2009-2017 Stanford University

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

 private:
	double kabsch(
			gsl_matrix *X,     /* the points to be moved */
			gsl_matrix *Y,     /* the points to move to */
			gsl_matrix *U,     /* the rotation matrix */
			gsl_vector *t      /* the translation vector */
	) ; // /* returns the rmsd between the optimized vector sets */

  double alignedrmsd(gsl_matrix *mat1, gsl_matrix *mat2, int N);
};

typedef struct
{
  double m[4][4];
} MATRIX;

#define vdiff2(a,b) ( ((a)[0]-(b)[0]) * ((a)[0]-(b)[0]) +	\
		((a)[1]-(b)[1]) * ((a)[1]-(b)[1]) + \
		((a)[2]-(b)[2]) * ((a)[2]-(b)[2]) )

double alignedrmsd(double *v1, double *v2, int N);
//void centroid(double *ret, double *v, int N);
//int getalignmtx(double *v1, double *v2, int N, MATRIX *mtx);
//void crossproduct(double *ans, double *pt1, double *pt2);
//void mtx_root(MATRIX *mtx);
//int almostequal(MATRIX *a, MATRIX *b);
//void mulpt(MATRIX *mtx, double *pt);
//void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y);
//void mtx_identity(MATRIX *mtx);
//bool isRightHanded(MATRIX *mtx);
//void mtx_trans(MATRIX *mtx, double x, double y, double z);
//int mtx_invert(double *mtx, int N);
//double absmaxv(double *v, int N);

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
//double rmsd(double *v1, double *v2, int N, double *mtx);


}
#endif
