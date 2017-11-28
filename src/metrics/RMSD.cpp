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

#include <cmath>
#include "metrics/RMSD.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include "math/gsl_helpers.h"

using namespace std;

namespace metrics {

RMSD::RMSD(Selection &selection) :
    Metric(selection) {}


double RMSD::distance(Configuration *c1, Configuration *c2) {

  const std::vector<Atom *> &atomsRMSD1 = m_selection.getSelectedAtoms(c1->getMolecule());

  Molecule* m1 = c1->getMolecule();
  Molecule* m2 = c2->getMolecule();

  // Build two vectors of atoms where each entry in one has matching chain name, residue id, and atom name in the other
  std::vector<Atom*> *matchedAtoms1;
  std::vector<Atom*> *matchedAtoms2;
  if (m1 == m2){
    matchedAtoms1 = new std::vector<Atom*>(atomsRMSD1);
    matchedAtoms2 = new std::vector<Atom*>(atomsRMSD1);
  }else {
    matchedAtoms1 = new std::vector<Atom*>();
    matchedAtoms2 = new std::vector<Atom*>();

    for (auto const &aIt : atomsRMSD1) {
      const std::string &name = aIt->getName();
      const std::string &chainName = aIt->getResidue()->getChain()->getName();
      int resId = aIt->getResidue()->getId();
      Atom *a1 = &(*aIt);
      Atom *a2 = m2->getAtom(chainName, resId, name);
      if (a2 != nullptr) {
        matchedAtoms1->push_back(a1);
        matchedAtoms2->push_back(a2);
      }
    }
  }

  if (atomsRMSD1.empty()) {
    cerr << "RMSD::distance - Atom-selection given to RMSD metric contained no atoms: " << m_selection << endl;
    exit(-1);
  }

  int atom_num = matchedAtoms1->size();
  assert(atom_num > 3);

  gsl_matrix *static_matrix = gsl_matrix_alloc(atom_num, 3);
  gsl_matrix *moving_matrix = gsl_matrix_alloc(atom_num, 3);

  c1->updateMolecule();
  unsigned int i = 0;
  for (auto const &aIt : *matchedAtoms1) {
    const Coordinate &c1 = (*aIt).m_position;
    gsl_matrix_set(static_matrix, i, 0, c1.x);
    gsl_matrix_set(static_matrix, i, 1, c1.y);
    gsl_matrix_set(static_matrix, i, 2, c1.z);

    assert(!std::isnan(c1.x));
    assert(!std::isnan(c1.y));
    assert(!std::isnan(c1.z));
    assert(c1.x < 10000.0);
    assert(c1.y < 10000.0);
    assert(c1.z < 10000.0);

    i++;
  }

  c2->updateMolecule();
  i = 0;
  for (auto const &aIt: *matchedAtoms2) {
    const Coordinate &c2 = (*aIt).m_position;
    gsl_matrix_set(moving_matrix, i, 0, c2.x);
    gsl_matrix_set(moving_matrix, i, 1, c2.y);
    gsl_matrix_set(moving_matrix, i, 2, c2.z);

    assert(!std::isnan(c2.x));
    assert(!std::isnan(c2.y));
    assert(!std::isnan(c2.z));
    assert(c2.x < 10000.0);
    assert(c2.y < 10000.0);
    assert(c2.z < 10000.0);

    i++;
  }

  gsl_matrix *rotate = gsl_matrix_alloc(3, 3);
  gsl_vector *transl = gsl_vector_alloc(3);

  delete matchedAtoms1;
  delete matchedAtoms2;

  /// Compute the optimal rotation and translation, returns remaining rmsd
  double rmsd = kabsch(moving_matrix, static_matrix, rotate, transl);

  gsl_matrix_free(moving_matrix);
  gsl_matrix_free(static_matrix);
  gsl_matrix_free(rotate);
  gsl_vector_free(transl);

  return rmsd;
}

double RMSD::distance_noOptimization(Configuration *c1, Configuration *c2) {

  vector<Atom *> &atomsRMSD1 = m_selection.getSelectedAtoms(c1->getMolecule());

  if (atomsRMSD1.empty()) {
    cout << "RMSD::distance_noOptimization - Found no atoms to align .. using all" << endl;
    exit(-1);
  }

  vector<Coordinate> p1_atoms;
  vector<Coordinate> p2_atoms;

  int resId;
  std::string name;
  std::string chainName;

  Molecule *p1 = c1->updatedMolecule();
  for (auto const& aIt : atomsRMSD1) {
    name = aIt->getName();
    chainName = aIt->getResidue()->getChain()->getName();
    resId = aIt->getResidue()->getId();
    Atom *a2 = c2->getMolecule()->getAtom(chainName, resId, name);
    if (a2 == nullptr) continue;
    p1_atoms.push_back(aIt->m_position);
  }

  Molecule *p2 = c2->updatedMolecule();
  for (auto const &aIt : atomsRMSD1) {
    name = aIt->getName();
    chainName = aIt->getResidue()->getChain()->getName();
    resId = aIt->getResidue()->getId();
    Atom *a2 = c2->getMolecule()->getAtom(chainName, resId, name);
    if (a2 == nullptr) continue;
    p2_atoms.push_back(a2->m_position);
  }

  double sum = 0.0;
  for (int i = 0; i != atomsRMSD1.size(); i++) {
    sum += p1_atoms[i].distanceSquared(p2_atoms[i]);
  }

  return sqrt(sum / atomsRMSD1.size());
}

double RMSD::align(Molecule *moving_mol, Molecule *static_mol) {
  const std::vector<Atom *> &atomsRMSD1 = m_selection.getSelectedAtoms(static_mol);

  Molecule* m1 = static_mol;
  Molecule* m2 = moving_mol;

  // Build two vectors of atoms where each entry in one has matching chain name, residue id, and atom name in the other
  std::vector<Atom*> *matchedAtoms1;
  std::vector<Atom*> *matchedAtoms2;
  if (m1 == m2){
    matchedAtoms1 = new std::vector<Atom*>(atomsRMSD1);
    matchedAtoms2 = new std::vector<Atom*>(atomsRMSD1);
  }else {
    matchedAtoms1 = new std::vector<Atom*>();
    matchedAtoms2 = new std::vector<Atom*>();

    for (auto const &aIt : atomsRMSD1) {
      const std::string &name = aIt->getName();
      const std::string &chainName = aIt->getResidue()->getChain()->getName();
      int resId = aIt->getResidue()->getId();
      Atom *a1 = &(*aIt);
      Atom *a2 = m2->getAtom(chainName, resId, name);
      if (a2 != nullptr) {
        matchedAtoms1->push_back(a1);
        matchedAtoms2->push_back(a2);
      }
    }
  }

  if (atomsRMSD1.empty()) {
    cerr << "RMSD::distance - Atom-selection given to RMSD metric contained no atoms: " << m_selection << endl;
    exit(-1);
  }

  int atom_num = matchedAtoms1->size();
  assert(atom_num > 3);

  gsl_matrix *static_matrix = gsl_matrix_alloc(atom_num, 3);
  gsl_matrix *moving_matrix = gsl_matrix_alloc(atom_num, 3);

  unsigned int i = 0;
  for (auto const &aIt : *matchedAtoms1) {
    const Coordinate &c1 = (*aIt).m_position;
    gsl_matrix_set(static_matrix, i, 0, c1.x);
    gsl_matrix_set(static_matrix, i, 1, c1.y);
    gsl_matrix_set(static_matrix, i, 2, c1.z);

    assert(!std::isnan(c1.x));
    assert(!std::isnan(c1.y));
    assert(!std::isnan(c1.z));
    assert(c1.x < 10000.0);
    assert(c1.y < 10000.0);
    assert(c1.z < 10000.0);

    i++;
  }

  i = 0;
  for (auto const &aIt: *matchedAtoms2) {
    const Coordinate &c2 = (*aIt).m_position;
    gsl_matrix_set(moving_matrix, i, 0, c2.x);
    gsl_matrix_set(moving_matrix, i, 1, c2.y);
    gsl_matrix_set(moving_matrix, i, 2, c2.z);

    assert(!std::isnan(c2.x));
    assert(!std::isnan(c2.y));
    assert(!std::isnan(c2.z));
    assert(c2.x < 10000.0);
    assert(c2.y < 10000.0);
    assert(c2.z < 10000.0);

    i++;
  }

  gsl_matrix *rotate = gsl_matrix_alloc(3, 3);
  gsl_vector *transl = gsl_vector_alloc(3);

  /// Compute the optimal rotation and translation, returns remaining rmsd
  double rmsd = kabsch(moving_matrix, static_matrix, rotate, transl);

  // Transform position in `moving_mol`
  gsl_matrix *moving_coords = gsl_matrix_alloc(3, moving_mol->getAtoms().size());
  gsl_matrix *moving_coords2 = gsl_matrix_alloc(3, moving_mol->getAtoms().size());
  gsl_matrix_set_zero(moving_coords);
  gsl_matrix_set_zero(moving_coords2);

  i = 0;
  /// First copy the positions
  for (auto const &aIt : moving_mol->getAtoms()) {
    Coordinate &pos = aIt->m_position;
    gsl_matrix_set(moving_coords, 0, i, pos.x);
    gsl_matrix_set(moving_coords, 1, i, pos.y);
    gsl_matrix_set(moving_coords, 2, i, pos.z);
    i++;
  }
  /// Do the rotation
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rotate, moving_coords, 0.0, moving_coords2);

  double dx = gsl_vector_get(transl, 0);
  double dy = gsl_vector_get(transl, 1);
  double dz = gsl_vector_get(transl, 2);

  ///Update molecule positions
  //Overall Px,new = R*(Px,old - cx)+cy = R*Px,old - R*cx + cy = R*Px,old + t
  i = 0;
  for (auto const &aIt : moving_mol->getAtoms()) {
    Coordinate &pos = aIt->m_position;
    Coordinate &refPos = aIt->m_referencePosition;
    pos.x = gsl_matrix_get(moving_coords2, 0, i) + dx;
    pos.y = gsl_matrix_get(moving_coords2, 1, i) + dy;
    pos.z = gsl_matrix_get(moving_coords2, 2, i) + dz;
    i++;
  }

  gsl_matrix_free(static_matrix);
  gsl_matrix_free(moving_matrix);
  gsl_matrix_free(rotate);
  gsl_vector_free(transl);
  gsl_matrix_free(moving_coords);
  gsl_matrix_free(moving_coords2);

  return rmsd;
}

/* gsl does not provide it */
static inline void gsl_vector_cross(
    const gsl_vector *a,
    const gsl_vector *b,
    gsl_vector *c
) {
  double a0 = gsl_vector_get(a, 0);
  double a1 = gsl_vector_get(a, 1);
  double a2 = gsl_vector_get(a, 2);
  double b0 = gsl_vector_get(b, 0);
  double b1 = gsl_vector_get(b, 1);
  double b2 = gsl_vector_get(b, 2);
  gsl_vector_set(c, 0, a1 * b2 - b1 * a2);
  gsl_vector_set(c, 1, a2 * b0 - b2 * a0);
  gsl_vector_set(c, 2, a0 * b1 - b0 * a1);
}

#define NORM_EPS 0.00000001


double RMSD::kabsch(
    gsl_matrix *X,     /* the points to be moved */
    gsl_matrix *Y,     /* the points to move to */
    gsl_matrix *U,     /* the rotation matrix */
    gsl_vector *t     /* the translation vector */
) {
  unsigned int i, j, k;
  int U_ok = 1;
  unsigned int size = X->size1;
  double n = 1.0 / size;
  gsl_vector *cx = gsl_vector_alloc(3);     /* centroid of X */
  gsl_vector *cy = gsl_vector_alloc(3);     /* centroid of Y */
  gsl_matrix *R = gsl_matrix_alloc(3, 3);    /* Kabsch's R */
  gsl_matrix *RTR = gsl_matrix_alloc(3, 3);  /* R_trans * R (and Kabsch's bk) */
  gsl_eigen_symmv_workspace *espace = gsl_eigen_symmv_alloc(3);
  gsl_matrix *evec = gsl_matrix_alloc(3, 3); /* eigenvectors (and Kabsch's ak) */
  gsl_vector *eval = gsl_vector_alloc(3);   /* vector of eigenvalues */
  gsl_matrix *Xrot = gsl_matrix_alloc(size,3);     /* the optimized points of X, overlaid on moved X */

  /* compute centroid of X */
  gsl_vector_set_zero(cx);
//  gsl_vector *row = gsl_vector_alloc(3);
  for (i = size; i > 0;) {
    gsl_vector_const_view row = gsl_matrix_const_row(X, --i);
    gsl_vector_add(cx, &row.vector);
//    gsl_matrix_get_row(row,X, --i);
//    gsl_vector_scale(row,n);
//    gsl_vector_add(cx, row); //scaling each entry for numerical reasons
  }
  gsl_vector_scale(cx, n);

  /* compute centroid of Y */
  gsl_vector_set_zero(cy);
  for (i = size; i > 0;) {
    gsl_vector_const_view row = gsl_matrix_const_row(Y, --i);
    gsl_vector_add(cy, &row.vector);//scaling each entry for numerical reasons
//    gsl_matrix_get_row(row,Y, --i);
//    gsl_vector_scale(row,n);
//    gsl_vector_add(cy, row); //scaling each entry for numerical reasons
  }
  gsl_vector_scale(cy, n);

  /* move X to origin */
  for (i = size; i > 0;) {
    gsl_vector_view row = gsl_matrix_row(X, --i);
    gsl_vector_sub(&row.vector, cx);
  }
  /* move Y to origin */
  for (i = size; i > 0;) {
    gsl_vector_view row = gsl_matrix_row(Y, --i);
    gsl_vector_sub(&row.vector, cy);
  }

  if (size == 1) {
    /* just one point, so U is trival */
    gsl_matrix_set_identity(U);
  } else {
    /* compute R */
    gsl_matrix_set_zero(R);
    for (k = size; k > 0;) {
      --k;
      for (i = 3; i > 0;) {
        --i;
        for (j = 3; j > 0;) {
          --j;
          gsl_matrix_set(R, i, j,
                         gsl_matrix_get(R, i, j) +
                         gsl_matrix_get(Y, k, i) * gsl_matrix_get(X, k, j)
          );
        }
      }
    }

    /* compute RTR = R_trans * R */
    gsl_matrix_set_zero(RTR);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, R, R, 0.0, RTR);

    /* compute orthonormal eigenvectors */
    gsl_eigen_symmv(RTR, eval, evec, espace);  /* RTR will be modified! */
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
    if (gsl_vector_get(eval, 1) > NORM_EPS) {
      /* compute ak's (as columns of evec) and bk's (as columns of RTR) */
      double norm_b0, norm_b1, norm_b2;
      gsl_vector_const_view a0 = gsl_matrix_const_column(evec, 0);
      gsl_vector_const_view a1 = gsl_matrix_const_column(evec, 1);
      gsl_vector_view a2 = gsl_matrix_column(evec, 2);
      gsl_vector_view b0 = gsl_matrix_column(RTR, 0);
      gsl_vector_view b1 = gsl_matrix_column(RTR, 1);
      gsl_vector_view b2 = gsl_matrix_column(RTR, 2);
      gsl_vector_cross(&a0.vector, &a1.vector, &a2.vector); /* a2 = a0 x a1 */
      gsl_blas_dgemv(CblasNoTrans, 1.0, R, &a0.vector, 0.0, &b0.vector);
      norm_b0 = gsl_blas_dnrm2(&b0.vector);
      gsl_blas_dgemv(CblasNoTrans, 1.0, R, &a1.vector, 0.0, &b1.vector);
      norm_b1 = gsl_blas_dnrm2(&b1.vector);
      if (norm_b0 > NORM_EPS && norm_b1 > NORM_EPS) {
        gsl_vector_scale(&b0.vector, 1.0 / norm_b0);         /* b0 = ||R * a0|| */
        gsl_vector_scale(&b1.vector, 1.0 / norm_b1);         /* b1 = ||R * a1|| */
        gsl_vector_cross(&b0.vector, &b1.vector, &b2.vector);  /* b2 = b0 x b1 */

        norm_b2 = gsl_blas_dnrm2(&b2.vector);
        if (norm_b2 > NORM_EPS) {
          /* we reach this point only if all bk different from 0 */
          /* compute U = B * A_trans (use RTR as B and evec as A) */
          gsl_matrix_set_zero(U); /* to avoid nan */
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, RTR, evec, 0.0, U);
        } else {
          U_ok = 0;
          gsl_matrix_set_identity(U);
        }
      } else {
        U_ok = 0;
        gsl_matrix_set_identity(U);
      }
    } else {
      U_ok = 0;
      gsl_matrix_set_identity(U);
    }
  }

  /* compute t = cy - U * cx  */
  gsl_vector_memcpy(t, cy);
  gsl_blas_dgemv(CblasNoTrans, -1.0, U, cx, 1.0, t);

  /* compute the rmsd (cheapest at this point) */
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,X,U,0.0,Xrot);
  double rmsd = alignedrmsd(Y,Xrot,size);

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  gsl_eigen_symmv_free(espace);
  gsl_matrix_free(RTR);
  gsl_matrix_free(R);
  gsl_vector_free(cy);
  gsl_vector_free(cx);
  gsl_matrix_free(Xrot);

  if(U_ok == 1)
    return rmsd;
  else
    return -1.0;
}

/*
   calculate rmsd between two aligned structures (trivial)
Params: mat1 - first structure
mat2 - second structure
N - number of points
Returns: rmsd
*/
double RMSD::alignedrmsd(gsl_matrix *mat1, gsl_matrix *mat2, int N)
{
  double answer =0.0;
  int i,j;
  for(i=0;i<N;i++) {
      answer += (gsl_matrix_get(mat1,i,0)-gsl_matrix_get(mat2,i,0))*(gsl_matrix_get(mat1,i,0)-gsl_matrix_get(mat2,i,0))
                + (gsl_matrix_get(mat1,i,1)-gsl_matrix_get(mat2,i,1))*(gsl_matrix_get(mat1,i,1)-gsl_matrix_get(mat2,i,1))
                + (gsl_matrix_get(mat1,i,2)-gsl_matrix_get(mat2,i,2))*(gsl_matrix_get(mat1,i,2)-gsl_matrix_get(mat2,i,2));
  }
  return sqrt(answer/N);
}


//double rmsd(double *v1, double *v2, int N, double *mtx) {
//  double cent1[3];
//  double cent2[3];
//  MATRIX tmtx;
//  MATRIX tempmtx;
//  MATRIX move1;
//  MATRIX move2;
//  int i;
//  double answer;
//  int err;
//
//  assert(N > 3);
//
//  double* temp1 = (double*) malloc(N * 3 * sizeof(double));
//  double* temp2 = (double*) malloc(N * 3 * sizeof(double));
//  if(!temp1 || !temp2)
//    goto error_exit;
//
//  centroid(cent1, v1, N);
//  centroid(cent2, v2, N);
//  for(i=0;i<N;i++)
//  {
//    temp1[i*3+0] = v1[i*3+0] - cent1[0];
//    temp1[i*3+1] = v1[i*3+1] - cent1[1];
//    temp1[i*3+2] = v1[i*3+2] - cent1[2];
//
//    temp2[i*3+0] = v2[i*3+0] - cent2[0];
//    temp2[i*3+1] = v2[i*3+1] - cent2[1];
//    temp2[i*3+2] = v2[i*3+2] - cent2[2];
//  }
//
//  err = getalignmtx(temp1, temp2, N, &tmtx);
//  if(err == -1)
//    goto error_exit;
//
//  mtx_trans(&move1, -cent2[0], -cent2[1], -cent2[2]);
//  mtx_mul(&tempmtx, &move1, &tmtx);
//  mtx_trans(&move2, cent1[0], cent1[1], cent1[2]);
//  mtx_mul(&tmtx, &tempmtx, &move2);
//  memcpy(temp2, v2, N * sizeof(double) * 3);
//  for(i=0;i<N;i++)
//    mulpt(&tmtx, temp2 + i * 3);
//  answer = alignedrmsd(v1, temp2, N);
//  free(temp1);
//  free(temp2);
//  if(mtx)
//    memcpy(mtx, &tmtx.m, 16 * sizeof(double));
//
//  return answer;
//  error_exit:
//  free(temp1);
//  free(temp2);
//  if(mtx)
//  {
//    for(i=0;i<16;i++)
//      mtx[i] = 0;
//  }
//  return sqrt(-1.0);
//}

/*
   calculate rmsd between two aligned structures (trivial)
Params: v1 - first structure
v2 - second structure
N - number of points
Returns: rmsd
*/
//double alignedrmsd(double *v1, double *v2, int N)
//{
//  double answer =0;
//  int i;
//
//  for(i=0;i<N;i++)
//    answer += vdiff2(v1 + i *3, v2 + i * 3);
//  return sqrt(answer/N);
//}

/* compute the centroid */
//void centroid(double *ret, double *v, int N)
//{
//  int i;
//
//  ret[0] = 0;
//  ret[1] = 0;
//  ret[2] = 0;
//  for(i=0;i<N;i++){
//    ret[0] += v[i*3+0]/N;
//    ret[1] += v[i*3+1]/N;
//    ret[2] += v[i*3+2]/N;
//    assert(ret[0]<10000);
//    assert(ret[1]<10000);
//    assert(ret[2]<10000);
//  }
//}

/*
   get the matrix needed to align two structures
Params: v1 - reference structure
v2 - structure to align
N - number of points
mtx - return for rigid body alignment matrix
Notes: only calculates rotation part of matrix.
assumes input has been aligned to centroids
*/
//int getalignmtx(double *v1, double *v2, int N, MATRIX *mtx)
//{
//  MATRIX A = { {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}} };
//  MATRIX At;
//  MATRIX Ainv;
//  MATRIX temp;
//  double tv[3];
//  double tw[3];
//  double tv2[3];
//  double tw2[3];
//  int k, i, j;
//  int flag = 0;
//  double correction;
//
//  correction = absmaxv(v1, N * 3) * absmaxv(v2, N * 3);
//
//  for(k=0;k<N;k++)
//    for(i=0;i<3;i++)
//      for(j=0;j<3;j++)
//        A.m[i][j] += (v1[k*3+i] * v2[k*3+j])/correction;
//
//  while(flag < 3)
//  {
//    for(i=0;i<4;i++)
//      for(j=0;j<4;j++)
//        At.m[i][j] = A.m[j][i];
//
//    memcpy(&Ainv, &A, sizeof(MATRIX));
//    /* this will happen if all points are in a plane */
//    if( mtx_invert((double *) &Ainv, 4) == -1)
//    {
//      if(flag == 0)
//      {
//        crossproduct(tv, v1, v1+3);
//        crossproduct(tw, v2, v2+3);
//      }
//      else
//      {
//        crossproduct(tv2, tv, v1);
//        crossproduct(tw2, tw, v2);
//        memcpy(tv, tv2, 3 * sizeof(double));
//        memcpy(tw, tw2, 3 * sizeof(double));
//      }
//      for(i=0;i<3;i++)
//        for(j=0;j<3;j++)
//          A.m[i][j] += tv[i] * tw[j];
//
//      flag++;
//    }
//    else
//      flag = 5;
//  }
//  if(flag != 5)
//    return -1;
//
//  mtx_mul(&temp, &At, &A);
//  mtx_root(&temp);
//  mtx_mul(mtx, &temp, &Ainv);
//  ////////NEW FOR CORRECTION
//  /// check right-handedness of rotation matrix
//  if (!isRightHanded(mtx)){
//    mtx->m[0][0] = - mtx->m[0][0];
//    mtx->m[1][0] = - mtx->m[1][0];
//    mtx->m[2][0] = - mtx->m[2][0];
//    cerr<<"Detected unallowed reflection of the molecule"<<endl;
////      exit(-1);
//  }
//  ////////DONE WITH CORRECTION
//  return 0;
//}

/*
   get the crossproduct of two vectors.
Params: ans - return pinter for answer.
pt1 - first vector
pt2 - second vector.
Notes: crossproduct is at right angles to the two vectors.
*/
//void crossproduct(double *ans, double *pt1, double *pt2)
//{
//  ans[0] = pt1[1] * pt2[2] - pt1[2] * pt2[1];
//  ans[1] = pt1[0] * pt2[2] - pt1[2] * pt2[0];
//  ans[2] = pt1[0] * pt2[1] - pt1[1] * pt2[0];
//}

/*
   Denman-Beavers square m_root iteration
   */
//void mtx_root(MATRIX *mtx)
//{
//  MATRIX Y = *mtx;
//  MATRIX Z;
//  MATRIX Y1;
//  MATRIX Z1;
//  MATRIX invY;
//  MATRIX invZ;
//  MATRIX Y2;
//  int iter = 0;
//  int i, ii;
//
//  mtx_identity(&Z);
//
//  do
//  {
//    invY = Y;
//    invZ = Z;
//    if( mtx_invert((double *) &invY, 4) == -1)
//      return;
//    if( mtx_invert((double *) &invZ, 4) == -1)
//      return;
//    for(i=0;i<4;i++)
//      for(ii=0;ii<4;ii++)
//      {
//        Y1.m[i][ii] = 0.5 * (Y.m[i][ii] + invZ.m[i][ii]);
//        Z1.m[i][ii] = 0.5 * (Z.m[i][ii] + invY.m[i][ii]);
//      }
//    Y = Y1;
//    Z = Z1;
//
//    mtx_mul(&Y2, &Y, &Y);
//  }
//  while(!almostequal(&Y2, mtx) && iter++ < 20 );
//
//  *mtx = Y;
//}

/*
   Check two matrices for near-enough equality
Params: a - first matrix
b - second matrix
Returns: 1 if almost equal, else 0, epsilon 0.0001f.
*/
//int almostequal(MATRIX *a, MATRIX *b)
//{
//  int i, ii;
//  double epsilon = 0.001f;
//
//  for(i=0;i<4;i++)
//    for(ii=0;ii<4;ii++)
//      if(fabs(a->m[i][ii] - b->m[i][ii]) > epsilon)
//        return 0;
//  return 1;
//}

/*
   multiply a point by a matrix.
Params: mtx - matrix
pt - the point (transformed)
*/
//void mulpt(MATRIX *mtx, double *pt)
//{
//  double ans[4] = {0};
//  int i;
//  int ii;
//
//  for(i=0;i<4;i++)
//  {
//    for(ii=0;ii<3;ii++)
//    {
//      ans[i] += pt[ii] * mtx->m[ii][i];
//    }
//    ans[i] += mtx->m[3][i];
//  }
//  pt[0] = ans[0];
//  pt[1] = ans[1];
//  pt[2] = ans[2];
//}

/*
   multiply two matrices.
Params: ans - return pointer for answer.
x - first matrix
y - second matrix.
Notes: ans may not be equal to x or y.
*/
//void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y)
//{
//  int i;
//  int ii;
//  int iii;
//
//  for(i=0;i<4;i++)
//    for(ii=0;ii<4;ii++)
//    {
//      ans->m[i][ii] = 0;
//      for(iii=0;iii<4;iii++)
//        ans->m[i][ii] += x->m[i][iii] * y->m[iii][ii];
//    }
//}


/*
   create an identity matrix.
Params: mtx - return pointer.
*/
//void mtx_identity(MATRIX *mtx)
//{
//  int i;
//  int ii;
//
//  for(i=0;i<4;i++)
//    for(ii=0;ii<4;ii++)
//    {
//      if(i==ii)
//        mtx->m[i][ii] = 1.0f;
//      else
//        mtx->m[i][ii] = 0;
//    }
//}

/*
 check the determinant of the (rotation) matrix
if it is right-handed, or left-handed
*/
//bool isRightHanded(MATRIX *mtx)
//{
//  //compute determinant
//  float det1 = mtx->m[1][1]*(mtx->m[2][2]*mtx->m[3][3] - mtx->m[2][3]*mtx->m[3][2]);
//  float det2 = - mtx->m[1][2]*(mtx->m[2][1]*mtx->m[3][3] - mtx->m[2][3]*mtx->m[3][1]);
//  float det3 = mtx->m[1][3]*(mtx->m[2][1]*mtx->m[3][2] - mtx->m[2][2]*mtx->m[3][1]);
//
//  if( (det1 + det2 +det3) > 0)
//    return 1;
//  else
//    return 0;
//}

/*
   create a translation matrix.
Params: mtx - return pointer for matrix.
x - x translation.
y - y translation.
z - z translation
*/
//void mtx_trans(MATRIX *mtx, double x, double y, double z)
//{
//  mtx->m[0][0] = 1;
//  mtx->m[0][1] = 0;
//  mtx->m[0][2] = 0;
//  mtx->m[0][3] = 0;
//
//  mtx->m[1][0] = 0;
//  mtx->m[1][1] = 1;
//  mtx->m[1][2] = 0;
//  mtx->m[1][3] = 0;
//
//  mtx->m[2][0] = 0;
//  mtx->m[2][1] = 0;
//  mtx->m[2][2] = 1;
//  mtx->m[2][3] = 0;
//
//  mtx->m[3][0] = x;
//  mtx->m[3][1] = y;
//  mtx->m[3][2] = z;
//  mtx->m[3][3] = 1;
//}

/*
   matrix invert routine
Params: mtx - the matrix in raw format, in/out
N - width and height
Returns: 0 on success, -1 on fail
*/
//int mtx_invert(double *mtx, int N)
//{
//  int indxc[100]; /* these 100s are the only restriction on matrix size */
//  int indxr[100];
//  int ipiv[100];
//  int i, j, k;
//  int irow, icol;
//  double big;
//  double pinv;
//  int l, ll;
//  double dum;
//  double temp;
//
//  assert(N <= 100);
//
//  for(i=0;i<N;i++)
//    ipiv[i] = 0;
//
//  for(i=0;i<N;i++)
//  {
//    big = 0.0;
//
//    /* find biggest element */
//    for(j=0;j<N;j++)
//      if(ipiv[j] != 1)
//        for(k=0;k<N;k++)
//          if(ipiv[k] == 0)
//            if(fabs(mtx[j*N+k]) >= big)
//            {
//              big = fabs(mtx[j*N+k]);
//              irow = j;
//              icol = k;
//            }
//
//    ipiv[icol]=1;
//
//    if(irow != icol)
//      for(l=0;l<N;l++)
//      {
//        temp = mtx[irow * N + l];
//        mtx[irow * N + l] = mtx[icol * N + l];
//        mtx[icol * N + l] = temp;
//      }
//
//    indxr[i] = irow;
//    indxc[i] = icol;
//
//
//    /* if biggest element is zero matrix is singular, bail */
//    if(mtx[icol* N + icol] == 0)
//      goto error_exit;
//
//    pinv = 1.0/mtx[icol * N + icol];
//
//    mtx[icol * N + icol] = 1.0;
//
//    for(l=0;l<N;l++)
//      mtx[icol * N + l] *= pinv;
//
//    for(ll=0;ll<N;ll++)
//      if(ll != icol)
//      {
//        dum = mtx[ll * N + icol];
//        mtx[ll * N + icol] = 0.0;
//        for(l=0;l<N;l++)
//          mtx[ll * N + l] -= mtx[icol * N + l]*dum;
//      }
//  }
//
//
//  /* unscramble matrix */
//  for (l=N-1;l>=0;l--)
//  {
//    if (indxr[l] != indxc[l])
//      for (k=0;k<N;k++)
//      {
//        temp = mtx[k * N + indxr[l]];
//        mtx[k * N + indxr[l]] = mtx[k * N + indxc[l]];
//        mtx[k * N + indxc[l]] = temp;
//      }
//  }
//
//  return 0;
//
//  error_exit:
//  return -1;
//}

/*
   get the asolute maximum of an array
   */
//double absmaxv(double *v, int N)
//{
//  double answer=0;
//  int i;
//
//  for(i=0;i<N;i++)
//    if(answer < fabs(v[i]))
//      answer = fabs(v[i]);
//  return answer;
//}

#include <stdio.h>

/*
   debug utlitiy
   */
void printmtx(FILE *fp, MATRIX *mtx)
{
  int i, ii;

  for(i=0;i<4;i++)
  {
    for(ii=0;ii<4;ii++)
      fprintf(fp, "%f, ", mtx->m[i][ii]);
    fprintf(fp, "\n");
  }
}

}//end namespace
