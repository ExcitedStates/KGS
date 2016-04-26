#include "metrics/RMSD.h"
#include "../SamplingOptions.h"
#include "core/Chain.h"

using namespace std;

namespace metrics{

RMSD::RMSD()
{}

double RMSD::distance(Configuration* c1, Configuration* c2)
{
  const std::vector<Atom*>* atomsRMSD = SamplingOptions::getOptions()->getAtomsAlign();

  // If atomsAlign is nullptr, align the entire m_protein
  if( atomsRMSD->empty() ) {
    //cout<<"Found no atoms to align .. using all."<<endl;
    atomsRMSD = &(c1->getMolecule()->atoms);//choose all atoms
  }

  Molecule * protein = c1->updatedMolecule();
  int atom_num = atomsRMSD->size();
  assert(atom_num>3);
  float* v1 = new float[atom_num*3];
  float* v2 = new float[atom_num*3];

  int resId;
  string name, chainName;
  unsigned int i=0;
  for (vector<Atom*>::const_iterator it=atomsRMSD->begin(); it!=atomsRMSD->end(); ++it) {
    name = (*it)->getName();
    chainName = (*it)->getResidue()->getChain()->getName();
    resId = (*it)->getResidue()->getId();
    Atom* atom = protein->getAtom(chainName,resId, name);
    v1[i]=atom->m_Position.x;
    v1[i+1]=atom->m_Position.y;
    v1[i+2]=atom->m_Position.z;
    i+=3;
  }

  protein= c2->updatedMolecule();
  i=0;
  for (vector<Atom*>::const_iterator it=atomsRMSD->begin(); it!=atomsRMSD->end(); ++it) {
    name = (*it)->getName();
    chainName = (*it)->getResidue()->getChain()->getName();
    resId = (*it)->getResidue()->getId();
    Atom* atom = protein->getAtom(chainName,resId, name);
    v2[i]=atom->m_Position.x;
    v2[i+1]=atom->m_Position.y;
    v2[i+2]=atom->m_Position.z;
    i+=3;
  }

  // Prepare the m_transformation matrix
  MATRIX mtx;
  // Call the rmsd procedure
  double diff = rmsd(v1,v2,atom_num,(float *) mtx.m);

  // clear memory
  delete[] v1;
  delete[] v2;

  return diff;
}

double RMSD::distance_noOptimization (Configuration *c1, Configuration *c2) {

  const std::vector<Atom*>* atomsRMSD = SamplingOptions::getOptions()->getAtomsMoving();

  // If atomsAlign is nullptr, align the entire m_protein
  if (atomsRMSD==nullptr) {
    cout<<"Found no atoms to align .. using all"<<endl;
    atomsRMSD = &(c1->getMolecule()->atoms);//choose all atoms
  }

  vector<Coordinate> p1_atoms;
  vector<Coordinate> p2_atoms;

  double sum=0;
  int atom_size = atomsRMSD->size();
  int resId;
  std::string name;
  std::string chainName;

  Molecule * p1 = c1->updatedMolecule();
  for (auto const& aIt : *atomsRMSD ) {
    name = aIt->getName();
    chainName = aIt->getResidue()->getChain()->getName();
    resId = aIt->getResidue()->getId();
    p1_atoms.push_back(p1->getAtom(chainName,resId, name)->m_Position);
  }

  Molecule * p2 = c2->updatedMolecule();
  for (auto const& aIt : *atomsRMSD ) {
    name = aIt->getName();
    chainName = aIt->getResidue()->getChain()->getName();
    resId = aIt->getResidue()->getId();
    p2_atoms.push_back(p2->getAtom(chainName,resId, name)->m_Position);
  }

  for(int i=0;i!=atom_size;i++){
    sum += p1_atoms[i].distanceSquared(p2_atoms[i]);
  }

  return sqrt(sum/atom_size);
}

double RMSD::align(Molecule * other, Molecule * base) {

  const std::vector<Atom*>* atomsAlign = SamplingOptions::getOptions()->getAtomsAlign();

  // If atomsAlign is nullptr, align the entire m_protein
  if (atomsAlign == nullptr) {
    cout<<"Found no atoms to align"<<endl;
    atomsAlign = &(base->atoms);//choose all atoms
  }

  int atom_size = atomsAlign->size();
  int coord_num = atom_size * 3;
  float* v1 = new float[coord_num];
  float* v2 = new float[coord_num];
  unsigned int i=0;
  Coordinate c1, c2;
  Atom* a2;
  std::string name;
  std::string chainName;

  for (auto const& aIt : *atomsAlign ) {
    name = aIt->getName();
    chainName = aIt->getResidue()->getChain()->getName();
    int resId = aIt->getResidue()->getId();
    c1 = base->getAtom(chainName, resId, name)->m_Position;
    a2 = other->getAtom(chainName, resId, name);
    if (a2) {
      c2 = a2->m_Position;
      v1[i] = c1.x;
      v1[i + 1] = c1.y;
      v1[i + 2] = c1.z;

      v2[i] = c2.x;
      v2[i + 1] = c2.y;
      v2[i + 2] = c2.z;
      i += 3;
    }
  }

  // Prepare the m_transformation matrix
  MATRIX mtx;
  // Call the rmsd procedure, v1 is base, v2 is mobile
  double diff = rmsd(v1,v2,atom_size,(float *) mtx.m);

  // Transform position in m_molecule other
  float* v3 = new float[3];
  for (vector<Atom*>::iterator it=other->atoms.begin(); it != other->atoms.end(); ++it) {
    Coordinate pos = (*it)->m_Position;
    v3[0] = pos.x;
    v3[1] = pos.y;
    v3[2] = pos.z;
    mulpt(&mtx,v3);
    (*it)->m_Position.x = v3[0];
    (*it)->m_Position.y = v3[1];
    (*it)->m_Position.z = v3[2];
  }

  // clear memory
  delete[] v1;
  delete[] v2;
  delete[] v3;

  return diff;
}


double rmsd(float *v1, float *v2, int N, float *mtx) {
  float cent1[3];
  float cent2[3];
  MATRIX tmtx;
  MATRIX tempmtx;
  MATRIX move1;
  MATRIX move2;
  int i;
  double answer;
  float *temp1 = 0;
  float *temp2 = 0;
  int err;

  assert(N > 3);

  temp1 = (float*) malloc(N * 3 * sizeof(float));
  temp2 = (float*) malloc(N * 3 * sizeof(float));
  if(!temp1 || !temp2)
    goto error_exit;

  centroid(cent1, v1, N);
  centroid(cent2, v2, N);
  for(i=0;i<N;i++)
  {
    temp1[i*3+0] = v1[i*3+0] - cent1[0];
    temp1[i*3+1] = v1[i*3+1] - cent1[1];
    temp1[i*3+2] = v1[i*3+2] - cent1[2];

    temp2[i*3+0] = v2[i*3+0] - cent2[0];
    temp2[i*3+1] = v2[i*3+1] - cent2[1];
    temp2[i*3+2] = v2[i*3+2] - cent2[2];
  }

  err = getalignmtx(temp1, temp2, N, &tmtx);
  if(err == -1)
    goto error_exit;

  mtx_trans(&move1, -cent2[0], -cent2[1], -cent2[2]);
  mtx_mul(&tempmtx, &move1, &tmtx);
  mtx_trans(&move2, cent1[0], cent1[1], cent1[2]);
  mtx_mul(&tmtx, &tempmtx, &move2);
  memcpy(temp2, v2, N * sizeof(float) * 3);
  for(i=0;i<N;i++)
    mulpt(&tmtx, temp2 + i * 3);
  answer = alignedrmsd(v1, temp2, N);
  free(temp1);
  free(temp2);
  if(mtx)
    memcpy(mtx, &tmtx.m, 16 * sizeof(float));

  return answer;
  error_exit:
  free(temp1);
  free(temp2);
  if(mtx)
  {
    for(i=0;i<16;i++)
      mtx[i] = 0;
  }
  return sqrt(-1.0);
}

/*
   calculate rmsd between two aligned structures (trivial)
Params: v1 - first structure
v2 - second structure
N - number of points
Returns: rmsd
*/
  double alignedrmsd(float *v1, float *v2, int N)
  {
    double answer =0;
    int i;

    for(i=0;i<N;i++)
      answer += vdiff2(v1 + i *3, v2 + i * 3);
    return sqrt(answer/N);
  }

  /*
     compute the centroid
     */
  void centroid(float *ret, float *v, int N)
  {
    int i;

    ret[0] = 0;
    ret[1] = 0;
    ret[2] = 0;
    for(i=0;i<N;i++)
    {
      ret[0] += v[i*3+0];
      ret[1] += v[i*3+1];
      ret[2] += v[i*3+2];
    }
    ret[0] /= N;
    ret[1] /= N;
    ret[2] /= N;
  }

  /*
     get the matrix needed to align two structures
Params: v1 - reference structure
v2 - structure to align
N - number of points
mtx - return for rigid body alignment matrix
Notes: only calculates rotation part of matrix.
assumes input has been aligned to centroids
*/
  int getalignmtx(float *v1, float *v2, int N, MATRIX *mtx)
  {
    MATRIX A = { {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}} };
    MATRIX At;
    MATRIX Ainv;
    MATRIX temp;
    float tv[3];
    float tw[3];
    float tv2[3];
    float tw2[3];
    int k, i, j;
    int flag = 0;
    float correction;

    correction = absmaxv(v1, N * 3) * absmaxv(v2, N * 3);

    for(k=0;k<N;k++)
      for(i=0;i<3;i++)
        for(j=0;j<3;j++)
          A.m[i][j] += (v1[k*3+i] * v2[k*3+j])/correction;

    while(flag < 3)
    {
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          At.m[i][j] = A.m[j][i];

      memcpy(&Ainv, &A, sizeof(MATRIX));
      /* this will happen if all points are in a plane */
      if( mtx_invert((float *) &Ainv, 4) == -1)
      {
        if(flag == 0)
        {
          crossproduct(tv, v1, v1+3);
          crossproduct(tw, v2, v2+3);
        }
        else
        {
          crossproduct(tv2, tv, v1);
          crossproduct(tw2, tw, v2);
          memcpy(tv, tv2, 3 * sizeof(float));
          memcpy(tw, tw2, 3 * sizeof(float));
        }
        for(i=0;i<3;i++)
          for(j=0;j<3;j++)
            A.m[i][j] += tv[i] * tw[j];

        flag++;
      }
      else
        flag = 5;
    }
    if(flag != 5)
      return -1;

    mtx_mul(&temp, &At, &A);
    mtx_root(&temp);
    mtx_mul(mtx, &temp, &Ainv);
    return 0;
  }

  /*
     get the crossproduct of two vectors.
Params: ans - return pinter for answer.
pt1 - first vector
pt2 - second vector.
Notes: crossproduct is at right angles to the two vectors.
*/
  void crossproduct(float *ans, float *pt1, float *pt2)
  {
    ans[0] = pt1[1] * pt2[2] - pt1[2] * pt2[1];
    ans[1] = pt1[0] * pt2[2] - pt1[2] * pt2[0];
    ans[2] = pt1[0] * pt2[1] - pt1[1] * pt2[0];
  }

  /*
     Denman-Beavers square root iteration
     */
  void mtx_root(MATRIX *mtx)
  {
    MATRIX Y = *mtx;
    MATRIX Z;
    MATRIX Y1;
    MATRIX Z1;
    MATRIX invY;
    MATRIX invZ;
    MATRIX Y2;
    int iter = 0;
    int i, ii;

    mtx_identity(&Z);

    do
    {
      invY = Y;
      invZ = Z;
      if( mtx_invert((float *) &invY, 4) == -1)
        return;
      if( mtx_invert((float *) &invZ, 4) == -1)
        return;
      for(i=0;i<4;i++)
        for(ii=0;ii<4;ii++)
        {
          Y1.m[i][ii] = 0.5 * (Y.m[i][ii] + invZ.m[i][ii]);
          Z1.m[i][ii] = 0.5 * (Z.m[i][ii] + invY.m[i][ii]);
        }
      Y = Y1;
      Z = Z1;

      mtx_mul(&Y2, &Y, &Y);
    }
    while(!almostequal(&Y2, mtx) && iter++ < 20 );

    *mtx = Y;
  }

  /*
     Check two matrices for near-enough equality
Params: a - first matrix
b - second matrix
Returns: 1 if almost equal, else 0, epsilon 0.0001f.
*/
  int almostequal(MATRIX *a, MATRIX *b)
  {
    int i, ii;
    float epsilon = 0.001f;

    for(i=0;i<4;i++)
      for(ii=0;ii<4;ii++)
        if(fabs(a->m[i][ii] - b->m[i][ii]) > epsilon)
          return 0;
    return 1;
  }

  /*
     multiply a point by a matrix.
Params: mtx - matrix
pt - the point (transformed)
*/
  void mulpt(MATRIX *mtx, float *pt)
  {
    float ans[4] = {0};
    int i;
    int ii;

    for(i=0;i<4;i++)
    {
      for(ii=0;ii<3;ii++)
      {
        ans[i] += pt[ii] * mtx->m[ii][i];
      }
      ans[i] += mtx->m[3][i];
    }
    pt[0] = ans[0];
    pt[1] = ans[1];
    pt[2] = ans[2];
  }

  /*
     multiply two matrices.
Params: ans - return pointer for answer.
x - first matrix
y - second matrix.
Notes: ans may not be equal to x or y.
*/
  void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y)
  {
    int i;
    int ii;
    int iii;

    for(i=0;i<4;i++)
      for(ii=0;ii<4;ii++)
      {
        ans->m[i][ii] = 0;
        for(iii=0;iii<4;iii++)
          ans->m[i][ii] += x->m[i][iii] * y->m[iii][ii];
      }
  }


  /*
     create an identity matrix.
Params: mtx - return pointer.
*/
  void mtx_identity(MATRIX *mtx)
  {
    int i;
    int ii;

    for(i=0;i<4;i++)
      for(ii=0;ii<4;ii++)
      {
        if(i==ii)
          mtx->m[i][ii] = 1.0f;
        else
          mtx->m[i][ii] = 0;
      }
  }

  /*
     create a translation matrix.
Params: mtx - return pointer for matrix.
x - x translation.
y - y translation.
z - z translation
*/
  void mtx_trans(MATRIX *mtx, float x, float y, float z)
  {
    mtx->m[0][0] = 1;
    mtx->m[0][1] = 0;
    mtx->m[0][2] = 0;
    mtx->m[0][3] = 0;

    mtx->m[1][0] = 0;
    mtx->m[1][1] = 1;
    mtx->m[1][2] = 0;
    mtx->m[1][3] = 0;

    mtx->m[2][0] = 0;
    mtx->m[2][1] = 0;
    mtx->m[2][2] = 1;
    mtx->m[2][3] = 0;

    mtx->m[3][0] = x;
    mtx->m[3][1] = y;
    mtx->m[3][2] = z;
    mtx->m[3][3] = 1;
  }

  /*
     matrix invert routine
Params: mtx - the matrix in raw format, in/out
N - width and height
Returns: 0 on success, -1 on fail
*/
  int mtx_invert(float *mtx, int N)
  {
    int indxc[100]; /* these 100s are the only restriction on matrix size */
    int indxr[100];
    int ipiv[100];
    int i, j, k;
    int irow, icol;
    double big;
    double pinv;
    int l, ll;
    double dum;
    double temp;

    assert(N <= 100);

    for(i=0;i<N;i++)
      ipiv[i] = 0;

    for(i=0;i<N;i++)
    {
      big = 0.0;

      /* find biggest element */
      for(j=0;j<N;j++)
        if(ipiv[j] != 1)
          for(k=0;k<N;k++)
            if(ipiv[k] == 0)
            if(fabs(mtx[j*N+k]) >= big)
            {
              big = fabs(mtx[j*N+k]);
              irow = j;
              icol = k;
            }

      ipiv[icol]=1;

      if(irow != icol)
        for(l=0;l<N;l++)
        {
          temp = mtx[irow * N + l];
          mtx[irow * N + l] = mtx[icol * N + l];
          mtx[icol * N + l] = temp;
        }

      indxr[i] = irow;
      indxc[i] = icol;


      /* if biggest element is zero matrix is singular, bail */
      if(mtx[icol* N + icol] == 0)
        goto error_exit;

      pinv = 1.0/mtx[icol * N + icol];

      mtx[icol * N + icol] = 1.0;

      for(l=0;l<N;l++)
        mtx[icol * N + l] *= pinv;

      for(ll=0;ll<N;ll++)
        if(ll != icol)
        {
          dum = mtx[ll * N + icol];
          mtx[ll * N + icol] = 0.0;
          for(l=0;l<N;l++)
            mtx[ll * N + l] -= mtx[icol * N + l]*dum;
        }
    }


    /* unscramble matrix */
    for (l=N-1;l>=0;l--)
    {
      if (indxr[l] != indxc[l])
        for (k=0;k<N;k++)
        {
          temp = mtx[k * N + indxr[l]];
          mtx[k * N + indxr[l]] = mtx[k * N + indxc[l]];
          mtx[k * N + indxc[l]] = temp;
        }
    }

    return 0;

    error_exit:
    return -1;
  }

  /*
     get the asolute maximum of an array
     */
  float absmaxv(float *v, int N)
  {
    float answer=0;
    int i;

    for(i=0;i<N;i++)
      if(answer < fabs(v[i]))
        answer = fabs(v[i]);
    return answer;
  }

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
}
