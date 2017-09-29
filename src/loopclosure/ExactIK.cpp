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

//!----------------------------------------------------------------------
//! This code was adapted from code which is protected by the following copyright:
//!	Copyright (C) 2003
//!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill
//!      UCSF and University of New Mexico
//! 	Witten by Chaok Seok 2003.
//!----------------------------------------------------------------------------

#include "ExactIK.h"

#include <cassert>
#include <cmath>

#include "core/Molecule.h"
#include "core/Chain.h"

using namespace std;

std::vector<Configuration*> ExactIK::rebuildLoop(
    const Residue* res1,
    const Residue* res2,
    const Residue* res3
)
{

//  cout<<"ExactIK:rebuildLoop("<<res1->getId()<<","<<res2->getId()<<","<<res3->getId()<<")"<<endl;
  double blen[6];
  double bang[7];
  double tang[2];

  initializeIKParams(res1, res2, res3, blen, bang, tang);
  initialize_loop_closure(blen, bang, tang);

  double r_n1[3], r_a1[3], r_a3[3], r_c3[3];

  for (int i=0; i<3; i++){
    r_n1[i] = res1->getAtom( "N")->m_position[i];
    r_a1[i] = res1->getAtom("CA")->m_position[i];
    r_a3[i] = res3->getAtom("CA")->m_position[i];
    r_c3[i] = res3->getAtom( "C")->m_position[i];
  }

  double r_soln_n[16][3][3], r_soln_a[16][3][3], r_soln_c[16][3][3];
  int n_soln;

  solve_3pep_poly(r_n1, r_a1, r_a3, r_c3,
                  r_soln_n,  r_soln_a,
                  r_soln_c, &n_soln);


//  for(int i=0;i<n_soln;i++){
//    cout<<"MODEL "<<(i+1)<<endl;
//    for(int r=0;r<3;r++) {
//      const Residue* res = res1;
//      if(r==1) res=res2;
//      if(r==2) res=res3;
//      printf("ATOM      1  N   ALA A %3d    %8.3f%8.3f%8.3f  1.00  1.00           N\n", res->getId(),
//             r_soln_n[i][r][0], r_soln_n[i][r][1], r_soln_n[i][r][2]);
//      printf("ATOM      2  CA  ALA A %3d    %8.3f%8.3f%8.3f  1.00  1.00           C\n", res->getId(),
//             r_soln_a[i][r][0], r_soln_a[i][r][1], r_soln_a[i][r][2]);
//      printf("ATOM      3  C   ALA A %3d    %8.3f%8.3f%8.3f  1.00  1.00           C\n", res->getId(),
//             r_soln_c[i][r][0], r_soln_c[i][r][1], r_soln_c[i][r][2]);
//    }
//    cout<<"ENDMDL"<<endl;
//  }

  Configuration* parent = res1->getChain()->getMolecule()->m_conf;
  if(parent==nullptr) parent = new Configuration(res1->getChain()->getMolecule());
  vector<Configuration*> ret;

  for(int i=0;i<n_soln;i++){
    Configuration* child = new Configuration(parent);
    for(int d=0;d<parent->getNumDOFs();d++)
      child->m_dofs[d] = parent->m_dofs[d];

    Atom* C0 =res1->getPrevResidue()->getAtom("C");
    Atom* N1 = res1->getAtom( "N");
    Atom* A1 = res1->getAtom("CA");
    Atom* C1 = res1->getAtom( "C");
    Atom* N2 = res2->getAtom( "N");
    Atom* A2 = res2->getAtom("CA");
    Atom* C2 = res2->getAtom( "C");
    Atom* N3 = res3->getAtom( "N");
    Atom* A3 = res3->getAtom("CA");
    Atom* C3 = res3->getAtom( "C");
    Atom* N4 = res3->getNextResidue()->getAtom("N");
    Coordinate _C1(r_soln_c[i][0][0],r_soln_c[i][0][1],r_soln_c[i][0][2]);
    Coordinate _N2(r_soln_n[i][1][0],r_soln_n[i][1][1],r_soln_n[i][1][2]);
    Coordinate _A2(r_soln_a[i][1][0],r_soln_a[i][1][1],r_soln_a[i][1][2]);
    Coordinate _C2(r_soln_c[i][1][0],r_soln_c[i][1][1],r_soln_c[i][1][2]);
    Coordinate _N3(r_soln_n[i][2][0],r_soln_n[i][2][1],r_soln_n[i][2][2]);

    double oldPhi1 = TorsionalAngle(C0->m_position, N1->m_position, A1->m_position, C1->m_position);
    double newPhi1 = TorsionalAngle(C0->m_position, N1->m_position, A1->m_position, _C1);
    double delPhi1 = newPhi1 - oldPhi1;

    double oldPsi1 = TorsionalAngle(N1->m_position, A1->m_position, C1->m_position, N2->m_position);
    double newPsi1 = TorsionalAngle(N1->m_position, A1->m_position, _C1, _N2);
    double delPsi1 = newPsi1 - oldPsi1;

    double oldPhi2 = TorsionalAngle(C1->m_position, N2->m_position, A2->m_position, C2->m_position);
    double newPhi2 = TorsionalAngle(_C1, _N2, _A2, _C2);
    double delPhi2 = newPhi2 - oldPhi2;

    double oldPsi2 = TorsionalAngle(N2->m_position, A2->m_position, C2->m_position, N3->m_position);
    double newPsi2 = TorsionalAngle(_N2, _A2, _C2, _N3);
    double delPsi2 = newPsi2 - oldPsi2;

    double oldPhi3 = TorsionalAngle(C2->m_position, N3->m_position, A3->m_position, C3->m_position);
    double newPhi3 = TorsionalAngle(_C2, _N3, A3->m_position, C3->m_position);
    double delPhi3 = newPhi3 - oldPhi3;

    double oldPsi3 = TorsionalAngle(N3->m_position, A3->m_position, C3->m_position, N4->m_position);
    double newPsi3 = TorsionalAngle(_N3, A3->m_position, C3->m_position, N4->m_position);
    double delPsi3 = newPsi3 - oldPsi3;

    double e = 0.001;
    if(   fabs(delPhi1)<e && fabs(delPsi1)<e &&
          fabs(delPhi2)<e && fabs(delPsi2)<e &&
          fabs(delPhi3)<e && fabs(delPsi3)<e ) {
      delete child;
      continue;
    }
//    cout<<"Setting dofs"<<endl;

    for(auto const& edge: parent->getMolecule()->m_spanningTree->Edges){
      Bond* b = edge->getBond();
      if(b==nullptr) continue;
      if(b->m_atom1==N1 && b->m_atom2==A1) { child->m_dofs[edge->getDOF()->getIndex()] += delPhi1; }
      if(b->m_atom1==A1 && b->m_atom2==C1) { child->m_dofs[edge->getDOF()->getIndex()] += delPsi1; }
      if(b->m_atom1==N2 && b->m_atom2==A2) { child->m_dofs[edge->getDOF()->getIndex()] += delPhi2; }
      if(b->m_atom1==A2 && b->m_atom2==C2) { child->m_dofs[edge->getDOF()->getIndex()] += delPsi2; }
      if(b->m_atom1==N3 && b->m_atom2==A3) { child->m_dofs[edge->getDOF()->getIndex()] += delPhi3; }
      if(b->m_atom1==A3 && b->m_atom2==C3) { child->m_dofs[edge->getDOF()->getIndex()] += delPsi3; }
    }

    ret.push_back(child);
  }

  return ret;
}


bool ExactIK::validRebuildLoop(const Residue* res1, const Residue* res2, const Residue* res3 ) const
{
  //Ensure that res1 < res2 < res3
  if(res1->getId()>=res2->getId()) return false;
  if(res2->getId()>=res3->getId()) return false;

  //Prolines have fixed phi-angle so wont work with rebuild
  if(res1->getName()=="PRO") return false;
  if(res2->getName()=="PRO") return false;
  if(res3->getName()=="PRO") return false;

  //Check that all residues are in the same chain
  if(res1->getChain()!=res2->getChain()) return false;
  if(res1->getChain()!=res2->getChain()) return false;

  //Check that all residues are amino acids
  if(res1->getAtom("C")==nullptr) return false;
  if(res1->getAtom("CA")==nullptr) return false;
  if(res1->getAtom("N")==nullptr) return false;
  if(res2->getAtom("C")==nullptr) return false;
  if(res2->getAtom("CA")==nullptr) return false;
  if(res2->getAtom("N")==nullptr) return false;
  if(res3->getAtom("C")==nullptr) return false;
  if(res3->getAtom("CA")==nullptr) return false;
  if(res3->getAtom("N")==nullptr) return false;

  //Check for constraints between res1 and res3
//  for(auto const& edge: res1->getChain()->getMolecule()->m_spanningTree->Edges){
//    Bond* bond = edge->getBond();
//    TODO: check for constraints
//  }

  return true;
}

vector< tuple<Residue*> > ExactIK::findAllValidTriples(Configuration* conf) const
{
  vector< tuple<Residue*> > ret;

  Molecule* mol = conf->getMolecule();
  for(auto const& chain: mol->m_chains){

  }

  return ret;
}

void ExactIK::initializeIKParams(
    const Residue* res1,
    const Residue* res2,
    const Residue* res3,
    double blen[6],
    double bang[7],
    double tang[2]
)
{
  Coordinate& N1 = res1->getAtom( "N")->m_position;
  Coordinate& A1 = res1->getAtom("CA")->m_position;
  Coordinate& C1 = res1->getAtom( "C")->m_position;
  Coordinate& N2 = res2->getAtom( "N")->m_position;
  Coordinate& A2 = res2->getAtom("CA")->m_position;
  Coordinate& C2 = res2->getAtom( "C")->m_position;
  Coordinate& N3 = res3->getAtom( "N")->m_position;
  Coordinate& A3 = res3->getAtom("CA")->m_position;
  Coordinate& C3 = res3->getAtom( "C")->m_position;

  //Initialize lengths
  blen[0] = A1.distanceTo(C1);
  blen[1] = C1.distanceTo(N2);
  blen[2] = N2.distanceTo(A2);
  blen[3] = A2.distanceTo(C2);
  blen[4] = C2.distanceTo(N3);
  blen[5] = N3.distanceTo(A3);

  //Initialize angles
  bang[0] = Angle(N1,A1,C1);
  bang[1] = Angle(A1,C1,N2);
  bang[2] = Angle(C1,N2,A2);
  bang[3] = Angle(N2,A2,C2);
  bang[4] = Angle(A2,C2,N3);
  bang[5] = Angle(C2,N3,A3);
  bang[6] = Angle(N3,A3,C3);

  //Initialize torsions
  tang[0] = TorsionalAngle( A1, C1, N2, A2 );
  tang[1] = TorsionalAngle( A2, C2, N3, A3 );
}

void ExactIK::initialize_loop_closure(double b_len[6], double b_ang[7], double t_ang[2])
{
  double len1, len2, a_min, a_max;
  double axis[3], rr_a1[3], rr_c1[3], rr_n2[3], rr_a2[3], rr_n2a2_ref[3], rr_c1a1[3];
  double rr_a1a2[3], dr[3], bb_c1a1[3], bb_a1a2[3], bb_a2n2[3];
  double p[4], Us[3][3];
  double mulpro[3];
  double tmp_val[3];
  double tol_secant = 1.0e-15;
  int max_iter_sturm = 100;
  int max_iter_secant = 20;
  int i, j;

  initialize_sturm(&tol_secant, &max_iter_sturm, &max_iter_secant);

  for(i=0;i<6;i++)
    len0[i] = b_len[i];
  for(i=0;i<7;i++)
    b_ang0[i] = b_ang[i];
  for(i=0;i<2;i++)
    t_ang0[i] = t_ang[i];

  for(i=0;i<3;i++)
    rr_c1[i] = 0.;
  axis[0] = 1.;
  axis[1] = 0.;
  axis[2] = 0.;

  for(i=0;i<2;i++)
  {
    rr_a1[0] = cos(b_ang0[3*i+1])*len0[3*i];
    rr_a1[1] = sin(b_ang0[3*i+1])*len0[3*i];
    rr_a1[2] = 0.0e0;
    rr_n2[0] = len0[3*i+1];
    rr_n2[1] = 0.0e0;
    rr_n2[2] = 0.0e0;
    for(j=0;j<3;j++)
      rr_c1a1[j] = rr_a1[j] - rr_c1[j];
    rr_n2a2_ref[0] = -cos(b_ang0[3*i+2])*len0[3*i+2];
    rr_n2a2_ref[1] = sin(b_ang0[3*i+2])*len0[3*i+2];
    rr_n2a2_ref[2] = 0.0e0;
    quaternion(axis, t_ang0[i]*0.25e0, p);
    rotation_matrix(p, Us);
    matmul(Us, rr_n2a2_ref, mulpro);
    for(j=0;j<3;j++)
    {
      rr_a2[j] =  mulpro[j] + rr_n2[j];
      rr_a1a2[j] = rr_a2[j] - rr_a1[j];
      dr[j] = rr_a1a2[j];
    }
    len2 = dot_product(dr, dr);
    len1 = sqrt(len2);
    len_aa[i+1] = len1;
    for(j=0;j<3;j++)
    {
      bb_c1a1[j] = rr_c1a1[j]/len0[3*i];
      bb_a1a2[j] = rr_a1a2[j]/len1;
      bb_a2n2[j] = (rr_n2[j] - rr_a2[j])/len0[3*i+2];
    }
    for(j=0;j<3;j++)
      tmp_val[j] = -bb_a1a2[j];
    calc_bnd_ang(tmp_val, bb_a2n2, &xi[i+1]);
    for(j=0;j<3;j++)
      tmp_val[j] = -bb_c1a1[j];
    calc_bnd_ang(bb_a1a2, tmp_val, &eta[i]);
    calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2, &delta[i+1]);
    delta[i+1] = M_PI - delta[i+1];
  }

  a_min = b_ang[3] - (xi[1] + eta[1]);
  a_max = std::min((b_ang[3] + (xi[1] + eta[1])), M_PI);

  aa13_min_sqr = pow(len_aa[1],2) + pow(len_aa[2],2) - 2.0e0*len_aa[1]*len_aa[2]*cos(a_min);
  aa13_max_sqr = pow(len_aa[1],2) + pow(len_aa[2],2) - 2.0e0*len_aa[1]*len_aa[2]*cos(a_max);
}

double ExactIK::dot_product(double va[3], double vb[3])
{
  return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
}

void ExactIK::matmul(double ma[3][3], double mb[3], double mc[3])
{
  int i, j;

  for(i=0;i<3;i++)
  {
    mc[i] = 0.;
    for(j=0;j<3;j++)
      mc[i] += ma[i][j]*mb[j];
  }
  return;
}

double ExactIK::sign(double a, double b)
{
  if (b>=0.)
    return fabs(a);
  else
    return -fabs(a);
}


void ExactIK::solve_3pep_poly(
    double r_n1[3],
    double r_a1[3],
    double r_a3[3],
    double r_c3[3],
    double r_soln_n[16][3][3],
    double r_soln_a[16][3][3],
    double r_soln_c[16][3][3],
    int *n_soln
)
{
  double poly_coeff[deg_pol+1], roots[max_soln];
  get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3);
  if (*n_soln == 0)
    return;

  get_poly_coeff(poly_coeff);
  solve_sturm(&deg_pol, n_soln, poly_coeff, roots);

  if (*n_soln == 0)
    return;

  coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c);

  return;
}


//! calculate quaternion, given rotation axis and angle.
void ExactIK::quaternion(double axis[3], double quarter_ang, double p[4]) const
{
  double tan_w = tan(quarter_ang);
  double tan_sqr = tan_w * tan_w;
  double tan1 = 1.0e0 + tan_sqr;
  double cosine = (1.0e0 - tan_sqr)/tan1;
  double sine = 2.0e0*tan_w/tan1;
  p[0] = cosine;
  p[1] = axis[0] * sine;
  p[2] = axis[1] * sine;
  p[3] = axis[2] * sine;
}


void ExactIK::rotation_matrix(double q[4], double U[3][3]) const
{
//! constructs rotation matrix U from quaternion q.
  double q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33;

  q0 = q[0];
  q1 = q[1];
  q2 = q[2];
  q3 = q[3];
  b0 = 2.0e0*q0;
  b1 = 2.0e0*q1;
  q00 = b0*q0-1.0e0;
  q02 = b0*q2;
  q03 = b0*q3;
  q11 = b1*q1;
  q12 = b1*q2;
  q13 = b1*q3;
  b2 = 2.0e0*q2;
  b3 = 2.0e0*q3;
  q01 = b0*q1;
  q22 = b2*q2;
  q23 = b2*q3;
  q33 = b3*q3;
  U[0][0] = q00+q11;
  U[0][1] = q12-q03;
  U[0][2] = q13+q02;
  U[1][0] = q12+q03;
  U[1][1] = q00+q22;
  U[1][2] = q23-q01;
  U[2][0] = q13-q02;
  U[2][1] = q23+q01;
  U[2][2] = q00+q33;

  return;
}

void ExactIK::get_input_angles(int *n_soln, double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3])
{
//! Input angles and vectors (later used in coordinates)
  double dr_sqr;
  double tmp_val[3];
  int i;
  char cone_type[2];

  *n_soln = max_soln;

  //  ! vertual bond

  for(i=0;i<3;i++)
    r_a1a3[i] = r_a3[i] - r_a1[i];
  dr_sqr = dot_product(r_a1a3,r_a1a3);
  len_aa[0] = sqrt(dr_sqr);

  if ((dr_sqr < aa13_min_sqr) || (dr_sqr > aa13_max_sqr))
  {
    *n_soln = 0;
    return;
  }

//  ! bond lengths
  for(i=0;i<3;i++)
    r_a1n1[i] = r_n1[i] - r_a1[i];
  len_na[0] = sqrt(dot_product(r_a1n1,r_a1n1));
  len_na[1] = len0[2];
  len_na[2] = len0[5];
  for(i=0;i<3;i++)
    r_a3c3[i] = r_c3[i] - r_a3[i];
  len_ac[0] = len0[0];
  len_ac[1] = len0[3];
  len_ac[2] = sqrt(dot_product(r_a3c3,r_a3c3));

//  ! unit vectors
  for(i=0;i<3;i++)
  {
    b_a1n1[i] = r_a1n1[i]/len_na[0];
    b_a3c3[i] = r_a3c3[i]/len_ac[2];
    b_a1a3[i] = r_a1a3[i]/len_aa[0];
  }

//  ! delta(3): dih of N(1)CA(1)CA(3)C(3)
  for(i=0;i<3;i++)
    tmp_val[i] = -b_a1n1[i];
  calc_dih_ang(tmp_val, b_a1a3, b_a3c3, &delta[3]);
  delta[0] = delta[3];

//  ! xi(1)
  for(i=0;i<3;i++)
    tmp_val[i] = -b_a1a3[i];
  calc_bnd_ang(tmp_val, b_a1n1, &xi[0]);

//  ! eta(3)
  calc_bnd_ang(b_a1a3, b_a3c3, &eta[2]);

  for(i=0;i<3;i++)
  {
    cos_delta[i+1] = cos(delta[i+1]);
    sin_delta[i+1] = sin(delta[i+1]);
    cos_xi[i] = cos(xi[i]);
    sin_xi[i] = sin(xi[i]);
    sin_xi[i] = sin(xi[i]);
    cos_eta[i] = cos(eta[i]);
    cos_eta[i] = cos(eta[i]);
    sin_eta[i] = sin(eta[i]);
    sin_eta[i] = sin(eta[i]);
  }
  cos_delta[0] = cos_delta[3];
  sin_delta[0] = sin_delta[3];

//  ! theta (N, CA, C) bond angle
  theta[0] = b_ang0[0];
  theta[1] = b_ang0[3];
  theta[2] = b_ang0[6];
  for(i=0;i<3;i++)
    cos_theta[i] = cos(theta[i]);

//  ! alpha
  cos_alpha[0] = -(pow(len_aa[0],2) + pow(len_aa[1],2) - pow(len_aa[2],2))/(2.0e0*len_aa[0]*len_aa[1]);
  alpha[0] = acos(cos_alpha[0]);
  sin_alpha[0] = sin(alpha[0]);
  cos_alpha[1] = (pow(len_aa[1],2) + pow(len_aa[2],2) - pow(len_aa[0],2))/(2.0e0*len_aa[1]*len_aa[2]);
  alpha[1] = acos(cos_alpha[1]);
  sin_alpha[1] = sin(alpha[1]);
  alpha[2] = M_PI - alpha[0] + alpha[1];
  cos_alpha[2] = cos(alpha[2]);
  sin_alpha[2] = sin(alpha[2]);

//  ! check for existence of soln
  for(i=0;i<3;i++)
  {
    test_two_cone_existence_soln(theta[i], xi[i], eta[i], alpha[i], n_soln, cone_type);
    if (*n_soln == 0)
      return;
  }

  return;
}


void ExactIK::test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln, char cone_type[2])
{
  double at, ex, abs_at, ap1, kx1, et1;
  double cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2;
  int s1, s2, t1, t2;
  int complicated = 0;

  *n_soln = max_soln;

  ap1 = ap;
  kx1 = kx;
  et1 = et;

  at = ap1 - tt;
  ex = kx1 + et1;
  abs_at = fabs(at);

//  ! case of no soln
  if (abs_at > ex)
  {
    *n_soln = 0;
    return;
  }

  if (complicated)
  {
    cos_tx1 = cos(tt+kx1);
    cos_tx2 = cos(tt-kx1);
    cos_te1 = cos(tt+et1);
    cos_te2 = cos(tt-et1);
    cos_ea1 = cos(et1+ap1);
    cos_ea2 = cos(et1-ap1);
    cos_xa1 = cos(kx1+ap1);
    cos_xa2 = cos(kx1-ap1);
    s1 = 0;
    s2 = 0;
    t1 = 0;
    t2 = 0;
    if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0e0)
      s1 = 0;
    if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0e0)
      s2 = 0;
    if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0e0)
      t1 = 0;
    if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0e0)
      t2 = 0;
  }
}

void ExactIK::get_poly_coeff(double poly_coeff[16+1])
{
  int i, j;
  double A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42;
  double B0[3], B1[3], B2[3], B3[3], B4[3], B5[3], B6[3], B7[3], B8[3];
  double u11[5][5], u12[5][5], u13[5][5], u31[5][5], u32[5][5], u33[5][5];
  double um1[5][5], um2[5][5], um3[5][5], um4[5][5], um5[5][5], um6[5][5], q_tmp[5][5];
  int p1 [2], p3[2], p_um1[2], p_um2[2], p_um3[2], p_um4[2], p_um5[2], p_um6[2], p_Q[2];
  int p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, p_f8, p_f9;
  int p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, p_f18;
  int p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26;
  int p_final;
  double f1[17], f2[17], f3[17], f4[17], f5[17], f6[17], f7[17], f8[17], f9[17];
  double f10[17], f11[17], f12[17], f13[17], f14[17], f15[17], f16[17], f17[17], f18[17];
  double f19[17], f20[17], f21[17], f22[17], f23[17], f24[17], f25[17], f26[17];

//  ! A0, B0
  for(i=0;i<3;i++)
  {
    A0 = cos_alpha[i]*cos_xi[i]*cos_eta[i] - cos_theta[i];
    A1 = -sin_alpha[i]*cos_xi[i]*sin_eta[i];
    A2 = sin_alpha[i]*sin_xi[i]*cos_eta[i];
    A3 = sin_xi[i]*sin_eta[i];
    A4 = A3*cos_alpha[i];
    j = i;
    A21 = A2*cos_delta[j];
    A22 = A2*sin_delta[j];
    A31 = A3*cos_delta[j];
    A32 = A3*sin_delta[j];
    A41 = A4*cos_delta[j];
    A42 = A4*sin_delta[j];
    B0[i] = A0 + A22 + A31;
    B1[i] = 2.0e0*(A1 + A42);
    B2[i] = 2.0e0*(A32 - A21);
    B3[i] = -4.0e0*A41;
    B4[i] = A0 + A22 - A31;
    B5[i] = A0 - A22 - A31;
    B6[i] = -2.0e0*(A21 + A32);
    B7[i] = 2.0e0*(A1 - A42);
    B8[i] = A0 - A22 + A31;
  }

//  ! C0i
  i = 0;
  C0[i][0] = B0[i];
  C0[i][1] = B2[i];
  C0[i][2] = B5[i];
  C1[i][0] = B1[i];
  C1[i][1] = B3[i];
  C1[i][2] = B7[i];
  C2[i][0] = B4[i];
  C2[i][1] = B6[i];
  C2[i][2] = B8[i];
  for(i=1;i<3;i++)
  {
    C0[i][0] = B0[i];
    C0[i][1] = B1[i];
    C0[i][2] = B4[i];
    C1[i][0] = B2[i];
    C1[i][1] = B3[i];
    C1[i][2] = B6[i];
    C2[i][0] = B5[i];
    C2[i][1] = B7[i];
    C2[i][2] = B8[i];
  }

//  ! first determinant
  for(i=0;i<3;i++)
  {
    u11[0][i] = C0[0][i];
    u12[0][i] = C1[0][i];
    u13[0][i] = C2[0][i];
    u31[i][0] = C0[1][i];
    u32[i][0] = C1[1][i];
    u33[i][0] = C2[1][i];
  }

  p1[0] = 2;
  p1[1] = 0;
  p3[0] = 0;
  p3[1] = 2;

  poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1);
  poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2);
  poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3);
  poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4);
  poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5);
  poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6);
  poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q);

  for(i=0;i<5;i++)
    for(j=0;j<5;j++)
      Q[i][j] = q_tmp[i][j];

//  ! second determinant
  for(i=0;i<3;i++)
    for(j=0;j<17;j++)
      R[i][j] = 0.;
  for(i=0;i<3;i++)
  {
    R[0][i] = C0[2][i];
    R[1][i] = C1[2][i];
    R[2][i] = C2[2][i];
  }
  p2 = 2;
  p4 = 4;

  poly_mul_sub1(R[1], R[1], R[0], R[2], p2, p2, p2, p2, f1, &p_f1);
  poly_mul1(R[1], R[2], p2, p2, f2, &p_f2);
  poly_mul_sub1(R[1], f1, R[0], f2, p2, p_f1, p2, p_f2, f3, &p_f3);
  poly_mul1(R[2], f1, p2, p_f1, f4, &p_f4);
  poly_mul_sub1(R[1], f3, R[0], f4, p2, p_f3, p2, p_f4, f5, &p_f5);

  poly_mul_sub1(Q[1], R[1], Q[0], R[2], p4, p2, p4, p2, f6, &p_f6);
  poly_mul_sub1(Q[2], f1, R[2], f6, p4, p_f1, p2, p_f6, f7, &p_f7);
  poly_mul_sub1(Q[3], f3, R[2], f7, p4, p_f3, p2, p_f7, f8, &p_f8);
  poly_mul_sub1(Q[4], f5, R[2], f8, p4, p_f5, p2, p_f8, f9, &p_f9);

  poly_mul_sub1(Q[3], R[1], Q[4], R[0], p4, p2, p4, p2, f10, &p_f10);
  poly_mul_sub1(Q[2], f1, R[0], f10, p4, p_f1, p2, p_f10, f11, &p_f11);
  poly_mul_sub1(Q[1], f3, R[0], f11, p4, p_f3, p2, p_f11, f12, &p_f12);

  poly_mul_sub1(Q[2], R[1], Q[1], R[2], p4, p2, p4, p2, f13, &p_f13);
  poly_mul_sub1(Q[3], f1, R[2], f13, p4, p_f1, p2, p_f13, f14, &p_f14);
  poly_mul_sub1(Q[3], R[1], Q[2], R[2], p4, p2, p4, p2, f15, &p_f15);
  poly_mul_sub1(Q[4], f1, R[2], f15, p4, p_f1, p2, p_f15, f16, &p_f16);
  poly_mul_sub1(Q[1], f14, Q[0], f16, p4, p_f14, p4, p_f16, f17, &p_f17);

  poly_mul_sub1(Q[2], R[2], Q[3], R[1], p4, p2, p4, p2, f18, &p_f18);
  poly_mul_sub1(Q[1], R[2], Q[3], R[0], p4, p2, p4, p2, f19, &p_f19);
  poly_mul_sub1(Q[3], f19, Q[2], f18, p4, p_f19, p4, p_f18, f20, &p_f20);
  poly_mul_sub1(Q[1], R[1], Q[2], R[0], p4, p2, p4, p2, f21, &p_f21);
  poly_mul1(Q[4], f21, p4, p_f21, f22, &p_f22);
  poly_sub1(f20, f22, p_f20, p_f22, f23, &p_f23);
  poly_mul1(R[0], f23, p2, p_f23, f24, &p_f24);
  poly_sub1(f17, f24, p_f17, p_f24, f25, &p_f25);
  poly_mul_sub1(Q[4], f12, R[2], f25, p4, p_f12, p2, p_f25, f26, &p_f26);
  poly_mul_sub1(Q[0], f9, R[0], f26, p4, p_f9, p2, p_f26, poly_coeff, &p_final);

  if (p_final != 16)
  {
    printf("Error. Degree of polynomial is not 16!\n");
    exit(1);
  }

  if (poly_coeff[16] < 0.0e0)
    for(i=0;i<17;i++)
      poly_coeff[i] *= -1.0;

}


void ExactIK::poly_mul_sub2(double u1[5][5], double u2[5][5], double u3[5][5], double u4[5][5], int p1[2], int p2[2], int p3[2], int p4[2], double u5[5][5], int p5[2])
{
  double d1[5][5], d2[5][5];
//  integer, dimension(2) :: pd1, pd2
  int pd1[2], pd2[2];

//  call poly_mul2(u1, u2, p1, p2, d1, pd1)
  poly_mul2(u1, u2, p1, p2, d1, pd1);
//  call poly_mul2(u3, u4, p3, p4, d2, pd2)
  poly_mul2(u3, u4, p3, p4, d2, pd2);
//  call poly_sub2(d1, d2, pd1, pd2, u5, p5)
  poly_sub2(d1, d2, pd1, pd2, u5, p5);
}

void ExactIK::poly_mul2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2])
{
//  implicit none
//  real(dp), dimension(0:4,0:4), intent(in) :: u1, u2
//  integer, dimension(2), intent(in) :: p1, p2
//  real(dp), dimension(0:4,0:4), intent(out) :: u3
//  integer, intent(out) :: p3(2)
//  integer :: i1, j1, i2, j2, i3, j3, p11, p12, p21, p22
  int i1, j1, i2, j2, i3, j3, p11, p12, p21, p22;
  int i, j;
//  real(dp) :: u1ij
  double u1ij;

//  p3(:) = p1(:) + p2(:)
  for(i=0;i<2;i++)
    p3[i] = p1[i] + p2[i];
  for(i=0;i<5;i++)
    for(j=0;j<5;j++)
      u3[i][j] = 0.0e0;

//  p11 = p1(1)
  p11 = p1[0];
//  p12 = p1(2)
  p12 = p1[1];
//  p21 = p2(1)
  p21 = p2[0];
//  p22 = p2(2)
  p22 = p2[1];

//  do i1 = 0, p12
  for(i1=0;i1<=p12;i1++)
  {
//     do j1 = 0, p11
    for(j1=0;j1<=p11;j1++)
    {
//        u1ij = u1(j1,i1)
      u1ij = u1[i1][j1];
//        do i2 = 0, p22
      for(i2=0;i2<=p22;i2++)
      {
//           i3 = i1 + i2
        i3 = i1 + i2;
//           do j2 = 0, p21
        for (j2=0;j2<=p21;j2++)
        {
//              j3 = j1 + j2
          j3 = j1 + j2;
//              u3(j3,i3) = u3(j3,i3) + u1ij*u2(j2,i2)
          u3[i3][j3] = u3[i3][j3] + u1ij*u2[i2][j2];
//           end do
        }
//        end do
      }
//     end do
    }
//  end do
  }
}

void ExactIK::poly_sub2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2])
{
//  implicit none
//  real(dp), dimension(0:4,0:4), intent(in) :: u1, u2
//  integer, intent(in) :: p1(2), p2(2)
//  real(dp), dimension(0:4,0:4), intent(out) :: u3
//  integer, intent(out) :: p3(2)
//  integer :: i, j, p11, p12, p21, p22
  int i, j, p11, p12, p21, p22;
//  logical :: i1_ok, i2_ok
  int i1_ok, i2_ok;

//  p11 = p1(1)
  p11 = p1[0];
//  p12 = p1(2)
  p12 = p1[1];
//  p21 = p2(1)
  p21 = p2[0];
//  p22 = p2(2)
  p22 = p2[1];
//  p3(1) = max(p11,p21)
  p3[0] = std::max(p11,p21);
//  p3(2) = max(p12,p22)
  p3[1] = std::max(p12,p22);

//  do i = 0, p3(2)
  for(i=0;i<=p3[1];i++)
  {
//     i1_ok = (i > p12)
    i1_ok = (i > p12);
//     i2_ok = (i > p22)
    i2_ok = (i > p22);
//     do j = 0, p3(1)
    for(j=0;j<=p3[0];j++)
    {
//        if (i2_ok .or. (j > p21)) then
//           u3(j,i) = u1(j,i)
      if (i2_ok || (j > p21))
        u3[i][j] = u1[i][j];
//        else if (i1_ok .or. (j > p11)) then
//           u3(j,i) = -u2(j,i)
      else if (i1_ok || (j > p11))
        u3[i][j] = -u2[i][j];
//        else
//           u3(j,i) = u1(j,i) - u2(j,i)
      else
        u3[i][j] = u1[i][j] - u2[i][j];
//        end if
//     end do
    }
//  end do
  }
}



void ExactIK::poly_mul_sub1(double u1[17], double u2[17], double u3[17], double u4[17], int p1, int p2, int p3, int p4, double u5[17], int *p5)
{
//  implicit none
//  real(dp), dimension(0:16), intent(in) :: u1, u2, u3, u4
//  integer, intent(in) :: p1, p2, p3, p4
//  real(dp), dimension(0:16), intent(out) :: u5
//  integer, intent(out) :: p5
//  real(dp), dimension(0:16) :: d1, d2
  double d1[17], d2[17];
//  integer :: pd1, pd2
  int pd1, pd2;

//  call poly_mul1(u1, u2, p1, p2, d1, pd1)
  poly_mul1(u1, u2, p1, p2, d1, &pd1);
//  call poly_mul1(u3, u4, p3, p4, d2, pd2)
  poly_mul1(u3, u4, p3, p4, d2, &pd2);
//  call poly_sub1(d1, d2, pd1, pd2, u5, p5)
  poly_sub1(d1, d2, pd1, pd2, u5, p5);
}


void ExactIK::poly_mul1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3)
{
//  implicit none
//  real(dp), dimension(0:16), intent(in) :: u1, u2
//  integer, intent(in) :: p1, p2
//  real(dp), dimension(0:16), intent(out) :: u3
//  integer, intent(out) :: p3
//  integer :: i1, i2, i3
  int i, i1, i2, i3;
//  real(dp) :: u1i
  double u1i;

//  p3 = p1 + p2
  *p3 = p1 + p2;
//  u3(:) = 0.0d0
  for(i=0;i<17;i++)
    u3[i] = 0.;

//  do i1 = 0, p1
  for(i1=0;i1<=p1;i1++)
  {
//     u1i = u1(i1)
    u1i = u1[i1];
//     do i2 = 0, p2
    for(i2=0;i2<=p2;i2++)
    {
//        i3 = i1 + i2
      i3 = i1 + i2;
//        u3(i3) = u3(i3) + u1i*u2(i2)
      u3[i3] = u3[i3] + u1i*u2[i2];
//     end do
    }
//  end do
  }
}


void ExactIK::poly_sub1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3)
{
//  implicit none
//  real(dp), dimension(0:16), intent(in) :: u1, u2
//  integer, intent(in) :: p1, p2
//  real(dp), dimension(0:16), intent(out) :: u3
//  integer, intent(out) :: p3
//  integer :: i
  int i;

//  p3 = max(p1, p2)
  *p3 = std::max(p1, p2);

//  do i = 0, p3
  for(i=0;i<=*p3;i++)
  {
//     if (i > p2) then
//        u3(i) = u1(i)
    if (i > p2)
      u3[i] = u1[i];
//     else if (i > p1) then
//        u3(i) = -u2(i)
    else if (i > p1)
      u3[i] = -u2[i];
//     else
//        u3(i) = u1(i) - u2(i)
    else
      u3[i] = u1[i] - u2[i];
//     end if
//  end do
  }

}


void ExactIK::coord_from_poly_roots(int *n_soln, double roots[16], double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3], double r_soln_n[16][3][3], double r_soln_a[16][3][3], double r_soln_c[16][3][3])
{
  double ex[3], ey[3], ez[3], b_a1a2[3], b_a3a2[3], r_tmp[3];
  double p_s[3][3], s1[3][3], s2[3][3], p_t[3][3], t1[3][3], t2[3][3];
  double p_s_c[3][3], s1_s[3][3], s2_s[3][3], p_t_c[3][3], t1_s[3][3], t2_s[3][3];
  double angle, sig1_init, half_tan[3];
  double cos_tau[4], sin_tau[4], cos_sig[3], sin_sig[3], ht, tmp, sig1;
  double r_s[3], r_t[3], r0[3], r_n[3][3], r_a[3][3], r_c[3][3], p[4], Us[3][3];
  int i_soln, i, j;
  double a1c1, c1n2, n2a2, a2c2, c2n3, n3a3, a1a2, a2a3;
  double rr_a1c1[3], rr_c1n2[3], rr_n2a2[3], rr_a2c2[3], rr_c2n3[3], rr_n3a3[3], rr_a1a2[3], rr_a2a3[3];
  double a3a1a2, a2a3a1, n1a1c1, n2a2c2, n3a3c3, a1c1n2a2, a2c2n3a3;
  double tmp_value, ex_tmp[3];
  double tmp_array[3], tmp_array1[3], tmp_array2[3], tmp_array3[3];
  double mat1[3], mat2[3], mat3[3], mat4[3], mat5[3];
  double mat11[3], mat22[3], mat33[3], mat44[3], mat55[3];

  if (*n_soln == 0)
    return;

//  ! Define body frame (ex, ey, ez)
  for(i=0;i<3;i++)
    ex[i] = b_a1a3[i];
  cross(r_a1n1, ex, ez);
//  ez(:) = ez(:)/sqrt(dot_product(ez,ez))
  tmp_value = sqrt(dot_product(ez,ez));
  for(i=0;i<3;i++)
    ez[i] = ez[i]/tmp_value;
  cross(ez, ex, ey);
//  ! vertual bond vectors in the reference plane
  for(i=0;i<3;i++)
  {
    b_a1a2[i] = -cos_alpha[0]*ex[i] + sin_alpha[0]*ey[i];
    b_a3a2[i] = cos_alpha[2]*ex[i] + sin_alpha[2]*ey[i];
  }
//  !! Define cone coordinates for each angle joint.
//  ! (p_s,s1,s2) and (p_t,t1,t2):  Right Orthonormal systems
//  ! residue 1
  for(i=0;i<3;i++)
  {
    p_s[0][i] = -ex[i];
    s1[0][i]  = ez[i];
    s2[0][i]  = ey[i];
    p_t[0][i] = b_a1a2[i];
    t1[0][i]  = ez[i];
    t2[0][i]  = sin_alpha[0]*ex[i] + cos_alpha[0]*ey[i];
  }
//  ! residue 2
  for(i=0;i<3;i++)
  {
    p_s[1][i] = -b_a1a2[i];
    s1[1][i]  = -ez[i];
    s2[1][i]  = t2[0][i];
    p_t[1][i] = -b_a3a2[i];
    t1[1][i]  = -ez[i];
    t2[1][i]  = sin_alpha[2]*ex[i] - cos_alpha[2]*ey[i];
  }
//  ! residue 3
  for(i=0;i<3;i++)
  {
    p_s[2][i] = b_a3a2[i];
    s2[2][i]  = t2[1][i];
    s1[2][i]  = ez[i];
    p_t[2][i] = ex[i];
    t1[2][i] =  ez[i];
    t2[2][i] = -ey[i];
  }
//  ! scale vectors
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      p_s_c[i][j] = p_s[i][j]*cos_xi[i];
      s1_s[i][j]  = s1[i][j]*sin_xi[i];
      s2_s[i][j]  = s2[i][j]*sin_xi[i];
      p_t_c[i][j] = p_t[i][j]*cos_eta[i];
      t1_s[i][j]  = t1[i][j]*sin_eta[i];
      t2_s[i][j]  = t2[i][j]*sin_eta[i];
    }

//  ! initial sig(1)
  for(i=0;i<3;i++)
    r_tmp[i] = (r_a1n1[i]/len_na[0] - p_s_c[0][i])/sin_xi[0];
  calc_bnd_ang(s1[0], r_tmp, &angle);
  sig1_init = sign(angle, dot_product(r_tmp,s2[0]));

//  ! CA
  for(i=0;i<3;i++)
  {
    r_a[0][i] = r_a1[i];
    r_a[1][i] = r_a1[i] + len_aa[1]*b_a1a2[i];
    r_a[2][i] = r_a3[i];
    r0[i] = r_a1[i];
  }

  for(i_soln=0;i_soln<*n_soln;i_soln++)
  {
    half_tan[2] = roots[i_soln];
    half_tan[1] = calc_t2(half_tan[2]);
    half_tan[0] = calc_t1(half_tan[2], half_tan[1]);
    for(i=1;i<=3;i++)
    {
      ht = half_tan[i-1];
      tmp = 1.0e0 + ht*ht;
      cos_tau[i] = (1.0e0 - ht*ht)/tmp;
      sin_tau[i] = 2.0e0*ht/tmp;
    }
    cos_tau[0] = cos_tau[3];
    sin_tau[0] = sin_tau[3];
    for(i=0;i<3;i++)
    {
      cos_sig[i] = cos_delta[i]*cos_tau[i] + sin_delta[i]*sin_tau[i];
      sin_sig[i] = sin_delta[i]*cos_tau[i] - cos_delta[i]*sin_tau[i];
    }
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
      {
        r_s[j] = p_s_c[i][j] + cos_sig[i]*s1_s[i][j] + sin_sig[i]*s2_s[i][j];
        r_t[j] = p_t_c[i][j] + cos_tau[i+1]*t1_s[i][j] + sin_tau[i+1]*t2_s[i][j];
        r_n[i][j] = r_s[j]*len_na[i] + r_a[i][j];
        r_c[i][j] = r_t[j]*len_ac[i] + r_a[i][j];
      }

//     ! rotate back atoms by -(sig(1) - sig1_init) around -ex
    sig1 = atan2(sin_sig[0], cos_sig[0]);
    ex_tmp[0] = -ex[0];
    ex_tmp[1] = -ex[1];
    ex_tmp[2] = -ex[2];
    tmp_value = -(sig1-sig1_init)*0.25;
    quaternion(ex_tmp, tmp_value, p);
    rotation_matrix(p, Us);

    for(i=0;i<3;i++)
    {
      mat11[i] = r_c[0][i]-r0[i];
      mat22[i] = r_n[1][i]-r0[i];
      mat33[i] = r_a[1][i]-r0[i];
      mat44[i] = r_c[1][i]-r0[i];
      mat55[i] = r_n[2][i]-r0[i];
    }
    matmul(Us,mat11,mat1);
    matmul(Us,mat22,mat2);
    matmul(Us,mat33,mat3);
    matmul(Us,mat44,mat4);
    matmul(Us,mat55,mat5);
    for(i=0;i<3;i++)
    {
      r_soln_n[i_soln][0][i] = r_n1[i];
      r_soln_a[i_soln][0][i] = r_a1[i];
      r_soln_c[i_soln][0][i] = mat1[i] + r0[i];
      r_soln_n[i_soln][1][i] = mat2[i] + r0[i];
      r_soln_a[i_soln][1][i] = mat3[i] + r0[i];
      r_soln_c[i_soln][1][i] = mat4[i] + r0[i];
      r_soln_n[i_soln][2][i] = mat5[i] + r0[i];
      r_soln_a[i_soln][2][i] = r_a3[i];
      r_soln_c[i_soln][2][i] = r_c3[i];
    }


    if (print_level > 0)
    {
      for(i=0;i<3;i++)
      {
        rr_a1c1[i] = r_soln_c[i_soln][0][i] - r_soln_a[i_soln][0][i];
        rr_c1n2[i] = r_soln_n[i_soln][1][i] - r_soln_c[i_soln][0][i];
        rr_n2a2[i] = r_soln_a[i_soln][1][i] - r_soln_n[i_soln][1][i];
        rr_a2c2[i] = r_soln_c[i_soln][1][i] - r_soln_a[i_soln][1][i];
        rr_c2n3[i] = r_soln_n[i_soln][2][i] - r_soln_c[i_soln][1][i];
        rr_n3a3[i] = r_soln_a[i_soln][2][i] - r_soln_n[i_soln][2][i];
        rr_a1a2[i] = r_soln_a[i_soln][1][i] - r_soln_a[i_soln][0][i];
        rr_a2a3[i] = r_soln_a[i_soln][2][i] - r_soln_a[i_soln][1][i];
      }

      a1c1 = sqrt(dot_product(rr_a1c1, rr_a1c1));
      c1n2 = sqrt(dot_product(rr_c1n2, rr_c1n2));
      n2a2 = sqrt(dot_product(rr_n2a2, rr_n2a2));
      a2c2 = sqrt(dot_product(rr_a2c2, rr_a2c2));
      c2n3 = sqrt(dot_product(rr_c2n3, rr_c2n3));
      n3a3 = sqrt(dot_product(rr_n3a3, rr_n3a3));
      a1a2 = sqrt(dot_product(rr_a1a2, rr_a1a2));
      a2a3 = sqrt(dot_product(rr_a2a3, rr_a2a3));

      for(i=0;i<3;i++)
        tmp_array[i] = rr_a1a2[i]/a1a2;
      calc_bnd_ang(b_a1a3, tmp_array, &a3a1a2);
      for(i=0;i<3;i++)
        tmp_array[i] = rr_a2a3[i]/a2a3;
      calc_bnd_ang(tmp_array, b_a1a3, &a2a3a1);
      for(i=0;i<3;i++)
        tmp_array[i] = rr_a1c1[i]/a1c1;
      calc_bnd_ang(b_a1n1, tmp_array, &n1a1c1);
      for(i=0;i<3;i++)
        tmp_array[i] = -rr_n3a3[i]/n3a3;
      calc_bnd_ang(b_a3c3, tmp_array, &n3a3c3);
      for(i=0;i<3;i++)
      {
        tmp_array1[i] = rr_a2c2[i]/a2c2;
        tmp_array2[i] = -rr_n2a2[i]/n2a2;
      }
      calc_bnd_ang(tmp_array1, tmp_array2, &n2a2c2);
      for(i=0;i<3;i++)
      {
        tmp_array1[i] = rr_a1c1[i]/a1c1;
        tmp_array2[i] = rr_c1n2[i]/c1n2;
        tmp_array3[i] = rr_n2a2[i]/n2a2;
      }
      calc_dih_ang(tmp_array1, tmp_array2, tmp_array3, &a1c1n2a2);
      for(i=0;i<3;i++)
      {
        tmp_array1[i] = rr_a2c2[i]/a2c2;
        tmp_array2[i] = rr_c2n3[i]/c2n3;
        tmp_array3[i] = rr_n3a3[i]/n3a3;
      }
      calc_dih_ang(tmp_array1, tmp_array2, tmp_array3, &a2c2n3a3);
    }
  }
}


double ExactIK::calc_t2(double t0)
{
//  implicit none
//  real(dp), intent(in) :: t0
//  real(dp) :: calc_t2
  double tmp_value;
//  real(dp) :: B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3
  double B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3;
//  real(dp) :: K0, K1, K2, K3, t0_2, t0_3, t0_4
  double K0, K1, K2, K3, t0_2, t0_3, t0_4;

//  t0_2 = t0*t0
  t0_2 = t0*t0;
//  t0_3 = t0_2*t0
  t0_3 = t0_2*t0;
//  t0_4 = t0_3*t0
  t0_4 = t0_3*t0;

//  A0 = Q(0,0) + Q(1,0)*t0 + Q(2,0)*t0_2 + Q(3,0)*t0_3 + Q(4,0)*t0_4
  A0 = Q[0][0] + Q[0][1]*t0 + Q[0][2]*t0_2 + Q[0][3]*t0_3 + Q[0][4]*t0_4;
//  A1 = Q(0,1) + Q(1,1)*t0 + Q(2,1)*t0_2 + Q(3,1)*t0_3 + Q(4,1)*t0_4
  A1 = Q[1][0] + Q[1][1]*t0 + Q[1][2]*t0_2 + Q[1][3]*t0_3 + Q[1][4]*t0_4;
//  A2 = Q(0,2) + Q(1,2)*t0 + Q(2,2)*t0_2 + Q(3,2)*t0_3 + Q(4,2)*t0_4
  A2 = Q[2][0] + Q[2][1]*t0 + Q[2][2]*t0_2 + Q[2][3]*t0_3 + Q[2][4]*t0_4;
//  A3 = Q(0,3) + Q(1,3)*t0 + Q(2,3)*t0_2 + Q(3,3)*t0_3 + Q(4,3)*t0_4
  A3 = Q[3][0] + Q[3][1]*t0 + Q[3][2]*t0_2 + Q[3][3]*t0_3 + Q[3][4]*t0_4;
//  A4 = Q(0,4) + Q(1,4)*t0 + Q(2,4)*t0_2 + Q(3,4)*t0_3 + Q(4,4)*t0_4
  A4 = Q[4][0] + Q[4][1]*t0 + Q[4][2]*t0_2 + Q[4][3]*t0_3 + Q[4][4]*t0_4;

//  B0 = R(0,0) + R(1,0)*t0 + R(2,0)*t0_2
  B0 = R[0][0] + R[0][1]*t0 + R[0][2]*t0_2;
//  B1 = R(0,1) + R(1,1)*t0 + R(2,1)*t0_2
  B1 = R[1][0] + R[1][1]*t0 + R[1][2]*t0_2;
//  B2 = R(0,2) + R(1,2)*t0 + R(2,2)*t0_2
  B2 = R[2][0] + R[2][1]*t0 + R[2][2]*t0_2;

//  B2_2 = B2*B2
  B2_2 = B2*B2;
//  B2_3 = B2_2*B2
  B2_3 = B2_2*B2;

//  K0 = A2*B2 - A4*B0
  K0 = A2*B2 - A4*B0;
//  K1 = A3*B2 - A4*B1
  K1 = A3*B2 - A4*B1;
//  K2 = A1*B2_2 - K1*B0
  K2 = A1*B2_2 - K1*B0;
//  K3 = K0*B2 - K1*B1
  K3 = K0*B2 - K1*B1;

//  calc_t2 = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1)
  tmp_value = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1);

  return tmp_value;
}


double ExactIK::calc_t1(double t0, double t2)
{
//  implicit none
//  real(dp), intent(in) :: t0, t2
//  real(dp) :: calc_t1
  double tmp_value;
//  real(dp) :: U11, U12, U13, U31, U32, U33
  double U11, U12, U13, U31, U32, U33;
//  real(dp) :: t0_2, t2_2
  double t0_2, t2_2;

//  t0_2 = t0*t0
  t0_2 = t0*t0;
//  t2_2 = t2*t2
  t2_2 = t2*t2;

//  U11 = C0(0,1) + C0(1,1)*t0 + C0(2,1)*t0_2
  U11 = C0[0][0] + C0[0][1]*t0 + C0[0][2]*t0_2;
//  U12 = C1(0,1) + C1(1,1)*t0 + C1(2,1)*t0_2
  U12 = C1[0][0] + C1[0][1]*t0 + C1[0][2]*t0_2;
//  U13 = C2(0,1) + C2(1,1)*t0 + C2(2,1)*t0_2
  U13 = C2[0][0] + C2[0][1]*t0 + C2[0][2]*t0_2;
//  U31 = C0(0,2) + C0(1,2)*t2 + C0(2,2)*t2_2
  U31 = C0[1][0] + C0[1][1]*t2 + C0[1][2]*t2_2;
//  U32 = C1(0,2) + C1(1,2)*t2 + C1(2,2)*t2_2
  U32 = C1[1][0] + C1[1][1]*t2 + C1[1][2]*t2_2;
//  U33 = C2(0,2) + C2(1,2)*t2 + C2(2,2)*t2_2
  U33 = C2[1][0] + C2[1][1]*t2 + C2[1][2]*t2_2;

//  calc_t1 = (U31*U13-U11*U33)/(U12*U33-U13*U32)
  tmp_value = (U31*U13-U11*U33)/(U12*U33-U13*U32);

  return tmp_value;
}

void ExactIK::calc_dih_ang(double r1[3], double r2[3], double r3[3], double *angle)
{
//!-----------------------------------------------------------------------
//! r1=Rab, r2=Rbc, r3=Rcd : angle between planes abc and bcd
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), intent(in) :: r1(3), r2(3), r3(3)
//  real(dp), intent(out) :: angle
//  real(dp), dimension(3) :: p, q, s
  double p[3], q[3], s[3];
//  real(dp) :: arg
  double arg;

//  call cross(r1, r2, p)
  cross(r1, r2, p);
//  call cross(r2, r3, q)
  cross(r2, r3, q);
//  call cross(r3, r1, s)
  cross(r3, r1, s);
//  arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q))
  arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q));
//  arg = sign(min(abs(arg),1.0d0),arg) ! to be sure abs(arg)<=1
  arg = sign(std::min(fabs(arg),1.0e0),arg);
//  angle = sign(acos(arg), dot_product(s,r2))
  *angle = sign(acos(arg), dot_product(s,r2));

  return;
}

void ExactIK::calc_bnd_ang(double r1[3], double r2[3], double *angle)
{
//!-----------------------------------------------------------------------
//! assume that each vector is normalized
//! r1=Rba, r2=Rbc: angle between r1 and r2
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), intent(in) :: r1(3), r2(3)
//  real(dp), intent(out) :: angle
//  real(dp) :: arg
  double arg;

//  arg = dot_product(r1, r2)
  arg = dot_product(r1, r2);
//  arg = sign(min(abs(arg),1.0d0),arg) ! to be sure abs(arg)<=1
  arg = sign(std::min(fabs(arg),1.0e0),arg);
//  angle = acos(arg)
  *angle = acos(arg);

  return;
}

void ExactIK::cross(double p[3], double q[3], double s[3])
{
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), dimension(:), intent(in) :: p, q
//  real(dp), dimension(:), intent(out) :: s

//  s(1) = p(2)*q(3) - p(3)*q(2)
  s[0] = p[1]*q[2] - p[2]*q[1];
//  s(2) = p(3)*q(1) - p(1)*q(3)
  s[1] = p[2]*q[0] - p[0]*q[2];
//  s(3) = p(1)*q(2) - p(2)*q(1)
  s[2] = p[0]*q[1] - p[1]*q[0];

  return;
}


// ----------- Sturmian functions --------------

void ExactIK::initialize_sturm(double *tol_secant, int *max_iter_sturm, int *max_iter_secant)
{
  RELERROR = *tol_secant;
  MAXIT = *max_iter_sturm;
  MAX_ITER_SECANT = *max_iter_secant;
}

void ExactIK::solve_sturm(int *p_order, int *n_root, double *poly_coeffs, double *roots)
{
  poly sseq[MAX_ORDER*2];
  double min, max;
  int order, i, j, nroots, nchanges, np, atmin, atmax;

  order = *p_order;

  for (i = order; i >= 0; i--)
  {
    sseq[0].coef[i] = poly_coeffs[i];
  }

  if (PRINT_LEVEL > 0)
  {
    for (i = order; i >= 0; i--)
    {
      printf("coefficients in Sturm solver\n");
      printf("%d %lf\n", i, sseq[0].coef[i]);
    }
  }

  /*
   * build the Sturm sequence
   */
  np = buildsturm(order, sseq);

  if (PRINT_LEVEL > 0)
  {
    printf("Sturm sequence for:\n");
    for (i = order; i >= 0; i--)
      printf("%lf ", sseq[0].coef[i]);
    printf("\n\n");
    for (i = 0; i <= np; i++)
    {
      for (j = sseq[i].ord; j >= 0; j--)
        printf("%lf ", sseq[i].coef[j]);
      printf("\n");
    }
    printf("\n");
  }

  /*
   * get the number of real roots
   */

  nroots = numroots(np, sseq, &atmin, &atmax);


  if (nroots == 0)
  {
    // printf("solve: no real roots\n");
    *n_root = nroots;
    return;
  }

  if (PRINT_LEVEL > 0)
    printf("Number of real roots: %d\n", nroots);

  /*
   * calculate the bracket that the roots live in
   */
  min = -1.0;
  nchanges = numchanges(np, sseq, min);

  for (i = 0; nchanges != atmin && i != MAXPOW; i++) {
    min *= 10.0;
    nchanges = numchanges(np, sseq, min);
  }

  if (nchanges != atmin) {
    printf("solve: unable to bracket all negative roots\n");
    atmin = nchanges;
  }

  max = 1.0;
  nchanges = numchanges(np, sseq, max);
  for (i = 0; nchanges != atmax && i != MAXPOW; i++) {
    max *= 10.0;
    nchanges = numchanges(np, sseq, max);
  }

  if (nchanges != atmax) {
    printf("solve: unable to bracket all positive roots\n");
    atmax = nchanges;
  }

  nroots = atmin - atmax;

  /*
   * perform the bisection.
   */

  sbisect(np, sseq, min, max, atmin, atmax, roots);

  *n_root = nroots;

  /*
   * write out the roots...
   */
  if (PRINT_LEVEL > 0)
  {
    if (nroots == 1) {
      printf("\n1 distinct real m_root at x = %f\n", roots[0]);
    } else {
      printf("\n%d distinct real roots for x: \n", nroots);

      for (i = 0; i != nroots; i++)
      {
        printf("%f\n", roots[i]);
      }
    }
  }

  return;
}

/*
 * modp
 *
 *	calculates the modulus of u(x) / v(x) leaving it in r, it
 *  returns 0 if r(x) is a constant.
 *  note: this function assumes the leading coefficient of v
 *	is 1 or -1
 */

double ExactIK::hyper_tan(double a, double x)
{
  double exp_x1, exp_x2, ax;

  ax = a*x;
  if (ax > 100.0)
  {
    return(1.0);
  }
  else if (ax < -100.0)
  {
    return(-1.0);
  }
  else
  {
    exp_x1 = exp(ax);
    exp_x2 = exp(-ax);
    return (exp_x1 - exp_x2)/(exp_x1 + exp_x2);
  }
}

int ExactIK::modp(poly *u, poly *v, poly *r)
//	poly *u, *v, *z;
{
  int		k, j;
  double	*nr, *end, *uc;

  nr = r->coef;
  end = &u->coef[u->ord];

  uc = u->coef;
  while (uc <= end)
    *nr++ = *uc++;

  if (v->coef[v->ord] < 0.0) {


    for (k = u->ord - v->ord - 1; k >= 0; k -= 2)
      r->coef[k] = -r->coef[k];

    for (k = u->ord - v->ord; k >= 0; k--)
      for (j = v->ord + k - 1; j >= k; j--)
        r->coef[j] = -r->coef[j] - r->coef[v->ord + k]
                                   * v->coef[j - k];
  } else {
    for (k = u->ord - v->ord; k >= 0; k--)
      for (j = v->ord + k - 1; j >= k; j--)
        r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
  }

  k = v->ord - 1;
  while (k >= 0 && fabs(r->coef[k]) < SMALL_ENOUGH) {
    r->coef[k] = 0.0;
    k--;
  }

  r->ord = (k < 0) ? 0 : k;

  return(r->ord);
}

/*
 * buildsturm
 *
 *	build up a sturm sequence for a polynomial in smat, returning
 * the number of polynomials in the sequence
 */
int ExactIK::buildsturm(int ord, poly *sseq)
{
  int		i;
  double	f, *fp, *fc;
  poly	*sp;

  sseq[0].ord = ord;
  sseq[1].ord = ord - 1;

  /*
   * calculate the derivative and normalise the leading
   * coefficient.
   */
  f = fabs(sseq[0].coef[ord]*ord);

  fp = sseq[1].coef;
  fc = sseq[0].coef + 1;

  for (i = 1; i <= ord; i++)
  {
    *fp++ = *fc++ * i / f;
  }


  /*
   * construct the rest of the Sturm sequence
   */
  //	for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++) {
  for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++) {

    /*
     * reverse the sign and normalise
     */

    f = -fabs(sp->coef[sp->ord]);
    for (fp = &sp->coef[sp->ord]; fp >= sp->coef; fp--)
      *fp /= f;
  }

  sp->coef[0] = -sp->coef[0];	/* reverse the sign */

  return(sp - sseq);
}

/*
 * numroots
 *
 *	return the number of distinct real roots of the polynomial
 * described in sseq.
 */
int ExactIK::numroots(int np, poly *sseq, int *atneg, int *atpos)
//		int		np;
//		poly	*sseq;
//		int		*atneg, *atpos;
{
  int		atposinf, atneginf;
  poly	*s;
  double	f, lf;

  atposinf = atneginf = 0;


  /*
   * changes at positive infinity
   */
  lf = sseq[0].coef[sseq[0].ord];

  for (s = sseq + 1; s <= sseq + np; s++) {
    f = s->coef[s->ord];
    if (lf == 0.0 || lf * f < 0)
      atposinf++;
    lf = f;
  }

  /*
   * changes at negative infinity
   */
  if (sseq[0].ord & 1)
    lf = -sseq[0].coef[sseq[0].ord];
  else
    lf = sseq[0].coef[sseq[0].ord];

  for (s = sseq + 1; s <= sseq + np; s++) {
    if (s->ord & 1)
      f = -s->coef[s->ord];
    else
      f = s->coef[s->ord];
    if (lf == 0.0 || lf * f < 0)
      atneginf++;
    lf = f;
  }

  *atneg = atneginf;
  *atpos = atposinf;
  //	printf("atneginf, atposinf = %d %d\n", atneginf, atposinf);
  return(atneginf - atposinf);
}

/*
 * numchanges
 *
 *	return the number of sign changes in the Sturm sequence in
 * sseq at the value a.
 */
int ExactIK::numchanges(int np, poly *sseq, double a)
//	int		np;
//	poly	*sseq;
//	double	a;

{
  int		changes;
  double	f, lf;
  poly	*s;

  changes = 0;

  lf = evalpoly(sseq[0].ord, sseq[0].coef, a);

  for (s = sseq + 1; s <= sseq + np; s++) {
    f = evalpoly(s->ord, s->coef, a);
    if (lf == 0.0 || lf * f < 0)
      changes++;
    lf = f;
//			printf("lf %lf %d \n", f, changes);
  }

  //	printf("%d \n", changes);
  return(changes);
}

/*
 * sbisect
 *
 *	uses a bisection based on the sturm sequence for the polynomial
 * described in sseq to isolate intervals in which roots occur,
 * the roots are returned in the roots array in order of magnitude.
 */
void ExactIK::sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots)
//	int	np;
//	poly	*sseq;
//	double	min, max;
//	int		atmin, atmax;
//	double	*roots;
{
  double	mid=0.;
  int		n1 = 0, n2 = 0, its, atmid, nroot;

  if ((nroot = atmin - atmax) == 1)
  {

    /*
     * first try a less expensive technique.
     */
    //      printf("min max %lf, %lf \n", min, max);
    if (modrf(sseq->ord, sseq->coef, min, max, &roots[0]))
      return;

    //      printf("try hard way\n");
    /*
     * if we get here we have to evaluate the m_root the hard
     * way by using the Sturm sequence.
     */
    for (its = 0; its < MAXIT; its++)
    {
      mid = (min + max) / 2;

      atmid = numchanges(np, sseq, mid);

      if (fabs(mid) > RELERROR)
      {
        if (fabs((max - min) / mid) < RELERROR)
        {
          roots[0] = mid;
          return;
        }
      }
      else if (fabs(max - min) < RELERROR)
      {
        roots[0] = mid;
        return;
      }

      if ((atmin - atmid) == 0)
        min = mid;
      else
        max = mid;
    }

    if (its == MAXIT)
    {
      fprintf(stderr, "sbisect: overflow min %f max %f\
					diff %e nroot %d n1 %d n2 %d\n",
              min, max, max - min, nroot, n1, n2);
      roots[0] = mid;
    }
    return;
  }

  /*
   * more than one m_root in the interval, we have to bisect...
   */

  for (its = 0; its < MAXIT; its++)
  {

    mid = (min + max) / 2;

    atmid = numchanges(np, sseq, mid);

    n1 = atmin - atmid;
    n2 = atmid - atmax;

    if (n1 != 0 && n2 != 0)
    {
      sbisect(np, sseq, min, mid, atmin, atmid, roots);
      sbisect(np, sseq, mid, max, atmid, atmax, &roots[n1]);
      break;
    }

    if (n1 == 0)
      min = mid;
    else
      max = mid;
  }

  if (its == MAXIT)
  {
    for (n1 = atmax; n1 < atmin; n1++)
      roots[n1 - atmax] = mid;
  }
}

/*
 * evalpoly
 *
 *	evaluate polynomial defined in coef returning its value.
 */
double ExactIK::evalpoly(int ord, double *coef, double x)
//	int		ord;
//	double	*coef, x;
{
  double	*fp, f;

  fp = &coef[ord];
  f = *fp;

  for (fp--; fp >= coef; fp--)
    f = x * f + *fp;

  return(f);
}


/*
 * modrf
 *
 *	uses the modified regula-falsi method to evaluate the m_root
 * in interval [a,b] of the polynomial described in coef. The
 * m_root is returned is returned in *val. The routine returns zero
 * if it can't converge.
 */
int ExactIK::modrf(int ord, double *coef, double	a, double b, double *val)
{
  int its;
  double fa, fb, x, fx, lfx;
  double *fp, *scoef, *ecoef;

  scoef = coef;
  ecoef = &coef[ord];

  fb = fa = *ecoef;
  for (fp = ecoef - 1; fp >= scoef; fp--) {
    fa = a * fa + *fp;
    fb = b * fb + *fp;
  }

  /*
   * if there is no sign difference the method won't work
   */
  if (fa * fb > 0.0)
    return(0);


  lfx = fa;

  for (its = 0; its < MAX_ITER_SECANT; its++)
  {
    x = (fb * a - fa * b) / (fb - fa);

    // constrain that x stays in the bounds
    if (x < a || x > b)
      x = 0.5 * (a+b);

    fx = *ecoef;
    for (fp = ecoef - 1; fp >= scoef; fp--)
      fx = x * fx + *fp;

    if (fabs(x) > RELERROR)
    {
      if (fabs(fx / x) < RELERROR)
      {
        *val = x;
        //	      printf(" x, fx %lf %lf\n", x, fx);
        return(1);
      }
    }
    else if (fabs(fx) < RELERROR)
    {
      *val = x;
      //	  printf(" x, fx %lf %lf\n", x, fx);
      return(1);
    }

    if ((fa * fx) < 0)
    {
      b = x;
      fb = fx;
      if ((lfx * fx) > 0)
        fa /= 2;
    }
    else
    {
      a = x;
      fa = fx;
      if ((lfx * fx) > 0)
        fb /= 2;
    }

    lfx = fx;
  }
  return(0);
}


