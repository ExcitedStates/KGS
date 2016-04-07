//
// Created by Rasmus Fonseca on 29/03/16.
//

#ifndef KGS_EXACTIK_H
#define KGS_EXACTIK_H

#include <vector>

#include <core/Configuration.h>

// Used in Sturmian solver
#define PRINT_LEVEL 0
#define	MAX_ORDER 16
#define MAXPOW 32
#define	SMALL_ENOUGH 1.0e-18


class ExactIK {
 public:
  std::vector<Configuration*> rebuildLoop( const Residue* res1, const Residue* res2, const Residue* res3 );

  bool validRebuildLoop(const Residue* res1, const Residue* res2, const Residue* res3 ) const;
 private:

  void initializeIKParams(
      const Residue* res1,
      const Residue* res2,
      const Residue* res3,
      double blen[6],
      double bang[7],
      double tang[2]
  );

  int max_soln = 16;
  int deg_pol = 16;
  int print_level = 1;
  double len0[6], b_ang0[7], t_ang0[2];
  double aa13_min_sqr, aa13_max_sqr;
  double delta[4], xi[3], eta[3], alpha[3], theta[3];
  double cos_alpha[3], sin_alpha[3], cos_theta[3], sin_theta[3];
  double cos_delta[4], sin_delta[4];
  double cos_xi[3], cos_eta[3], sin_xi[3], sin_eta[3];
  double r_a1a3[3], r_a1n1[3], r_a3c3[3];
  double b_a1a3[3], b_a1n1[3], b_a3c3[3];
  double len_na[3], len_ac[3], len_aa[3];
  double C0[3][3], C1[3][3], C2[3][3];
  double Q[5][17], R[3][17];

  void solve_3pep_poly(double r_n1[3],
                       double r_a1[3],
                       double r_a3[3],
                       double r_c3[3],
                       double r_soln_n[16][3][3],
                       double r_soln_a[16][3][3],
                       double r_soln_c[16][3][3],
                       int *n_soln);
  void initialize_loop_closure(double b_len[6], double b_ang[7], double t_ang[2]);

  double dot_product(double va[3], double vb[3]);
  void matmul(double ma[3][3], double mb[3], double mc[3]);
  double sign(double a, double b);

  /** Calculate quaternion, given rotation axis and angle. */
  void quaternion(double axis[3], double quarter_ang, double p[4]) const;

  /** Constructs rotation matrix U from quaternion q.*/
  void rotation_matrix(double q[4], double U[3][3]) const;

  void get_input_angles(int *n_soln, double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3]);
  void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln, char cone_type[2]);
  void get_poly_coeff(double poly_coeff[16+1]);
  void poly_mul_sub2(double u1[5][5], double u2[5][5], double u3[5][5], double u4[5][5], int p1[2], int p2[2], int p3[2], int p4[2], double u5[5][5], int p5[2]);
  void poly_mul2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2]);
  void poly_sub2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2]);
  void poly_mul_sub1(double u1[17], double u2[17], double u3[17], double u4[17], int p1, int p2, int p3, int p4, double u5[17], int *p5);
  void poly_mul1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3);
  void poly_sub1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3);
  void coord_from_poly_roots(int *n_soln, double roots[16], double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3], double r_soln_n[16][3][3], double r_soln_a[16][3][3], double r_soln_c[16][3][3]);
  double calc_t2(double t0);
  double calc_t1(double t0, double t2);
  void calc_dih_ang(double r1[3], double r2[3], double r3[3], double *angle);
  void calc_bnd_ang(double r1[3], double r2[3], double *angle);
  void cross(double p[3], double q[3], double s[3]);

  // -- Sturmian calculations --

  typedef  	struct	p {
    int	ord;
    double	coef[MAX_ORDER+1];
  } poly;

  double RELERROR;
  int MAXIT;
  int MAX_ITER_SECANT;
  void initialize_sturm(double *tol_secant, int *max_iter_sturm, int *max_iter_secant);
  void solve_sturm(int *p_order, int *n_root, double *poly_coeffs, double *roots);
  double hyper_tan(double a, double x);
  int modp(poly *u, poly *v, poly *r);
  int buildsturm(int ord, poly *sseq);
  int numroots(int np, poly *sseq, int *atneg, int *atpos);
  int numchanges(int np, poly *sseq, double a);
  void sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots);
  double evalpoly(int ord, double *coef, double x);
  int modrf(int ord, double *coef, double	a, double b, double *val);


};


#endif //KGS_EXACTIK_H
