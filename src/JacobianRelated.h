/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

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
#ifndef JACOBIANRELATED_H
#define JACOBIANRELATED_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

///These values have to be chosen according to the numerical analysis
static const double SINGVAL_TOL = 0.000000000001; // only generic 10^-12
//static const double SINGVAL_TOL = 0.000001; // include non-generic 10^-10
static const double RIGID_TOL = 0.0000000001; //depends on molecule, but 10^-10 seems a good fit!
//static const double RIGID_TOL = 0.000001; //depends on molecule ToDo: Adapt to obtain on-the-fly parameter selection

class NullSpaceRet
{
public:
//	static gsl_matrix* Ut;//The Ut from the SVD
	static gsl_matrix* V;//The V from the SVD
	static gsl_vector* singularValues;
	int nullspaceSize;
	int m, n;

//	gsl_matrix* P;//Null-space projection matrix
	gsl_matrix* m_nullspaceBasis; //basis of the nullspace
//	gsl_matrix* Jinv;//Jacobian pseudo-inverse

    gsl_vector* m_rigidAngles;
    gsl_vector* m_rigidHBonds;

    int m_numCoordinated;
    int m_numRigid;
    int m_numRigidHBonds;

	NullSpaceRet (int input_m, int input_n);
	~NullSpaceRet ();
	void print();

	// A function to make the projection operation more efficient;
	void NullSpacePostCompute();

	/// A function that identifies rigidified dihedral angles
	void RigidityAnalysis(gsl_matrix* HBondJacobian);
	
	//A function to project a vector on the nullspace
	void ProjectOnNullSpace (gsl_vector *to_project, gsl_vector *after_project);

	// To manipulate owner configations of NullSpaceRet instances
	void numOwners(unsigned int numOwners);
	unsigned int numOwners() const;
private:
	unsigned int numOwners_;
};

void nr_svd(gsl_matrix* a_g, int m, int n, double w[], gsl_matrix* v_out);
void ComputeNullSpace(gsl_matrix* Jac, NullSpaceRet* Ret);
//void ProjectOnNullSpace (NullSpaceRet* nullspace, gsl_vector *e, gsl_vector *to_project, gsl_vector *after_project);
void TorsionUpdate(NullSpaceRet* nullspace, gsl_vector *e, double lambda, gsl_vector *dTheta);

#endif

