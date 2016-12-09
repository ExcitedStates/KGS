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

#include "VDWDirection.h"

#include <vector>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#include "core/Molecule.h"
#include "core/Grid.h"
#include "core/Chain.h"

void VDWDirection::computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret)
{
  gsl_vector_set_all(ret, 0.0);

  Molecule * protein = conf->updatedMolecule();
  gsl_matrix* atomJacobian1 = gsl_matrix_calloc(3,protein->totalDofNum());
  gsl_matrix* atomJacobian2 = gsl_matrix_calloc(3,protein->totalDofNum());
  gsl_vector* p12 = gsl_vector_calloc(3);
  gsl_vector* p_temp = gsl_vector_calloc(protein->totalDofNum());

  for (auto const& atom1: protein->getAtoms()) {
    std::vector<Atom*> neighbors = protein->getGrid()->getNeighboringAtomsVDW(atom1,      //atom
                                                                              true,       //neighborWithLargerId
                                                                              true,       //noCovBondNeighbor
                                                                              true,       //noSecondCovBondNeighbor
                                                                              true,       //noHbondNeighbor
                                                                              VDW_R_MAX); //radius

    computeAtomJacobian(atom1,atomJacobian1);

    for(auto const& atom2: neighbors){
      double r_12 = atom1->distanceTo(atom2);
      //Dimitars version
      //double atomContribution = (-12)*VDW_SIGMA*(pow(VDW_R0,6)*pow(r_12,-8)-pow(VDW_R0,12)*pow(r_12,-14));
      double vdw_r12 =      atom1->getRadius()  + atom2->getRadius()  ; // from CHARMM: arithmetic mean
      double eps_r12 = sqrt(atom1->getEpsilon() * atom2->getEpsilon()); // from CHARMM: geometric mean
      double ratio = vdw_r12/r_12;
      //double atomContribution = 4 * eps_r12 * (pow(ratio,12)-2*pow(ratio,6));
      double atomContribution = 12 * 4 * eps_r12 * (pow(ratio,6)-pow(ratio,12))/r_12;

      computeAtomJacobian(atom2,atomJacobian2);
      gsl_matrix_sub(atomJacobian1,atomJacobian2); // atomJacobian1 = atomJacobian1 - atomJacobian2
      Math3D::Vector3 p12_v3 = atom1->m_position - atom2->m_position;//TODO: Subtract  directly into gsl_vector
      Coordinate::copyToGslVector(p12_v3, p12);
      gsl_blas_dgemv(CblasTrans,1,atomJacobian1,p12,0,p_temp);
      //std::cout<<"VDWDirection::computeGradient - pair-gradient norm: "<<gsl_blas_dnrm2(p_temp)<<std::endl;
      gsl_vector_scale(p_temp, atomContribution/gsl_blas_dnrm2(p_temp));
      //std::cout<<"VDWDirection::computeGradient - after scaling: "<<gsl_blas_dnrm2(p_temp)<<std::endl;

      gsl_vector_add(ret,p_temp);
    }

  }

  gsl_vector_scale(ret,0.001);
  //std::cout<<"VDWDirection::computeGradient - total gradient norm: "<<gsl_blas_dnrm2(ret)<<std::endl;

  //TODO: Improve to minimize reallocations
  gsl_matrix_free(atomJacobian1);
  gsl_matrix_free(atomJacobian2);
  gsl_vector_free(p12);
  gsl_vector_free(p_temp);


}
void VDWDirection::computeAtomJacobian (Atom* atom, gsl_matrix* jacobian) {
  Molecule * protein = atom->getResidue()->getChain()->getMolecule();
  //KinVertex *vertex = protein->getRigidbodyGraphVertex(atom);
  KinVertex *vertex = atom->getRigidbody()->getVertex();
  while (vertex->m_parent!=nullptr) {
    KinEdge* edge = vertex->m_parent->findEdge(vertex);
    int dof_id = edge->getDOF()->getIndex();
//    Bond * bond_ptr = edge->getBond();
//    Coordinate bp1 = bond_ptr->Atom1->m_position;
//    Coordinate bp2 = bond_ptr->Atom2->m_position;
//    Math3D::Vector3 jacobian_entry = ComputeJacobianEntry(bp1,bp2,atom->m_position);
    Math3D::Vector3 jacobian_entry = edge->getDOF()->getDerivative(atom->m_position);
    gsl_matrix_set(jacobian,0,dof_id,jacobian_entry.x);
    gsl_matrix_set(jacobian,1,dof_id,jacobian_entry.y);
    gsl_matrix_set(jacobian,2,dof_id,jacobian_entry.z);
    vertex = vertex->m_parent;
  }

  //Todo: This should be optimizable using the sorted vertices and the Abé implementation of the MSD gradient
}
