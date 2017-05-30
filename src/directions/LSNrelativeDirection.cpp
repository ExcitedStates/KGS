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
#include "LSNrelativeDirection.h"

#include <cassert>
#include <cmath>

#include "math/gsl_helpers.h"


using namespace std;



LSNrelativeDirection::LSNrelativeDirection(
    Selection& atomsMoving,
    const std::vector< std::tuple<Atom*, Atom*, double> >& goal_distances
) :
    m_atomsMovingSelection(atomsMoving),
    goal_distances(goal_distances)
{
//  int numMovingAtoms = ExploreOptions::getOptions()->getAtomsMoving()->size();
//  m_TargetJacobian = gsl_matrix_calloc( numMovingAtoms*3, 1 );
}


void LSNrelativeDirection::computeGradient(Configuration* conf, Configuration* conf2, gsl_vector* ret)
{
  Molecule* protein = conf->getMolecule();

  vector<Atom*> atomList = m_atomsMovingSelection.getSelectedAtoms(protein);
  //current_q->ComputeCycleJacobianAndNullSpace();
  gsl_matrix* N = conf->getNullspace()->getBasis();
  int dof=protein->totalDofNum();
  gsl_matrix* fullN = gsl_matrix_calloc(dof,N->size2+dof-N->size1);
  gsl_matrix_set_zero(fullN);
  int fulldof=0;
  for (auto const& edge: protein->m_spanningTree->Edges) {
    int dof_id = edge->getDOF()->getIndex();
    int cycle_dof_id = edge->getDOF()->getCycleIndex();
    if ( cycle_dof_id!=-1 ) {
      for( int i=0; i!=N->size2; i++){
        gsl_matrix_set(fullN, dof_id,i, gsl_matrix_get(N,cycle_dof_id,i));
      }
    }
    else{
      gsl_matrix_set(fullN,dof_id, N->size2+fulldof, 1);
      fulldof++;
    }
  }

//  const vector<Atom*>* atomList = ExploreOptions::getOptions()->getAtomsMoving();
  gsl_matrix* TargetJacobian = gsl_matrix_calloc(goal_distances.size()*6,fullN->size1);
  gsl_matrix* TargetDirection = gsl_matrix_calloc(goal_distances.size()*6,fullN->size1);
//  gsl_matrix_set_zero(m_TargetJacobian); //Already done by calloc
//  gsl_matrix_set_zero(m_TargetPosition);
  fillmatrices(conf, TargetJacobian, TargetDirection);

  gsl_matrix* bestmove = determineBestMove(fullN,TargetJacobian,TargetDirection);

  for (int i=0;i<fullN->size1;i++){
    double val = gsl_matrix_get(bestmove, i, 0);
    assert( !std::isnan(val) );
    gsl_vector_set(ret,i,val); //TODO: Speed up with memcpy
  }
  gsl_matrix_free(bestmove);
  gsl_matrix_free(fullN);
  gsl_matrix_free(TargetJacobian);
  gsl_matrix_free(TargetDirection);

}


void LSNrelativeDirection::fillmatrices(Configuration* current_q,
                                        gsl_matrix* targetJacobian,
                                        gsl_matrix* targetPosition)
{

  int i=0;
  Molecule* protein=current_q->getMolecule();

//  const vector<Atom*>* atomList = ExploreOptions::getOptions()->getAtomsMoving();
//  const vector<Atom*>& atomList1 = m_atomsMovingSelection.getSelectedAtoms(protein);
//  const vector<Atom*>& atomList2 = m_atomsMovingSelection.getSelectedAtoms(target);
  const vector<Atom*>& atomList = m_atomsMovingSelection.getSelectedAtoms(protein);

//  if(atomList1.size()!=atomList2.size()){
//    cerr<<"LSNullspaceDirection::fillmatrices - Molecules contain different number of atoms"<<endl;
//    exit(-1);
//  }
  int nb_couple = goal_distances.size();
  //for(auto const& atom : atomList){//only use to extract res IDs and names
  gsl_vector* u = gsl_vector_alloc(3);//contains the m_direction between each couple of atoms.
  for (int i=0;i<nb_couple;i++){
    Atom* atom1 = get<0>(goal_distances[i]);
    Atom* atom2 = get<1>(goal_distances[i]);
    double dist_goal = get<2>(goal_distances[i]);
    /*string name = atom->getName();
    int resId = atom->getResidue()->getId();
    string chainName = atom->getResidue()->getChain()->getName();
//    Atom* atom = protein->getAtom(chainName,resId,name);
    Coordinate p = atom->m_position;
    Atom* aTarget = target->getAtom(chainName,resId,name);
    if(aTarget == NULL ) continue;//skip the non-existing atom*/
    gsl_vector_set(u,0,atom1->m_position.x - atom2->m_position.x);
    gsl_vector_set(u,1,atom1->m_position.y - atom2->m_position.y);
    gsl_vector_set(u,2,atom1->m_position.z - atom2->m_position.z);
    double d = gsl_vector_length(u);
    //cout<<"Distance: "<<d<<" "<<dist_goal<<endl;
    gsl_matrix_set(targetPosition,2*i*3+0,0, (gsl_vector_get(u,0)/d)*(dist_goal-d)/2);
    //cout<<"Direction : "<<(gsl_vector_get(u,0)/d)*(goal_distances[i]-d)/2<<endl;
    gsl_matrix_set(targetPosition,2*i*3+1,0, (gsl_vector_get(u,1)/d)*(dist_goal-d)/2);
    //cout<<"Direction : "<<(gsl_vector_get(u,1)/d)*(goal_distances[i]-d)/2<<endl;
    gsl_matrix_set(targetPosition,2*i*3+2,0, (gsl_vector_get(u,2)/d)*(dist_goal-d)/2);
    // cout<<"Direction : "<<(gsl_vector_get(u,2)/d)*(goal_distances[i]-d)/2<<endl;
    gsl_matrix_set(targetPosition,2*i*3+3,0, -(gsl_vector_get(u,0)/d)*(dist_goal-d)/2);
    // cout<<"Direction : "<<-(gsl_vector_get(u,0)/d)*(goal_distances[i]-d)/2<<endl;
    gsl_matrix_set(targetPosition,2*i*3+4,0, -(gsl_vector_get(u,1)/d)*(dist_goal-d)/2);
    //cout<<"Direction : "<<-(gsl_vector_get(u,1)/d)*(goal_distances[i]-d)/2<<endl;
    gsl_matrix_set(targetPosition,2*i*3+5,0, -(gsl_vector_get(u,2)/d)*(dist_goal-d)/2);
    //cout<<"Direction : "<<-(gsl_vector_get(u,2)/d)*(goal_distances[i]-d)/2<<endl;

    KinVertex* currVertex1 = atom1->getRigidbody()->getVertex();
    KinVertex* currVertex2 = atom2->getRigidbody()->getVertex();

    //Trace back until the m_root from currVertex
    while ( currVertex1->m_parent != NULL){
      KinEdge* p_edge = currVertex1->m_parent->findEdge(currVertex1);

      int dof_id = p_edge->getDOF()->getIndex();
      if (dof_id!=-1) { // this edge is a DOF
//        Atom* ea1 = p_edge->getBond()->Atom1;
//        Atom* ea2 = p_edge->getBond()->Atom2;
//        Math3D::Vector3 derivativeP = ComputeJacobianEntry(ea1->m_position,ea2->m_position,p);
        Math3D::Vector3 derivativeP = p_edge->getDOF()->getDerivative(atom1->m_position);

        gsl_matrix_set(targetJacobian,2*i*3+0,dof_id,derivativeP.x);
        gsl_matrix_set(targetJacobian,2*i*3+1,dof_id,derivativeP.y);
        gsl_matrix_set(targetJacobian,2*i*3+2,dof_id,derivativeP.z);
      }
      currVertex1 = currVertex1->m_parent;
    }
    //Trace back until the m_root from currVertex
    while ( currVertex2->m_parent != NULL){
      KinEdge* p_edge = currVertex2->m_parent->findEdge(currVertex2);

      int dof_id = p_edge->getDOF()->getIndex();
      if (dof_id!=-1) { // this edge is a DOF
//        Atom* ea1 = p_edge->getBond()->Atom1;
//        Atom* ea2 = p_edge->getBond()->Atom2;
//        Math3D::Vector3 derivativeP = ComputeJacobianEntry(ea1->m_position,ea2->m_position,p);
        Math3D::Vector3 derivativeP = p_edge->getDOF()->getDerivative(atom2->m_position);

        gsl_matrix_set(targetJacobian,2*i*3+3,dof_id,derivativeP.x);
        gsl_matrix_set(targetJacobian,2*i*3+4,dof_id,derivativeP.y);
        gsl_matrix_set(targetJacobian,2*i*3+5,dof_id,derivativeP.z);
      }
      currVertex2 = currVertex2->m_parent;
    }
  }
  gsl_vector_free(u);
}

gsl_matrix* LSNrelativeDirection::determineBestMove(gsl_matrix* N, gsl_matrix* targetJacobian, gsl_matrix* TargetDirection) {

  gsl_matrix *Mn = gsl_matrix_mul(targetJacobian, N);
  /*for (int i=0;i<Mn->size1;i++){
      cout<<i<<" MN ";
      for (int j=0;j<Mn->size2;j++)
          cout<<gsl_matrix_get(Mn,i,j)<<" ";
      cout<<endl;
  }*/
  SVD *svd2 = SVD::createSVD(Mn);//new SVDMKL(Mn);
  // gsl_matrix* inv=svd2->pseudoInverse();
  svd2->UpdateFromMatrix();
  gsl_vector *S = svd2->S;
  /*for (int i=0;i<S->size;i++){
      cout<<i<<" S "<<gsl_vector_get(S,i)<<endl;
  }*/
  gsl_matrix *V = svd2->V;
  gsl_matrix *U = svd2->U;
  gsl_matrix *Ut = gsl_matrix_trans(U);
  gsl_matrix *Tbis = gsl_matrix_mul(Ut, TargetDirection);
  gsl_matrix_free(Ut);
  //cout<<V->size1<<" "<<U->size1<<endl;
  gsl_matrix *moveb = gsl_matrix_calloc(V->size1,
                                        1); //TODO: Was S->size2. Make sure its correct. Amelie S->size2 is V->size1 now I guess
  gsl_matrix_set_zero(moveb);
  for (int i = 0; i < V->size1; i++) { //TODO: Same here Amelie S->size2 is V->size1 now I guess
    if (i < U->size1) { //TODO: Same here Amelie S->size1 is U->size1 now I guess
      if (gsl_vector_get(S, i) * gsl_vector_get(S, i) > 0.0000000000001) {
        //cout<<gsl_vector_get(S, i)<<endl;
        gsl_matrix_set(moveb, i, 0, gsl_matrix_get(Tbis, i, 0) / gsl_vector_get(S, i));
      }
    }
  }
  gsl_matrix_free(Tbis);
  gsl_matrix *move = gsl_matrix_mul(V, moveb);
  gsl_matrix_free(moveb);


  gsl_matrix *bestmove = gsl_matrix_mul(N, move);

  gsl_matrix_free(Mn);
  gsl_matrix_free(move);
  delete svd2;
  /*for (int i=0;i<bestmove->size1;i++){
      cout<<i<<" "<<gsl_matrix_get(bestmove,i,0)<<endl;
  }*/
  return bestmove;
}
