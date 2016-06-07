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
#include "LSNullspaceDirection.h"

#include <map>
#include <cassert>

#include "SamplingOptions.h"
#include "math/SVD.h"
#include "math/gsl_helpers.h"

#include "core/Molecule.h"
#include "core/Chain.h"


using namespace std;

gsl_matrix* LSNullspaceDirection::m_TargetJacobian = NULL;
gsl_matrix* LSNullspaceDirection::m_TargetPosition = NULL;


LSNullspaceDirection::LSNullspaceDirection()
{
  int numMovingAtoms = SamplingOptions::getOptions()->getAtomsMoving()->size();
  m_TargetJacobian = gsl_matrix_calloc( numMovingAtoms*3, 1 );
}


void LSNullspaceDirection::computeGradient(Configuration* conf, Configuration* targetConf, gsl_vector* ret)
{
  Molecule* protein = conf->getMolecule();
  Molecule* target = targetConf->getMolecule();

  //current_q->ComputeCycleJacobianAndNullSpace();
  gsl_matrix* N = conf->getNullspace()->getBasis();
  int dof=protein->totalDofNum();
  gsl_matrix* fullN = gsl_matrix_calloc(dof,N->size2+dof-N->size1);
  gsl_matrix_set_zero(fullN);
  int fulldof=0;
  for (auto const& edge: protein->m_spanning_tree->Edges) {
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
  const vector<Atom*>* atomList = SamplingOptions::getOptions()->getAtomsMoving();
  m_TargetJacobian = gsl_matrix_calloc(atomList->size()*3,fullN->size1);
  gsl_matrix_set_zero(m_TargetJacobian);
  gsl_matrix_set_zero(m_TargetPosition);
  fillmatrices(conf, targetConf);

  gsl_matrix* bestmove = determineBestMove(fullN,m_TargetJacobian,m_TargetPosition);

  for (int i=0;i<fullN->size1;i++){
    gsl_vector_set(ret,i,gsl_matrix_get(bestmove,i,0));
  }
  gsl_matrix_free(bestmove);
  gsl_matrix_free(fullN);

}


void LSNullspaceDirection::fillmatrices(Configuration* current_q, Configuration* pTarget){

  int i=0;
  Molecule* protein=current_q->getMolecule();
  Molecule* target=pTarget->getMolecule();//Todo: Here we have to be careful! This does not work if both configurations refer to the same protein!
  assert(protein!=target);

  const vector<Atom*>* atomList = SamplingOptions::getOptions()->getAtomsMoving();
  vector<Atom*>::const_iterator ait;
  for(ait=atomList->begin(); ait!= atomList->end(); ait++){//only use to extract res IDs and names
    string name = (*ait)->getName();
    int resId = (*ait)->getResidue()->getId();
    string chainName = (*ait)->getResidue()->getChain()->getName();
    Atom* atom = protein->getAtom(chainName,resId,name);
    Coordinate p = atom->m_Position;
    Atom* aTarget = target->getAtom(chainName,resId,name);
    if(aTarget == NULL ){
      //log("dominik")<<"Specified target atom does not exist!"<<endl;
      continue;//skip the non-existing atom
    }
    KinVertex* currVertex = atom->getRigidbody()->getVertex();
    gsl_matrix_set(m_TargetPosition,i*3+0,0,aTarget->m_Position.x-p.x);
    gsl_matrix_set(m_TargetPosition,i*3+1,0,aTarget->m_Position.y-p.y);
    gsl_matrix_set(m_TargetPosition,i*3+2,0,aTarget->m_Position.z-p.z);
    //       cout<<gsl_matrix_get(m_TargetPosition,i*3+0,0)<<" px "<<p.x<<" Tx "<<aTarget->Position.x<<gsl_matrix_get(m_TargetPosition,i*3+1,0)<<" "<<gsl_matrix_get(m_TargetPosition,i*3+2,0)<<" ";
    // trace back until the m_root from currVertex
    //while ( currVertex != protein->m_spanning_tree->Root ) {
    while ( currVertex->m_parent != NULL){
//      Edge* p_edge = currVertex->Parent->Edges.find(currVertex->VertexId)->second;
      KinEdge* p_edge = currVertex->m_parent->findEdge(currVertex);

      //int dof_id = p_edge->DOF_id;
      int dof_id = p_edge->getDOF()->getIndex();
      if (dof_id!=-1) { // this edge is a DO
        Atom* ea1 = p_edge->getBond()->Atom1;
        Atom* ea2 = p_edge->getBond()->Atom2;
        Math3D::Vector3 derivativeP = ComputeJacobianEntry(ea1->m_Position,ea2->m_Position,p);
        gsl_matrix_set(m_TargetJacobian,i*3+0,dof_id,derivativeP.x);
        gsl_matrix_set(m_TargetJacobian,i*3+1,dof_id,derivativeP.y);
        gsl_matrix_set(m_TargetJacobian,i*3+2,dof_id,derivativeP.z);
      }
      currVertex = currVertex->m_parent;
    }
    i++;
  }
  //cout<<endl;
}

gsl_matrix* LSNullspaceDirection::determineBestMove(gsl_matrix* N, gsl_matrix* targetJacobian, gsl_matrix* TargetPosition){

  gsl_matrix* Mn = gsl_matrix_mul(targetJacobian,N);
  SVD* svd2 = SVD::createSVD(Mn);//new MKLSVD(Mn);
  // gsl_matrix* inv=svd2->pseudoInverse();

  gsl_vector* S=svd2->S;
  gsl_matrix* V=svd2->V;
  gsl_matrix* U=svd2->U;
  gsl_matrix* Ut=gsl_matrix_trans(U);
  gsl_matrix* Tbis=gsl_matrix_mul(Ut,TargetPosition);
  gsl_matrix_free(Ut);
  gsl_matrix* moveb = gsl_matrix_calloc(S->size,1); //TODO: Was S->size2. Make sure its correct
  gsl_matrix_set_zero(moveb);
  for( int i=0;i<S->size;i++){ //TODO: Same here
    if (i<S->size){ //TODO: Same here
      if(gsl_vector_get(S,i)*gsl_vector_get(S,i)>0.0000000000001){
        gsl_matrix_set(moveb,i,0,gsl_matrix_get(Tbis,i,0)/gsl_vector_get(S,i));
      }
    }
  }
  gsl_matrix_free(Tbis);
  gsl_matrix* move=gsl_matrix_mul(V,moveb);
  gsl_matrix_free(moveb);


  gsl_matrix* bestmove=gsl_matrix_mul(N,move);

  gsl_matrix_free(Mn);
  gsl_matrix_free(move);
  delete svd2;

  return bestmove;
}

//void LSNullspaceDirection::clashFreeGradient(gsl_vector* gradient, gsl_vector* admissible_gradient, Molecule* protein){
//
////  gsl_matrix* N = Configuration::PreventClashNullSpace->m_nullspaceBasis;
//  gsl_matrix* N = Configuration::ClashAvoidingNullSpace->getBasis();
//  /*gsl_matrix* fullTargetJacobian = gsl_matrix_calloc(m_TargetJacobian->size1,N->size1);
//  gsl_matrix_set_zero(fullTargetJacobian);
//
//if(m_TargetJacobian->size2!= N->size1){
//      map<unsigned int, RigidbodyGraphVertex*>::iterator vit;
//for (vit=protein->m_spanning_tree->Vertex_map.begin(); vit!=protein->m_spanning_tree->Vertex_map.end(); vit++){
//  if( (*vit).second->isRibose ){
//    SugarVertex* v = reinterpret_cast<SugarVertex*>((*vit).second);
//    int dof_id = v->DOF_id;
//    int cycle_dof_id = v->Cycle_DOF_id;
//    if ( cycle_dof_id!=-1 ) {
//      for( int i=0; i!=m_TargetJacobian->size1; i++){
//        gsl_matrix_set(fullTargetJacobian, i, dof_id, gsl_matrix_get(m_TargetJacobian,i,cycle_dof_id));
//      }
//    }
//  }
//}
//for (vector<Edge*>::iterator eit=protein->m_spanning_tree->Edges.begin(); eit!=protein->m_spanning_tree->Edges.end(); ++eit) {
//  int dof_id = (*eit)->DOF_id;
//  int cycle_dof_id = (*eit)->Cycle_DOF_id;
//  if ( cycle_dof_id!=-1 ) {
//    for( int i=0; i!=m_TargetJacobian->size1; i++){
//      gsl_matrix_set(fullTargetJacobian, i, dof_id, gsl_matrix_get(m_TargetJacobian,i,cycle_dof_id));
//    }
//  }
//}
//      }
//
//  for (int j=0;j!=fullTargetJacobian->size2;j++){
//      for( int i=0; i!=m_TargetJacobian->size1; i++){
//              gsl_matrix_set(fullTargetJacobian, i, j, gsl_matrix_get(m_TargetJacobian,i,j));
//      }
//  }*/
//
//  gsl_matrix* bestmove = determineBestMove(N,m_TargetJacobian,m_TargetPosition);
//  for (int i=0;i<bestmove->size1;i++){
//    gsl_vector_set(admissible_gradient,i,gsl_matrix_get(bestmove,i,0));
//  }
//
//  gsl_matrix_free(bestmove);
//  //gsl_matrix_free(fullTargetJacobian);
//  //protein->ProjectOnClashFreeNullSpace(gradient,admissible_gradient);
//}

