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


#include "core/Configuration.h"
#include "core/graph/KinGraph.h"
#include "math3d/primitives.h"
#include "Transformation.h"
#include "Logger.h"



using namespace Math3D;
using namespace std;


void Confvec2MatrixGlobal(KinTree *pTree, Configuration *q, RigidTransform *ms)
{
  int n = q->getNumDOFs();
  int i;

  for(i=0; i<n; ++i)
    ms[i].setIdentity();

  ///	KinVertex *m_root = pTree->Vertex_map.begin()->second;
  KinVertex *root = pTree->m_root;
  KinVertex *node, *newNode;
  KinEdge *pEdge;
  Vector3 vec;

  list<KinVertex *>queue;
  RigidTransform localMat, m1, m2, m3;

  queue.push_back(root);
  root->m_transformation.setIdentity();

  // ms[0].setIdentity(); // the global matrix for the m_root
  // log("debug") << m_root->id << endl;
  while(queue.size()>0)
  {
    node = queue.front();
    queue.pop_front();

    auto children = node->m_edges;

//    for (auto eit=children.begin(); eit!=children.end(); ++eit){
    for (auto const& pEdge: children){
      int dof_id = pEdge->getDOF()->getIndex();
      m1.setIdentity();
      m2.setIdentity();
      m3.setIdentity();
      newNode = pEdge->EndVertex;

      vec = pEdge->getBond()->Atom2->m_Position - pEdge->getBond()->Atom1->m_Position;
      m1.setTranslate(pEdge->getBond()->Atom1->m_Position);
      m2.setRotate(FindRotationMatrix(vec, -q->m_dofs[dof_id])); // !!! Since the FindRotationMatrix is for left hand, choose the negative of the angle
      m3.setTranslate(-1.0 * (pEdge->getBond()->Atom1->m_Position) );

      //TODO: By deriving the closed form for localMat, the computation can be optimized
      localMat = m1 * m2 * m3;

      //log("debugRebuild")<<"RB:"<<pEdge->EndVertex->m_rigidbody->id()<<" "<<pEdge->getDOF()->getIndex()<<" had m_transformation changed by Transformation"<<endl;
      ms[dof_id] = node->m_transformation * localMat; // the global matrix for this node

      newNode->m_transformation = ms[dof_id];
      if(!newNode->m_transformation.isValid(0.0001)){
        cerr<<"Node m_transformation invalid: "<<dof_id<<endl;
      }

      queue.push_back(newNode);
    }

  }

  //Hack for sugars
  queue.push_back(root);
  while(queue.size()>0) {
    node = queue.front();
    queue.pop_front();

    auto children = node->m_edges;
    for (auto eit=node->m_edges.begin(); eit!=node->m_edges.end(); ++eit)
    {
      pEdge = *eit;
      newNode = pEdge->EndVertex;
      queue.push_back(newNode);
    }

  }


}

/**
  Update a local set of transformations defined by subVerts and starting at the specified m_root vertex
  */
void Confvec2MatrixLocal (KinVertex *root, Configuration *q, RigidTransform *ms, vector<KinVertex*> subVerts)
{

  //for(int i=0; i<q->m_numDOFs; ++i)
  //    ms[i].setIdentity();
  for(vector<KinVertex*>::iterator it = subVerts.begin(); it!=subVerts.end(); it++){
    KinVertex* v = *it;
    ms[v->m_rigidbody->id()].setIdentity();
  }

  KinVertex *node, *newNode;
  map<unsigned int,KinEdge*> children;
  map<unsigned int,KinEdge*>::iterator edge_itr;
  KinEdge *pEdge;
  Vector3 vec;

  list<KinVertex *>queue;
  RigidTransform localMat, m1, m2, m3;

  queue.push_back(root);

  while(queue.size()>0)
  {
    node = queue.front();
    queue.pop_front();

    for (auto const& edge: node->m_edges)
    {
      m1.setIdentity();
      m2.setIdentity();
      m3.setIdentity();
      //pEdge = edge_itr->second;
      pEdge = edge;
      newNode = pEdge->EndVertex;

      vec = pEdge->getBond()->Atom2->m_Position - pEdge->getBond()->Atom1->m_Position;
      m1.setTranslate(pEdge->getBond()->Atom1->m_Position);
      m2.setRotate(FindRotationMatrix(vec, -q->m_dofs[pEdge->getDOF()->getIndex()])); // !!! Since the FindRotationMatrix is for left hand, choose the negative of the angle
      m3.setTranslate(-1*pEdge->getBond()->Atom1->m_Position);

      //TODO: By deriving the closed form for localMat, the computation can be optimized
      localMat = m1 * m2 * m3;

      ms[pEdge->getDOF()->getIndex()] = node->m_transformation * localMat; // the global matrix for this node

      newNode->m_transformation = ms[pEdge->getDOF()->getIndex()];

      if(find(subVerts.begin(), subVerts.end(), newNode)!=subVerts.end())
        queue.push_back(newNode);
    }

  }

}

void Confvec2MatrixIndividual(Configuration *q, KinVertex *node, double* iniRef, RigidTransform *ms){

  map<unsigned int,KinEdge*> children;
  map<unsigned int,KinEdge*>::iterator edge_itr;
  KinEdge *pEdge;
  Vector3 vec;

  KinVertex *newNode;
  RigidTransform localMat, m1, m2, m3;


  // log("debug") << "m_parent " << node->id << ": ";
  //for (edge_itr=m_children.begin(); edge_itr != m_children.end(); ++edge_itr)
  for (auto const& pEdge: node->m_edges)
  {
    m1.setIdentity();
    m2.setIdentity();
    m3.setIdentity();
    ///			newNode = pEdge->StartVertex == node ? pEdge->EndVertex : pEdge->StartVertex;
    newNode = pEdge->EndVertex;
    //Measure desired global value with respect to the set value
//		cout<<"Initial ref: "<<iniRef[pEdge->getDOF()->getIndex()]<<", newVal: "<<pEdge->getBond()->getTorsion()<<endl;
//		cout<<"Now checking!"<<endl;

    Atom* atom1 = pEdge->getBond()->Atom1;
    Atom* atom2 = pEdge->getBond()->Atom2;
    int atom_id1 = atom1->getId(); // due to the assertion of Bond, atom_id1 must be smaller than atom_id2
    int atom_id2 = atom2->getId();
    Atom* atom3 = nullptr;
    for (vector<Atom*>::iterator aitr=atom1->Cov_neighbor_list.begin(); aitr!=atom1->Cov_neighbor_list.end(); ++aitr) {
      if ( (*aitr)->getId()==atom_id2 ) continue;
      if ( atom3==nullptr || (*aitr)->getId()<atom3->getId() ) {
        atom3 = *aitr;
      }
    }
    for (vector<Atom*>::iterator aitr=atom1->Hbond_neighbor_list.begin(); aitr!=atom1->Hbond_neighbor_list.end(); ++aitr) {
      if ( (*aitr)->getId()==atom_id2 ) continue;
      if ( atom3==nullptr ) { // || (*aitr)->getId()<atom3->getId()
        //			log("planner")<<"Using hbond neighbor for bond between "<<this->Atom1<<" and "<<this->Atom2<<endl;
        atom3 = *aitr;
      }
    }
    Atom* atom4 = nullptr;
    for (vector<Atom*>::iterator aitr=atom2->Cov_neighbor_list.begin(); aitr!=atom2->Cov_neighbor_list.end(); ++aitr) {
      if ( (*aitr)->getId()==atom_id1 ) continue;
      if ( atom4==nullptr || (*aitr)->getId()<atom4->getId() ) {
        atom4 = *aitr;
      }
    }
    for (vector<Atom*>::iterator aitr=atom2->Hbond_neighbor_list.begin(); aitr!=atom2->Hbond_neighbor_list.end(); ++aitr) {
      if ( (*aitr)->getId()==atom_id1 ) continue;
      if ( atom4==nullptr ) { // || (*aitr)->getId()<atom4->getId()
        atom4 = *aitr;
      }
    }
    Vector3 newPos;
    Vector3 oldPos = atom4->m_Position;
    //if( !atom4->m_bPositionModified ){
    newPos = node->m_transformation.R * oldPos;

    newPos.x += node->m_transformation.t.x;
    newPos.y += node->m_transformation.t.y;
    newPos.z += node->m_transformation.t.z;
    //}
    //else{
    //	newPos=oldPos;
    //}

    double currentTorsion = TorsionalAngle(atom3->m_Position,atom1->m_Position,atom2->m_Position,newPos);

    double rot = q->getGlobalTorsions(pEdge->getDOF()->getIndex()) - currentTorsion;
//		cout<<"Updated torsion from "<< pEdge->getBond()->getTorsion()<<" to "<<currentTorsion<<endl;
    rot = formatRangeRadian(rot);
//		log("planner")<<"Transformation at edge "<<pEdge->getDOF()->getIndex()<<" for node "<<newNode->id<<" with rot "<<q->getGlobalTorsions(pEdge->DOF_id)<<" - "<<currentTorsion<<" = "<<rot<<endl;

    vec = pEdge->getBond()->Atom2->m_Position - pEdge->getBond()->Atom1->m_Position;

    m1.setTranslate(pEdge->getBond()->Atom1->m_Position);
    m2.setRotate(FindRotationMatrix(vec, -rot)); // !!! Since the FindRotationMatrix is for left hand, choose the negative of the angle
    m3.setTranslate(-1*pEdge->getBond()->Atom1->m_Position);

    //TODO: By deriving the closed form for localMat, the computation can be optimized
    localMat = m1 * m2 * m3;
    //log("debugRebuild")<<"RB:"<<pEdge->EndVertex->m_rigidbody->id()<<" "<<pEdge->getDOF()->getIndex()<<" had m_transformation changed by Transformation"<<endl;
    ms[pEdge->getDOF()->getIndex()] = node->m_transformation * localMat; // the global matrix for this node

    newNode->m_transformation = ms[pEdge->getDOF()->getIndex()];

    if(!newNode->m_transformation.isValid(0.0001)){
      cerr<<"Node m_transformation invalid: "<<pEdge->getDOF()->getIndex()<<endl;
    }

  }
}

