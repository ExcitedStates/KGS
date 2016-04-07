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
#include <iostream>
#include <fstream>
#include <list>
#include <gsl/gsl_blas.h>
#include <queue>

#include "Molecule.h"
#include "core/Chain.h"
#include "Residue.h"
#include "core/Bond.h"
#include "ProteinHBond.h"
#include "Grid.h"
#include "CTKTimer.h"
#include "metrics/RMSD.h"
#include "Transformation.h"
#include "math/SVD.h"
#include "IO.h"
#include "Logger.h"
#include <iostream>
#include <iomanip> 		// std::setprecision
#include "SamplingOptions.h"
#include <fstream>
#include <sstream>
#include <math/MKLSVD.h>
#include <set>

const double VDW_SIGMA = 0.2; // sigma = 0.2 kcal/mol
const double VDW_R0 = 3.5; // r_0 = 3.5A
const double VDW_R_MAX = 8; // only compute vdw energy if R_ab <= 8A
//const double STD2 = 1;

using namespace std;


Molecule::Molecule() {
  name_ = "UNKNOWN";
  Atom_pos_index = NULL;
  backup_Atom_pos_index = NULL;
  m_spanning_tree = NULL;
  m_conf = NULL;
  m_Transformation = NULL;
  AtomJacobian1 = NULL;
  AtomJacobian2 = NULL;
  AtomJacobian3 = NULL;
}

Molecule::~Molecule() {
  // delete all atoms
  for (vector<Atom*>::iterator it=atoms.begin(); it != atoms.end(); ++it) {
    delete (*it);
  }

  // delete all chains
  for (auto const& chain: chains){
    delete chain;
  }

  // delete Atom_pos_index and backup_Atom_pos_index
  if (Atom_pos_index!=NULL)
    delete Atom_pos_index;
  if (backup_Atom_pos_index!=NULL)
    delete backup_Atom_pos_index;
  // delete all rigid bodies
  for (map<unsigned int,Rigidbody*>::iterator it=Rigidbody_map_by_id.begin(); it!=Rigidbody_map_by_id.end(); ++it) {
    if ( it->second != NULL )
      delete it->second;
  }

  for (list<Bond *>::iterator it=Cov_bonds.begin(); it != Cov_bonds.end(); ++it) {
    delete (*it);
  }
  for (list<Hbond *>::iterator it=H_bonds.begin(); it != H_bonds.end(); ++it) {
    delete (*it);
  }

  // delete m_spanning_tree
  if (m_spanning_tree!=NULL)
    delete m_spanning_tree;
  if (m_Transformation != NULL)
    delete [] m_Transformation;

  // delete atom jacobians
  if (AtomJacobian1!=NULL)
    gsl_matrix_free(AtomJacobian1);
  if (AtomJacobian2!=NULL)
    gsl_matrix_free(AtomJacobian2);
  if (AtomJacobian3!=NULL)
    gsl_matrix_free(AtomJacobian3);
}

void Molecule::setName (string& name) {
  name_ = name;
}

string Molecule::getName () const {
  return name_;
}

Atom*Molecule::addAtom(
    const std::string& chainName,
    const std::string& resName,
    const int& resId,
    const std::string& atomName,
    const int& atomId,
    const Coordinate& position )
{
  Chain* chain = getChain(chainName);
  if (chain==NULL) { // this is a new chain
    chain = addChain(chainName);
  }
  Atom* ret = chain->addAtom(resName,resId, atomName, atomId, position);
  atoms.push_back(ret);
  return ret;
}

Chain*Molecule::getChain (const string& chain_name) {
  for(auto const& chain: chains)
    if(chain->getName()==chain_name) return chain;

  return NULL;
}

void Molecule::printSummaryInfo() const {
  log() << "Molecule " << getName() << endl;
  for (auto const& chain: chains)
    chain->printSummaryInfo();
}

Chain*Molecule::addChain (const string& chainName) {
  //First check that chain name is not already used
  if(getChain(chainName) != NULL) {
    cerr << "Molecule::addChain - chain " << chainName << " already exists." << endl;
    exit(-1);
  }

  Chain* chain = new Chain(chainName, this);
  chains.push_back(chain);
  return chain;
}

int Molecule::getMaxResidueNumber(){
  int ret= -1000000;
  for(auto const& a: atoms){
    if( a->getResidue()->getId() > ret )
      ret = a->getResidue()->getId();
  }
  return ret;

}

int Molecule::getMinResidueNumber(){
  int ret= 1000000;
  for(vector<Atom*>::iterator ait = atoms.begin(); ait != atoms.end(); ait++){
    Atom* a = *ait;
    if( a->getResidue()->getId() < ret )
      ret = a->getResidue()->getId();
  }
  return ret;

}

/** Gets the atom specified by a residue number and a name. */
Atom* Molecule::getAtom(const string& chainName, const int& resNum, const string& name){
  Chain* chain = this->getChain(chainName);
  if(chain!=NULL){
    Residue* resi = chain->getResidue(resNum);
    if(resi!=NULL){
      return resi->getAtom(name);
    }
  }
  return NULL;
}

Atom* Molecule::getAtom (int atom_id) {
  for (auto const& a: atoms) {
    if(a->getId()==atom_id)
      return a;
  }
  return NULL;
}

void Molecule::printAllAtoms () const {
  for (vector<Atom*>::const_iterator it= atoms.begin(); it != atoms.end(); ++it) {
    (*it)->printSummaryInfo();
  }
}

int Molecule::size() const {
  return atoms.size();
}

void Molecule::updateAtom (int atom_id, Coordinate new_pos) {
  Atom* affected_atom = getAtom(atom_id);
  Atom_pos_index->removeAtom(affected_atom);
  affected_atom->m_Position = new_pos;
  Atom_pos_index->addAtom(affected_atom);
}

void Molecule::indexAtoms () {
  // Atom_pos_index is the current indexing
  if (Atom_pos_index!=NULL)
    delete Atom_pos_index;
  Atom_pos_index = new Grid(this);
}

void Molecule::backupAtomIndex () {
  if (backup_Atom_pos_index!=NULL)
    delete backup_Atom_pos_index;
  backup_Atom_pos_index = Atom_pos_index->deepClone();
}

void Molecule::restoreAtomIndex () {
  if (Atom_pos_index!=NULL)
    delete Atom_pos_index;
  Atom_pos_index = backup_Atom_pos_index->deepClone();
}

//---------------------------------------------------------
bool Molecule::inCollision (string collisionCheckAtoms ) const {

  for (vector<Atom*>::const_iterator itr= atoms.begin(); itr != atoms.end(); ++itr)
    if( (*itr)->isCollisionCheckAtom(collisionCheckAtoms ) )
    if ( Atom_pos_index->inCollision(*itr, Initial_collisions, collisionCheckAtoms)  )
      return true;
  return false;
}
//---------------------------------------------------------
double Molecule::minCollisionFactor (string collisionCheckAtoms) const {
  double minCollFactor = 10000;
  for (vector<Atom*>::const_iterator itr=atoms.begin(); itr!=atoms.end(); ++itr){
    if( (*itr)->isCollisionCheckAtom(collisionCheckAtoms ) ){
      double factor = Atom_pos_index->minFactorWithoutCollision(*itr, Initial_collisions, collisionCheckAtoms);
      if(factor < minCollFactor){
        minCollFactor = factor;
      }
    }
  }
  return minCollFactor;
}

//---------------------------------------------------------
/*
 * Get a list of colliding atoms in the current m_protein configuration.
 */
std::set< pair<Atom*,Atom*>> Molecule::getAllCollisions (std::string collisionCheckAtoms ) const {
  //log() << "InfogetAllCollisions)\t Obtaining a list of atom pair collisions in current m_protein structure..." << endl;
  set< pair<Atom*,Atom*>> collisions;
  for (vector<Atom*>::const_iterator itr= atoms.begin(); itr != atoms.end(); ++itr) {
    if( (*itr)->isCollisionCheckAtom( collisionCheckAtoms ) ) {
      Atom *atom = *itr;
      vector<Atom *> colliding_atoms = Atom_pos_index->getAllCollisions(*itr, Initial_collisions, collisionCheckAtoms);
      for (vector<Atom *>::const_iterator citr = colliding_atoms.begin(); citr != colliding_atoms.end(); ++citr) {
        Atom *colliding_atom = *citr;
        pair<Atom *, Atom *> collision_pair = make_pair(atom, colliding_atom);
        collisions.insert(collision_pair);
      }
    }
  }
  return collisions;
}
//---------------------------------------------------------

void Molecule::printAllCollisions () const {
  for (auto const& atom_pair: getAllCollisions()){
    log() << atom_pair.first->getId() << " " << atom_pair.second->getId() << endl;
  }
}

void Molecule::addCovBond (Bond * bond) {
  //log()<<"Molecule::addCovBond("<<bond<<")"<<endl;
  Cov_bonds.push_back(bond);
  bond->Atom1->addCovBond(bond);
  bond->Atom2->addCovBond(bond);
}

void Molecule::addHbond (Hbond * hb) {
  H_bonds.push_back(hb);
  hb->Atom1->addHbond(hb);
  hb->Atom2->addHbond(hb);
}

unsigned int Molecule::findBestRigidBodyMatch(int rootRBId, Molecule * target){

  unsigned int bestId = 0;
  unsigned int maxId = Rigidbody_map_by_id.size();
  if(rootRBId >= 0){//user-specified, or standard choice of rb id == 0
    if(rootRBId > maxId){
      cout<<"User-specified root id "<<bestId<<" out of bounds. Choosing standard Id."<<endl;
      return bestId;
    }
    cout<<"Choosing user-specified rigid body id "<<rootRBId<<" as root."<<endl;
    return (unsigned int)rootRBId;
  }
  else {
    if( target == NULL ){
      cout<<"No target to determine best root, choosing standard root id 0"<<endl;
      return bestId;
    }
    //Check the rmsd between individual vertices and choose the closest pair as root
    double bestSum = 9999;
    for (map<unsigned int,Rigidbody*>::iterator rbit=Rigidbody_map_by_id.begin(); rbit!=Rigidbody_map_by_id.end(); ++rbit) {
      vector<Atom*> *atomsRMSD = &(rbit->second->Atoms);
      double sum=0;
      bool allAtoms = true;
      Coordinate c1, c2;
      for (vector<Atom*>::iterator ait=atomsRMSD->begin(); ait!=atomsRMSD->end(); ++ait) {
        string name = (*ait)->getName();
        string chainName = (*ait)->getResidue()->getChain()->getName();
        int resId = (*ait)->getResidue()->getId();
        Atom* a1=this->getAtom(chainName,resId, name);
        Atom* a2=target->getAtom(chainName,resId, name);
        if(a1 && a2){
          c1 = a1->m_Position;
          c2 = a2->m_Position;
          sum += c1.distanceSquared(c2);
        } else{//only allow the rb's where all atoms are available in both proteins
          allAtoms = false;
        }
      }
      if(allAtoms){
        sum = sqrt(sum/(atomsRMSD->size()));
        if(sum<=bestSum){
          bestSum=sum;
          bestId=rbit->first;
        }
      }
    }
    cout<<"Choosing rigid body id "<<bestId<<" as root. Root RB rmsd: "<<bestSum<<endl;
    return bestId;
  }
}

void Molecule::buildRigidbodyTree (unsigned int rootRBId, bool flexibleSugar ) {
  //log() << "In buildRigidbodyTree" << endl;
  m_spanning_tree = new RigidbodyTree();

  m_spanning_tree->root = m_spanning_tree->addVertex(rootRBId, Rigidbody_map_by_id[rootRBId], flexibleSugar);


  // add all rigid bodies as vertices into Rigidbody_graph
  int dof_id = 0;
  for (map<unsigned int,Rigidbody*>::iterator it=Rigidbody_map_by_id.begin(); it!=Rigidbody_map_by_id.end(); ++it) {
    RigidbodyGraphVertex* vtx;

    //if (it->first == rootRBId ){
    if (it->second->id() == rootRBId ){
      vtx = m_spanning_tree->root;
    }else{
      vtx = m_spanning_tree->addVertex(it->second->id(),it->second, flexibleSugar);
    }

    if(vtx->isRibose){//RFonseca
      SugarVertex* v = reinterpret_cast<SugarVertex*>(vtx);
      v->DOF_id = dof_id++;
      m_spanning_tree->num_DOFs++;
    }

  }

  // Used to ensure that edges are only created once for any pair of rigid bodies
  bool visited[Rigidbody_map_by_id.size()+1];
  int count = 0;
  for (unsigned int i=0; i<=Rigidbody_map_by_id.size(); ++i) {
    visited[i] = false;
    count++;
  }

  list<Edge*> cycle_edges;
  list<RigidbodyGraphVertex*> queue;
  queue.push_back(m_spanning_tree->root);
  while (!queue.empty()) {
    RigidbodyGraphVertex* current_vertex = queue.front();
    visited[current_vertex->id] = true;
    Rigidbody* rb1 = current_vertex->Rb_ptr;

    for (vector<Bond *>::iterator bit1=rb1->Bonds.begin(); bit1 != rb1->Bonds.end(); ++bit1) {

      // for each unvisited vertex, look for common bonds
      for (map<unsigned int,Rigidbody*>::iterator it=Rigidbody_map_by_id.begin(); it!=Rigidbody_map_by_id.end(); ++it) {
        unsigned int i = it->first;
        if ( visited[i]==false ) {
          RigidbodyGraphVertex *vtx = m_spanning_tree->Vertex_map.find(i)->second;
          Rigidbody* rb2 = vtx->Rb_ptr;
          for (vector<Bond *>::iterator bit2=rb2->Bonds.begin(); bit2 != rb2->Bonds.end(); ++bit2) {
            if ( (*bit1)==(*bit2) ) {
              // if it's a H-bond, it closes a cycle. Add it in CycleAnchorEdges.
              if ( (*bit1)->isHbond() ) {
                Hbond *hb = (Hbond *)(*bit1);
                Edge *h_edge = new Edge(current_vertex,vtx,hb);
                cycle_edges.push_back(h_edge);
              }
                // if it's a covalent bond, add it into the tree edges
              else {
                queue.push_back(vtx);
                m_spanning_tree->addEdgeDirected(current_vertex,vtx,*bit1, dof_id++); // the edge and the bond associated with the edge are pointing from rb1 to rb2
                m_spanning_tree->num_DOFs++;
              }
            }
          }
        }
      }
    }
    queue.pop_front();
  } // end while

  // For each hbond Edge, find the common ancestor.
  // Assign the Cycle_DOF_id for all the edges from the hbond ends to the common ancestor
  int cycle_dof_id = 0;
  for ( list<Edge*>::iterator it=cycle_edges.begin(); it!=cycle_edges.end(); ++it) {
    Edge* h_edge = *it;
    RigidbodyGraphVertex *startv = h_edge->StartVertex;
    RigidbodyGraphVertex *endv = h_edge->EndVertex;
    RigidbodyGraphVertex *ancestor = m_spanning_tree->findCommonAncestor(startv,endv);
    m_spanning_tree->CycleAnchorEdges.push_back( make_pair(h_edge,ancestor) );
    RigidbodyGraphVertex *vertex, *parent;
    for (vertex=startv; vertex!=ancestor; vertex=parent) {
      parent = vertex->Parent;
      if(parent->isRibose){//RFonseca
        SugarVertex* v = reinterpret_cast<SugarVertex*>(parent);
        if(v->Cycle_DOF_id==-1)	v->Cycle_DOF_id = cycle_dof_id++;
      }

      Edge *edge = parent->findEdge(vertex);

      if (edge->Cycle_DOF_id==-1) {
        edge->Cycle_DOF_id = cycle_dof_id;
        ++cycle_dof_id;
      }
    }
    for (vertex=endv; vertex!=ancestor; vertex=parent) {
      parent = vertex->Parent;
      if(parent->isRibose){//RFonseca
        SugarVertex* v = reinterpret_cast<SugarVertex*>(parent);
        if(v->Cycle_DOF_id==-1)	v->Cycle_DOF_id = cycle_dof_id++;
      }
      Edge *edge = parent->findEdge(vertex);
      if (edge->Cycle_DOF_id==-1) {
        edge->Cycle_DOF_id = cycle_dof_id;
        ++cycle_dof_id;
      }
    }
  }
  m_spanning_tree->Cycle_DOF_num = cycle_dof_id;
  log("dimitar") << "InfoProtein)\t Molecule has " << m_spanning_tree->Cycle_DOF_num << " cycle DOFs." << endl;
  log("dimitar") << "InfoProtein)\t Molecule has " << m_spanning_tree->Edges.size() << " total edges." << endl;
  log("dimitar") << "InfoProtein)\t Molecule has m_numDOFs " << m_spanning_tree->num_DOFs << " total DOFs." << endl;

//	m_spanning_tree = Rigidbody_graph->ComputeSpanningTree();

  m_Transformation = new RigidTransform [m_spanning_tree->num_DOFs];//m_spanning_tree->Edges.size()];
/*	log() << "*** Spanning Tree ***" << endl;
	m_spanning_tree->print();
	Rigidbody_graph->findCycleClusters();
    log() << "*** Rigidbody Graph ***" << endl;
	Rigidbody_graph->print();
*/
  //m_spanning_tree->printHTML();
  //delete[] visited;

  //Create a sorted list of vertices for use in the euclideanGradient
  for (unsigned int i=0; i<=Rigidbody_map_by_id.size(); ++i) {
    visited[i] = false;
  }
  Edge* currEdge = m_spanning_tree->Edges.back();
  RigidbodyGraphVertex *currVertex = currEdge->EndVertex;
  int currId = currVertex->id;

  m_spanning_tree->m_sortedVertices.push_back(make_pair( currId, currVertex ));
  visited[currId] = true;

  while( currVertex != m_spanning_tree->root){
//		cout<<"Current vertex id "<<currId<<endl;fflush(stdout);
    RigidbodyGraphVertex *parent = currVertex->Parent;

    if(parent->edges.size() == 1){
//    		cout<<"Single edged m_parent vertex "<<m_parent->id<<endl;
      currVertex = parent;
      currId = currVertex->id;
//    		cout<<"Adding vertex with id "<<currId <<endl;fflush(stdout);
      m_spanning_tree->m_sortedVertices.push_back(make_pair( currId,currVertex));
      visited[currId] = true;

    }
    else{//multiple branches open
      currVertex = parent;
      int numEdgesLeft = currVertex->edges.size();
      while( numEdgesLeft != 0){
        //map<unsigned int,Edge*>::iterator eit;
        for( auto eit = currVertex->edges.begin(); eit != currVertex->edges.end(); eit++){
          if(visited[(*eit)->EndVertex->id] == true){
            numEdgesLeft--;
            continue;
          }else{

            currVertex = (*eit)->EndVertex;
            //currVertex = eit->second->EndVertex;
            numEdgesLeft = currVertex->edges.size();
            break;
          }
        }
      }
      currId = currVertex->id;
      m_spanning_tree->m_sortedVertices.push_back(make_pair(currId, currVertex) );
      visited[currId] = true;
    }
  }
}

RigidbodyGraphVertex*Molecule::getRigidbodyGraphVertex (Atom* atom) const {
  //TODO: Not sure whats happening here. Clean up
  int smallest_dof_id = m_spanning_tree->num_DOFs - 1;//Edges.size()-1; // the possible maximum DOF id
  RigidbodyGraphVertex* vertex_with_smallest_dof_id = m_spanning_tree->getVertex(atom->getRigidbody()->id());
  RigidbodyGraphVertex* vertex = m_spanning_tree->getVertex(atom->getRigidbody()->id());
  RigidbodyGraphVertex* parent = vertex->Parent;
  if (parent==NULL)
    return vertex;
  Edge* edge_to_parent = parent->findEdge(vertex);
  if ( edge_to_parent->DOF_id < smallest_dof_id ) {
    vertex_with_smallest_dof_id = vertex;
  }
  return vertex_with_smallest_dof_id;
}

void Molecule::computeAtomJacobian (Atom* atom, gsl_matrix **j_addr) {
  log("debug") << "in Molecule::computeAtomJacobian" << endl;
  if (*j_addr==NULL) {
    log("debug") << "in Molecule::computeAtomJacobian - AtomJacobian1 is NULL" << endl;
    *j_addr = gsl_matrix_calloc(3,totalDofNum());
  } else {
    log("debug") << "in Molecule::computeAtomJacobian - AtomJacobian1 is not NULL" << endl;
    gsl_matrix_set_all(*j_addr,0);
  }
  log("debug") << "in Molecule::computeAtomJacobian - after IF" << endl;
  gsl_matrix* jacobian = *j_addr;
  RigidbodyGraphVertex *vertex = getRigidbodyGraphVertex(atom);
  //RigidbodyGraphVertex *vertex1 = getRigidbodyGraphVertex( getAtom( 35 ) );
  atom->printSummaryInfo();
  //getAtom( 35 )->printSummaryInfo();
  log("debug") << "in Molecule::computeAtomJacobian - vertex->id= " << vertex->id << endl;
  //log("debug") << "in Molecule::computeAtomJacobian - vertex1->id= " << vertex1->id << endl;
  if ( vertex->Parent==NULL )
    log("debug") << " vertex->Parent==NULL " << endl;
  //if ( vertex1->Parent==NULL )
  //	log("debug") << " vertex1->Parent==NULL " << endl;
  //Edge* edge1 = vertex1->Parent->Edges.find(vertex1->id)->second;
  //edge1->printShort();
  log("debug") << "in Molecule::computeAtomJacobian - before while loop" << endl;
  while (vertex!=m_spanning_tree->root) {
    log("debug") << "in Molecule::computeAtomJacobian - in while loop" << endl;
    RigidbodyGraphVertex *parent;
    if ( vertex->Parent!=NULL )
      parent = vertex->Parent;
    else
      parent = vertex;
    log("debug") << " in Molecule::computeAtomJacobian - m_parent->Edges.size() " << parent->edges.size() << endl;
    Edge* edge = parent->findEdge(vertex);
    log("debug") << "in Molecule::computeAtomJacobian - in while loop - after m_parent->Edges.find(vertex->id)->second" << endl;
    int dof_id = edge->DOF_id;
    Bond * bond_ptr = edge->getBond();
    Coordinate bp1 = bond_ptr->Atom1->m_Position;
    Coordinate bp2 = bond_ptr->Atom2->m_Position;
    log("debug") << "in Molecule::computeAtomJacobian - in while loop - before jacobian_entry " << endl;
    Vector3 jacobian_entry = ComputeJacobianEntry(bp1,bp2,atom->m_Position);
    gsl_matrix_set(jacobian,0,dof_id,jacobian_entry.x);
    gsl_matrix_set(jacobian,1,dof_id,jacobian_entry.y);
    gsl_matrix_set(jacobian,2,dof_id,jacobian_entry.z);
    vertex = parent;
  }
  log("debug") << "in Molecule::computeAtomJacobian - after while loop" << endl;
}


gsl_vector*Molecule::getEndEffectors(){
  gsl_vector* ret = gsl_vector_alloc(  (m_spanning_tree->CycleAnchorEdges).size()*6  );

  int i=0;
  for (vector< pair<Edge*,RigidbodyGraphVertex*> >::iterator it=m_spanning_tree->CycleAnchorEdges.begin(); it!=m_spanning_tree->CycleAnchorEdges.end(); ++it) {
    // get end-effectors
    Edge* edge_ptr = it->first;
    Hbond * bond_ptr = (Hbond *)(edge_ptr->getBond());
    Vector3 p1 = bond_ptr->Atom1->m_Position;
    Vector3 p2 = bond_ptr->Atom2->m_Position;
    Vector3 p1Diff = (p1-bond_ptr->getIdealHPoint());
    Vector3 p2Diff = (p2-bond_ptr->getIdealAcceptorPoint());
    gsl_vector_set(ret, i+0, p1Diff.x);
    gsl_vector_set(ret, i+1, p1Diff.y);
    gsl_vector_set(ret, i+2, p1Diff.z);
    gsl_vector_set(ret, i+3, p2Diff.x);
    gsl_vector_set(ret, i+4, p2Diff.y);
    gsl_vector_set(ret, i+5, p2Diff.z);
    i+=6;
  }
  return ret;
}


void Molecule::ProjectOnCycleNullSpace (gsl_vector *to_project, gsl_vector *after_project) {

  // No cycle
  //if(!m_conf->CycleNullSpace)
  if(!m_conf->getNullspace())
  {
    gsl_vector_memcpy(after_project, to_project);
    return;
  }

  //if (to_project->size > m_conf->CycleNullSpace->n) {
  if( to_project->size > m_conf->getNullspace()->NumDOFs() ) {
    // The input vectors contain all DOFs, however, the null space only contains DOFs in cycles.
    // Convert the DOFs in the input vectors to DOFs in cycles.
    //gsl_vector *to_proj_short = gsl_vector_calloc(m_conf->CycleNullSpace->n);
    //gsl_vector *after_proj_short = gsl_vector_calloc(m_conf->CycleNullSpace->n);
    gsl_vector *to_proj_short = gsl_vector_calloc(m_conf->getNullspace()->NumDOFs());
    gsl_vector *after_proj_short = gsl_vector_calloc(m_conf->getNullspace()->NumDOFs());
    map<unsigned int, RigidbodyGraphVertex*>::iterator vit;
    for (vit=m_spanning_tree->Vertex_map.begin(); vit!=m_spanning_tree->Vertex_map.end(); vit++){
      if( (*vit).second->isRibose ){
        SugarVertex* v = reinterpret_cast<SugarVertex*>((*vit).second);
        int dof_id = v->DOF_id;
        int cycle_dof_id = v->Cycle_DOF_id;
        if ( cycle_dof_id!=-1 ) {
          gsl_vector_set(to_proj_short,cycle_dof_id,gsl_vector_get(to_project,dof_id));
        }
      }
    }
    for (vector<Edge*>::iterator eit=m_spanning_tree->Edges.begin(); eit!=m_spanning_tree->Edges.end(); ++eit) {
      int dof_id = (*eit)->DOF_id;
      int cycle_dof_id = (*eit)->Cycle_DOF_id;
      if ( cycle_dof_id!=-1 ) {
        gsl_vector_set(to_proj_short,cycle_dof_id,gsl_vector_get(to_project,dof_id));
      }
    }
    // Project onto the null space
    //cout<<"Molecule::ProjectOnCycleNullSpace(..) .. "<<to_proj_short->size<<endl;;

    //m_conf->CycleNullSpace->ProjectOnNullSpace(to_proj_short,after_proj_short);
    //New version
    m_conf->getNullspace()->ProjectOnNullSpace(to_proj_short, after_proj_short);

    // Convert back to full length DOFs vector
    for (vit=m_spanning_tree->Vertex_map.begin(); vit!=m_spanning_tree->Vertex_map.end(); vit++){
      if( (*vit).second->isRibose ){
        SugarVertex* v = reinterpret_cast<SugarVertex*>((*vit).second);
        int dof_id = v->DOF_id;
        int cycle_dof_id = v->Cycle_DOF_id;
        if ( cycle_dof_id!=-1 ) {
          gsl_vector_set(after_project,dof_id,gsl_vector_get(after_proj_short,cycle_dof_id));
        }
        else if ( dof_id!=-1 ) {
          gsl_vector_set(after_project,dof_id,gsl_vector_get(to_project,dof_id));
        }
      }
    }
    for (vector<Edge*>::iterator eit=m_spanning_tree->Edges.begin(); eit!=m_spanning_tree->Edges.end(); ++eit) {
      int dof_id = (*eit)->DOF_id;
      int cycle_dof_id = (*eit)->Cycle_DOF_id;
      if ( cycle_dof_id!=-1 ) {
        gsl_vector_set(after_project,dof_id,gsl_vector_get(after_proj_short,cycle_dof_id));
      }
      else if ( dof_id!=-1 ) {
        gsl_vector_set(after_project,dof_id,gsl_vector_get(to_project,dof_id));
      }
    }
    gsl_vector_free(to_proj_short);
    gsl_vector_free(after_proj_short);
  }
  else {
    m_conf->getNullspace()->ProjectOnNullSpace(to_project, after_project);
  }
}

void Molecule::alignReferencePositionsTo(Molecule * base){
  this->RestoreAtomPos();
  //Align conformations
  metrics::RMSD::align(this,base);

  for (vector<Atom*>::const_iterator it=atoms.begin(); it!=atoms.end(); ++it) {
    (*it)->m_referencePosition = (*it)->m_Position;
  }
}

void Molecule::RestoreAtomPos(){
  for (auto const& a: atoms)
    a->m_Position = a->m_referencePosition;

  m_conf = NULL;
  restoreAtomIndex();

}

void Molecule::SetConfiguration(Configuration *q){
  if( m_conf==q) return;

  RestoreAtomPos();

  m_conf = q;
  _SetConfiguration(q);

  if(q->getGlobalTorsions() == NULL){
    log("dominik")<<"Now updating global torsions"<<endl;
    q->updateGlobalTorsions();
  }
}



// set the positions of atoms at configuration q (according to the spanning tree)
void Molecule::_SetConfiguration(Configuration *q ){
  // assume the base vector is 0, relative rotation to original position
  Confvec2MatrixGlobal(m_spanning_tree, q, m_Transformation);

  //Initialize queue
  list<RigidbodyGraphVertex *> queue;
  RigidbodyGraphVertex *root = m_spanning_tree->root;
  queue.push_back(root);

  while(queue.size()>0){
    RigidbodyGraphVertex* node = queue.front();
    queue.pop_front();

    for(auto const& pEdge: node->edges){
      pEdge->EndVertex->TransformAtomPosition(m_Transformation+pEdge->DOF_id);
      queue.push_back(pEdge->EndVertex);
    }
  }

  indexAtoms();
}


/**
  Set the positions of atoms only within subVerts so they correspond to configuration q.
  Assumes root is member of subVerts.
  */
void Molecule::_SetConfiguration(Configuration *q, RigidbodyGraphVertex* root, vector<RigidbodyGraphVertex*>& subVerts){//, bool usePosition2){
  m_conf = q;

  // assume the base vector is 0
  Confvec2MatrixLocal(root, q, m_Transformation, subVerts);

  list<RigidbodyGraphVertex *>queue;

  queue.push_back(root);
  while(queue.size()>0)
  {
    RigidbodyGraphVertex* node = queue.front();
    queue.pop_front();

    //map<unsigned int,Edge*> m_children = node->Edges;

    //for (map<unsigned int,Edge*>::iterator edge_itr=m_children.begin(); edge_itr != m_children.end(); ++edge_itr){
    for (auto const& pEdge: node->edges){
      //Edge* pEdge = edge_itr->second;
      RigidbodyGraphVertex* newNode = pEdge->EndVertex;

      newNode->TransformAtomPosition(m_Transformation+pEdge->DOF_id);//,usePosition2);

      if(find(subVerts.begin(), subVerts.end(), newNode)!=subVerts.end())
        queue.push_back(newNode);
    }
  }
  //if (!usePosition2)
  indexAtoms();
}

int Molecule::totalDofNum () const {
  if (m_spanning_tree==NULL) {
    cerr << "Error: to get the total number of DOFs in the m_protein, you have to call Molecule::buildRigidbodyTree() first." << endl;
    exit(1);
  }
  return m_spanning_tree->num_DOFs;
}

gsl_vector*Molecule::vdwGradient () { // minimize the L-J potential
  log("debug") << "in Molecule::vdwGradient" << endl;
  gsl_vector* gradient = gsl_vector_calloc(totalDofNum());
  gsl_vector* p12 = gsl_vector_calloc(3);
  gsl_vector* p_temp = gsl_vector_calloc(totalDofNum());

  // For each atom, look for it's neighbors which are not in the same rigid body and whose ID is bigger (to avoid double count). For each such neighbor, compute dU(R_ab)/d(delta_q).
  // 1. For each atom, find it's neighbors which are not in the same rigid body
  for (vector<Atom*>::iterator ait= atoms.begin(); ait != atoms.end(); ++ait) {
    log("debug") << "in Molecule::vdwGradient - start 1st loop" << endl;
    Atom* atom1 = *ait;
    computeAtomJacobian(atom1,&AtomJacobian1);
    log("debug") << "in Molecule::vdwGradient - after computeAtomJacobian(atom1,&AtomJacobian1);" << endl;
    vector<Atom*> neighbors = Atom_pos_index->getNeighboringAtoms(atom1,true,true,true,VDW_R_MAX);
    // 2. For each such neighbor, compute dU(R_ab)/d(q)
    // The formula is -12*VDW_SIGMA*(VDW_R0^6*r_12^(-8)-VDW_R0^12*r_12^(-14))*(J1-J2)'(P1-P2)),
    // where ' is transpose.
    log("debug") << "in Molecule::vdwGradient - before 2nd loop" << endl;
    for (vector<Atom*>::iterator ait2=neighbors.begin(); ait2!=neighbors.end(); ++ait2) {
      if (atom1->inSameRigidbody(*ait2)) continue;
      Atom* atom2 = *ait2;
      double r_12 = atom1->distanceTo(atom2);
      double front_constant_part = (-12)*VDW_SIGMA*(pow(VDW_R0,6)*pow(r_12,-8)-pow(VDW_R0,12)*pow(r_12,-14));
      computeAtomJacobian(atom2,&AtomJacobian2);
      if (AtomJacobian3==NULL)
        AtomJacobian3 = gsl_matrix_calloc(3,totalDofNum());
      gsl_matrix_memcpy(AtomJacobian3,AtomJacobian1);
      gsl_matrix_sub(AtomJacobian3,AtomJacobian2); // AtomJacobian3 holds the result of substraction
      Vector3 p12_v3 = atom1->m_Position - atom2->m_Position;
      Coordinate::copyToGslVector(p12_v3, p12);
      gsl_blas_dgemv(CblasTrans,1,AtomJacobian3,p12,0,p_temp);
      gsl_vector_scale(p_temp,front_constant_part);
      gsl_vector_add(gradient,p_temp);
    }
  }
  log("debug") << "in Molecule::vdwGradient - after loop" << endl;

  gsl_vector_free(p12);
  gsl_vector_free(p_temp);
  return gradient;
}

pair<double,double> Molecule::vdwEnergy (set< pair<Atom*,Atom*> >* allCollisions, string collisionCheck) { // compute the total vdw energy, excluding the covalent bonds and atoms in the same rigid body

  double energy=0, collFreeEnergy=0,d_12, ratio, vdw_r1, vdw_d12, epsilon1, epsilon_12;
  set< pair<Atom*,Atom*> >::iterator mit;
  // for each atom, look for it's neighbors.
  // OLD: For each such neighbor, compute U(R_ab)=VDW_SIGMA*((VDW_R0/r_12)^12-2*(VDW_R0/r_12)^6)
  // Corrected: For each such neighbor, compute U(R_ab)=epsilon_ij*(vdw_r12/r_12)^12-2*(VDW_R0/r_12)^6) and sum up everything.
  // CHARMM: http://www.charmmtutorial.org/index.php/The_Energy_Function#Energy_calculation
  for (vector<Atom*>::const_iterator ait=atoms.begin(); ait!=atoms.end(); ++ait) {
    Atom* atom1 = *ait;
    if(!(atom1->isCollisionCheckAtom(collisionCheck)) ){//we only use atoms that are also used for clash detection
      continue;
    }
    vdw_r1 = atom1->getRadius();
    epsilon1 = atom1->getEpsilon();
    vector<Atom*> neighbors = Atom_pos_index->getNeighboringAtomsVDW(atom1,true,true,true,true,VDW_R_MAX);
    for (vector<Atom*>::const_iterator ait2=neighbors.begin(); ait2!=neighbors.end(); ++ait2) {
      Atom* atom2 = *ait2;

      if (!(atom2->isCollisionCheckAtom(collisionCheck)) || atom1->inSameRigidbody(atom2))
        continue;//we only use atoms that are also used for clash detection (otherwise they can be too close)

      //Check initial collisions --> always excluded
      pair<Atom*,Atom*> collision_pair = make_pair(atom1,atom2);
//      set< pair<Atom*,Atom*>,int >::const_iterator mit=Initial_collisions.find(collision_pair);
      auto mit = Initial_collisions.find(collision_pair);
      if ( mit!=Initial_collisions.end() )//ignore initial collision atoms
        continue;

      d_12 = atom1->distanceTo(atom2);
      vdw_d12 = vdw_r1 + atom2->getRadius(); // from CHARMM: Todo: do we need the arithmetic mean or the sum?
      ratio = vdw_d12/d_12;
      epsilon_12 = sqrt(epsilon1 * (atom2->getEpsilon())); //from CHARMM: geometric mean
      double atomContribution = 4 * epsilon_12 * (pow(ratio,12)-2*pow(ratio,6));

      //Full enthalpy including atoms in clash constraints
      energy += atomContribution;

      //Clash-constraints: excluded for special 'clash-constraint-free' enthalpy
      mit = allCollisions->find(collision_pair);
      if ( mit!=allCollisions->end() )
        continue;

      collFreeEnergy += atomContribution;
    }
  }
  return make_pair(energy, collFreeEnergy);
}

///Create a set of hbonds from the hbond list of another protein
void Molecule::setToHbondIntersection (Molecule * p2) {
  Hbond * hBond;
  Atom *hatom, *acceptor, *donor, *AAatom;
  list<Hbond *> intersection;
  int count1=0, count2=0;
  for (list<Hbond *>::iterator itr2=p2->H_bonds.begin(); itr2 != p2->H_bonds.end(); ++itr2) {
    hBond= (*itr2);

    hatom = this->getAtom(hBond->Hatom->getResidue()->getChain()->getName(),hBond->Hatom->getResidue()->getId(), hBond->Hatom->getName());
    acceptor = this->getAtom(hBond->Acceptor->getResidue()->getChain()->getName(),hBond->Acceptor->getResidue()->getId(), hBond->Acceptor->getName());
    donor = this->getAtom(hBond->Donor->getResidue()->getChain()->getName(),hBond->Donor->getResidue()->getId(), hBond->Donor->getName());
    AAatom = this->getAtom(hBond->AA->getResidue()->getChain()->getName(),hBond->AA->getResidue()->getId(), hBond->AA->getName());

    if( hatom && acceptor){
      Hbond * new_hb = new Hbond(hatom, acceptor, donor, AAatom);
      this->addHbond(new_hb);
      intersection.push_back(hBond);
      count1++;
    }
    else{
      cout<<"Could not find specified hbond in other protein: ";
      cout<<hBond->Hatom->getResidue()->getId()<<" "<<hBond->Hatom->getName()<<", "<<hBond->Acceptor->getResidue()->getId()<<" "<<hBond->Acceptor->getName();
      cout<<" Deleting to make the same set!"<<endl;
      Atom* hAtom = hBond->Hatom;
      Atom* acceptor = hBond->Acceptor;
      hAtom->removeHbond(hBond);
      acceptor->removeHbond(hBond);
//			p2->H_bonds.erase(itr2);
      count2++;
    }
  }
  p2->H_bonds=intersection;
}

bool Molecule::hasCycle() const {
  return (m_spanning_tree->CycleAnchorEdges.size() != 0);
}

void Molecule::printBackboneAngleAndLength (string length_file, string angle_file) const {
  ofstream length_output, angle_output;
  if (length_file!="") {
    length_output.open(length_file.c_str());
  }
  if (angle_file!="") {
    angle_output.open(angle_file.c_str());
  }

  for(auto const& chain: chains){
    for (auto const& res: chain->getResidues()) {
      log() << "Chain " << chain->getName() << " Residue " << res->getId() << " " << res->getName() << endl;
      Atom* curN = res->name_to_atom_map.find("N")->second;
      Atom* curCA = res->name_to_atom_map.find("CA")->second;
      Atom* curC = res->name_to_atom_map.find("C")->second;
      double angNCAC = VectorAngle(curN->m_Position -curCA->m_Position,curC->m_Position -curCA->m_Position);
      // Angle(N,CA,C)
      if (angle_file!="")
        angle_output << toDegree(angNCAC) << endl;
      else
        log() << toDegree(angNCAC) << endl;
      // Length(N,CA) and Length(CA,C)
      if (length_file!="")
        length_output << curN->distanceTo(curCA) << endl << curCA->distanceTo(curC) << endl;
      else
        log() << curN->distanceTo(curCA) << endl << curCA->distanceTo(curC) << endl;
      if (res!=chain->getResidues().back()) {
        Atom* nextN = res->getNextResidue()->name_to_atom_map.find("N")->second;
        Atom* nextCA = res->getNextResidue()->name_to_atom_map.find("CA")->second;
        double angCACN = VectorAngle(curCA->m_Position -curC->m_Position,nextN->m_Position -curC->m_Position);
        double angCNCA = VectorAngle(curC->m_Position -nextN->m_Position,nextCA->m_Position -nextN->m_Position);
        // Angle(CA,C,+N) and Angle(C,+N,+CA)
        if (angle_file!="")
          angle_output << toDegree(angCACN) << endl << toDegree(angCNCA) << endl;
        else
          log() << toDegree(angCACN) << endl << toDegree(angCNCA) << endl;
        // Length(C,+N)
        if (length_file!="")
          length_output << curC->distanceTo(nextN) << endl;
        else
          log() << curC->distanceTo(nextN) << endl;
      }
    }
  }

  if (length_file!="")
    length_output.close();
  if (angle_file!="")
    angle_output.close();
}

int Molecule::countOriginalDofs () const {
  int num = 0;
  for (list<Bond *>::const_iterator it=Cov_bonds.begin(); it != Cov_bonds.end(); ++it) {
    Atom* a1 = (*it)->Atom1;
    Atom* a2 = (*it)->Atom2;
    if ( a1->Cov_neighbor_list.size()==1 || a2->Cov_neighbor_list.size()==1 )
      continue;
    ++num;
  }
  num += H_bonds.size();
  return num;
}

Coordinate Molecule::centerOfMass () const {
  Coordinate cur_position, center_of_mass;
  double mass, sum_of_mass=0;
  for (vector<Atom*>::const_iterator it= atoms.begin(); it != atoms.end(); ++it) {
    cur_position = (*it)->m_Position;
    mass = (*it)->getMass();
    center_of_mass += cur_position * mass;
    sum_of_mass += mass;
  }
  center_of_mass /= sum_of_mass;
  return center_of_mass;
}

Coordinate Molecule::centerOfGeometry () const {
  Coordinate center;
  for (vector<Atom*>::const_iterator it= atoms.begin(); it != atoms.end(); ++it) {
    center += (*it)->m_Position;
  }
  center /= atoms.size();
  return center;
}

void Molecule::checkCycleClosure(Configuration *q){
  SetConfiguration(q);
  //Todo: Use intervals for hydrogen bond angles and lengths
  vector< pair<Edge*,RigidbodyGraphVertex*> >::iterator pair_it;
  Edge *pEdge;
  int id=1;
  double maxViolation = 0.0;
  for (pair_it=m_spanning_tree->CycleAnchorEdges.begin(); pair_it!=m_spanning_tree->CycleAnchorEdges.end(); ++pair_it) {
    pEdge = pair_it->first;
    Atom* atom1 = pEdge->getBond()->Atom1;
    Atom* atom2 = pEdge->getBond()->Atom2;
    Hbond * hBond = reinterpret_cast<Hbond *>(pEdge->getBond());
    float distanceChange = atom1->m_Position.distanceTo(atom2->m_Position)-hBond->iniLength;
    float rightAngleChange = hBond->getRightAngle() - hBond->iniOrientRight;
    float leftAngleChange = hBond->getLeftAngle() - hBond->iniOrientLeft;

    double distanceViolation = distanceChange / hBond->iniLength * 100;
    double absViolation = Abs(distanceViolation);

    log("report")<<"hBond strain at "<<id<<" between "<<atom1->getId()<<" and "<<atom2->getId()<<" in res "<<atom1->getResidue()->getName()<<atom1->getResidue()->getId()<<": " << distanceViolation <<" %"<<endl;
    if(absViolation > maxViolation){
      q->m_maxConstraintViolation = absViolation;
      maxViolation = absViolation;
    }
//		if(Abs(distanceViolation) > 10){//10 % change of length
//			log("report") <<"Distance violation at "<<id<<" between "<<atom1->getId()<<" and "<<atom2->getId()<<" in res "<<atom1->getResidue()->getName()<<atom1->getResidue()->getId()<<": " << distanceViolation <<" %"<<endl;
//		}
//		float rightAngleChange = hBond->getRightAngle() - hBond->iniOrientRight;
//		if(Abs(rightAngleChange) > 0.00001){
//			log("report") <<"Right angle violation "<<id<<" between "<<atom1->getId()<<" and "<<atom2->getId()<<" in res "<<atom1->getResidue()->getName()<<atom1->getResidue()->getId()<<": " << rightAngleChange<<endl;
//		}
//		float leftAngleChange = hBond->getLeftAngle() - hBond->iniOrientLeft;
//		if(Abs(leftAngleChange) > 0.00001){
//			log("report") <<"Left angle violation "<<id<<" between "<<atom1->getId()<<" and "<<atom2->getId()<<" in res "<<atom1->getResidue()->getName()<<atom1->getResidue()->getId()<<": " << leftAngleChange<<endl;
//		}

    id++;
  }
  log("report")<<endl;fflush(stdout);

}



/* Resample all sugar conformations in a segment of an RNA/DNA chain specified by startRes (included) and
 * endRes (not included) and reclose using the localRebuild method.
 * Returns the configuration and leaves the m_protein with that configuration set.
 */
Configuration*Molecule::resampleSugars(int startRes, int endRes, Configuration* cur, int aggression){
  //log("debugRas")<<"resampleSugars("<<startRes<<", "<<endRes<<" ... "<<aggression<<")"<<endl;
  vector<int> ignoreDOFs;
  vector<int> resetDOFs;
  vector<double> resetValues;
  vector<int> recloseDOFs;
  for(vector<Edge*>::iterator eit = m_spanning_tree->Edges.begin(); eit!=m_spanning_tree->Edges.end(); eit++){
    Edge* e = *eit;
    int res1 = e->getBond()->Atom1->getResidue()->getId();
    int res2 = e->getBond()->Atom2->getResidue()->getId();
    if( res1==(startRes-1) && res2==(startRes-1) && e->getBond()->Atom1->getName()=="C3'" && e->getBond()->Atom2->getName()=="O3'" ){
      recloseDOFs.push_back(e->DOF_id);
      if(aggression>=2){
        resetDOFs.push_back(e->DOF_id);
        resetValues.push_back(RandomAngleUniform(3.1415));
      }
    }
    if( res1==endRes && res2==endRes && (
        e->getBond()->Atom1->getName()=="P" ||
        e->getBond()->Atom1->getName()=="O5'"
    ) ){
      recloseDOFs.push_back(e->DOF_id);
      if(aggression>=2){
        resetDOFs.push_back(e->DOF_id);
        resetValues.push_back(RandomAngleUniform(3.1415));
      }
    }
    if( (res1>=startRes && res1<endRes) || (res2>=startRes && res2<endRes) ){
      if(	e->getBond()->Atom1->getName()=="P" || e->getBond()->Atom2->getName()=="P" ||
           e->getBond()->Atom1->getName()=="C5'" || e->getBond()->Atom2->getName()=="C5'" ||
           e->getBond()->Atom1->getName()=="C4'" || e->getBond()->Atom2->getName()=="C4'" ||
           e->getBond()->Atom1->getName()=="C3'" || e->getBond()->Atom2->getName()=="C3'" ||
           e->getBond()->Atom1->getName()=="O3'" || e->getBond()->Atom2->getName()=="O3'" ||
           e->getBond()->Atom1->getName()=="O5'" || e->getBond()->Atom2->getName()=="O5'" ){
        recloseDOFs.push_back(e->DOF_id);
        if(aggression>=2){
          resetDOFs.push_back(e->DOF_id);
          resetValues.push_back(RandomAngleUniform(3.1415));
        }
      }else if(
          e->getBond()->Atom1->getName()=="C2'" || e->getBond()->Atom2->getName()=="C2'" ||
          e->getBond()->Atom1->getName()=="C1'" || e->getBond()->Atom2->getName()=="C1'" ){
        if(aggression>=1){
          resetDOFs.push_back(e->DOF_id);
          resetValues.push_back(RandomAngleUniform(3.1415));
        }else{
          ignoreDOFs.push_back(e->DOF_id);
        }
      }
    }
  }

  for(map<unsigned int, RigidbodyGraphVertex*>::iterator vit = m_spanning_tree->Vertex_map.begin();vit!=m_spanning_tree->Vertex_map.end();vit++){
    RigidbodyGraphVertex* vtx = vit->second;
    if(vtx->isRibose){
      SugarVertex* v = reinterpret_cast<SugarVertex*>(vtx);
      int res = v->Rb_ptr->Atoms[0]->getResidue()->getId();
      //log("debug")<<"Res "<<res<<endl;;
      if( res>=startRes && res<endRes ){
        resetDOFs.push_back(v->DOF_id);

        double sampledTau = 0;
        //sampledTau = RandomAngleUniform(3.1415);
        if( Random01()<0.6 )	sampledTau = RandomNormalNPiPPi(-2.68, 0.2);
        else					sampledTau = RandomNormalNPiPPi(0.76, 0.3);
        double iniTau = acos(v->initTorsion/v->Amplitude)*v->initUp;
        resetValues.push_back(sampledTau-iniTau);
        //resetValues.push_back(sampledTau);
      }

    }
  }
  if(resetDOFs.empty()) {
    log("verbose")<<"Molecule::resampleSugars("<<startRes<<","<<endRes<<"..) called on non-nucleic acid chain. Exiting.";
    return cur;
  }

  return localRebuild(resetDOFs, resetValues, recloseDOFs, ignoreDOFs, cur);
}


/* Rebuilds a local a part of the molecule by resetting the DOFs indicated by resetDOFs to values indicated 
 * by resetValues and reclosing the region using the DOFs indicated by recloseDOFs. The same DOF can be
 * member of both resetDOFs and recloseDOFs. If the DOFs indicated by resetDOFs and recloseDOFs do not form
 * a connected subgraph the behavior of this method is unspecified.
 * Returns the configuration and leaves the m_protein with that configuration set.
 */
Configuration*Molecule::localRebuild(vector<int>& resetDOFs, vector<double>& resetValues, vector<int>& recloseDOFs, vector<int>& ignoreDOFs, Configuration* cur){
  //enableLogger("debugRebuild");

  double start_time, end_time, total_time;
  CTKTimer timer;
  start_time = timer.getTimeNow();

  //Configuration* ret = new Configuration(m_spanning_tree->m_numDOFs);

  //Configuration* ret = new Configuration(this);
  //ret->Copy(cur);
  Configuration* ret = cur->clone();

  //RestoreAtomPos(false);
  SetConfiguration(ret);
  ret->computeCycleJacobianAndNullSpace();

  //Find smallest connected subgraph that contains both resetDOFS and recloseDOFs (TODO: Approximate Steiner tree)
  vector<RigidbodyGraphVertex*> subVerts;
  vector<Edge*> subEdges;
  Edge* entry = NULL;
  for(vector<Edge*>::iterator eit = m_spanning_tree->Edges.begin(); eit!=m_spanning_tree->Edges.end(); eit++){
    Edge* e = *eit;
    if(	find(resetDOFs.begin(), 	resetDOFs.end(), 	e->DOF_id)!=resetDOFs.end() ||
         find(recloseDOFs.begin(), 	recloseDOFs.end(), 	e->DOF_id)!=recloseDOFs.end() ||
         find(ignoreDOFs.begin(), 	ignoreDOFs.end(), 	e->DOF_id)!=ignoreDOFs.end() ) {
      subEdges.push_back(e);
      if(find(subVerts.begin(), subVerts.end(), e->StartVertex)==subVerts.end()) 	subVerts.push_back(e->StartVertex);
      if(find(subVerts.begin(), subVerts.end(), e->EndVertex	)==subVerts.end()) 	subVerts.push_back(e->EndVertex);

      if(entry==NULL || entry->DOF_id>e->DOF_id) entry = e;
    }
  }
  for(map<unsigned int, RigidbodyGraphVertex*>::iterator vit = m_spanning_tree->Vertex_map.begin();vit!=m_spanning_tree->Vertex_map.end();vit++){
    RigidbodyGraphVertex* vtx = vit->second;
    if(vtx->isRibose){
      SugarVertex* v = reinterpret_cast<SugarVertex*>(vtx);
      if(	find(resetDOFs.begin(), resetDOFs.end(), v->DOF_id)!=resetDOFs.end() ||
           find(recloseDOFs.begin(), recloseDOFs.end(), v->DOF_id)!=recloseDOFs.end() ||
           find(ignoreDOFs.begin(), ignoreDOFs.end(), v->DOF_id)!=ignoreDOFs.end())
        subVerts.push_back(vtx);
    }
  }
  //log("debugRas")<<"SubEdges:"<<endl;
  //for(vector<Edge*>::iterator eit = subEdges.begin(); eit!=subEdges.end(); eit++){
  //    log("debugRas")<<"> "<<(*eit)<<endl;
  //}

  //int resi = entry->getBond()->Atom1->getResidue()->getId();

  //Collect edges with endpoints in subgraph and choose the covalent edge nearest to the root
  vector<Edge*> boundary;
  for(vector<Edge*>::iterator eit = m_spanning_tree->Edges.begin(); eit!=m_spanning_tree->Edges.end(); eit++){
    Edge* e = *eit;
    bool firstInSub = find(subVerts.begin(), subVerts.end(), e->StartVertex)!=subVerts.end();
    bool lastInSub = find(subVerts.begin(), subVerts.end(), e->EndVertex)!=subVerts.end();
    //if(!firstInSub &&  lastInSub) entry = e;

    if( firstInSub && !lastInSub && e->StartVertex!=entry->StartVertex) boundary.push_back(e);
  }
  for(vector<pair<Edge*,RigidbodyGraphVertex*> >::iterator it = m_spanning_tree->CycleAnchorEdges.begin(); it!=m_spanning_tree->CycleAnchorEdges.end(); it++){
    Edge* e = it->first;
    bool firstInSub = find(subVerts.begin(), subVerts.end(), e->StartVertex)!=subVerts.end();
    bool lastInSub = find(subVerts.begin(), subVerts.end(), e->EndVertex)!=subVerts.end();
    if(!firstInSub &&  lastInSub) boundary.push_back(e);
    if( firstInSub && !lastInSub) boundary.push_back(e);
  }
  //TODO: Also add hydrogen bonds to boundary
  //log("debugRas")<<"Boundary:"<<endl;
  //for(vector<Edge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
  //    log("debugRas")<<"> "<<(*eit)<<endl;
  //}
  //log("debugRas")<<"Entry: "<<entry<<endl;


  //Store the positions of endpoints and torsions of edges at boundary of subgraph
  vector<Coordinate> storedPositions;
  vector<double> storedTorsions;
  for(vector<Edge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
    Edge* e = *eit;
    if(e->getBond()!=NULL) {
      storedPositions.push_back(e->getBond()->Atom1->m_Position);
      storedPositions.push_back(e->getBond()->Atom2->m_Position);
      //log("debugRas")<<"Cylinder["<<e->getBond()->Atom1->m_Position[0]<<", "<<e->getBond()->Atom1->m_Position[1]<<", "<<e->getBond()->Atom1->m_Position[2]<<", ";
      //log("debugRas")<<e->getBond()->Atom2->m_Position[0]<<", "<<e->getBond()->Atom2->m_Position[1]<<", "<<e->getBond()->Atom2->m_Position[2]<<", 0.1, 0.9,0.9,0.2]"<<endl;

      //log("debugRebuild")<<"1: storedTorsion["<<storedTorsions.size()<<"] .. ";
      storedTorsions.push_back(e->getBond()->getTorsion());
    }
  }


  //Make sure only endpoints within subgraph and boundary can move
  for(vector<RigidbodyGraphVertex*>::iterator vit = subVerts.begin(); vit!=subVerts.end(); vit++){
    for(vector<Atom*>::iterator ait=(*vit)->Rb_ptr->Atoms.begin();ait!=(*vit)->Rb_ptr->Atoms.end();ait++){
      (*ait)->m_Position = (*ait)->m_referencePosition;
    }
  }
  for(vector<Edge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
    //(*eit)->getBond()->Atom1->m_bPositionModified = false;
    //(*eit)->getBond()->Atom2->m_bPositionModified = false;
    (*eit)->getBond()->Atom1->m_Position = entry->getBond()->Atom1->m_referencePosition;
    (*eit)->getBond()->Atom2->m_Position = entry->getBond()->Atom2->m_referencePosition;
  }
  //entry->getBond()->Atom1->m_bPositionModified = false;
  //entry->getBond()->Atom2->m_bPositionModified = false;
  entry->getBond()->Atom1->m_Position = entry->getBond()->Atom1->m_referencePosition;
  entry->getBond()->Atom2->m_Position = entry->getBond()->Atom2->m_referencePosition;
  //TODO: Also sugars

  //Change DOFs in resetDOFs to values indicated in resetValues
  for(int i=0;i<resetDOFs.size();i++){
    ret->m_dofs[resetDOFs[i]] = resetValues[i];
  }
  _SetConfiguration(ret, entry->StartVertex, subVerts);//, false);

  //stringstream fname;
  //fname<<"test_"<<resi<<"_0.pdb";
  //IO::writePdb(this, fname.str());
  double last_len = 1e10;
  double e_scaling = 0.1;
  int it;

  for(it=0;;it++){
    //log("debugRas")<<"Iteration "<<it<<endl;

    //Build end-effector vector: e
    gsl_vector* e = gsl_vector_alloc(storedPositions.size()*3);
    int c = 0;
    for(vector<Edge*>::iterator bit = boundary.begin(); bit!=boundary.end(); bit++){
      Edge* eBoundary = *bit;
      Vector3 v1 = (storedPositions[c*2+0])-(eBoundary->getBond()->Atom1->m_Position);
      Vector3 v2 = (storedPositions[c*2+1])-(eBoundary->getBond()->Atom2->m_Position);
      //log("debugRas")<<"Stored["<<c*2+0<<"]"
      gsl_vector_set(e, c*6+0, v1[0]*e_scaling);
      gsl_vector_set(e, c*6+1, v1[1]*e_scaling);
      gsl_vector_set(e, c*6+2, v1[2]*e_scaling);
      gsl_vector_set(e, c*6+3, v2[0]*e_scaling);
      gsl_vector_set(e, c*6+4, v2[1]*e_scaling);
      gsl_vector_set(e, c*6+5, v2[2]*e_scaling);

      c++;
    }
    //log("debugRas")<<"e"<<endl;
    //gsl_vector_log("debugRas")(e);
    double len = 0;
    for(int i=0;i<e->size;i++) len+=gsl_vector_get(e,i)*gsl_vector_get(e,i);
    len = sqrt(len)/e_scaling;
    //log("debugRas")<<"e length "<<len<<endl;

    //Stop if the end-effectors are sufficiently close to target
    if(len < 0.01) break;
    //Use smaller scaling if largest component increases
    if(len > 1.5*last_len){
      last_len = len;
      e_scaling *= 0.1;
      continue;
    }
    if( len>0.98*last_len ) {
      delete ret;
      return NULL;//break;
    }
    last_len = len;


    //Build Jacobian matrix
    gsl_matrix* J = gsl_matrix_calloc( storedPositions.size()*3, recloseDOFs.size() );
    for(int i=0;i<storedPositions.size()*3;i++)
      for(int j=0;j<recloseDOFs.size();j++)
        gsl_matrix_set(J, i,j, 0);
    c = 0;
    for(vector<Edge*>::iterator bit = boundary.begin(); bit!=boundary.end(); bit++){
      Edge* eBoundary = *bit;

      RigidbodyGraphVertex* v = eBoundary->StartVertex;
      Edge* e;

      do{ //find(subVerts.begin(), subVerts.end(), v)!=subVerts.end() ){
        if(v->Parent==NULL) break;
        e = v->Parent->findEdge(v);
        v = v->Parent;

        int dof = find(recloseDOFs.begin(), recloseDOFs.end(), e->DOF_id)-recloseDOFs.begin();
        //log("debugRas")<<"DOF: "<<dof<<endl;
        if(dof<recloseDOFs.size()){ // Continue if the edge is not in the recloseDOFs
          Vector3 v1 = ComputeJacobianEntry(e->getBond()->Atom1->m_Position, e->getBond()->Atom2->m_Position, eBoundary->getBond()->Atom1->m_Position);
          Vector3 v2 = ComputeJacobianEntry(e->getBond()->Atom1->m_Position, e->getBond()->Atom2->m_Position, eBoundary->getBond()->Atom2->m_Position);
          gsl_matrix_set(J, c*6+0, dof, v1[0]);
          gsl_matrix_set(J, c*6+1, dof, v1[1]);
          gsl_matrix_set(J, c*6+2, dof, v1[2]);
          gsl_matrix_set(J, c*6+3, dof, v2[0]);
          gsl_matrix_set(J, c*6+4, dof, v2[1]);
          gsl_matrix_set(J, c*6+5, dof, v2[2]);
        }

        SugarVertex* sv;
        if(v->isRibose) sv = reinterpret_cast<SugarVertex*>(v); else continue;
        dof = find(recloseDOFs.begin(), recloseDOFs.end(), sv->DOF_id)-recloseDOFs.begin();
        if(dof<recloseDOFs.size()){ // Continue if the sugar is not in the recloseDOFs
          Vector3 v1 = sv->computeJacobianEntry(e, ret->m_dofs, eBoundary->getBond()->Atom1->m_Position);
          Vector3 v2 = sv->computeJacobianEntry(e, ret->m_dofs, eBoundary->getBond()->Atom2->m_Position);
          gsl_matrix_set(J, c*6+0, dof, v1[0]);
          gsl_matrix_set(J, c*6+1, dof, v1[1]);
          gsl_matrix_set(J, c*6+2, dof, v1[2]);
          gsl_matrix_set(J, c*6+3, dof, v2[0]);
          gsl_matrix_set(J, c*6+4, dof, v2[1]);
          gsl_matrix_set(J, c*6+5, dof, v2[2]);
        }

      }while(e!=entry);

      c++;
    }
    //log("debugRas")<<"J"<<endl;
    //gsl_matrix_log("debugRas")(J);


    //Find pseudo-inverse
    MKLSVD svd(J);
    //svd.print();
    gsl_matrix* J_dag = svd.PseudoInverse();
    //log("debugRas")<<"J_dag"<<endl;
    //gsl_matrix_cout(J_dag);
    //gsl_matrix* J_dag = gsl_matrix_trans(J);

    //Apply it to e
    gsl_vector* dT = gsl_vector_alloc(J_dag->size1);
    gsl_blas_dgemv (CblasNoTrans, 1.0, J_dag, e, 0.0, dT);
    //gsl_vector* dT = gsl_matrix_vector_mul(J_dag, e);


    //log("debugRas")<<"dT"<<endl;
    //gsl_vector_cout(dT);
    for(int i=0;i<recloseDOFs.size();i++){
      ret->m_dofs[recloseDOFs[i]]+=gsl_vector_get(dT, i);
    }

    gsl_matrix_free(J_dag);
    gsl_vector_free(e);
    gsl_vector_free(dT);

    //Make sure only endpoints within subgraph and boundary can move
    //TODO: This needs to be updated for efficiency.
    //for(vector<Atom*>::iterator ait = atoms.begin(); ait!=atoms.end(); ait++){
    //	Atom* a = *ait;
    //	//a->m_bPositionModified = true;
    //}
    for(vector<RigidbodyGraphVertex*>::iterator vit = subVerts.begin(); vit!=subVerts.end(); vit++){
      for(vector<Atom*>::iterator ait=(*vit)->Rb_ptr->Atoms.begin();ait!=(*vit)->Rb_ptr->Atoms.end();ait++){
        //(*ait)->m_bPositionModified = false;
        (*ait)->m_Position = (*ait)->m_referencePosition;

      }
    }
    for(vector<Edge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
      //j(*eit)->getBond()->Atom1->m_bPositionModified = false;
      //(*eit)->getBond()->Atom2->m_bPositionModified = false;
      (*eit)->getBond()->Atom1->m_Position = (*eit)->getBond()->Atom1->m_referencePosition;
      (*eit)->getBond()->Atom2->m_Position = (*eit)->getBond()->Atom2->m_referencePosition;
    }
    //entry->getBond()->Atom1->m_bPositionModified = false;
    //entry->getBond()->Atom2->m_bPositionModified = false;
    entry->getBond()->Atom1->m_Position = entry->getBond()->Atom1->m_referencePosition;
    entry->getBond()->Atom2->m_Position = entry->getBond()->Atom2->m_referencePosition;

    //Move according to J_dag·e
    _SetConfiguration(ret, entry->StartVertex, subVerts);//, false);

    //stringstream ss; ss<<"test_"<<resi<<"_"<<(it+1)<<".pdb";
    //IO::writePdb(this, ss.str());
  }

  for(int i=0;i<boundary.size();i++){
    double oldTorsion = storedTorsions[i];
    double newTorsion = boundary[i]->getBond()->getTorsion();
    ret->m_dofs[boundary[i]->DOF_id] += (newTorsion - oldTorsion);

    //{
    //	SetConfiguration(cur);
    //	stringstream ss; ss<<"test_"<<resi<<"_init.pdb";
    //	IO::writePdb(this, ss.str());
    //}
    //{
    //	SetConfiguration(ret);
    //	stringstream ss; ss<<"test_"<<resi<<"_final2.pdb";
    //	IO::writePdb(this, ss.str());
    //}

  }


  end_time = timer.getTimeNow();
  total_time = end_time - start_time;
  log("debugRebuild")<<"Local rebuild took "<<total_time<<" secs. Stopped after "<<it<<" iterations."<<endl;
  //if(max_comp<0.01) 	log("debugRas")<<" Rebuild was succesful (no more than 0.05Å from target).";
  //else				log("debugRas")<<" Rebuild was NOT succesful.";
  ////if(max_comp>(0.98*last_max_comp)) log("debugRas")<<" Rebuild converged (changed less than 2%).";
  //log("debugRas")<<endl;

  return ret;
}
