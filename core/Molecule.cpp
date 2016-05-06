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

#include "DisjointSets.h"
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
#include <core/dofs/GlobalRotateDOF.h>
#include <core/dofs/GlobalTranslateDOF.h>
#include <core/dofs/TorsionDOF.h>
#include <cmath>

const double VDW_SIGMA = 0.2; // sigma = 0.2 kcal/mol
const double VDW_R0 = 3.5; // r_0 = 3.5A
const double VDW_R_MAX = 8; // only compute vdw energy if R_ab <= 8A
//const double STD2 = 1;

using namespace std;


Molecule::Molecule() {
  name_ = "UNKNOWN";
  Atom_pos_index = nullptr;
  backup_Atom_pos_index = nullptr;
  m_spanning_tree = nullptr;
  m_conf = nullptr;
  AtomJacobian1 = nullptr;
  AtomJacobian2 = nullptr;
  AtomJacobian3 = nullptr;
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
  if (Atom_pos_index!=nullptr)
    delete Atom_pos_index;
  if (backup_Atom_pos_index!=nullptr)
    delete backup_Atom_pos_index;
  // delete all rigid bodies
  for (map<unsigned int,Rigidbody*>::iterator it=Rigidbody_map_by_id.begin(); it!=Rigidbody_map_by_id.end(); ++it) {
    if ( it->second != nullptr )
      delete it->second;
  }

  for (list<Bond *>::iterator it=Cov_bonds.begin(); it != Cov_bonds.end(); ++it) {
    delete (*it);
  }
  for (list<Hbond *>::iterator it=H_bonds.begin(); it != H_bonds.end(); ++it) {
    delete (*it);
  }

  // delete m_spanning_tree
  if (m_spanning_tree!=nullptr)
    delete m_spanning_tree;

  // delete atom jacobians
  if (AtomJacobian1!=nullptr)
    gsl_matrix_free(AtomJacobian1);
  if (AtomJacobian2!=nullptr)
    gsl_matrix_free(AtomJacobian2);
  if (AtomJacobian3!=nullptr)
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
  if (chain==nullptr) { // this is a new chain
    chain = addChain(chainName);
  }
  Atom* ret = chain->addAtom(resName,resId, atomName, atomId, position);
  atoms.push_back(ret);
  return ret;
}

Chain*Molecule::getChain (const string& chain_name) {
  for(auto const& chain: chains)
    if(chain->getName()==chain_name) return chain;

  return nullptr;
}

void Molecule::printSummaryInfo() const {
  log() << "Molecule " << getName() << endl;
  for (auto const& chain: chains)
    chain->printSummaryInfo();
}

Chain*Molecule::addChain (const string& chainName) {
  //First check that chain name is not already used
  if(getChain(chainName) != nullptr) {
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
  if(chain!=nullptr){
    Residue* resi = chain->getResidue(resNum);
    if(resi!=nullptr){
      return resi->getAtom(name);
    }
  }
  return nullptr;
}

Atom* Molecule::getAtom (int atom_id) {
  for (auto const& a: atoms) {
    if(a->getId()==atom_id)
      return a;
  }
  return nullptr;
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
  if (Atom_pos_index!=nullptr)
    delete Atom_pos_index;

  Atom_pos_index = new Grid(this, SamplingOptions::getOptions()->collisionFactor);
}

void Molecule::backupAtomIndex () {
  if (backup_Atom_pos_index!=nullptr)
    delete backup_Atom_pos_index;
  backup_Atom_pos_index = Atom_pos_index->deepClone();
}

void Molecule::restoreAtomIndex () {
  if (Atom_pos_index!=nullptr)
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
 * Get a list of colliding atoms in the current protein configuration.
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
      cout<<"User-specified m_root id "<<bestId<<" out of bounds. Choosing standard Id."<<endl;
      return bestId;
    }
    cout<<"Choosing user-specified rigid body id "<<rootRBId<<" as m_root."<<endl;
    return (unsigned int)rootRBId;
  }
  else {
    if( target == nullptr ){
      cout<<"No target to determine best m_root, choosing standard m_root id 0"<<endl;
      return bestId;
    }
    //Check the rmsd between individual vertices and choose the closest pair as m_root
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
    cout<<"Choosing rigid body id "<<bestId<<" as m_root. Root RB rmsd: "<<bestSum<<endl;
    return bestId;
  }
}

void Molecule::buildSpanningTree() {
  unsigned int rootRBId = 0;
  //log() << "In buildSpanningTree" << endl;
  m_spanning_tree = new KinTree();

  // add all rigid bodies as vertices into Rigidbody_graph
  for (auto const& id_rb_pair: Rigidbody_map_by_id) {
    Rigidbody* rb = id_rb_pair.second;
    m_spanning_tree->addVertex(rb);
  }

  m_spanning_tree->m_root = m_spanning_tree->addVertex(nullptr);
  size_t numVertices = Rigidbody_map_by_id.size();

  list<KinEdge*> cycleEdges;
  int dofId = 0;
  std::set<KinVertex*> visitedVertices;

  //Initialize chain-roots but adding them to the queue and setting up edges from the super-m_root
  list<KinVertex*> queue;
  for(auto const& chain: chains) {
    //Add first vertex in chain to queue
    Residue* firstRes = chain->getResidues()[0];
    Atom* firstAtom = firstRes->getAtoms().front();
    KinVertex* firstVertex = firstAtom->getRigidbody()->getVertex();
    queue.push_back(firstVertex);

    //Connect to super-m_root
    KinVertex* v2 = m_spanning_tree->addVertex(nullptr);
    KinVertex* v3 = m_spanning_tree->addVertex(nullptr);
    KinVertex* v4 = m_spanning_tree->addVertex(nullptr);
    KinVertex* v5 = m_spanning_tree->addVertex(nullptr);
    KinVertex* v6 = m_spanning_tree->addVertex(nullptr);

    KinEdge* e1 = m_spanning_tree->addEdgeDirected(m_spanning_tree->m_root, v2, nullptr);
    KinEdge* e2 = m_spanning_tree->addEdgeDirected(v2, v3, nullptr);
    KinEdge* e3 = m_spanning_tree->addEdgeDirected(v3, v4, nullptr);
    KinEdge* e4 = m_spanning_tree->addEdgeDirected(v4, v5, nullptr);
    KinEdge* e5 = m_spanning_tree->addEdgeDirected(v5, v6, nullptr);
    KinEdge* e6 = m_spanning_tree->addEdgeDirected(v6, firstVertex, nullptr);

    e1->setDOF(new GlobalRotateDOF(e1,0));
    e2->setDOF(new GlobalRotateDOF(e2,1));
    e3->setDOF(new GlobalRotateDOF(e3,2));
    e4->setDOF(new GlobalTranslateDOF(e4,0));
    e5->setDOF(new GlobalTranslateDOF(e5,1));
    e6->setDOF(new GlobalTranslateDOF(e6,2));
    //e4->setDOF(new GlobalTranslateDOF(e4, 0));
    //e5->setDOF(new GlobalTranslateDOF(e5, 1));
    //e6->setDOF(new GlobalTranslateDOF(e6, 2));
  }

  //Perform breadth-first-search from all queue vertices and construct KinEdges
  //A variant of Prims algorithm is used to get the orientation of the tree correct in the first go
  while (!queue.empty()) {
    KinVertex* current_vertex = queue.front();
    queue.pop_front();
    visitedVertices.insert(current_vertex);

    for (auto const& bond: current_vertex->m_rigidbody->Bonds) {
      // Determine which other rb bond1 is connected to
      KinVertex* bonded_vertex = bond->Atom2->getRigidbody()->getVertex();
      if(bonded_vertex==current_vertex)
        bonded_vertex = bond->Atom1->getRigidbody()->getVertex();

      if(current_vertex==bonded_vertex || visitedVertices.count(bonded_vertex)>0)
        continue;

      if ( bond->isHbond() ) {
        // If it's a H-bond, it closes a cycle. Add it in CycleAnchorEdges.
        KinEdge *h_edge = new KinEdge(current_vertex,bonded_vertex,(Hbond*)bond);//TODO: Might be a problem. Idx changed from -1 to 0
        cycleEdges.push_back(h_edge);

      }else {
        // If it's a covalent bond, add it into the tree m_edges
        queue.push_back(bonded_vertex);
        KinEdge* edge = m_spanning_tree->addEdgeDirected(current_vertex,bonded_vertex,bond);
        edge->setDOF(new TorsionDOF(edge));
      }
    }
  } // end while

  m_spanning_tree->collectDOFs();

  // For each hbond KinEdge, find the lowest common ancestor (LCA) of its end-vertices and put all DOFs from the
  // end-points to the LCA into m_spanning_tree->m_cycleDOFs.
  for ( auto const& h_edge: cycleEdges) {
    KinVertex* lca = m_spanning_tree->findCommonAncestor(h_edge->StartVertex, h_edge->EndVertex);
    m_spanning_tree->CycleAnchorEdges.push_back( make_pair(h_edge,lca) );

    for(KinVertex* v=h_edge->StartVertex; v!=lca; v=v->m_parent){
      KinEdge* edge = v->m_parent->findEdge(v);
      m_spanning_tree->addCycleDOF(edge->getDOF());
    }

    for(KinVertex* v=h_edge->EndVertex; v!=lca; v=v->m_parent){
      KinEdge* edge = v->m_parent->findEdge(v);
      m_spanning_tree->addCycleDOF(edge->getDOF());
    }
  }


  /*
  int cycle_dof_id = 0;
  for ( auto const& h_edge: cycleEdges) {
    KinVertex *startv = h_edge->StartVertex;
    KinVertex *endv = h_edge->EndVertex;
    KinVertex *ancestor = m_spanning_tree->findCommonAncestor(startv,endv);
    m_spanning_tree->CycleAnchorEdges.push_back( make_pair(h_edge,ancestor) );
    KinVertex *vertex, *parent;
    for (vertex=startv; vertex!=ancestor; vertex=parent) {
      parent = vertex->m_parent;

      KinEdge *edge = parent->findEdge(vertex);

      if (edge->Cycle_DOF_id==-1) {
        edge->Cycle_DOF_id = cycle_dof_id;
        ++cycle_dof_id;
      }
    }
    for (vertex=endv; vertex!=ancestor; vertex=parent) {
      parent = vertex->m_parent;
      KinEdge *edge = parent->findEdge(vertex);
      if (edge->Cycle_DOF_id==-1) {
        edge->Cycle_DOF_id = cycle_dof_id;
        ++cycle_dof_id;
      }
    }
  }
  m_spanning_tree->m_numCycleDOFs = cycle_dof_id;

  m_Transformation = new RigidTransform [m_spanning_tree->getNumDOFs()];
  */

  /*
  //Create a sorted list of vertices for use in the euclideanGradient
  bool visited[Rigidbody_map_by_id.size()+1];
  for (unsigned int i=0; i<=Rigidbody_map_by_id.size(); ++i) {
    visited[i] = false;
  }
  KinEdge* currEdge = m_spanning_tree->Edges.back();
  KinVertex *currVertex = currEdge->EndVertex;
  int currId = currVertex->m_rigidbody->id();
  //int currId = currVertex->id;

  m_spanning_tree->m_sortedVertices.push_back(make_pair( currId, currVertex ));
  visited[currId] = true;

//  while( currVertex != m_spanning_tree->m_root){
  while( currVertex->m_rigidbody != nullptr){
    KinVertex *parent = currVertex->m_parent;

    if(parent->m_edges.size() == 1){
      currVertex = parent;
      currId = currVertex->m_rigidbody->id();
//      currId = currVertex->id;
      m_spanning_tree->m_sortedVertices.push_back(make_pair( currId,currVertex));
      visited[currId] = true;

    }else{ //multiple branches open
      currVertex = parent;
      int numEdgesLeft = currVertex->m_edges.size();
      while( numEdgesLeft != 0){
        for( KinEdge*& edge: currVertex->m_edges){
          if(visited[edge->EndVertex->m_rigidbody->id()] == true){
            numEdgesLeft--;
            continue;
          }else{
            currVertex = edge->EndVertex;
            //currVertex = eit->second->EndVertex;
            numEdgesLeft = currVertex->m_edges.size();
            break;
          }
        }
      }
      //currId = currVertex->id;
      currId = currVertex->m_rigidbody->id();
      m_spanning_tree->m_sortedVertices.push_back(make_pair(currId, currVertex) );
      visited[currId] = true;
    }
  }
   */

}

void Molecule::computeAtomJacobian (Atom* atom, gsl_matrix **j_addr) {
  if (*j_addr==nullptr) {
    *j_addr = gsl_matrix_calloc(3,totalDofNum());
  } else {
    gsl_matrix_set_all(*j_addr,0);
  }
  gsl_matrix* jacobian = *j_addr;
  KinVertex *vertex = atom->getRigidbody()->getVertex();
  atom->printSummaryInfo();
  while (vertex!=m_spanning_tree->m_root) {
    KinVertex *parent;
    if ( vertex->m_parent!=nullptr )
      parent = vertex->m_parent;
    else
      parent = vertex;
    KinEdge* edge = parent->findEdge(vertex);
    int dof_id = edge->getDOF()->getIndex();
    //Bond * bond_ptr = edge->getBond();
    //Coordinate bp1 = bond_ptr->Atom1->m_Position;
    //Coordinate bp2 = bond_ptr->Atom2->m_Position;
    //Math3D::Vector3 jacobian_entry = ComputeJacobianEntry(bp1,bp2,atom->m_Position);
    Math3D::Vector3 jacobian_entry = edge->getDOF()->getDerivative(atom->m_Position);
    gsl_matrix_set(jacobian,0,dof_id,jacobian_entry.x);
    gsl_matrix_set(jacobian,1,dof_id,jacobian_entry.y);
    gsl_matrix_set(jacobian,2,dof_id,jacobian_entry.z);
    vertex = parent;
  }
}


gsl_vector*Molecule::getEndEffectors(){
  gsl_vector* ret = gsl_vector_alloc(  (m_spanning_tree->CycleAnchorEdges).size()*6  );

  int i=0;
  for (vector< pair<KinEdge*,KinVertex*> >::iterator it=m_spanning_tree->CycleAnchorEdges.begin(); it!=m_spanning_tree->CycleAnchorEdges.end(); ++it) {
    // get end-effectors
    KinEdge* edge_ptr = it->first;
    Hbond * bond_ptr = (Hbond *)(edge_ptr->getBond());
    Math3D::Vector3 p1 = bond_ptr->Atom1->m_Position;
    Math3D::Vector3 p2 = bond_ptr->Atom2->m_Position;
    Math3D::Vector3 p1Diff = (p1-bond_ptr->getIdealHPoint());
    Math3D::Vector3 p2Diff = (p2-bond_ptr->getIdealAcceptorPoint());
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
    for (auto const& edge: m_spanning_tree->Edges){
      int dof_id = edge->getDOF()->getIndex();
      int cycle_dof_id = edge->getDOF()->getCycleIndex();
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
    for( auto const& edge:m_spanning_tree->Edges){
      int dof_id = edge->getDOF()->getIndex();
      int cycle_dof_id = edge->getDOF()->getCycleIndex();
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

void Molecule::translateReferencePositionsToRoot(Molecule * base)
{
  Coordinate& thisRoot = m_spanning_tree->m_root->m_rigidbody->Atoms[0]->m_referencePosition;
  Coordinate& thatRoot = base->m_spanning_tree->m_root->m_rigidbody->Atoms[0]->m_referencePosition;
  Math3D::Vector3 diff = thatRoot-thisRoot;

  for(auto const& atom: atoms){
    atom->m_referencePosition+=diff;
  }
}

void Molecule::RestoreAtomPos(){
  for (auto const& a: atoms)
    a->m_Position = a->m_referencePosition;

  m_conf = nullptr;
  restoreAtomIndex();
}

void Molecule::SetConfiguration(Configuration *q){
  if( m_conf==q) return;

  RestoreAtomPos();

  m_conf = q;
  _SetConfiguration(q);

  if(q->getGlobalTorsions() == nullptr){
    log("dominik")<<"Now updating global torsions"<<endl;
    q->updateGlobalTorsions();
  }
}



// set the positions of atoms at configuration q (according to the spanning tree)
void Molecule::_SetConfiguration(Configuration *q ){
  assert(this==q->getMolecule());

  for(size_t id=0 ; id<m_spanning_tree->getNumDOFs(); ++id){
    m_spanning_tree->getDOF(id)->setValue(q->m_dofs[id]);
  }

  KinVertex *root = m_spanning_tree->m_root;
  root->forwardPropagate();

  indexAtoms();
}


/**
  Set the positions of atoms only within subVerts so they correspond to configuration q.
  Assumes m_root is member of subVerts.
  */
void Molecule::_SetConfiguration(Configuration *q, KinVertex* root, vector<KinVertex*>& subVerts){//, bool usePosition2){
  m_conf = q;

  // assume the base vector is 0
//  Confvec2MatrixLocal(m_root, q, m_Transformation, subVerts);
  cerr<<"Molecule::_SetConfiguration(Configuration*,KinVertex*,vector<KinVertex*>) - This function is outdated. To be removed.";
  exit(-1);

  list<KinVertex *>queue;

  queue.push_back(root);
  while(queue.size()>0)
  {
    KinVertex* node = queue.front();
    queue.pop_front();

    //map<unsigned int,KinEdge*> m_children = node->Edges;

    //for (map<unsigned int,KinEdge*>::iterator edge_itr=m_children.begin(); edge_itr != m_children.end(); ++edge_itr){
    for (auto const& pEdge: node->m_edges){
      //KinEdge* pEdge = edge_itr->second;
      KinVertex* newNode = pEdge->EndVertex;

      //newNode->TransformAtomPosition(m_Transformation+pEdge->getDOF()->getIndex());//,usePosition2);

      if(find(subVerts.begin(), subVerts.end(), newNode)!=subVerts.end())
        queue.push_back(newNode);
    }
  }

  indexAtoms();
}

int Molecule::totalDofNum () const {
  if (m_spanning_tree==nullptr) {
    cerr << "Error: to get the total number of DOFs in the m_protein, you have to call Molecule::buildSpanningTree() first." << endl;
    exit(1);
  }
  return m_spanning_tree->getNumDOFs();
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
      if (AtomJacobian3==nullptr)
        AtomJacobian3 = gsl_matrix_calloc(3,totalDofNum());
      gsl_matrix_memcpy(AtomJacobian3,AtomJacobian1);
      gsl_matrix_sub(AtomJacobian3,AtomJacobian2); // AtomJacobian3 holds the result of substraction
      Math3D::Vector3 p12_v3 = atom1->m_Position - atom2->m_Position;
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

double Molecule::vdwEnergy (string collisionCheck) {// compute the total vdw energy, excluding covalently bonded atoms,
  // atoms in the same rigid body, and atoms we don't check clashes for

  double energy=0, d_12, ratio, vdw_r1, vdw_d12, epsilon1, epsilon_12;
  set< pair<Atom*,Atom*> >::iterator mit;
  // for each atom, look for it's neighbors.
  //For each such neighbor, compute U(R_ab)=epsilon_ij*(vdw_r12/r_12)^12-2*(VDW_R0/r_12)^6) and sum up everything.
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
    }
  }
  return energy;
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
      cout<<"Could not find specified hbond atoms in other protein: ";
      cout<<hBond->Hatom->getResidue()->getId()<<" "<<hBond->Hatom->getName()<<", "<<hBond->Acceptor->getResidue()->getId()<<" "<<hBond->Acceptor->getName();
      cout<<" Deleting to make hbond sets match!"<<endl;
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
  vector< pair<KinEdge*,KinVertex*> >::iterator pair_it;
  KinEdge *pEdge;
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
    double absViolation = std::fabs(distanceViolation);

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
  for(vector<KinEdge*>::iterator eit = m_spanning_tree->Edges.begin(); eit!=m_spanning_tree->Edges.end(); eit++){
    KinEdge* e = *eit;
    int res1 = e->getBond()->Atom1->getResidue()->getId();
    int res2 = e->getBond()->Atom2->getResidue()->getId();
    if( res1==(startRes-1) && res2==(startRes-1) && e->getBond()->Atom1->getName()=="C3'" && e->getBond()->Atom2->getName()=="O3'" ){
      recloseDOFs.push_back(e->getDOF()->getIndex());
      if(aggression>=2){
        resetDOFs.push_back(e->getDOF()->getIndex());
        resetValues.push_back(RandomAngleUniform(3.1415));
      }
    }
    if( res1==endRes && res2==endRes && (
        e->getBond()->Atom1->getName()=="P" ||
        e->getBond()->Atom1->getName()=="O5'"
    ) ){
      recloseDOFs.push_back(e->getDOF()->getIndex());
      if(aggression>=2){
        resetDOFs.push_back(e->getDOF()->getIndex());
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
        recloseDOFs.push_back(e->getDOF()->getIndex());
        if(aggression>=2){
          resetDOFs.push_back(e->getDOF()->getIndex());
          resetValues.push_back(RandomAngleUniform(3.1415));
        }
      }else if(
          e->getBond()->Atom1->getName()=="C2'" || e->getBond()->Atom2->getName()=="C2'" ||
          e->getBond()->Atom1->getName()=="C1'" || e->getBond()->Atom2->getName()=="C1'" ){
        if(aggression>=1){
          resetDOFs.push_back(e->getDOF()->getIndex());
          resetValues.push_back(RandomAngleUniform(3.1415));
        }else{
          ignoreDOFs.push_back(e->getDOF()->getIndex());
        }
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

  //Configuration* ret = new Configuration(m_spanning_tree->getNumDOFs());

  //Configuration* ret = new Configuration(this);
  //ret->Copy(cur);
  Configuration* ret = cur->clone();

  //RestoreAtomPos(false);
  SetConfiguration(ret);
  ret->computeCycleJacobianAndNullSpace();

  //Find smallest connected subgraph that contains both resetDOFS and recloseDOFs (TODO: Approximate Steiner tree)
  vector<KinVertex*> subVerts;
  vector<KinEdge*> subEdges;
  KinEdge* entry = nullptr;
  for(vector<KinEdge*>::iterator eit = m_spanning_tree->Edges.begin(); eit!=m_spanning_tree->Edges.end(); eit++){
    KinEdge* e = *eit;
    if(	find(resetDOFs.begin(), 	resetDOFs.end(), 	e->getDOF()->getIndex())!=resetDOFs.end() ||
         find(recloseDOFs.begin(), 	recloseDOFs.end(), 	e->getDOF()->getIndex())!=recloseDOFs.end() ||
         find(ignoreDOFs.begin(), 	ignoreDOFs.end(), 	e->getDOF()->getIndex())!=ignoreDOFs.end() ) {
      subEdges.push_back(e);
      if(find(subVerts.begin(), subVerts.end(), e->StartVertex)==subVerts.end()) 	subVerts.push_back(e->StartVertex);
      if(find(subVerts.begin(), subVerts.end(), e->EndVertex	)==subVerts.end()) 	subVerts.push_back(e->EndVertex);

      if(entry==nullptr || entry->getDOF()->getIndex()>e->getDOF()->getIndex()) entry = e;
    }
  }
  //log("debugRas")<<"SubEdges:"<<endl;
  //for(vector<KinEdge*>::iterator eit = subEdges.begin(); eit!=subEdges.end(); eit++){
  //    log("debugRas")<<"> "<<(*eit)<<endl;
  //}

  //int resi = entry->getBond()->Atom1->getResidue()->getId();

  //Collect m_edges with endpoints in subgraph and choose the covalent edge nearest to the m_root
  vector<KinEdge*> boundary;
  for(vector<KinEdge*>::iterator eit = m_spanning_tree->Edges.begin(); eit!=m_spanning_tree->Edges.end(); eit++){
    KinEdge* e = *eit;
    bool firstInSub = find(subVerts.begin(), subVerts.end(), e->StartVertex)!=subVerts.end();
    bool lastInSub = find(subVerts.begin(), subVerts.end(), e->EndVertex)!=subVerts.end();
    //if(!firstInSub &&  lastInSub) entry = e;

    if( firstInSub && !lastInSub && e->StartVertex!=entry->StartVertex) boundary.push_back(e);
  }
  for(vector<pair<KinEdge*,KinVertex*> >::iterator it = m_spanning_tree->CycleAnchorEdges.begin(); it!=m_spanning_tree->CycleAnchorEdges.end(); it++){
    KinEdge* e = it->first;
    bool firstInSub = find(subVerts.begin(), subVerts.end(), e->StartVertex)!=subVerts.end();
    bool lastInSub = find(subVerts.begin(), subVerts.end(), e->EndVertex)!=subVerts.end();
    if(!firstInSub &&  lastInSub) boundary.push_back(e);
    if( firstInSub && !lastInSub) boundary.push_back(e);
  }
  //TODO: Also add hydrogen bonds to boundary
  //log("debugRas")<<"Boundary:"<<endl;
  //for(vector<KinEdge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
  //    log("debugRas")<<"> "<<(*eit)<<endl;
  //}
  //log("debugRas")<<"Entry: "<<entry<<endl;


  //Store the positions of endpoints and torsions of m_edges at boundary of subgraph
  vector<Coordinate> storedPositions;
  vector<double> storedTorsions;
  for(vector<KinEdge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
    KinEdge* e = *eit;
    if(e->getBond()!=nullptr) {
      storedPositions.push_back(e->getBond()->Atom1->m_Position);
      storedPositions.push_back(e->getBond()->Atom2->m_Position);
      //log("debugRas")<<"Cylinder["<<e->getBond()->Atom1->m_Position[0]<<", "<<e->getBond()->Atom1->m_Position[1]<<", "<<e->getBond()->Atom1->m_Position[2]<<", ";
      //log("debugRas")<<e->getBond()->Atom2->m_Position[0]<<", "<<e->getBond()->Atom2->m_Position[1]<<", "<<e->getBond()->Atom2->m_Position[2]<<", 0.1, 0.9,0.9,0.2]"<<endl;

      //log("debugRebuild")<<"1: storedTorsion["<<storedTorsions.size()<<"] .. ";
      storedTorsions.push_back(e->getBond()->getTorsion());
    }
  }


  //Make sure only endpoints within subgraph and boundary can move
  for(vector<KinVertex*>::iterator vit = subVerts.begin(); vit!=subVerts.end(); vit++){
    for(vector<Atom*>::iterator ait=(*vit)->m_rigidbody->Atoms.begin();ait!=(*vit)->m_rigidbody->Atoms.end();ait++){
      (*ait)->m_Position = (*ait)->m_referencePosition;
    }
  }
  for(vector<KinEdge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
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
    for(vector<KinEdge*>::iterator bit = boundary.begin(); bit!=boundary.end(); bit++){
      KinEdge* eBoundary = *bit;
      Math3D::Vector3 v1 = (storedPositions[c*2+0])-(eBoundary->getBond()->Atom1->m_Position);
      Math3D::Vector3 v2 = (storedPositions[c*2+1])-(eBoundary->getBond()->Atom2->m_Position);
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
      return nullptr;//break;
    }
    last_len = len;


    //Build Jacobian matrix
    gsl_matrix* J = gsl_matrix_calloc( storedPositions.size()*3, recloseDOFs.size() );
    for(int i=0;i<storedPositions.size()*3;i++)
      for(int j=0;j<recloseDOFs.size();j++)
        gsl_matrix_set(J, i,j, 0);
    c = 0;
    for(vector<KinEdge*>::iterator bit = boundary.begin(); bit!=boundary.end(); bit++){
      KinEdge* eBoundary = *bit;

      KinVertex* v = eBoundary->StartVertex;
      KinEdge* e;

      do{ //find(subVerts.begin(), subVerts.end(), v)!=subVerts.end() ){
        if(v->m_parent==nullptr) break;
        e = v->m_parent->findEdge(v);
        v = v->m_parent;

        int dof = find(recloseDOFs.begin(), recloseDOFs.end(), e->getDOF()->getIndex())-recloseDOFs.begin();
        //log("debugRas")<<"DOF: "<<dof<<endl;
        if(dof<recloseDOFs.size()){ // Continue if the edge is not in the recloseDOFs
          //Math3D::Vector3 v1 = ComputeJacobianEntry(e->getBond()->Atom1->m_Position, e->getBond()->Atom2->m_Position, eBoundary->getBond()->Atom1->m_Position);
          //Math3D::Vector3 v2 = ComputeJacobianEntry(e->getBond()->Atom1->m_Position, e->getBond()->Atom2->m_Position, eBoundary->getBond()->Atom2->m_Position);
          Math3D::Vector3 v1 = e->getDOF()->getDerivative(eBoundary->getBond()->Atom1->m_Position);
          Math3D::Vector3 v2 = e->getDOF()->getDerivative(eBoundary->getBond()->Atom2->m_Position);
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
    for(vector<KinVertex*>::iterator vit = subVerts.begin(); vit!=subVerts.end(); vit++){
      for(vector<Atom*>::iterator ait=(*vit)->m_rigidbody->Atoms.begin();ait!=(*vit)->m_rigidbody->Atoms.end();ait++){
        //(*ait)->m_bPositionModified = false;
        (*ait)->m_Position = (*ait)->m_referencePosition;

      }
    }
    for(vector<KinEdge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
      (*eit)->getBond()->Atom1->m_Position = (*eit)->getBond()->Atom1->m_referencePosition;
      (*eit)->getBond()->Atom2->m_Position = (*eit)->getBond()->Atom2->m_referencePosition;
    }
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
    ret->m_dofs[boundary[i]->getDOF()->getIndex()] += (newTorsion - oldTorsion);

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
