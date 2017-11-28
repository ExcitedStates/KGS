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
#include "HBond.h"
#include "Grid.h"
#include "CTKTimer.h"
#include "metrics/RMSD.h"
#include "Transformation.h"
#include "math/SVD.h"
#include "IO.h"
#include "Logger.h"
#include <iostream>
#include <iomanip> 		// std::setprecision
#include <fstream>
#include <sstream>
#include <math/SVDMKL.h>
#include <set>
#include <core/dofs/GlobalRotateDOF.h>
#include <core/dofs/GlobalTranslateDOF.h>
#include <core/dofs/TorsionDOF.h>
#include <cmath>
#include <math/gsl_helpers.h>
#include <math/MathUtility.h>

const double VDW_SIGMA = 0.2; // sigma = 0.2 kcal/mol
const double VDW_R0 = 3.5; // r_0 = 3.5A
const double VDW_R_MAX = 8; // only compute vdw energy if R_ab <= 8A
//const double STD2 = 1;

using namespace std;


Molecule::Molecule():
  m_name("UNKNOWN"),
  m_grid(nullptr),
  m_spanningTree(nullptr),
  m_conf(nullptr),
  m_collisionFactor(1.0)
{
}

Molecule::~Molecule() {
  // delete all atoms
  for (vector<Atom *>::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it) {
    delete (*it);
  }

  // delete all chains
  for (auto const &chain: m_chains) {
    delete chain;
  }

  // delete m_grid and m_backupGrid
  if (m_grid != nullptr){
    delete m_grid;
    m_grid = nullptr;
  }

  // delete all rigid bodies
  for (map<unsigned int,Rigidbody*>::iterator it=m_rigidBodyMap.begin(); it!=m_rigidBodyMap.end(); ++it) {
    if ( it->second != nullptr )
      delete it->second;
  }

  for (list<Bond *>::iterator it=m_covBonds.begin(); it != m_covBonds.end(); ++it) {
    delete (*it);
  }
  for (list<Hbond *>::iterator it=m_hBonds.begin(); it != m_hBonds.end(); ++it) {
    delete (*it);
  }
  for (list<DBond *>::iterator it=m_dBonds.begin(); it != m_dBonds.end(); ++it) {
    delete (*it);
  }

  // delete m_spanningTree
  if (m_spanningTree!=nullptr)
    delete m_spanningTree;

}

void Molecule::setName (const string& name) {
  m_name = name;
}

string Molecule::getName () const {
  return m_name;
}

Atom* Molecule::addAtom(
    const bool& hetatm,
    const std::string& chainName,
    const std::string& resName,
    const int& resId,
    const std::string& atomName,
    const int& atomId,
    const Coordinate& position )
{
  Chain* chain = getChain(chainName);
  if (chain==nullptr)  // this is a new chain
    chain = addChain(chainName);

  Atom* ret = chain->addAtom(hetatm, resName,resId, atomName, atomId, position);
  m_atoms.push_back(ret);

  return ret;
}

//Bond* Molecule::addCovBond (Residue* res1, Residue* res2, const string& atom_name1, const string& atom_name2) {
Bond* Molecule::addCovBond (Atom* atom1, Atom* atom2) {
//  Atom* atom1 = res1->getAtom(atom_name1);
//  Atom* atom2 = res2->getAtom(atom_name2);
  assert(atom1!=nullptr);
  assert(atom2!=nullptr);

//  if (atom1!=nullptr && atom2!=nullptr) {
    if(atom1->getId()>atom2->getId())
      std::swap(atom1,atom2);

    if(std::find(atom1->Cov_neighbor_list.begin(), atom1->Cov_neighbor_list.end(), atom2)==atom1->Cov_neighbor_list.end()) {
      Bond *new_cb = new Bond(atom1, atom2, "COV");
      addCovBond(new_cb);
      return new_cb;
    }
    return nullptr;
//  }else{
//    cerr << "Molecule::addCovBond["<< __LINE__<<"] - Atom not found so bond cant be created: ";
//    cerr << res1->getId() << "/" << atom_name1 << " - ";
//    cerr << res2->getId() << "/" << atom_name2 << endl;
//    exit(-1);
//  }
}

Chain*Molecule::getChain (const string& chain_name) const{
  for(auto const& chain: m_chains)
    if(chain->getName()==chain_name) return chain;

  return nullptr;
}

Chain*Molecule::addChain (const string& chainName) {
  //First check that chain name is not already used
  if(getChain(chainName) != nullptr) {
    cerr << "Molecule::addChain - chain " << chainName << " already exists." << endl;
    exit(-1);
  }

  Chain* chain = new Chain(chainName, this);
  m_chains.push_back(chain);
  return chain;
}

int Molecule::getMaxResidueNumber(){
  int ret= -1000000;
  for(auto const& a: m_atoms){
    if( a->getResidue()->getId() > ret )
      ret = a->getResidue()->getId();
  }
  return ret;

}

int Molecule::getMinResidueNumber(){
  int ret= 1000000;
  for(vector<Atom*>::iterator ait = m_atoms.begin(); ait != m_atoms.end(); ait++){
    Atom* a = *ait;
    if( a->getResidue()->getId() < ret )
      ret = a->getResidue()->getId();
  }
  return ret;

}

/** Gets the atom specified by a residue number and a name. */
Atom* Molecule::getAtom(const string& chainName, const int& resNum, const string& name) const{
  Chain* chain = this->getChain(chainName);
  if(chain!=nullptr){
    Residue* resi = chain->getResidue(resNum);
    if(resi!=nullptr){
      return resi->getAtom(name);
    }
  }
  return nullptr;
}

Atom* Molecule::getAtom (int atom_id) const{
  for (auto const& a: m_atoms) {
    if(a->getId()==atom_id)
      return a;
  }
  return nullptr;
}

const std::vector<Atom*>& Molecule::getAtoms() const {
  return m_atoms;
}

std::vector<Atom*>& Molecule::getAtoms(){
  return m_atoms;
}

const std::list<Bond*>& Molecule::getCovBonds() const {
  return m_covBonds;
}

std::list<Bond*>& Molecule::getCovBonds(){
  return m_covBonds;
}

const std::list<Hbond*>& Molecule::getHBonds() const {
  return m_hBonds;
}

const std::list<DBond*>& Molecule::getDBonds() const {
  return m_dBonds;
}

std::list<Hbond*>& Molecule::getHBonds(){
  return m_hBonds;
}

std::list<DBond*>& Molecule::getDBonds(){
  return m_dBonds;
}

int Molecule::size() const {
  return m_atoms.size();
}

Grid* Molecule::getGrid() {
  if(m_grid==nullptr)
    indexAtoms();
  return m_grid;
}

std::set< std::pair<Atom*,Atom*> >& Molecule::getInitialCollisions()
{
  return m_initialCollisions;
};

void Molecule::indexAtoms () {
  // m_grid is the current indexing
  if (m_grid != nullptr){
    delete m_grid;
    m_grid = nullptr;
  }

  if(m_conf==nullptr) {
    restoreAtomPos();
  }

  m_grid = new Grid(this, m_collisionFactor);
}

void Molecule::setCollisionFactor(double collisionFactor)
{
  m_collisionFactor = collisionFactor;

  getGrid()->setCollisionFactor(collisionFactor);

  //Recompute initial collisions
  m_initialCollisions.clear();
  m_initialCollisions = getAllCollisions();
  log("debug")<<"Molecule::setCollisionFactor("<<collisionFactor<<") - Initial collisions: "<<m_initialCollisions.size()<<endl;
}

//void Molecule::backupAtomIndex () {
//  if (m_backupGrid!=nullptr)
//    delete m_backupGrid;
//
//  m_backupGrid = m_grid->deepClone();
//}

//void Molecule::restoreAtomIndex () {
//  if (m_grid!=nullptr)
//    delete m_grid;
//  if (m_backupGrid==nullptr)
//    backupAtomIndex();
//
//  m_grid = m_backupGrid->deepClone();
//}


bool Molecule::inCollision (string collisionCheckAtoms ) {

  Grid* grid = getGrid();

  for (vector<Atom*>::const_iterator itr= m_atoms.begin(); itr != m_atoms.end(); ++itr)
    if( (*itr)->isCollisionCheckAtom(collisionCheckAtoms ) )
    if ( grid->inCollision(*itr, m_initialCollisions, collisionCheckAtoms)  )
      return true;
  return false;
}

double Molecule::minCollisionFactor (string collisionCheckAtoms) {
  double minCollFactor = 10000;
  Grid* grid = getGrid();
  for (vector<Atom*>::const_iterator itr=m_atoms.begin(); itr!=m_atoms.end(); ++itr){
    if( (*itr)->isCollisionCheckAtom(collisionCheckAtoms ) ){
      double factor = grid->minFactorWithoutCollision(*itr, m_initialCollisions, collisionCheckAtoms);
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
std::set< pair<Atom*,Atom*> > Molecule::getAllCollisions (std::string collisionCheckAtoms ){
  if(m_conf==nullptr) {
    cerr << "Molecule::getAllCollisions - No configuration set" << endl;
    throw "Molecule::getAllCollisions - No configuration set";
  }

  Grid* grid = getGrid();

  set< pair<Atom*,Atom*>> collisions;
  for (auto const& atom: m_atoms) {
    if( atom->isCollisionCheckAtom( collisionCheckAtoms ) ) {
      vector<Atom *> colliding_atoms = grid->getAllCollisions(atom, m_initialCollisions, collisionCheckAtoms);
      for (auto const& colliding_atom : colliding_atoms) {
        pair<Atom *, Atom *> collision_pair = make_pair(atom, colliding_atom);
        collisions.insert(collision_pair);
      }
    }
  }
  return collisions;
}
//---------------------------------------------------------

void Molecule::printAllCollisions () {
  for (auto const& atom_pair: getAllCollisions()){
    log() << atom_pair.first->getId() << " " << atom_pair.second->getId() << endl;
  }
}

void Molecule::addCovBond (Bond * bond) {
  //log()<<"Molecule::addCovBond("<<bond<<")"<<endl;
  m_covBonds.push_back(bond);
  bond->m_atom1->addCovBond(bond);
  bond->m_atom2->addCovBond(bond);
}

void Molecule::addHbond (Hbond * hb) {
  m_hBonds.push_back(hb);
  hb->m_atom1->addHbond(hb);
  hb->m_atom2->addHbond(hb);
}
void Molecule::addDBond (DBond * db) {
  m_dBonds.push_back(db);
//  db->m_atom1->addDBond(db);
//  db->m_atom2->addDBond(db);
}

std::vector<int> Molecule::findBestRigidBodyMatch(std::vector<int> chainRoots, Molecule * target){

  int standardId = 1; //standard choice for root rigid body
  int maxId = m_atoms.size();
  int chainCount = 0;
  int rootId = 1;

  std::vector<int> bestRootIDs;
  for (const auto chain: this->getChains() ){
    rootId = chainRoots[chainCount];
    if (rootId >= 1) {//user-specified, or standard choice of atom id == 0
      if (rootId > maxId) {
        log("samplingStatus") << "User-specified root id " << rootId << " for chain " << chainCount <<" out of bounds. Choosing standard id "<<standardId << endl;
        bestRootIDs.push_back(standardId);
      }
      else {
        log("samplingStatus") << "Choosing user-specified atom id " << rootId << " as root for chain " << chainCount << endl;
        bestRootIDs.push_back(rootId);
      }
    }
    else { // root specified as <=0 --> therefore align and pick best rigid body
      if (target == nullptr) {
        log("samplingStatus") << "No target to determine best root, choosing standard root id "<<standardId << endl;
        bestRootIDs.push_back(standardId);
      }
      //Check the rmsd between individual vertices and choose the closest pair as m_root
      double bestSum = 9999;
      int bestId = 1;
      for (map<unsigned int, Rigidbody *>::iterator rbit = m_rigidBodyMap.begin();
           rbit != m_rigidBodyMap.end(); ++rbit) {

        vector<Atom *> *atomsRMSD = &(rbit->second->Atoms);
        if (atomsRMSD->front()->getResidue()->getChain() != chain) {
          continue;//skip if not this chain
        }
        double sum = 0;
        bool allAtoms = true;
        Coordinate c1, c2;
        Atom* rootAtom = nullptr;

        //loop through the rigid body atoms
        for (vector<Atom *>::iterator ait = atomsRMSD->begin(); ait != atomsRMSD->end(); ++ait) {
          Atom* atom = (*ait);
          string name = atom->getName();
          string chainName = atom->getResidue()->getChain()->getName();
          int resId = atom->getResidue()->getId();
          if (!rootAtom && atom->isHeavyAtom() && atom->isBackboneAtom() ){
            rootAtom = atom; //choose first heavy backbone atom as root
          }
          Atom *a2 = target->getAtom(chainName, resId, name);
          if (atom && a2) {
            c1 = atom->m_position;
            c2 = a2->m_position;
            sum += c1.distanceSquared(c2);
          } else {//only allow the rb's where all atoms are available in both proteins
            allAtoms = false;
          }
        }

        if (allAtoms && rootAtom ) {
          sum = sqrt(sum / (atomsRMSD->size()));
          if (sum <= bestSum) {
            bestSum = sum;
            bestId = rootAtom->getId();
          }
        }
      }
      log("samplingStatus") << this->getName()<<": Choosing atom id " << bestId << " as root. Root RB rmsd: " << bestSum << endl;
      bestRootIDs.push_back(bestId);
    }
    chainCount++;
  }
  return bestRootIDs;
}

void Molecule::initializeTree(Selection& movingResidues,double collisionFactor, const std::vector<int>& roots,Molecule* target) {
  this->sortHbonds();
  this->buildRigidBodies(movingResidues); //Necessary to do before building spanning tree
  ///adapted to multi-chain root choices based on input structures
  vector<int> bestRoots = this->findBestRigidBodyMatch(roots,target);
  this->buildSpanningTree(bestRoots); //Necessary before conformations are defined
  this->setConfiguration(new Configuration(this));
  this->setCollisionFactor(collisionFactor); //Sets the initial collisions //ToDo: Do we really need this here? Better when we know collision factor
}


void Molecule::buildRigidBodies(Selection& movingResidues, int collapseLevel) {
  //Create disjoint set
  DisjointSets ds(getAtoms()[size() - 1]->getId() + 1); //Assumes the last atom has the highest id.

  //For all pairs of atoms in residues not in movingResidues call Union (rigidifies everything not in movingResidues)
  for(auto const& atom: m_atoms) {
    if(movingResidues.inSelection(atom)) continue;

    for(auto const& neighbor: atom->Cov_neighbor_list){
      if(movingResidues.inSelection(neighbor)) continue;

      //Both atom and neighbor are outside movingResidues .. rigidify them
      ds.Union(atom->getId(), neighbor->getId());
      log("debug") << "IO::buildRigidBodies["<< __LINE__<<"] - Rigidifying bond " << atom->getId() << " - ";
      log("debug") << neighbor->getId() << " as they're not in residueNetwork" << endl;
    }
  }
  //for(auto const& chain: chains) {
  //  Atom* lastAtom = nullptr;
  //  for (auto const& res: chain->getResidues()) {
  //    if(movingResidues.inSelection(res)) {
  //      log("debug") << "IO::buildRigidBodies["<< __LINE__<<"] - Not rigidifying residue " << res->getId() << endl;
  //      continue; //Skip residue if its in movingResidues
  //    }else {
  //      log("debug") << "IO::buildRigidBodies["<< __LINE__<<"] - Rigidifying residue " << res->getId() << endl;
  //    }

  //    for (Atom *const &res_atom: res->getAtoms()){
  //      if(lastAtom==nullptr) {
  //        lastAtom = res_atom;
  //        continue;
  //      }
  //      ds.Union(lastAtom->getId(), res_atom->getId());
  //      log("debug") << "IO::buildRigidBodies["<< __LINE__<<"] - Joining " << lastAtom->getId() << " - " << res_atom->getId() << endl;
  //    }
  //  }

  //For each atom, a1, with exactly one cov neighbor and not participating in an hbond, a2, call Union(a1,a2)
  for (int i=0;i<size();i++){
    Atom* atom = getAtoms()[i];
    if(atom->Cov_neighbor_list.size()==1 && atom->Hbond_neighbor_list.size()==0){
      ds.Union(atom->getId(), atom->Cov_neighbor_list[0]->getId());
      log("debug") << "IO::buildRigidBodies["<< __LINE__<<"] - Joining " << atom->getId() << " - " << atom->Cov_neighbor_list[0]->getId() << endl;
      //cout<<"Only one neighbor: "<<atom->getName()<<" "<<atom->getId()<<" - "<<atom->Cov_neighbor_list[0]->getName()<<" "<<atom->Cov_neighbor_list[0]->getId()<<endl;
    }
  }


  int count=0;
  //For each fixed bond (a1,a2) call Union(a1,a2)
  for (auto const& bond: getCovBonds()){
    if( bond->Bars == 6 || bond->rigidified){//This is fixed in the Bond -> isLocked and from rigidity analysis
      count++;
      ds.Union(bond->m_atom1->getId(), bond->m_atom2->getId());
      log("debug") << "IO::buildRigidBodies["<< __LINE__<<"] - Joining " << bond->m_atom1->getId() << " - " << bond->m_atom2->getId() << endl;
      continue;
    }
  }
  log("debug")<<"IO::buildRigidBodies["<< __LINE__<<"] - Rigidified "<<count<<" covalent bonds."<<endl;
  count=0;
  //For each fixed bond (a1,a2) call Union(a1,a2)
  if(collapseLevel == 2) {
    for( auto const &bond: getHBonds()) {
      if( bond->Bars == 6 || bond->rigidified ) {//This is fixed in the Bond -> isLocked and from rigidity analysis
        count++;
        ds.Union(bond->m_atom1->getId(), bond->m_atom2->getId());
        log("debug") << "IO::readRigidbody[" << __LINE__ << "] - Joining " << bond->m_atom1->getId() << " - "
                     << bond->m_atom2->getId() << endl;
        continue;
      }
    }
  }
  log("debug")<<"IO::buildRigidBodies["<< __LINE__<<"] - Rigidified "<<count<<" hydrogen bonds."<<endl;

  int c=0;
  map<int,int> idMap;//Maps atom id's to rigid body id's for use in the DS structure.

  //Map the set-ID's to RB-ID's and add bonded atoms to RBs.
  for (int i=0;i<size();i++){
    Atom* atom = getAtoms()[i];

    //Map the set-id to the RB-id
    int set_id = ds.FindSet(atom->getId());
    int body_id;
    if(idMap.find(set_id)!=idMap.end()) body_id = idMap.find(set_id)->second;
    else {
      body_id = c++;
      idMap.insert( make_pair(set_id, body_id) );
    }
    //If the set containing a1 is not a rigid body: create one
    if ( m_rigidBodyMap.find(body_id)==m_rigidBodyMap.end() ) {
      Rigidbody* new_rb = new Rigidbody(body_id);
      m_rigidBodyMap.insert( make_pair(body_id,new_rb) );
    }
    Rigidbody* rb = m_rigidBodyMap.find(body_id)->second;
    if (!rb->containsAtom(atom)) rb->addAtom(atom);

  }

  //Delete small RBs and sort atoms within each RB
  map<unsigned int, Rigidbody*>::iterator it = m_rigidBodyMap.begin();
  while ( it!=m_rigidBodyMap.end() ) {
    Rigidbody *rb = it->second;
    if ( rb->Atoms.size() <= 0 ) {
      cerr << "Error: rigid body " << rb->id() << " has no atoms." << endl;
      m_rigidBodyMap.erase(it++);
      delete rb;
    }
    else { // sort atom ids
#ifndef WIN32 // Liangjun Zhang's tmp code
      vector<Atom*>::iterator sit = rb->Atoms.begin();
      vector<Atom*>::iterator eit = rb->Atoms.end();

      sort(sit,eit,Atom::compare);
#endif
      // Determine if Rigidbody is on Mainchain
//      rb->setMainchainRb();
      ++it;
    }
  }

  //Store bonds in rigid bodies
  for (auto const& bond: getCovBonds()){
    int setId1 = ds.FindSet(bond->m_atom1->getId());
    m_rigidBodyMap.find( idMap.find(setId1)->second )->second->addBond(bond);
    int setId2 = ds.FindSet(bond->m_atom2->getId());
    if(setId1!=setId2)
      m_rigidBodyMap.find( idMap.find(setId2)->second )->second->addBond(bond);
  }
  for (auto const& bond: getHBonds()) {
    int setId1 = ds.FindSet(bond->m_atom1->getId());
    m_rigidBodyMap.find( idMap.find(setId1)->second )->second->addBond(bond);
    int setId2 = ds.FindSet(bond->m_atom2->getId());
    if(setId1!=setId2)
      m_rigidBodyMap.find( idMap.find(setId2)->second )->second->addBond(bond);
  }
  for (auto const& bond: getDBonds()) {
    int setId1 = ds.FindSet(bond->m_atom1->getId());
    m_rigidBodyMap.find( idMap.find(setId1)->second )->second->addBond(bond);
    int setId2 = ds.FindSet(bond->m_atom2->getId());
    if(setId1!=setId2)
      m_rigidBodyMap.find( idMap.find(setId2)->second )->second->addBond(bond);
  }
}


void Molecule::buildSpanningTree(const vector<int>& rootIds) {
  //Collect rigid bodies
  vector<Rigidbody*> rigidBodies;
  for(const auto& id_rb: m_rigidBodyMap)
    rigidBodies.push_back(id_rb.second);

  //Collect root atoms
  vector<Atom*> roots;
  for(const int& rootId: rootIds)
    roots.push_back( getAtom(rootId) );

  m_spanningTree = new KinTree(rigidBodies, roots);
}

double Molecule::alignReferencePositionsTo(Molecule * base, Selection &sel){
  this->restoreAtomPos();
  metrics::RMSD rmsd(sel);
  double rmsdVal = rmsd.align(this,base);

  for (auto const& atom: m_atoms)
    atom->m_referencePosition = atom->m_position;

  return rmsdVal;
}

void Molecule::translateReferencePositionsToRoot(Molecule * base)
{
  Coordinate& thisRoot = m_spanningTree->m_root->m_rigidbody->Atoms[0]->m_referencePosition;
  Coordinate& thatRoot = base->m_spanningTree->m_root->m_rigidbody->Atoms[0]->m_referencePosition;
  Math3D::Vector3 diff = thatRoot-thisRoot;

  for(auto const& atom: m_atoms){
    atom->m_referencePosition+=diff;
  }
}

void Molecule::restoreAtomPos(){
  for (auto const& a: m_atoms)
    a->m_position = a->m_referencePosition;

  m_conf = nullptr;

  //restoreAtomIndex();
  if(m_grid!=nullptr) {
    delete m_grid;
    m_grid = nullptr;
  }

}

void Molecule::forceUpdateConfiguration(Configuration *q){
  assert(m_spanningTree!=nullptr);

  restoreAtomPos();
  m_conf = q;
  if(q==nullptr) return;

  _SetConfiguration(q);

}

void Molecule::setConfiguration(Configuration *q){
  assert(m_spanningTree!=nullptr);

  if(m_conf==q) return;

  restoreAtomPos();
  m_conf = q;
  if(q==nullptr) return;

  _SetConfiguration(q);

//  if(q->getGlobalTorsions() == nullptr){
//    log("planner")<<"Now updating global torsions"<<endl;
//    q->updateGlobalTorsions();
//  }
}



// set the positions of atoms at configuration q (according to the spanning tree)
void Molecule::_SetConfiguration(Configuration *q ){
  assert(this==q->getMolecule());

  for(size_t id=0 ; id<m_spanningTree->getNumDOFs(); ++id){
    m_spanningTree->getDOF(id)->setValue(q->m_dofs[id]);
  }

  KinVertex *root = m_spanningTree->m_root;
  root->forwardPropagate();

//  indexAtoms();
  m_grid = nullptr;
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

    //map<unsigned int,KinEdge*> m_children = node->m_edges;

    //for (map<unsigned int,KinEdge*>::iterator edge_itr=m_children.begin(); edge_itr != m_children.end(); ++edge_itr){
    for (auto const& pEdge: node->m_edges){
      //KinEdge* pEdge = edge_itr->second;
      KinVertex* newNode = pEdge->EndVertex;

      //newNode->TransformAtomPosition(m_Transformation+pEdge->getDOF()->getIndex());//,usePosition2);

      if(find(subVerts.begin(), subVerts.end(), newNode)!=subVerts.end())
        queue.push_back(newNode);
    }
  }

//  indexAtoms();
  m_grid = nullptr;
}

int Molecule::totalDofNum () const {
  if (m_spanningTree==nullptr) {
    cerr << "Error: to get the total number of DOFs in the m_molecule, you have to call Molecule::buildSpanningTree() first." << endl;
    exit(1);
  }
  return m_spanningTree->getNumDOFs();
}


pair<double,double> Molecule::vdwEnergy (set< pair<Atom*,Atom*> >* allCollisions, string collisionCheck) { // compute the total vdw energy, excluding the covalent bonds and atoms in the same rigid body

  double energy=0, collFreeEnergy=0,d_12, ratio, vdw_r1, vdw_d12, epsilon1, epsilon_12;
  set< pair<Atom*,Atom*> >::iterator mit;
  // for each atom, look for it's neighbors.
  // OLD: For each such neighbor, compute U(R_ab)=VDW_SIGMA*((VDW_R0/r_12)^12-2*(VDW_R0/r_12)^6)
  // Corrected: For each such neighbor, compute U(R_ab)=epsilon_ij*(vdw_r12/r_12)^12-2*(VDW_R0/r_12)^6) and sum up everything.
  // CHARMM: http://www.charmmtutorial.org/index.php/The_Energy_Function#Energy_calculation
  for (vector<Atom*>::const_iterator ait=m_atoms.begin(); ait!=m_atoms.end(); ++ait) {
    Atom* atom1 = *ait;
    if(!(atom1->isCollisionCheckAtom(collisionCheck)) ){//we only use atoms that are also used for clash detection
      continue;
    }
    vdw_r1 = atom1->getRadius();
    epsilon1 = atom1->getEpsilon();
    vector<Atom*> neighbors = getGrid()->getNeighboringAtomsVDW(atom1,true,true,true,true,VDW_R_MAX);
    for (vector<Atom*>::const_iterator ait2=neighbors.begin(); ait2!=neighbors.end(); ++ait2) {
      Atom* atom2 = *ait2;

      if (!(atom2->isCollisionCheckAtom(collisionCheck)) || atom1->inSameRigidbody(atom2))
        continue;//we only use atoms that are also used for clash detection (otherwise they can be too close)

      //Check initial collisions --> always excluded
      pair<Atom*,Atom*> collision_pair = make_pair(atom1,atom2);
//      set< pair<Atom*,Atom*>,int >::const_iterator mit=m_initialCollisions.find(collision_pair);
      auto mit = m_initialCollisions.find(collision_pair);
      if ( mit!=m_initialCollisions.end() )//ignore initial collision atoms
        continue;

      d_12 = atom1->distanceTo(atom2);
      vdw_d12 = (vdw_r1 + atom2->getRadius())/2.0; // from CHARMM: arithmetic mean
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
  for (vector<Atom*>::const_iterator ait=m_atoms.begin(); ait!=m_atoms.end(); ++ait) {
    Atom* atom1 = *ait;
    if(!(atom1->isCollisionCheckAtom(collisionCheck)) ){//we only use atoms that are also used for clash detection
      continue;
    }
    vdw_r1 = atom1->getRadius();
    epsilon1 = atom1->getEpsilon();
    vector<Atom*> neighbors = getGrid()->getNeighboringAtomsVDW(atom1,true,true,true,true,VDW_R_MAX);
    for (vector<Atom*>::const_iterator ait2=neighbors.begin(); ait2!=neighbors.end(); ++ait2) {
      Atom* atom2 = *ait2;

      if (!(atom2->isCollisionCheckAtom(collisionCheck)) || atom1->inSameRigidbody(atom2))
        continue;//we only use atoms that are also used for clash detection (otherwise they can be too close)

      //Check initial collisions --> always excluded
      pair<Atom*,Atom*> collision_pair = make_pair(atom1,atom2);
//      set< pair<Atom*,Atom*>,int >::const_iterator mit=m_initialCollisions.find(collision_pair);
      auto mit = m_initialCollisions.find(collision_pair);
      if ( mit!=m_initialCollisions.end() )//ignore initial collision atoms
        continue;

      d_12 = atom1->distanceTo(atom2);
      vdw_d12 = (vdw_r1 + atom2->getRadius())/2.0; // from CHARMM: arithmetic mean
      ratio = vdw_d12/d_12;
      epsilon_12 = sqrt(epsilon1 * (atom2->getEpsilon())); //from CHARMM: geometric mean
      double atomContribution = 4 * epsilon_12 * (pow(ratio,12)-2*pow(ratio,6));

      //Full enthalpy including atoms in clash constraints
      energy += atomContribution;
    }
  }
  return energy;
}

/////Create a set of common hbonds from the hbond list of another protein
void Molecule::setToHbondIntersection (Molecule * p2) {

  Hbond *ownHbond, *otherHbond;
  Atom *ownH, *ownA, *otherH, *otherA;
  list<Hbond *> ownIntersection, otherIntersection;
  list<Hbond *> deleteOwn, deleteOther;
  bool found = false;

  for (list<Hbond *>::iterator itr1=this->m_hBonds.begin(); itr1 != this->m_hBonds.end(); ++itr1) {
    found=false;
    ownHbond = (*itr1);
    ownH = ownHbond->Hatom;
    ownA = ownHbond->Acceptor;

    for (list<Hbond *>::iterator itr2 = p2->m_hBonds.begin(); itr2 != p2->m_hBonds.end(); ++itr2) {

      otherHbond = (*itr2);
      otherH = otherHbond->Hatom;
      otherA = otherHbond->Acceptor;

      if (ownH->getName() == otherH->getName() &&
          ownA->getName() == otherA->getName() &&
          ownH->getResidue()->getId() == otherH->getResidue()->getId() &&
          ownA->getResidue()->getId() == otherA->getResidue()->getId() &&
          ownH->getResidue()->getChain()->getName() == otherH->getResidue()->getChain()->getName() &&
          ownA->getResidue()->getChain()->getName() == otherA->getResidue()->getChain()->getName()) {
        found = true;
        break;
      }
    }
    if (found) {
      ownIntersection.push_back(ownHbond);
      otherIntersection.push_back(otherHbond);
    }
  }
  ///Delete hbonds non-existing in both molecules
  for (list<Hbond *>::iterator itr1=this->m_hBonds.begin(); itr1 != this->m_hBonds.end(); ++itr1) {
    ownHbond = (*itr1);
    if(find(ownIntersection.begin(),ownIntersection.end(),ownHbond) == ownIntersection.end()) {
      ownHbond->Hatom->removeHbond(ownHbond); //moved to destructor
      ownHbond->Acceptor->removeHbond(ownHbond); //moved to destructor
      delete ownHbond;
      ownHbond = nullptr;
    }
  }
  for (list<Hbond *>::iterator itr2=p2->m_hBonds.begin(); itr2 != p2->m_hBonds.end(); ++itr2) {
    otherHbond = (*itr2);
    if (find(otherIntersection.begin(), otherIntersection.end(), otherHbond) == otherIntersection.end()) {
      otherHbond->Hatom->removeHbond(otherHbond); //moved to destructor
      otherHbond->Acceptor->removeHbond(otherHbond); //moved to destructor
      delete otherHbond;
      otherHbond = nullptr;
    }
  }
  //Reset the intersection pointers
  this->m_hBonds = ownIntersection;
  p2->m_hBonds = otherIntersection;
}

void Molecule::sortHbonds() {
  m_hBonds.sort(Bond::compareIDs);
}

int Molecule::countOriginalDofs () const {
  int num = 0;
  for (list<Bond *>::const_iterator it=m_covBonds.begin(); it != m_covBonds.end(); ++it) {
    Atom* a1 = (*it)->m_atom1;
    Atom* a2 = (*it)->m_atom2;
    if ( a1->Cov_neighbor_list.size()==1 || a2->Cov_neighbor_list.size()==1 )
      continue;
    ++num;
  }
  num += m_hBonds.size();
  return num;
}

Coordinate Molecule::centerOfMass () const {
  Coordinate cur_position, center_of_mass;
  double mass, sum_of_mass=0;
  for (vector<Atom*>::const_iterator it= m_atoms.begin(); it != m_atoms.end(); ++it) {
    cur_position = (*it)->m_position;
    mass = (*it)->getMass();
    center_of_mass += cur_position * mass;
    sum_of_mass += mass;
  }
  center_of_mass /= sum_of_mass;
  return center_of_mass;
}

Coordinate Molecule::centerOfGeometry () const {
  Coordinate center;
  for (vector<Atom*>::const_iterator it= m_atoms.begin(); it != m_atoms.end(); ++it) {
    center += (*it)->m_position;
  }
  center /= m_atoms.size();
  return center;
}

double Molecule::checkCycleClosure(Configuration *q){
  setConfiguration(q);
  //Todo: Use intervals for hydrogen bond angles and lengths
  vector< pair<KinEdge*,KinVertex*> >::iterator pair_it;
  KinEdge *cycleEdge;
  int id=1;
  double maxViolation = 0.0;
  double normOfViolation = 0.0;

  log("report")<<"Conformation "<<q->m_id<<endl;

  for (pair_it=m_spanningTree->m_cycleAnchorEdges.begin(); pair_it!=m_spanningTree->m_cycleAnchorEdges.end(); ++pair_it) {

    // get end-effectors
    cycleEdge = pair_it->first;
    if( !cycleEdge->getBond()->isHBond()) continue; // for now have to check for hbonds
    //Todo: Make this more general using a torsional constraint class which is the cycle edges

    KinVertex* common_ancestor = pair_it->second;
    Hbond * hBond = reinterpret_cast<Hbond *>(cycleEdge->getBond());

    //End-effectors and their positions, corresponds to a and b
    Atom* atom1 = hBond->m_atom1;
    Atom* atom2 = hBond->m_atom2;
    Coordinate p1 = atom1->m_position; //end-effector, position 1
    Coordinate p2 = atom2->m_position; //end-effector, position 2

    if(hBond->AA== nullptr || hBond->Donor == nullptr){
      cerr<<"Not enough neighbors for h-bond between "<<atom1->getId()<<" and "<<atom2->getId()<<endl;
      exit(-1);
    }

    KinVertex* vertex1 = cycleEdge->StartVertex;
    KinVertex* vertex2 = cycleEdge->EndVertex;
    if(find(vertex1->m_rigidbody->Atoms.begin(),vertex1->m_rigidbody->Atoms.end(),atom1) == vertex1->m_rigidbody->Atoms.end()){
      vertex1=cycleEdge->EndVertex;
      vertex2=cycleEdge->StartVertex;
    }


    Math3D::Vector3 p2_test = vertex1->m_transformation * atom2->m_referencePosition;
    Math3D::Vector3 p1_test = vertex2->m_transformation * atom1->m_referencePosition;
    Math3D::Vector3 translationCons = 0.5*( (p1+p2_test) - (p1_test+p2) );

    //3 translation constraints, violation is difference to 0
    double translationViol = translationCons.normSquared();

    // Only to validate the distance violation as well
    double distanceChange = hBond->getLength() - hBond->getIniLength();

    //Angular violations
    double rightAngleChange = formatRangeRadian( hBond->getAngle_H_A_AA() - hBond->getIniAngle_H_A_AA() );
    double leftAngleChange = formatRangeRadian( hBond->getAngle_D_H_A() - hBond->getIniAngle_D_H_A() );

    //Norm
    double absoluteViolation = translationViol + rightAngleChange * rightAngleChange + leftAngleChange * leftAngleChange; //sum of all squares
    normOfViolation += absoluteViolation;

    double distanceViolation = std::fabs(distanceChange / hBond->getIniLength() * 100 );
    rightAngleChange = Math::RtoD( rightAngleChange );
    leftAngleChange = Math::RtoD( leftAngleChange );

    log("report")<<"hBond strain at "<<id<<" between "<<atom1->getId()<<" and "<<atom2->getId()<<" in res "<<atom1->getResidue()->getName()<<atom1->getResidue()->getId()<<": " << distanceViolation <<" %";
    log("report")<<", "<<rightAngleChange<<" deg rangle, "<<leftAngleChange<<" deg langle"<<endl;

    if(absoluteViolation > maxViolation){
      q->m_maxConstraintViolation = absoluteViolation;
      maxViolation = absoluteViolation;
    }


//		if(distanceViolation > 10){//10 % change of length
//			log("report") <<"Distance violation at "<<id<<" between "<<atom1->getId()<<" and "<<atom2->getId()<<" in res "<<atom1->getResidue()->getName()<<atom1->getResidue()->getId()<<": " << distanceViolation <<" %"<<endl;
//		}
//		if(std::fabs(rightAngleChange) > 0.1){
//			log("report") <<"Right angle violation "<<id<<" between "<<atom1->getId()<<" and "<<atom2->getId()<<" in res "<<atom1->getResidue()->getName()<<atom1->getResidue()->getId()<<": " << rightAngleChange<<endl;
//		}
//		if(std::fabs(leftAngleChange) > 0.1){
//			log("report") <<"Left angle violation "<<id<<" between "<<atom1->getId()<<" and "<<atom2->getId()<<" in res "<<atom1->getResidue()->getName()<<atom1->getResidue()->getId()<<": " << leftAngleChange<<endl;
//		}

    id++;
  }
  normOfViolation = sqrt(normOfViolation);

  log("report")<<"Norm of violation: "<<normOfViolation<<endl;
  log("report")<<endl;fflush(stdout);

  return normOfViolation;
}

void Molecule::computeCycleViolation(Configuration *q, gsl_vector* currentViolation){
  setConfiguration(q);
  //Todo: Use intervals for hydrogen bond angles and lengths
  vector< pair<KinEdge*,KinVertex*> >::iterator pair_it;
  KinEdge *hBondEdge;
  int id=0;
  double maxViolation = 0.0;
  double normOfViolation = 0.0;

  log("report")<<"Conformation "<<q->m_id<<endl;

  for (pair_it=m_spanningTree->m_cycleAnchorEdges.begin(); pair_it!=m_spanningTree->m_cycleAnchorEdges.end(); ++pair_it) {

    // get end-effectors
    hBondEdge = pair_it->first;
    KinVertex* common_ancestor = pair_it->second;
    Hbond * hBond = reinterpret_cast<Hbond *>(hBondEdge->getBond());

    //End-effectors and their positions, corresponds to a and b
    Atom* atom1 = hBond->m_atom1;
    Atom* atom2 = hBond->m_atom2;
    Coordinate p1 = atom1->m_position; //end-effector, position 1
    Coordinate p2 = atom2->m_position; //end-effector, position 2

    KinVertex* vertex1 = hBondEdge->StartVertex;
    KinVertex* vertex2 = hBondEdge->EndVertex;
    if(find(vertex1->m_rigidbody->Atoms.begin(),vertex1->m_rigidbody->Atoms.end(),atom1) == vertex1->m_rigidbody->Atoms.end()){
      vertex1=hBondEdge->EndVertex;
      vertex2=hBondEdge->StartVertex;
    }

    Math3D::Vector3 p2_test = vertex1->m_transformation * atom2->m_referencePosition;
    Math3D::Vector3 p1_test = vertex2->m_transformation * atom1->m_referencePosition;
    Math3D::Vector3 translationCons = 0.5*( (p1+p2_test) - (p1_test+p2) );
    //Angular violations
    double rightAngleChange = formatRangeRadian( hBond->getAngle_H_A_AA() - hBond->getIniAngle_H_A_AA() );
    double leftAngleChange = formatRangeRadian( hBond->getAngle_D_H_A() - hBond->getIniAngle_D_H_A() );

    gsl_vector_set(currentViolation,id*5+0,translationCons.x);
    gsl_vector_set(currentViolation,id*5+1,translationCons.y);
    gsl_vector_set(currentViolation,id*5+2,translationCons.z);
    gsl_vector_set(currentViolation,id*5+3,rightAngleChange);
    gsl_vector_set(currentViolation,id*5+4,leftAngleChange);

    id++;
  }

}

/* Resample all sugar conformations in a segment of an RNA/DNA chain specified by startRes (included) and
 * endRes (not included) and reclose using the localRebuild method.
 * Returns the configuration and leaves the m_molecule with that configuration set.
 */
Configuration*Molecule::resampleSugars(int startRes, int endRes, Configuration* cur, int aggression){
  //log("debugRas")<<"resampleSugars("<<startRes<<", "<<endRes<<" ... "<<aggression<<")"<<endl;
  vector<int> ignoreDOFs;
  vector<int> resetDOFs;
  vector<double> resetValues;
  vector<int> recloseDOFs;
  for(vector<KinEdge*>::iterator eit = m_spanningTree->m_edges.begin(); eit!=m_spanningTree->m_edges.end(); eit++){
    KinEdge* e = *eit;
    int res1 = e->getBond()->m_atom1->getResidue()->getId();
    int res2 = e->getBond()->m_atom2->getResidue()->getId();
    if( res1==(startRes-1) && res2==(startRes-1) && e->getBond()->m_atom1->getName()=="C3'" && e->getBond()->m_atom2->getName()=="O3'" ){
      recloseDOFs.push_back(e->getDOF()->getIndex());
      if(aggression>=2){
        resetDOFs.push_back(e->getDOF()->getIndex());
        resetValues.push_back(RandomAngleUniform(3.1415));
      }
    }
    if( res1==endRes && res2==endRes && (
        e->getBond()->m_atom1->getName()=="P" ||
        e->getBond()->m_atom1->getName()=="O5'"
    ) ){
      recloseDOFs.push_back(e->getDOF()->getIndex());
      if(aggression>=2){
        resetDOFs.push_back(e->getDOF()->getIndex());
        resetValues.push_back(RandomAngleUniform(3.1415));
      }
    }
    if( (res1>=startRes && res1<endRes) || (res2>=startRes && res2<endRes) ){
      if(	e->getBond()->m_atom1->getName()=="P" || e->getBond()->m_atom2->getName()=="P" ||
           e->getBond()->m_atom1->getName()=="C5'" || e->getBond()->m_atom2->getName()=="C5'" ||
           e->getBond()->m_atom1->getName()=="C4'" || e->getBond()->m_atom2->getName()=="C4'" ||
           e->getBond()->m_atom1->getName()=="C3'" || e->getBond()->m_atom2->getName()=="C3'" ||
           e->getBond()->m_atom1->getName()=="O3'" || e->getBond()->m_atom2->getName()=="O3'" ||
           e->getBond()->m_atom1->getName()=="O5'" || e->getBond()->m_atom2->getName()=="O5'" ){
        recloseDOFs.push_back(e->getDOF()->getIndex());
        if(aggression>=2){
          resetDOFs.push_back(e->getDOF()->getIndex());
          resetValues.push_back(RandomAngleUniform(3.1415));
        }
      }else if(
          e->getBond()->m_atom1->getName()=="C2'" || e->getBond()->m_atom2->getName()=="C2'" ||
          e->getBond()->m_atom1->getName()=="C1'" || e->getBond()->m_atom2->getName()=="C1'" ){
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
 * Returns the configuration and leaves the m_molecule with that configuration set.
 */
Configuration* Molecule::localRebuild(vector<int>& resetDOFs, vector<double>& resetValues, vector<int>& recloseDOFs, vector<int>& ignoreDOFs, Configuration* cur){
  //enableLogger("debugRebuild");

  double start_time, end_time, total_time;
  CTKTimer timer;
  start_time = timer.getTimeNow();

  //Configuration* ret = new Configuration(m_spanningTree->getNumDOFs());

  //Configuration* ret = new Configuration(this);
  //ret->Copy(cur);
  Configuration* ret = cur->clone();

  //restoreAtomPos(false);
  setConfiguration(ret);
  //ret->computeCycleJacobianAndNullSpace();

  //Find smallest connected subgraph that contains both resetDOFS and recloseDOFs (TODO: Approximate Steiner tree)
  vector<KinVertex*> subVerts;
  vector<KinEdge*> subEdges;
  KinEdge* entry = nullptr;
  for(vector<KinEdge*>::iterator eit = m_spanningTree->m_edges.begin(); eit!=m_spanningTree->m_edges.end(); eit++){
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
  for(vector<KinEdge*>::iterator eit = m_spanningTree->m_edges.begin(); eit!=m_spanningTree->m_edges.end(); eit++){
    KinEdge* e = *eit;
    bool firstInSub = find(subVerts.begin(), subVerts.end(), e->StartVertex)!=subVerts.end();
    bool lastInSub = find(subVerts.begin(), subVerts.end(), e->EndVertex)!=subVerts.end();
    //if(!firstInSub &&  lastInSub) entry = e;

    if( firstInSub && !lastInSub && e->StartVertex!=entry->StartVertex) boundary.push_back(e);
  }
  for(vector<pair<KinEdge*,KinVertex*> >::iterator it = m_spanningTree->m_cycleAnchorEdges.begin(); it!=m_spanningTree->m_cycleAnchorEdges.end(); it++){
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
      storedPositions.push_back(e->getBond()->m_atom1->m_position);
      storedPositions.push_back(e->getBond()->m_atom2->m_position);
      //log("debugRas")<<"Cylinder["<<e->getBond()->Atom1->m_position[0]<<", "<<e->getBond()->Atom1->m_position[1]<<", "<<e->getBond()->Atom1->m_position[2]<<", ";
      //log("debugRas")<<e->getBond()->m_atom2->m_position[0]<<", "<<e->getBond()->m_atom2->m_position[1]<<", "<<e->getBond()->m_atom2->m_position[2]<<", 0.1, 0.9,0.9,0.2]"<<endl;

      //log("debugRebuild")<<"1: storedTorsion["<<storedTorsions.size()<<"] .. ";
      storedTorsions.push_back(e->getBond()->getTorsion());
    }
  }


  //Make sure only endpoints within subgraph and boundary can move
  for(vector<KinVertex*>::iterator vit = subVerts.begin(); vit!=subVerts.end(); vit++){
    for(vector<Atom*>::iterator ait=(*vit)->m_rigidbody->Atoms.begin();ait!=(*vit)->m_rigidbody->Atoms.end();ait++){
      (*ait)->m_position = (*ait)->m_referencePosition;
    }
  }
  for(vector<KinEdge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
    //(*eit)->getBond()->Atom1->m_bPositionModified = false;
    //(*eit)->getBond()->m_atom2->m_bPositionModified = false;
    (*eit)->getBond()->m_atom1->m_position = entry->getBond()->m_atom1->m_referencePosition;
    (*eit)->getBond()->m_atom2->m_position = entry->getBond()->m_atom2->m_referencePosition;
  }
  //entry->getBond()->Atom1->m_bPositionModified = false;
  //entry->getBond()->m_atom2->m_bPositionModified = false;
  entry->getBond()->m_atom1->m_position = entry->getBond()->m_atom1->m_referencePosition;
  entry->getBond()->m_atom2->m_position = entry->getBond()->m_atom2->m_referencePosition;
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
      Math3D::Vector3 v1 = (storedPositions[c*2+0])-(eBoundary->getBond()->m_atom1->m_position);
      Math3D::Vector3 v2 = (storedPositions[c*2+1])-(eBoundary->getBond()->m_atom2->m_position);
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
          //Math3D::Vector3 v1 = ComputeJacobianEntry(e->getBond()->Atom1->m_position, e->getBond()->m_atom2->m_position, eBoundary->getBond()->Atom1->m_position);
          //Math3D::Vector3 v2 = ComputeJacobianEntry(e->getBond()->Atom1->m_position, e->getBond()->m_atom2->m_position, eBoundary->getBond()->m_atom2->m_position);
          Math3D::Vector3 v1 = e->getDOF()->getDerivative(eBoundary->getBond()->m_atom1->m_position);
          Math3D::Vector3 v2 = e->getDOF()->getDerivative(eBoundary->getBond()->m_atom2->m_position);
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
    SVDMKL svd(J);
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
        (*ait)->m_position = (*ait)->m_referencePosition;

      }
    }
    for(vector<KinEdge*>::iterator eit = boundary.begin(); eit!=boundary.end(); eit++){
      (*eit)->getBond()->m_atom1->m_position = (*eit)->getBond()->m_atom1->m_referencePosition;
      (*eit)->getBond()->m_atom2->m_position = (*eit)->getBond()->m_atom2->m_referencePosition;
    }
    entry->getBond()->m_atom1->m_position = entry->getBond()->m_atom1->m_referencePosition;
    entry->getBond()->m_atom2->m_position = entry->getBond()->m_atom2->m_referencePosition;

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
    //	setConfiguration(cur);
    //	stringstream ss; ss<<"test_"<<resi<<"_init.pdb";
    //	IO::writePdb(this, ss.str());
    //}
    //{
    //	setConfiguration(ret);
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

Molecule* Molecule::deepClone() const{
  Molecule* ret = new Molecule();
  ret->setName(m_name);

  //Clone atoms
  for(auto const& atom: getAtoms()){
    ret->addAtom(
        atom->isHetatm(),
        atom->getResidue()->getChain()->getName(),
        atom->getResidue()->getName(),
        atom->getResidue()->getId(),
        atom->getName(),
        atom->getId(),
        atom->m_referencePosition
    );

//    if(atom->is)
//    atom->set
  }

  //Clone covalent bonds
  for(auto const& covBond: getCovBonds()){
    Atom* a1_new = ret->getAtom( covBond->m_atom1->getId() );
    Atom* a2_new = ret->getAtom( covBond->m_atom2->getId() );
//    Bond* newBond = ret->addCovBond( a1_new->getResidue(), a2_new->getResidue(), a1_new->getName(), a2_new->getName() );
    Bond* newBond = ret->addCovBond( a1_new, a2_new );
    if(newBond){
      newBond->Bars = covBond->Bars;
      newBond->rigidified = covBond->rigidified;
    }
  }

  //Clone hbonds
  for(auto const& hbond: getHBonds()) {
    Atom *hatom = ret->getAtom(hbond->Hatom->getId());
    Atom *donor = ret->getAtom(hbond->Donor->getId());
    Atom *acc = ret->getAtom(hbond->Acceptor->getId());
    Atom *AA = ret->getAtom(hbond->AA->getId());
    Hbond *new_hb = new Hbond(hatom, acc, donor, AA, hbond->getIniEnergy());
    new_hb->Bars = hbond->Bars;
    new_hb->rigidified = hbond->rigidified;
    ret->addHbond(new_hb);
  }

  //Clone dbonds
  for(auto const& dbond: getDBonds()) {
    Atom *a1 = ret->getAtom(dbond->m_atom1->getId());
    Atom *a2 = ret->getAtom(dbond->m_atom2->getId());
    DBond *new_db = new DBond(a1, a2);
    new_db->Bars = dbond->Bars;
    ret->addDBond(new_db);
  }

  // Fill in the second_cov_neighbor_list
  // Cannot do this step when the bond is still in creation because it won't know the neighbors of its neighbors yet
  for( auto const &atom: ret->getAtoms()) {
    for (auto const &n1: atom->Cov_neighbor_list) {
      for (auto const &n2: n1->Cov_neighbor_list) {
        // check if n2 is ait itself. if yes, ignore it.
        if (n2 == atom)
          continue;
        // check whether n2 is already in the second_cov_neighbor_list of ait
        bool got_already = false;
        for (auto const &n3: atom->Second_cov_neighbor_list) {
          if (n3 == n2) {
            got_already = true;
            break;
          }
        }
        if (!got_already)
          atom->Second_cov_neighbor_list.push_back(n2);
      }
    }
  }
  return ret;
}

Molecule* Molecule::collapseRigidBonds(int collapseLevel) {
  assert(collapseLevel == 1 || collapseLevel == 2);

  m_conf->rigidityAnalysis();//identifies rigidified bonds necessary for collapsing

  Molecule *ret=deepClone();

  int i=0; //indexing for hBonds
  //To collapse molecule, we turn rigid h-bonds into covalent bonds
  if(collapseLevel==2) {
    for( auto const &edge_nca_pair : m_spanningTree->m_cycleAnchorEdges ) {

      KinEdge *edge=edge_nca_pair.first;

      //Get corresponding rigidity information
      if( m_conf->getNullspace()->isHBondRigid(i++)) {
        //If its a rigid hbond convert it to a rigid covalent bond
        if( edge->getBond()->isHBond()) {
          Atom *a1_new=ret->getAtom(edge->getBond()->m_atom1->getId());
          Atom *a2_new=ret->getAtom(edge->getBond()->m_atom2->getId());
          Bond *newBond=ret->addCovBond(a1_new, a2_new);
          if( newBond ) {
            newBond->rigidified=true;
//          log("debug")<<"Molecule.cpp:"<<__LINE__<<" covalently connecting "<<a1_new<<"-"<<a2_new<<" with rigid bond"<<endl;
          }
        }
      }//end if
    }
  }

  //Recreate roots vector. For edge leaving the root, descend until a rigid body is found and then pick the first atom
  vector<int> roots;
  for(auto const& e: m_spanningTree->m_root->m_edges){
    KinVertex* v = e->EndVertex;
    while(v->m_rigidbody==nullptr)
      v = v->m_edges[0]->EndVertex;
    Atom* firstChainRootAtom = v->m_rigidbody->Atoms[0];
    roots.push_back(firstChainRootAtom->getId());
  }

  Selection movingResidues("all");
  ret->buildRigidBodies(movingResidues, collapseLevel); //Necessary to do before building spanning tree
  ret->buildSpanningTree(roots); //Necessary before conformations are defined
  ret->setConfiguration(new Configuration(ret));
  ret->setCollisionFactor(m_collisionFactor); //Sets the initial collisions

  return ret;
}

void Molecule::writeRigidbodyIDToBFactor()
{
  //We sort them so the biggest clusters always have the smallest IDs and the same colors
  std::vector< std::pair<int, unsigned int> > sortedRBs; //pair of int size, unsigned int ID

  unsigned int maxSize = 0;
  unsigned int maxIndex=0;

  for(auto const& it : m_rigidBodyMap ){
    int rbSize = it.second->size();
    sortedRBs.push_back(make_pair( rbSize , it.first ));
    if(rbSize > maxSize){
      maxSize = rbSize;
      maxIndex = it.first;
    }
  }

  vector< pair<int, unsigned int> >::iterator vsit = sortedRBs.begin();
  vector< pair<int, unsigned int> >::iterator veit = sortedRBs.end();

  sort(vsit, veit,Rigidbody::compareSize); //sorts them by size

  int outputID = 0;
  for(auto const& idPair : sortedRBs){
    ///get rb id sorted by size, access rb in protein, color all atoms to id
    Rigidbody* currentRB = m_rigidBodyMap[idPair.second];
//    log("debug")<<"WriteRBID: sorted number "<<outputID<<" rigid body ID "<<currentRB->id()<<" size "<<currentRB->size()<<endl;
    for(auto const& atom: currentRB->Atoms){
//      atom->setBFactor(float(outputID)/100);
      atom->setBFactor(float(currentRB->id())/100);
    }
    outputID++;
  }

  ///Store information in the molecule
  m_conf->m_numClusters = m_rigidBodyMap.size();
  m_conf->m_maxSize = maxSize;
  m_conf->m_maxIndex = maxIndex;
}

const std::vector<Rigidbody*> Molecule::getRigidbodies() const
{
  vector<Rigidbody*> ret;
  for(auto const& it: m_rigidBodyMap){
    ret.push_back(it.second);
  }
  return ret;
}
