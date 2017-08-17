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
#include <string>
#include <fstream>
#include <iomanip>
#include <math/math.h>
#include <assert.h>
#include <set>
#include <math/SVDGSL.h>
#include <math/SVDMKL.h>
#include <gsl/gsl_matrix_double.h>
#include <math/gsl_helpers.h>
#include <math/NullspaceSVD.h>

#include "Configuration.h"
#include "Molecule.h"
#include "CTKTimer.h"
#include "core/Bond.h"
#include "HBond.h"
#include "ResidueProfiles.h"
#include "DisjointSets.h"
#include "Logger.h"

using namespace std;

double jacobianAndNullspaceTime = 0;
double rigidityTime = 0;

gsl_matrix* Configuration::CycleJacobian = nullptr;
gsl_matrix* Configuration::HBondJacobian = nullptr;
Configuration* Configuration::CycleJacobianOwner = nullptr;
SVD* Configuration::JacobianSVD = nullptr;

//gsl_matrix* Configuration::ClashAvoidingJacobian = nullptr;
//Nullspace* Configuration::ClashAvoidingNullSpace = nullptr;

Configuration::Configuration(Molecule * mol):
  m_molecule(mol),
  nullspace(nullptr),
  m_parent(nullptr),
  m_dofs_global(nullptr),
  m_treeDepth(0)
{
  assert(m_molecule!=nullptr);

  m_id 										 = 0;
  m_vdwEnergy 						 = 0;
  m_deltaH                 = 0;
  m_distanceToTarget       = 99999;
  m_paretoFrontDistance 	 = 99999;
  m_distanceToParent       = 0;
  m_distanceToIni          = 0;
  m_maxIndex               = 0;
  m_maxSize                = 0;
  m_maxConstraintViolation = 99999;
  m_numClusters            = 0;
  m_minCollisionFactor     = 0;
  m_usedClashPrevention    = false;
  m_clashFreeDofs          = m_molecule->m_spanningTree->getNumDOFs();

  // Set up DOF-values and set them to 0
  m_dofs = new double[getNumDOFs()];
  //m_sumProjSteps = new double[getNumDOFs()];
  for(int i=0; i<getNumDOFs(); ++i){
    m_dofs[i] = 0;
    //m_sumProjSteps[i]=0;
  }
}

Configuration::Configuration(Configuration* parent_):
    m_molecule(parent_->m_molecule),
    m_parent(parent_),
    m_dofs_global(nullptr),
    nullspace(nullptr),
    m_treeDepth(parent_->m_treeDepth +1)
{
  assert(m_molecule!=nullptr);
  if(m_molecule==NULL){
    std::cerr<<"Configuration(..) - molecule is NULL"<<std::endl;
  }
  m_id                     = -1;//Setting id to -1 by default is important. Check out PoissonSampler2.cpp for example
  m_vdwEnergy              = 0;
  m_deltaH                 = 0;
  m_distanceToTarget       = 99999;
  m_paretoFrontDistance    = 99999;
  m_distanceToParent       = 0;
  m_distanceToIni          = 0;
  m_maxIndex               = 0;
  m_maxSize                = 0;
  m_maxConstraintViolation = 99999;
  m_numClusters            = 0;
  m_minCollisionFactor     = 0;
  m_usedClashPrevention    = false;
  m_clashFreeDofs          = m_molecule->m_spanningTree->getNumDOFs();

  parent_->m_children.push_back(this);

  // Set up DOF-values and set them to 0
  m_dofs = new double[getNumDOFs()];
  //m_sumProjSteps = new double[getNumDOFs()];
  for(int i=0; i<getNumDOFs(); ++i){
    m_dofs[i] = 0;
    //m_sumProjSteps[i]=0;
  }
}

Configuration::~Configuration(){
  if (m_molecule->m_conf==this)
    m_molecule->setConfiguration(nullptr);

  // Remove DOF-value arrays
  if (m_dofs != nullptr)
    delete[] m_dofs;
  if(m_dofs_global != nullptr)
    delete[] m_dofs_global;
  //if( m_sumProjSteps != nullptr)
  //  delete[] m_sumProjSteps;

//  //m_biggerRBMap.clear();
//  for(auto const& rbPair: m_biggerRBMap)
//    delete rbPair.second;

  //m_sortedRBs.clear();
  if(nullspace)
    delete nullspace;

  if( m_parent!=nullptr )
    m_parent->m_children.remove(this);
}

void Configuration::computeCycleJacobianAndNullSpace() {
  CTKTimer timer;
  timer.Reset();
  double old_time = timer.LastElapsedTime();

  //Compute the Jacobian matrix
  computeJacobians();

  if (JacobianSVD!=nullptr) {
    nullspace = new NullspaceSVD(JacobianSVD);
    nullspace->updateFromMatrix();
  }

  double new_time = timer.ElapsedTime();
  jacobianAndNullspaceTime += new_time - old_time;

//  if(CycleJacobian!=nullptr) {
//    nullspace->performRigidityAnalysis(HBondJacobian);
////    identifyBiggerRigidBodies();
//  }

  double new_time_2 = timer.ElapsedTime();
  rigidityTime += new_time_2 - new_time;
}

void Configuration::rigidityAnalysis() {

  //Checks if Jacobians need update
  computeJacobians();

  if(nullspace==nullptr){ //update nullspace if necessary
    computeCycleJacobianAndNullSpace();
  }
  CTKTimer timer;
  timer.Reset();
  double old_time = timer.LastElapsedTime();

  if(CycleJacobian!=nullptr) {///identifies rigid/rotatable hbonds and bonds based on set cut-off
    nullspace->performRigidityAnalysis(HBondJacobian);
  }

  int i=0; //indexing for hBonds

  //First, constrain cycle edge bonds
  for (auto const &edge_nca_pair : m_molecule->m_spanningTree->m_cycleAnchorEdges) {

    KinEdge *edge = edge_nca_pair.first;
    KinVertex *common_ancestor = edge_nca_pair.second;

    //Get corresponding rigidity information
    if (nullspace->isHBondRigid(i++)) {
      Bond* bond = edge->getBond();
      //If its a rigid hbond convert it to a rigid covalent bond
      if (bond->isHBond()) {
        bond->rigidified = true;
      }
    }//end if

    //Now, the dihedral angles
    KinVertex* vertex1 = edge->StartVertex;
    KinVertex* vertex2 = edge->EndVertex;

    //Trace back along dof m_edges for vertex 1
    while ( vertex1 != common_ancestor ) {
      KinVertex* parent = vertex1->m_parent;
      KinEdge* p_edge = parent->findEdge(vertex1);

      int dof_id = p_edge->getDOF()->getCycleIndex();
      if (dof_id!=-1) { // this edge is a cycle DOF, dof_id is the corresponding column!
        if( nullspace->isCovBondRigid(dof_id) ) {
          Bond* bond = p_edge->getBond();
          bond->rigidified = true;
        }
      }
      vertex1 = parent;
    }

    //Trace back along edges from vertex 2
    while ( vertex2 != common_ancestor ) {
      KinVertex* parent = vertex2->m_parent;
      KinEdge* p_edge = parent->findEdge(vertex2);

      int dof_id = p_edge->getDOF()->getCycleIndex();
      if (dof_id!=-1) { // this edge is a cycle DOF, dof_id is the corresponding column!
        if( nullspace->isCovBondRigid(dof_id) ) {
          Bond* bond = p_edge->getBond();
          bond->rigidified = true;
        }
      }
      vertex2 = parent;
    }
  }

  double new_time = timer.ElapsedTime();
  rigidityTime += new_time - old_time;
}

//void Configuration::identifyBiggerRigidBodies(){
//
//  /// Now, all dihedrals and hBonds that are fixed have the set flag constrained = true!
//  m_biggerRBMap.clear();
//  m_sortedRBs.clear();
//
//  readBiggerSet();
//
//  /// Now, we have the map of bigger rigid bodies with atoms
//  /// We will order the atoms by ID for a nicer display
//  m_numClusters = 0;
//  m_maxSize = 0;
//  m_maxIndex=0;
//
////	cout<<"Map of bigger rigid bodies:"<<endl;
//  auto it = m_biggerRBMap.begin();
//  while( it != m_biggerRBMap.end() ){
//    m_sortedRBs.push_back( make_pair( (it->second)->size(), it->first ));
//    m_numClusters++;
//    int num = 0;
//
//    ///Sorting
//    vector<Atom*>::iterator sit = it->second->Atoms.begin();
//    vector<Atom*>::iterator eit = it->second->Atoms.end();
//    sort(sit,eit,Atom::compare);
//
////		cout<<"Bigger rb with ID "<<it->first<<" contains atoms with IDs: "<<endl;
//    auto ait = it->second->Atoms.begin();
//    while( ait != it->second->Atoms.end() ){
//      //cout<<(*ait)->m_id<<endl;
//      ait++;
//      num++;
//    }
//
//    if( num > m_maxSize ){
//      m_maxSize = num;
//      m_maxIndex = it->first;
//    }
//    ++it;
//  }
//
//  vector< pair<int, unsigned int> >::iterator vsit = m_sortedRBs.begin();
//  vector< pair<int, unsigned int> >::iterator veit = m_sortedRBs.end();
//
//  sort(vsit, veit,compareSize);
//
//}

////---------------------------------------------------------
//bool Configuration::compareSize(pair<int, unsigned int> firstEntry, pair<int, unsigned int> secondEntry) {
//  if( firstEntry.first > secondEntry.first )
//    return true;
//  if( firstEntry.first < secondEntry.first )
//    return false;
//
//  return (firstEntry.second < secondEntry.second);
//}

//void Configuration::readBiggerSet(){
//
//  //Create disjoint set
//  DisjointSets ds(m_molecule->getAtoms()[m_molecule->size() - 1]->getId() + 1); //Assumes the last atom has the highest id.
//
//  //For each atom, a1, with exactly one cov neighbor, a2, call Union(a1,a2)
//  for (int i=0;i< m_molecule->size();i++){
//    Atom* atom = m_molecule->getAtoms()[i];
//    if(atom->Cov_neighbor_list.size()==1 && atom->Hbond_neighbor_list.size()==0){
//      ds.Union(atom->getId(), atom->Cov_neighbor_list[0]->getId());
//    }
//  }
//
//
//  //For each fixed or rigidified bond (a1,a2) call Union(a1,a2)
//  for (auto const& bond: m_molecule->getCovBonds()){
//    //First, simply check if bond is rigidified
//    if( bond->rigidified || bond->Bars == 6){
//      ds.Union(bond->Atom1->getId(), bond->m_atom2->getId());
//    }
//  }
//
//  //Also, do the same thing for the hydrogen bonds insert the h-bonds at the correct place
//  for(auto const& bond: m_molecule->getHBonds()){
//    if( bond->rigidified ){
//      ds.Union(bond->Atom1->getId(), bond->m_atom2->getId());
//    }
//  }
//
//  //All disjoint sets have been united if they are rigidified or bonded fix!
//  //We now just add the covalent neighbors to have the representation with atoms in multiple rbs!
//
//  int c=0;
//  map<int,int> idMap;//Maps atom id's to rigid body id's for use in the DS structure.
//
//  //Map the set-ID's (first map entry) to RB-ID's (second map entry) and add bonded atoms to RBs.
//  for (int i=0;i< m_molecule->size();i++){
//    Atom* atom = m_molecule->getAtoms()[i];
//
//    //Map the set-id to the RB-id
//    int set_id = ds.FindSet(atom->getId());
//    int body_id;
//    if(idMap.find(set_id)!=idMap.end())
//      body_id = idMap.find(set_id)->second;
//    else {
//      body_id = c++;
//      idMap[set_id] = body_id;
//    }
//    //If the set containing a1 is not a rigid body: create one
//    if ( m_biggerRBMap.find(body_id) == m_biggerRBMap.end() ) {
//      Rigidbody* new_rb = new Rigidbody(body_id);
//      m_biggerRBMap[body_id] = new_rb;
//    }
//    Rigidbody* rb = m_biggerRBMap[body_id];
//    if (!rb->containsAtom(atom)){
//      rb->Atoms.push_back(atom);
//      atom->setBiggerRigidbody(rb);
//      //For graphical display, we assign this rb-id to the rbColumn variable of the atom
//      atom->setBFactor( float( rb->id() ) );
//    }
//
//  }
//
//}

//------------------------------------------------------
void Configuration::writeQToBfactor(){
  for (auto const& edge: m_molecule->m_spanningTree->Edges) {
    if(edge->EndVertex->m_rigidbody == NULL)
      continue; //global dofs
    float value = m_dofs[edge->getDOF()->getIndex() ];
    for (auto const& atom : edge->EndVertex->m_rigidbody->Atoms ){
      atom->setBFactor(value);
    }
  }
}

//------------------------------------------------------
void Configuration::updateGlobalTorsions(){
  if(m_dofs_global != nullptr)
    return;

  updateMolecule();
  m_dofs_global = new double[getNumDOFs()];
  for(int i=0; i<getNumDOFs(); ++i){
    m_dofs_global[i] = m_molecule->m_spanningTree->getDOF(i)->getGlobalValue();
  }
  //for (vector<KinEdge*>::iterator itr= m_molecule->m_spanningTree->Edges.begin(); itr!= m_molecule->m_spanningTree->Edges.end(); ++itr) {
  //  KinEdge* pEdge = (*itr);
  //  int dof_id = pEdge->DOF_id;
  //  if (dof_id==-1)
  //    continue;
  //  m_dofs_global[dof_id] = pEdge->getBond()->getTorsion();
  //}
}

double Configuration::getGlobalTorsion( int i ) {
  assert(i>=0 && i<getNumDOFs());
  updateGlobalTorsions();
  return m_dofs_global[i];
}

double* Configuration::getGlobalTorsions() {
  updateGlobalTorsions();
  return m_dofs_global;
}

unsigned int Configuration::getNumDOFs() const {
  return m_molecule->m_spanningTree->getNumDOFs();
}

Configuration* Configuration::clone() const {

  Configuration* ret;

  if (m_parent) {
    ret = new Configuration(m_parent);
  }else{
    ret = new Configuration(m_molecule);
  }

  ret->m_id = m_id;
  ret->m_vdwEnergy = m_vdwEnergy;

  for(int i=0; i <getNumDOFs(); ++i){
    ret->m_dofs[i]        = m_dofs[i];
  }
  if(m_dofs_global!=nullptr){
    ret->m_dofs_global = new double[getNumDOFs()];
    for(int i=0; i<getNumDOFs(); ++i) {
      ret->m_dofs_global[i]  = m_dofs_global[i];
    }
  }else{
    ret->m_dofs_global = nullptr;
  }
  return ret;
}

void Configuration::Print () {
  for (int i=0; i < getNumDOFs(); ++i){
    cout << "Relative: " << m_dofs[i] << " ";
    cout << "Global: " << m_dofs_global[i] << " ";
  }
  cout << endl;
}


void Configuration::computeJacobians() {
//  log("debug")<<"computeJacobians"<<endl;
  if(CycleJacobianOwner==this) return;
  CycleJacobianOwner = this;

  updateMolecule();

  // No cycles
  if(m_molecule->m_spanningTree->m_cycleAnchorEdges.size() == 0) {
    CycleJacobian = nullptr; //TODO: Memory leak
    return;
  }

  int hBond_row_num = 0;
  int row_num = 0;
  for(auto const& edge: m_molecule->m_spanningTree->m_cycleAnchorEdges){
    if(edge.first->getBond()->isDBond()) {
      row_num += 3;
    }else{
      row_num += 5;
      hBond_row_num += 1;
    }
  }

//  int hBond_row_num = (m_molecule->m_spanningTree->m_cycleAnchorEdges).size();
//  int row_num = hBond_row_num * 5; // 5 times the number of cycles, non-redundant description
  int col_num = m_molecule->m_spanningTree->getNumCycleDOFs(); // number of DOFs in cycles

  if(CycleJacobian==nullptr){
    CycleJacobian = gsl_matrix_calloc(row_num,col_num);
    JacobianSVD = SVD::createSVD(CycleJacobian);//new SVDMKL(CycleJacobian);
  }else if(CycleJacobian->size1==row_num && CycleJacobian->size2==col_num){
    gsl_matrix_set_zero(CycleJacobian);
  }else{
    gsl_matrix_free(CycleJacobian);
    CycleJacobian = gsl_matrix_calloc(row_num,col_num);
    JacobianSVD = SVD::createSVD(CycleJacobian);//new SVDMKL(CycleJacobian);
  }

  ///HBond Jacobian
  if(HBondJacobian==nullptr){
      HBondJacobian = gsl_matrix_calloc(hBond_row_num,col_num);

  }else if(HBondJacobian->size1==hBond_row_num && HBondJacobian->size2==col_num){
      gsl_matrix_set_zero(HBondJacobian);

  }	else{
      gsl_matrix_free(HBondJacobian);
      HBondJacobian = gsl_matrix_calloc(hBond_row_num,col_num);
  }

  // for each cycle, fill in the Jacobian entries
  int i=0;
  int hbidx=0;
  for (std::pair<KinEdge*,KinVertex*>& edge_vertex_pair: m_molecule->m_spanningTree->m_cycleAnchorEdges)
  {
    // get end-effectors
    KinEdge* edge_ptr = edge_vertex_pair.first;
    KinVertex* common_ancestor = edge_vertex_pair.second;
    Bond * bond_ptr = edge_ptr->getBond();
//    int bondIndex = bond_ptr

    //End-effectors and their positions, corresponds to a and b
    Atom* atom1 = bond_ptr->m_atom1;
    Atom* atom2 = bond_ptr->m_atom2;
    Coordinate p1 = atom1->m_position; //end-effector, position 1
    Coordinate p2 = atom2->m_position; //end-effector, position 2
    log("debug")<<"Jacobian row "<<i<<", atom 1: "<<atom1->getId()<<", atom 2: "<<atom2->getId()<<endl;

    KinVertex* vertex1 = edge_ptr->StartVertex;
    KinVertex* vertex2 = edge_ptr->EndVertex;
    if(find(vertex1->m_rigidbody->Atoms.begin(),vertex1->m_rigidbody->Atoms.end(),atom1) == vertex1->m_rigidbody->Atoms.end()){
      vertex1=edge_ptr->EndVertex;
      vertex2=edge_ptr->StartVertex;
    }

    //Use the covalently bonded atoms to find A-1 and B-1
    Atom* atom1_prev = atom1->Cov_neighbor_list[0] == atom2 ? atom1->Cov_neighbor_list[1] : atom1->Cov_neighbor_list[0];
    Atom* atom2_prev = atom2->Cov_neighbor_list[0] == atom1 ? atom2->Cov_neighbor_list[1] : atom2->Cov_neighbor_list[0];

     //Make sure a1 is the covalent neighbor of Atom1 with lexicographically smallest name
    for(vector<Atom*>::iterator ait = atom1->Cov_neighbor_list.begin(); ait!=atom1->Cov_neighbor_list.end(); ait++){
        Atom* a = *ait;
        if(a!=atom2 && a->getName()<atom1_prev->getName())
          atom1_prev = a;
    }
    //Make sure a4 is the covalent neighbor of m_atom2 with lexicographically smallest name
    for(vector<Atom*>::iterator ait = atom2->Cov_neighbor_list.begin(); ait!=atom2->Cov_neighbor_list.end(); ait++){
        Atom* a = *ait;
        if(a!=atom1 && a->getName()<atom2_prev->getName())
          atom2_prev = a;
    }

    if(!atom1 || !atom1_prev || !atom2 || !atom2_prev){
      cerr<<"Not all neighboring atoms for non-redundant Jacobian defined, quitting!"<<endl;
      exit(-1);
    }

    Coordinate p1_prev = atom1_prev->m_position;
    Coordinate p2_prev = atom2_prev->m_position;
    //Now we have all atoms we need --> Calculate Jacobian


    // trace back until the common ancestor from vertex1
    while ( vertex1 != common_ancestor ) {
      KinVertex* parent = vertex1->m_parent;
      KinEdge* p_edge = parent->findEdge(vertex1);

      //int dof_id = p_edge->Cycle_DOF_id;
      int dof_id = p_edge->getDOF()->getCycleIndex();
      if (dof_id!=-1) { // this edge is a DOF

        Math3D::Vector3 derivativeP1      = p_edge->getDOF()->getDerivative(p1);
        Math3D::Vector3 derivativeP2      = p_edge->getDOF()->getDerivative(p2);

        Math3D::Vector3 jacobianEntryTrans=0.5*(derivativeP1 + derivativeP2);
        gsl_matrix_set(CycleJacobian,i + 0, dof_id, jacobianEntryTrans.x); //set: Matrix, row, column, what to set
        gsl_matrix_set(CycleJacobian,i + 1, dof_id, jacobianEntryTrans.y);
        gsl_matrix_set(CycleJacobian,i + 2, dof_id, jacobianEntryTrans.z);
        if(bond_ptr->isHBond()) {
          Math3D::Vector3 derivativeP1_prev = p_edge->getDOF()->getDerivative(p1_prev);
          double jacobianEntryRot1 = dot((p2 - p1),(derivativeP1-derivativeP1_prev));
          double jacobianEntryRot2 = dot((p2 - p2_prev), (derivativeP1 - derivativeP2));
          gsl_matrix_set(CycleJacobian, i + 3, dof_id, jacobianEntryRot1);
          gsl_matrix_set(CycleJacobian, i + 4, dof_id, jacobianEntryRot2);

          ///Matrix to check hBond Rotation
          double hBondEntry = dot((p1_prev - p2_prev), (derivativeP1_prev));
          gsl_matrix_set(HBondJacobian, hbidx, dof_id, hBondEntry);
        }

      }
      vertex1 = parent;
    }
    // trace back until the common ancestor from vertex2
    while ( vertex2 != common_ancestor ) {
      KinVertex* parent = vertex2->m_parent;
      KinEdge* p_edge = parent->findEdge(vertex2);

//			int dof_id = p_edge->Cycle_DOF_id;
      int dof_id = p_edge->getDOF()->getCycleIndex();
      if (dof_id!=-1) { // this edge is a DOF
        //if(dof_id==20)
        //	cout<<"KinEdge "<<p_edge<<endl;
        /*
        Atom* ea1 = p_edge->getBond()->Atom1;
        Atom* ea2 = p_edge->getBond()->m_atom2;

        // Jacobian_entry is the derivative of the vertices of the hydrogen bond
        // Now, we also have to calculate the derivatives of the other neighboring atoms!
        Math3D::Vector3 derivativeP1 = ComputeJacobianEntry(ea1->m_position,ea2->m_position,p1); //a
        Math3D::Vector3 derivativeP2 = ComputeJacobianEntry(ea1->m_position,ea2->m_position,p2); //b
        //Math3D::Vector3 derivativeP1_prev = ComputeJacobianEntry(ea1->m_position,ea2->m_position,p1_prev);
        Math3D::Vector3 derivativeP2_prev = ComputeJacobianEntry(ea1->m_position,ea2->m_position,p2_prev);
         */
        Math3D::Vector3 derivativeP1      = p_edge->getDOF()->getDerivative(p1);
        Math3D::Vector3 derivativeP2      = p_edge->getDOF()->getDerivative(p2);
        Math3D::Vector3 derivativeP2_prev = p_edge->getDOF()->getDerivative(p2_prev);

        Math3D::Vector3 jacobianEntryTrans= -0.5*(derivativeP1 + derivativeP2);
//				cout<<"Jac at i="<<dof_id<<", Trans: "<<jacobianEntryTrans<<", Rot1: "<<jacobianEntryRot1<<", Rot2: "<<jacobianEntryRot2<<endl;
        gsl_matrix_set(CycleJacobian,i + 0, dof_id, jacobianEntryTrans.x); //set: Matrix, row, column, what to set
        gsl_matrix_set(CycleJacobian,i + 1, dof_id, jacobianEntryTrans.y);
        gsl_matrix_set(CycleJacobian,i + 2, dof_id, jacobianEntryTrans.z);
        if(bond_ptr->isHBond()) {
          double jacobianEntryRot1 = dot((p1 - p1_prev), (derivativeP2 - derivativeP1));
          double jacobianEntryRot2 = dot((p1 - p2), (derivativeP2 - derivativeP2_prev));
          gsl_matrix_set(CycleJacobian, i + 3, dof_id, jacobianEntryRot1);
          gsl_matrix_set(CycleJacobian, i + 4, dof_id, jacobianEntryRot2);

          ///Matrix to check hBond Rotation
          double hBondEntry = dot((p1_prev - p2_prev), (- derivativeP2_prev));
          gsl_matrix_set(HBondJacobian, hbidx, dof_id, hBondEntry);
        }
      }
      vertex2 = parent;
    }
//    ++i;
    if(bond_ptr->isHBond()){
      i+=5;
      hbidx++;
    }
    else if(bond_ptr->isDBond()) i+=3;
  }

}

//------------------------------------------------------------
//We compute another Jacobian that also considers motion along clash normal directions as a constraint
//--> allowed motions will not move clashing atoms further into each other

/*
void Configuration::ComputeClashAvoidingJacobianAndNullSpace (std::map< std::pair<Atom*,Atom*>,int > allCollisions,bool firstTime, bool projectConstraints) {
//	CTKTimer timer;
//	timer.Reset();
//	double old_time = timer.LastElapsedTime();

  //First, recompute the Jacobian for the current configuration
  //As CycleJacobian is static, the current Jacobian might not correspond to the current configuration
  if(firstTime)//only have to do this in the first clash-prevention trial, after the Jacobian is the same
    computeJacobians();
  //Add the clash constraints
  computeClashAvoidingJacobian(allCollisions,projectConstraints);

  if (JacobianSVD != nullptr) {
    // CycleNullSpace needs to be deleted to free memory.
    if( ClashAvoidingNullSpace!=nullptr ) {
      delete ClashAvoidingNullSpace;
    }
    ClashAvoidingNullSpace = new Nullspace(JacobianSVD);
    ClashAvoidingNullSpace->updateFromMatrix();
  }

//	double new_time = timer.ElapsedTime();
//	clashJacobianTime += new_time - old_time;
}

//---------------------------------------------------------
void Configuration::computeClashAvoidingJacobian (std::map< std::pair<Atom*,Atom*>,int > allCollisions, bool projectConstraints) {

  updateMolecule();

  //The clash Jacobian is the regular Jacobian's constraints, plus one constraint per pair of clashing atoms
  int numCollisions = allCollisions.size();
  int rowNum, colNum;

  //Clashes can occur also for previously free dihedrals!
  //Therefore, we use the full set of dihedrals to determine this matrix!

  if( CycleJacobian != nullptr){
    rowNum = CycleJacobian->size1 + numCollisions;
    colNum = m_molecule->m_spanningTree->getNumDOFs();
  }
  else{
    rowNum = numCollisions;
    colNum = m_molecule->m_spanningTree->getNumDOFs();
  }

  if(ClashAvoidingJacobian==nullptr){
    ClashAvoidingJacobian = gsl_matrix_calloc(rowNum,colNum);
    JacobianSVD = SVD::createSVD(ClashAvoidingJacobian);//new SVDMKL(ClashAvoidingJacobian);
  }else if(ClashAvoidingJacobian->size1==rowNum && ClashAvoidingJacobian->size2==colNum){
    gsl_matrix_set_zero(ClashAvoidingJacobian);
  }else{
    gsl_matrix_free(ClashAvoidingJacobian);
    ClashAvoidingJacobian = gsl_matrix_calloc(rowNum,colNum);
    delete JacobianSVD;
    JacobianSVD = SVD::createSVD(ClashAvoidingJacobian);//new SVDMKL(ClashAvoidingJacobian);
  }

  //Convert the cycle Jacobian to a full Jacobian
  //Columns correspond to cycle_dof_ids
  if(projectConstraints){
//		for (vector<KinEdge*>::iterator eit=m_molecule->m_spanningTree->Edges.begin(); eit!=m_molecule->m_spanningTree->Edges.end(); ++eit) {
    for (auto const& edge: m_molecule->m_spanningTree->Edges){
      int dof_id = edge->getDOF()->getIndex();
      int cycle_dof_id = edge->getDOF()->getCycleIndex();
      if ( cycle_dof_id!=-1 ) {
        for( int i=0; i!=CycleJacobian->size1; i++){
          gsl_matrix_set(ClashAvoidingJacobian, i, dof_id, gsl_matrix_get(CycleJacobian,i,cycle_dof_id));
        }
      }
    }
  } // else: we maintain a zero-valued Jacobian to not consider the constraints (only for testing of constraint limitations)

  int i=CycleJacobian->size1;//start entries after all cycles

  //Quick trick for plotting the clash-cycle network
//	int consCounter = 1;
//	vector<Atom*>::iterator ait;
//	for(ait=protein->Atom_list.begin(); ait != protein->Atom_list.end(); ait++){
//		(*ait)->m_assignedBiggerRB_id = consCounter;
//	}

  for (std::map< std::pair<Atom*,Atom*>,int >::const_iterator mit=allCollisions.begin(); mit!=allCollisions.end(); ++mit) {
//		consCounter++;
    Atom* atom1 = mit->first.first;
    Atom* atom2 = mit->first.second;
    log("dominik") << "Using clash constraint for atoms: "<<atom1->getId() << " " << atom2->getId() << endl;

    Coordinate p1 = atom1->m_position; //end-effector, position 1
    Coordinate p2 = atom2->m_position; //end-effector, position 2

    Math3D::Vector3 clashNormal = p2-p1;
    clashNormal.getNormalized(clashNormal);

    //Vertices
    KinVertex* vertex1 = atom1->getRigidbody()->getVertex();
    KinVertex* vertex2 = atom2->getRigidbody()->getVertex();
    KinVertex* common_ancestor = m_molecule->m_spanningTree->findCommonAncestor(vertex1, vertex2);

    // trace back until the common ancestor from vertex1
    while ( vertex1 != common_ancestor ) {
      KinVertex* parent = vertex1->m_parent;
      KinEdge* p_edge = parent->findEdge(vertex1);

      int dof_id = p_edge->getDOF()->getIndex();
      if (dof_id!=-1) { // this edge is a DOF

//				Math3D::Vector3 derivativeP1 = ComputeJacobianEntry(p_edge->getBond()->Atom1->m_position,p_edge->getBond()->m_atom2->m_position,p1);
        Math3D::Vector3 derivativeP1 = p_edge->getDOF()->getDerivative(p1);
        double jacobianEntryClash = dot(clashNormal, derivativeP1);

        gsl_matrix_set(ClashAvoidingJacobian,i,dof_id,jacobianEntryClash); //set: Matrix, row, column, what to set
//				log("dominik")<<"Setting clash entry on left branch at "<<dof_id<<endl;
      }

      //Quick trick for plotting the clash cycles (don't use in real sampling)
//			for (ait=vertex1->m_rigidbody->Atoms.begin(); ait !=vertex1->m_rigidbody->Atoms.end(); ait++){
//				Atom* atom = (*ait);
//				atom->m_assignedBiggerRB_id = consCounter;
//			}
      vertex1 = parent;
    }

    // trace back until the common ancestor from vertex2
    while ( vertex2 != common_ancestor ) {
      KinVertex* parent = vertex2->m_parent;
      KinEdge* p_edge = parent->findEdge(vertex2);

      int dof_id = p_edge->getDOF()->getCycleIndex();
      if (dof_id!=-1) { // this edge is a DOF

//				Math3D::Vector3 derivativeP2 = ComputeJacobianEntry(
//            p_edge->getBond()->Atom1->m_position,
//            p_edge->getBond()->m_atom2->m_position,
//            p2); //b
        Math3D::Vector3 derivativeP2 = p_edge->getDOF()->getDerivative(p2);
        double jacobianEntryClash = - dot(clashNormal, derivativeP2);

        gsl_matrix_set(ClashAvoidingJacobian,i,dof_id,jacobianEntryClash); //set: Matrix, row, column, what to set
//				log("dominik")<<"Setting clash entry on right branch at "<<dof_id<<endl;
      }

//			//Quick trick for plotting the clash cycles (don't use in real sampling)
//			for (ait=vertex2->m_rigidbody->Atoms.begin(); ait !=vertex2->m_rigidbody->Atoms.end(); ait++){
//				Atom* atom = (*ait);
//				atom->m_assignedBiggerRB_id = consCounter;
//			}
      vertex2 = parent;
    }
    ++i;
  }
}

 */

/// This function converts a vector with all dof entries to a vector with cycle dofs only, by extracting correct entries
void Configuration::convertAllDofsToCycleDofs( gsl_vector *cycleDofs, gsl_vector *allDofs){

  Molecule* M = getMolecule();

  for (auto const& edge: M->m_spanningTree->Edges){
    int dof_id = edge->getDOF()->getIndex();
    int cycle_dof_id = edge->getDOF()->getCycleIndex();
    if ( cycle_dof_id!=-1 ) {
      gsl_vector_set(cycleDofs,cycle_dof_id,gsl_vector_get(allDofs,dof_id));
    }
  }

}

void Configuration::convertCycleDofsToAllDofs( gsl_vector *allDofsAfter, gsl_vector *cycleDofs, gsl_vector *allDofsBefore){

  Molecule* M = getMolecule();

  // Convert back to full length DOFs vector
  for( auto const& edge: M->m_spanningTree->Edges){
    int dof_id = edge->getDOF()->getIndex();
    int cycle_dof_id = edge->getDOF()->getCycleIndex();
    if ( cycle_dof_id!=-1 ) {
      gsl_vector_set(allDofsAfter,dof_id,gsl_vector_get(cycleDofs,cycle_dof_id));
    }
    else if ( dof_id!=-1 ) {
      if (allDofsBefore == nullptr)
        gsl_vector_set(allDofsAfter,dof_id,0);
      else
        gsl_vector_set(allDofsAfter,dof_id,gsl_vector_get(allDofsBefore,dof_id));
    }
  }
}

void Configuration::projectOnCycleNullSpace (gsl_vector *to_project, gsl_vector *after_project) {
  Nullspace* N = getNullspace();
  //Since we're only using this for converting Cycle-DOF ids the mol doesnt have to be updated
  Molecule* M = getMolecule();

  if(N==nullptr){
    gsl_vector_memcpy(after_project, to_project);
    return;
  }

  if( to_project->size > N->getNumDOFs() ) {
    // The input vectors contain all DOFs, however, the null space only contains DOFs in cycles.
    // Convert the DOFs in the input vectors to DOFs in cycles.

    gsl_vector *to_proj_short = gsl_vector_calloc(N->getNumDOFs());
    convertAllDofsToCycleDofs(to_proj_short, to_project);

    // Project onto the null space
    double normBefore = gsl_vector_length(to_proj_short);
    gsl_vector *after_proj_short = gsl_vector_calloc(N->getNumDOFs());
    N->projectOnNullSpace(to_proj_short, after_proj_short);
    double normAfter = gsl_vector_length(after_proj_short);

    //Scale projected gradient to same norm as unprojected
    if(normAfter>0.0000001)
      gsl_vector_scale(after_proj_short, normBefore/normAfter);

    // Convert back to full length DOFs vector
    convertCycleDofsToAllDofs(after_project,after_proj_short,to_project);

    gsl_vector_free(to_proj_short);
    gsl_vector_free(after_proj_short);
  }
  else {
    double normBefore = gsl_vector_length(to_project);
    N->projectOnNullSpace(to_project, after_project);
    double normAfter = gsl_vector_length(after_project);
    gsl_vector_scale(after_project, normBefore/normAfter);
  }
}

Molecule * Configuration::getMolecule() const
{
  return m_molecule;
}

Molecule * Configuration::updatedMolecule()
{
  m_molecule->setConfiguration(this);
  return m_molecule;
}

void Configuration::updateMolecule()
{
  m_molecule->setConfiguration(this);
}

gsl_matrix* Configuration::getCycleJacobian()
{
  computeJacobians();
  return CycleJacobian;
}

Nullspace* Configuration::getNullspace()
{
  if(nullspace==nullptr){
    computeCycleJacobianAndNullSpace();
  }

  return nullspace;
}

Configuration* Configuration::getParent()
{
  return m_parent;
}

std::list<Configuration*>& Configuration::getChildren()
{
  return m_children;
}
