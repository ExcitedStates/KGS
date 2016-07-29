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

#include "MSDDirection.h"
#include "metrics/RMSD.h"
#include "SamplingOptions.h"
#include "core/Molecule.h"
#include "core/Chain.h"

#include <stack>
#include <cassert>
#include <cmath>
#include <gsl/gsl_vector_double.h>

using namespace std;

MSDDirection::MSDDirection(Selection& resNetwork, bool alignAlways):
    m_preprocessedTree(nullptr),
    m_resNetwork(resNetwork),
    m_rmsd(resNetwork),
    m_alignAlways(alignAlways)
{
}

void MSDDirection::computeGradient(Configuration* conf, Configuration* c_target, gsl_vector* ret) {
  //Gradient is the vector to be filled, gradientMethod decides whether blending with randomized entries takes place
  //atom_selection defines the atoms considered for the distance gradient
  Molecule *protein=conf->updatedMolecule();
  Molecule *target=c_target->updatedMolecule();
  if( target == nullptr ) {
    std::cerr << "MSDDirection::computeGradient - Target not set" << std::endl;
    exit(-1);
  }
  if( target == protein ) {
    std::cerr << "MSDDirection::computeGradient - Source and target molecule must differ" << std::endl;
    exit(-1);
  }

  if( m_alignAlways ) {
    m_rmsd.align(target, protein);
  }

  vector<double> gradient(ret->size, 0.0);
  vector<int> counts(ret->size, 0);
//  int count=0;

  for( auto const &atom: m_resNetwork.getSelectedAtoms(protein) ) {


    //get corresponding target atom (if existing) - not based on ID
    Atom *aTarget=target->getAtom(atom->getResidue()->getChain()->getName(),
                                  atom->getResidue()->getId(),
                                  atom->getName());
    if( aTarget==nullptr ) continue;

    Math3D::Vector3 diff = aTarget->m_position - atom->m_position;

    //Compute gradient contribution all the way to the root
    KinVertex* v = atom->getRigidbody()->getVertex();
    while(v!=protein->m_spanning_tree->m_root){
      KinEdge* parentEdge = v->m_parent->findEdge(v);
      int dof_id = parentEdge->getDOF()->getIndex();
      Math3D::Vector3 deriv = parentEdge->getDOF()->getDerivative(atom->m_position);

      gradient[dof_id] += deriv.dot(diff);
      assert( !std::isnan(gradient[dof_id]) );
      counts[dof_id]++;

      v = v->m_parent;
    }
//    count++;
  }

  //Scaling with appropriate pre-factor from derivative of MSD
  for(int i=0;i<ret->size;i++){
    if(counts[i]>0)
      gradient[i] *= 2.0 / counts[i];
//    gradient[i] *= 2.0 / count;
//    gradient[i] = formatRangeRadian(gradient[i]);
  }
  std::copy(&gradient[0], &gradient[ret->size], ret->data);
  //double factor = 2.0 / count;
  //gsl_vector_scale(ret, factor);
}

/* //This all assumes dihedral DOFs only. With other types of DOFs we need a more straightforward approach
void MSDDirection::computeGradient(Configuration* conf, Configuration* c_target, gsl_vector* ret)
{
//Gradient is the vector to be filled, gradientMethod decides whether blending with randomized entries takes place
//atom_selection defines the atoms considered for the distance gradient
Molecule * protein = conf->updatedMolecule();
Molecule * target = c_target->updatedMolecule();
if(target==nullptr){
  std::cerr<<"MSDDirection::computeGradient - Target not set"<<std::endl;
  exit(-1);
}
if(target==protein){
  std::cerr<<"MSDDirection::computeGradient - Source and target m_molecule must differ"<<std::endl;
  exit(-1);
}

SamplingOptions *options;
options = SamplingOptions::getOptions();
if(options->alignAlways){
  metrics::RMSD::align(target, protein);
}

//TODO: change this to a solid selection using options "gradientSelection" and "selectAtoms"
const vector<int>& resNetwork = options->residueNetwork;
bool allResidues = resNetwork.size() == 0 ? true:false;
string atom_selection="HEAVY"; // "CA" "ALL" "BACKBONE" "RES" "RESCA" "HEAVY" "RESHEAVY"
if(!allResidues){
  atom_selection="RESHEAVY";
}

//Ensure that m_sortedVertices is filled with vertices from `protein`.
if(m_preprocessedTree==nullptr || protein->m_spanning_tree!=m_preprocessedTree){
  m_sortedVertices.clear();
  collectVerticesPostorder(protein->m_spanning_tree->m_root);
  m_preprocessedTree = protein->m_spanning_tree;
}

Math3D::Vector3 f, g, zeros(0.0);
stack<Math3D::Vector3> stackF, stackG;
int count=0;

//Now we compute the MSD gradient with a "fast" implementation based on Abé's paper from 1984
//  for ( vit = orderedVertices->begin(); vit != orderedVertices->end(); vit++ ){//traverse tree by vertices starting at the final leaves
for ( KinVertex* const& currVertex: m_sortedVertices ){
  //traverse tree by vertices in a bottom-up fashion (starting at leaves)

  //KinVertex *currVertex = vit->second;
  if(currVertex == protein->m_spanning_tree->m_root)
    break;

  KinEdge* currEdge = currVertex->m_parent->findEdge(currVertex);
  Bond * bond_ptr = currEdge->getBond();
  Atom* atom1 = bond_ptr->Atom1;
  Atom* atom2 = bond_ptr->Atom2;

  f.setZero();
  g.setZero();

  for (vector<Atom*>::iterator ait=currVertex->m_rigidbody->Atoms.begin(); ait!=currVertex->m_rigidbody->Atoms.end(); ait++ ){
    Atom* atom = *ait;
    //Don't include the bond atoms in the calculation, not necessary
    if(atom == atom1 || atom == atom2)
      continue;
    //Check selection
//      if (atom_selection=="CA" && ( atom->getName()!="CA" || atom->getName()!="C1'") )
//        continue;
//      else if (atom_selection=="BACKBONE" && ( atom->getName()!="N" || atom->getName()!="C" || atom->getName()!="CA" || atom->getName()!="O" ||
//                                               atom->getName()!="C1'" || atom->getName()!="C2'" || atom->getName()!="C3'" || atom->getName()!="C4'" || atom->getName()!="C5'" || atom->getName()!="O5'" || atom->getName()!="O3'" || atom->getName()!="P" )
//          )
//        continue;
    else if (atom_selection=="HEAVY" && !atom->isHeavyAtom())
      continue;
//      else if (atom_selection=="RES" && find(resNetwork.begin(), resNetwork.end(), atom->getResidue()->getId()) == resNetwork.end() )//Atom not in residue selection
//        continue;
//      else if (atom_selection=="RESCA" && ( find(resNetwork.begin(), resNetwork.end(), atom->getResidue()->getId()) == resNetwork.end() || atom->getName()!="CA" || atom->getName()!="C4'" ) )//C_alpha atom not in residue selection
//        continue;
    else if (atom_selection=="RESHEAVY" && ( find(resNetwork.begin(), resNetwork.end(), atom->getResidue()->getId()) == resNetwork.end() || !atom->isHeavyAtom() ))//HEAVY atom not in residue selection
      continue;
    else{//all atoms
      //get corresponding target atom (if existing), not based on ID's anymore
      Atom* aTarget = target->getAtom(atom->getResidue()->getChain()->getName(),
                                      atom->getResidue()->getId(),
                                      atom->getName() );
      if(aTarget == nullptr ){//|| aTarget->getResidue()->getProperName() != atom->getResidue()->getProperName()){
        continue;//skip the non-existing atom
      }
      Math3D::Vector3 distance = aTarget->m_position - atom->m_position;
      if(distance != zeros ){
        g += distance;
        Math3D::Vector3 crossP = cross(aTarget->m_position,atom->m_position);
        f += crossP;
        ++count;
      }
    }
  }
  if( currVertex->m_edges.size() == 0){//push on stack
    stackF.push(f);
    stackG.push(g);
  }
  else if(currVertex->m_edges.size() == 1){//add with correct existing stack entry
    stackF.top() += f;
    stackG.top() += g;
  }
  else{
    for(int i=2; i <= currVertex->m_edges.size(); i++){
      //add two and pop for each edge if equal or more than 2
      //this is necessary as we specify h-bonds as additional covalent bonds (stack is not limited to three as in Abé's paper)
      f += stackF.top();
      stackF.pop();
      stackF.top() += f;
      f.setZero();

      g += stackG.top();
      stackG.pop();
      stackG.top() += g;
      g.setZero();
    }
  }
  //Check that only "active" m_edges are used in the gradient, where a torsion can be defined properly
  //Inactive m_edges are currently disabled (see IO read rigid body), so not necessary
  Atom* atom3 = nullptr;
  vector<Atom*>::iterator ait;
  for (ait=atom1->Cov_neighbor_list.begin(); ait!=atom1->Cov_neighbor_list.end(); ++ait) {
    if ( (*ait)->getId()==atom2->getId() ) continue;
    if ( atom3==nullptr || (*ait)->getId()<atom3->getId() ) {
      atom3 = *ait;
    }
  }
  Atom* atom4 = nullptr;
  for (ait=atom2->Cov_neighbor_list.begin(); ait!=atom2->Cov_neighbor_list.end(); ++ait) {
    if ( (*ait)->getId()==atom2->getId() ) continue;
    if ( atom4==nullptr || (*ait)->getId()<atom4->getId() ) {
      atom4 = *ait;
    }
  }
  if(atom3 != nullptr && atom4 != nullptr){
    Math3D::Vector3 r_i = atom2->m_position - atom1->m_position;
    r_i.setNormalized(r_i);
    Math3D::Vector3 rp_cross = cross(r_i,atom1->m_position);

    double deltaQ_i = -r_i.dot(stackF.top() ) - rp_cross.dot( stackG.top() ); //final formula without pre-factor
    gsl_vector_set(ret, currEdge->getDOF()->getIndex(), deltaQ_i);
  }
  else{
    cerr<<"Found inactive torsion"<<endl;
  }
}

  double factor = 2.0 / count;
  gsl_vector_scale(ret, factor);//scaling with appropriate pre-factor from derivative of MSD
  for (int i=0; i<ret->size; ++i){
    gsl_vector_set(ret, i,formatRangeRadian(gsl_vector_get(ret, i)));//format range
  }

}
 */

void MSDDirection::collectVerticesPostorder(KinVertex* v)
{
  //First add children
  for(KinEdge* const& edge: v->m_edges){
    collectVerticesPostorder(v);
  }

  //Next, add v itself
  if(v->m_parent!=nullptr)
    m_sortedVertices.push_back(v);
}
