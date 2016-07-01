#include "metrics/RMSD.h"
#include "../SamplingOptions.h"
#include "core/Chain.h"

using namespace std;

namespace metrics{

RMSD::RMSD():
  m_atomsRMSD(nullptr)
{}

RMSD::RMSD(const vector<Atom*>* atomsRMSD):
  m_atomsRMSD(atomsRMSD)
{}

double RMSD::distance(Configuration* c1, Configuration* c2)
{
  const std::vector<Atom*>* atomsRMSD = SamplingOptions::getOptions()->getAtomsMoving();

  // If atomsAlign is nullptr, align the entire m_protein
  if (atomsRMSD==nullptr) {
    cout<<"Found no atoms to align .. using all"<<endl;
    atomsRMSD = &(c1->getMolecule()->atoms);//choose all atoms
  }

  vector<Coordinate> p1_atoms;
  vector<Coordinate> p2_atoms;

  double sum=0;
  int atom_size = atomsRMSD->size();
  int resId;
  std::string name;
  std::string chainName;

  Molecule * p1 = c1->updatedMolecule();
  for (auto const& aIt : *atomsRMSD ) {
    name = aIt->getName();
    chainName = aIt->getResidue()->getChain()->getName();
    resId = aIt->getResidue()->getId();
    p1_atoms.push_back(p1->getAtom(chainName,resId, name)->m_Position);
  }

  Molecule * p2 = c2->updatedMolecule();
  for (auto const& aIt : *atomsRMSD ) {
    name = aIt->getName();
    chainName = aIt->getResidue()->getChain()->getName();
    resId = aIt->getResidue()->getId();
    p2_atoms.push_back(p2->getAtom(chainName,resId, name)->m_Position);
  }

  for(int i=0;i!=atom_size;i++){
    sum += p1_atoms[i].distanceSquared(p2_atoms[i]);
  }

  return sqrt(sum/atom_size);
}
