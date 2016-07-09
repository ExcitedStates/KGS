#include "metrics/Dihedral.h"
#include "metrics/Metric.h"
#include "core/Molecule.h"
#include "math/gsl_helpers.h"
#include "Logger.h"
#include <math.h>

using namespace std;

namespace metrics{

Dihedral::Dihedral(Selection& selection):
Metric(selection)
{}


double Dihedral::distance(Configuration* c1, Configuration* c2)
{
  const std::vector<Bond*>& bonds1 = m_selection.getSelectedBonds(c1->getMolecule());
  const std::vector<Bond*>& bonds2 = m_selection.getSelectedBonds(c2->getMolecule());

  if( bonds1.empty() || bonds2.empty() ){
    cerr<<"Dihedral::distance - Selection given to Dihedral metric matched no bonds: "<<m_selection<<endl;
    exit(-1);
  }

  if(bonds1.size() != bonds2.size()){
    cerr<<"Dihedral::distance(..)";
    cerr<<" - Configurations have different number of bonds ("<<bonds1.size()<<" vs "<<bonds2.size()<<")"<<endl;
    exit(-1);
  }

  Molecule * m_protein = c1->getMolecule();
  bool useGlobals = m_protein!=c2->getMolecule();
  //TODO: Extract c1->getGlobalTorsions and c2->getGlobalTorsions to avoid complicated loop and multiple global guards

  int count = 0;
  double distance=0.0;
  for(KinEdge*& edge: m_protein->m_spanning_tree->Edges){
    if(edge->getBond()==nullptr || !m_selection.inSelection(edge->getBond())) continue;
    int dofId = edge->getDOF()->getIndex();
    double angle_diff;
    if(useGlobals) {
      angle_diff = fabs(c2->getGlobalTorsion(dofId) - c1->getGlobalTorsion(dofId));
      if(angle_diff>M_PI)
        angle_diff = 2*M_PI - angle_diff;
    }else{
      angle_diff = fabs(c2->m_dofs[dofId] - c1->m_dofs[dofId]);
      if(angle_diff>M_PI)
        angle_diff = 2*M_PI - angle_diff;
    }

    distance += angle_diff*angle_diff;
    count++;
  }

  return sqrt(distance/count);
}

}
