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
  int DOF = c1->getNumDOFs();

  if(DOF != c2->getNumDOFs()){
    cerr<<"Configurations to compare do not have the same number of DOFs!"<<endl;
    exit(-1);
  }
  //if(c1->getGlobalTorsions()== nullptr || c2->getGlobalTorsions()== nullptr){
  //	cerr<<"Need to set configuration and determine global torsions first!"<<endl;
  //	exit(-1);
  //}


  Molecule * m_protein = c1->getMolecule();
  bool useGlobals = m_protein!=c2->getMolecule();
  if(useGlobals){
    c1->updateGlobalTorsions();
    c2->updateGlobalTorsions();
  }


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

  double ret = sqrt(distance/count);
  return ret;
}

}
