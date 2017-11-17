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
  for(KinEdge*& edge: m_protein->m_spanningTree->m_edges){
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
