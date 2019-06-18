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



#include "DihedralDirection.h"

#include <applications/options/ExploreOptions.h>
#include <cassert>
#include <cmath>
#include "core/Molecule.h"

DihedralDirection::DihedralDirection(Selection& resNetwork):
    m_resNetwork(resNetwork)
{}


void DihedralDirection::computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret)
{
  assert(target!=nullptr);
  Molecule * protein = conf->updatedMolecule();
//  const std::vector<int>& resNetwork = ExploreOptions::getOptions()->residueNetwork;
//  bool allResidues = resNetwork.size() == 0 ? true:false;

  for (auto const& edge: protein->m_spanningTree->m_edges){
    if( edge->getBond()==nullptr ) continue; //This excludes global dofs from being computed here
    if( !m_resNetwork.inSelection(edge->getBond()) ) continue;

    int dofId = edge->getDOF()->getIndex();

    double angle_diff = target->getGlobalTorsion(dofId) - conf->getGlobalTorsion(dofId);
    angle_diff = formatRangeRadian(angle_diff);

    gsl_vector_set(ret,dofId,angle_diff);

    assert( !std::isnan(angle_diff) );

//    if(allResidues){//gradient for all residues
//      gsl_vector_set(ret,dofId,angle_diff);
//    }
//    else if( std::find(resNetwork.begin(), resNetwork.end(), resId)!= resNetwork.end() ){//only for specified residues
//      gsl_vector_set(ret,dofId,angle_diff);
//    }
  }
}
