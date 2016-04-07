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

#include "DihedralDirection.h"

#include <SamplingOptions.h>
#include "core/Molecule.h"

DihedralDirection::DihedralDirection():
    m_target(NULL)
{}

void DihedralDirection::setTarget(Configuration* conf)
{
  m_target = conf;
}

void DihedralDirection::computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret)
{
  Molecule * protein = conf->updatedProtein();
  const std::vector<int>& resNetwork = SamplingOptions::getOptions()->residueNetwork;
  bool allResidues = resNetwork.size() == 0 ? true:false;

  if(m_target->getGlobalTorsions() == NULL){
    std::cerr<<"DihedralDirection::computeGradient - No global torsions, please calculate them first!"<<std::endl;
    exit(-1);
  }

  for (auto const& edge: protein->m_spanning_tree->Edges){
    int dofId = edge->DOF_id;
    int resId = edge->getBond()->Atom1->getResidue()->getId();
    double angle_diff = m_target->getGlobalTorsions(dofId) - conf->getGlobalTorsions(dofId);
    angle_diff = formatRangeRadian(angle_diff);

    if(allResidues){//gradient for all residues
      gsl_vector_set(ret,dofId,angle_diff);
    }
    else if( std::find(resNetwork.begin(), resNetwork.end(), resId)!= resNetwork.end() ){//only for specified residues
      gsl_vector_set(ret,dofId,angle_diff);
    }
  }
}
