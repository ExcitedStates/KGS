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

#include "NullspaceMove.h"

#include <math/gsl_helpers.h>
#include <gsl/gsl_vector_double.h>
#include "Logger.h"

using namespace std;

NullspaceMove::NullspaceMove(double maxRotation):
    Move(maxRotation)
{}

Configuration* NullspaceMove::performMove(Configuration* current, gsl_vector* gradient) {
  double currNorm = gsl_vector_length(gradient);
  log("planner") << "Norm of gradient: " << currNorm << endl;

  // Project the gradient onto the null space of current
  gsl_vector *projected_gradient = gsl_vector_calloc(current->getNumDOFs());
  current->projectOnCycleNullSpace(gradient, projected_gradient);

//  double currProjNorm = gsl_vector_length(projected_gradient);
//  log("planner") << "Norm of projected gradient: " << currProjNorm << endl;

  //Scale max entry
  if(m_scale)
    gsl_vector_scale_max_component(projected_gradient,m_maxRotation);

  Configuration *new_q = new Configuration(current);
  for (int i = 0; i < new_q->getNumDOFs(); ++i)
    new_q->m_dofs[i] = formatRangeRadian(current->m_dofs[i] + gsl_vector_get(projected_gradient, i));
  gsl_vector_free(projected_gradient);

  return new_q;
}

//gsl_vector* NullspaceMove::projectOnNullspace(Configuration* conf, gsl_vector* gradient, gsl_vector* ret)
//{
//  if(ret==nullptr) ret = gsl_vector_alloc(gradient->size);
//  else assert(gradient->size==ret->size);
//
//  // No cycles
//  if(!conf->getNullspace()) {
//    gsl_vector_memcpy(ret, gradient);
//    return ret;
//  }
//
//  if( gradient->size >conf->getNullspace()->getNumDOFs() ) {
//    // The input vectors contain all DOFs, however, the null space only contains DOFs in cycles.
//    // Convert the DOFs in the input vectors to DOFs in cycles.
//
//    //Fill reduced gradient
//    gsl_vector *to_proj_short = gsl_vector_calloc(conf->getNullspace()->getNumDOFs());
//    for (auto const& edge: conf->getMolecule()->m_spanningTree->m_edges){
//      int dof_id = edge->getDOF()->getIndex();
//      int cycle_dof_id = edge->getDOF()->getCycleIndex();
//      if ( cycle_dof_id!=-1 ) {
//        gsl_vector_set(to_proj_short,cycle_dof_id,gsl_vector_get(gradient, dof_id));
//      }
//    }
//
//    //Project reduced gradient
//    double normBefore = gsl_vector_length(to_proj_short);
//    gsl_vector *after_proj_short = gsl_vector_calloc(conf->getNullspace()->getNumDOFs());
//    conf->getNullspace()->projectOnNullSpace(to_proj_short, after_proj_short);
//    double normAfter = gsl_vector_length(after_proj_short);
//
//    //Scale so the norm is the same as before projection
//    gsl_vector_scale(after_proj_short, normBefore/normAfter);
//
//    //Convert back to full length DOFs vector
//    for( auto const& edge: conf->getMolecule()->m_spanningTree->m_edges){
//      int dof_id = edge->getDOF()->getIndex();
//      int cycle_dof_id = edge->getDOF()->getCycleIndex();
//      if ( cycle_dof_id!=-1 ) {
//        gsl_vector_set(ret, dof_id,gsl_vector_get(after_proj_short,cycle_dof_id));
//      }
//      else if ( dof_id!=-1 ) {
//        gsl_vector_set(ret, dof_id,gsl_vector_get(gradient,dof_id));
//      }
//    }
//    gsl_vector_free(to_proj_short);
//    gsl_vector_free(after_proj_short);
//  } else {
//    double normBefore = gsl_vector_length(gradient);
//    conf->getNullspace()->projectOnNullSpace(gradient, ret);
//    double normAfter = gsl_vector_length(ret);
//    gsl_vector_scale(ret, normBefore/normAfter);
//  }
//  return ret;
//
//}
