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


#include "CompositeMove.h"
#include "math/MathUtility.h"
#include "moves/NullspaceMove.h"
#include "moves/RebuildMove.h"

using namespace std;

CompositeMove::CompositeMove():
	total_weight(0.0)
{
	double nullspace_move_frequency = 1.0;

	if(SamplingOptions::getOptions()->rebuild_frequency>0.0){
		nullspace_move_frequency -= SamplingOptions::getOptions()->rebuild_frequency;

		addMove( new RebuildMove(), SamplingOptions::getOptions()->rebuild_frequency );
	}

	addMove( new NullspaceMove(SamplingOptions::getOptions()->maxRotation), nullspace_move_frequency );
}
CompositeMove::~CompositeMove()
{
	for( pair< Move*, double > &child : child_moves ){
		delete child.first;
	}
}

void CompositeMove::addMove(Move* m, double weight)
{
	child_moves.push_back( std::make_pair(m, weight) );
	total_weight+=weight;
}

Configuration* CompositeMove::performMove(Configuration* current, gsl_vector* gradient)
{
	double sampled_val = Random01() * total_weight;

	for( pair< Move*, double > &child : child_moves ){
		sampled_val -= child.second;

		if(sampled_val<0){
			return child.first->move(current, gradient);
		}
	}

	throw "CompositeMove::performMove - Error: Shouldnt reach this point";
}

void CompositeMove::setScalingFlag(bool scale) {

  for( pair< Move*, double > &child : child_moves ) {
    child.first->setScalingFlag(scale);
  }
}
