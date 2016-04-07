
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

	addMove( new NullspaceMove(), nullspace_move_frequency );
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
