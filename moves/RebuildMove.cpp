#include "RebuildMove.h"

#include "core/Chain.h"
#include "core/Molecule.h"


RebuildMove::RebuildMove(int fragmentLength, int aggression):
    m_fragmentLength(fragmentLength),
    m_rebuildAggression(aggression)
//		m_fragmentLength(SamplingOptions::getOptions()->rebuild_fragment_length),
//		m_rebuildAggression(SamplingOptions::getOptions()->rebuildAggression)
{}

Configuration* RebuildMove::performMove(Configuration* current, gsl_vector*)
{
	Molecule * m_protein = current->updatedMolecule();

	//Find random free segment
	int freeCount = 0, i, j;
	Chain* chain = m_protein->chains[0];
	int residues = chain->getResidues()[chain->getResidues().size()-1]->getId()+1;
	for(i=0;i<residues;i++)
		if(m_protein->residueAnnotations[i]==0) freeCount++;

  int fragmentLength = m_fragmentLength;

	int randFree = (int)(Random01()*freeCount);
	//cout<<"randFree: "<<randFree<<endl;
	int freeIdx=0;//Reusing variable
	for(i=0;i<residues;i++){
		if(m_protein->residueAnnotations[i]==0) {
			if(randFree==freeIdx){// This is our fragment of interest
				int nextNonFree,prevNonFree;
				for(nextNonFree=i;;nextNonFree++)
					if(m_protein->residueAnnotations[nextNonFree]!=0) break;
				for(prevNonFree=i;;prevNonFree--)
					if(m_protein->residueAnnotations[prevNonFree]!=0) break;
				int dPrev = i-prevNonFree;
				int dNext = nextNonFree-i;
				if(nextNonFree-prevNonFree-1 < fragmentLength)
					fragmentLength = nextNonFree-prevNonFree-1;
				if(dNext>=dPrev){
					if(dNext<m_fragmentLength) i-=m_fragmentLength-dNext;
				}else{
					if(dPrev<m_fragmentLength) i=prevNonFree+1;
					else i-=m_fragmentLength-1;
				}

				break;
			}
			freeIdx++;
		}
	}
	//i is now the index of the first residue we want to rebuild
	j = i+fragmentLength;

	Configuration* new_q = nullptr;
	int count = 0;
	while(count++<20 && new_q==nullptr) {
		new_q = m_protein->resampleSugars(i, j, current, m_rebuildAggression);
	}

	return new_q;

}
