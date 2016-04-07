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
#include <iostream>
#include <stack>
#include <list>
#include <queue>

#include "KinGraph.h"
#include "SugarVertex.h"
#include "Logger.h"

using namespace std;


KinGraph::KinGraph () {
}

KinGraph::~KinGraph () {
	map<unsigned int, KinVertex*>::iterator it;
	m_sortedVertices.clear();
	for (it=Vertex_map.begin(); it!=Vertex_map.end(); ++it) {
		delete it->second;
	}

}

KinVertex* KinGraph::addVertex (unsigned int rb_id, Rigidbody* rb, bool flexibleRibose){
    //log("debugRas")<<"addVertex("<<rb_id<<" , "<<rb<<" )"<<endl;
	KinVertex* new_vertex;
	//if( (rb->getAtom("O4'") && rb->getAtom("C3'")) || (rb->getAtom("CB") && rb->getAtom("CB")->getResidue()->getName()=="PRO") ){
  if( flexibleRibose && (rb->getAtom("O4'") && rb->getAtom("C3'")) ){
		new_vertex = new SugarVertex(rb_id,rb);
	}else{
		new_vertex = new KinVertex(rb_id,rb);
	}
	//Vertex_map.insert(make_pair(rb_id,new_vertex));
	Vertex_map[rb_id] = new_vertex;
	return new_vertex;
}

Edge* KinVertex::findEdge(KinVertex* v) const
{
  for(auto const& edge: edges){
    if( edge->EndVertex==v )
      return edge;
  }
  return NULL;
}

// Add a directed edge from rb_id1 to rb_id2
Edge* KinGraph::addEdgeDirected (KinVertex *vertex1, KinVertex *vertex2, Bond * bond, int DOF_id)
{
    //log("debugRas")<<"KinGraph::addEdgeDirected("<<vertex1->Rb_ptr<<", "<<vertex2->Rb_ptr<<", "<<bond<<"..)"<<endl;
	Atom *atom2, *atom3, *atom4;
	Bond *bond_copy = new Bond(*bond);
	atom2 = bond_copy->Atom1;
	atom3 = bond_copy->Atom2;
	atom4 = NULL;

	// Find out the atom that covalently bonded to atom3 with smallest Id. It participates in the definition of the torsional angle.
	for (vector<Atom*>::iterator aitr=atom3->Cov_neighbor_list.begin(); aitr!=atom3->Cov_neighbor_list.end(); ++aitr) {
		if ( (*aitr)==atom2 ) continue;
		if ( atom4==NULL || (*aitr)->getId()<atom4->getId() ) {
			atom4 = *aitr;
		}
	}
	for (vector<Atom*>::iterator aitr=atom3->Hbond_neighbor_list.begin(); aitr!=atom3->Hbond_neighbor_list.end(); ++aitr) {
		if ( (*aitr)==atom2 ) continue;
		if ( atom4==NULL || (*aitr)->getId()<atom4->getId() ) {
			atom4 = *aitr;
		}
	}

	// If atom4 is in vertex1, then should flip atom2 and atom3 so that the bond is pointing from vertex1 to vertex2
	Atom *tmp_atom;
	for (vector<Atom*>::iterator svIT = vertex1->Rb_ptr->Atoms.begin(); svIT != vertex1->Rb_ptr->Atoms.end(); ++svIT)
	{
		if((*svIT) == atom4) 
		{
			tmp_atom = bond_copy->Atom1;
			bond_copy->Atom1 = bond_copy->Atom2;
			bond_copy->Atom2 = tmp_atom;
			break;
		}
	}

  //Old and weird
//  Edge *edge1 = new Edge(vertex1,vertex2,bond_copy);
//  edge1->Bond = bond;
	Edge *edge1 = new Edge(vertex1,vertex2,bond);
	vertex1->addEdge(vertex2->id, edge1);
	vertex2->setParent(vertex1);
	Edges.push_back(edge1);
	edge1->DOF_id = DOF_id;
	return edge1;
}

bool KinGraph::hasVertex (int rb_id) {
	return (Vertex_map.find(rb_id)!=Vertex_map.end());
}

void KinGraph::print () {
	map<unsigned int, KinVertex*>::iterator it;
	for (it=Vertex_map.begin(); it!=Vertex_map.end(); ++it) {
		KinVertex *vertex = it->second;
		log() << "Rigidbody " << vertex->id << ": " << vertex->Rb_ptr->Atoms.size() << " atoms and " << vertex->edges.size() << " edges" << endl;
		vertex->print();
		//for (map<unsigned int,Edge*>::iterator eit=vertex->Edges.begin(); eit!=vertex->Edges.end(); ++eit) {
		for (auto eit=vertex->edges.begin(); eit!=vertex->edges.end(); ++eit) {
			(*eit)->print();
		}
		log() << endl;
	}

    log() << "Total number of edges = " << Edges.size()/2 << endl;
    log() << "Total number of rigid bodies = " << Vertex_map.size() << endl;
}


