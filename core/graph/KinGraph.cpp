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
#include "Logger.h"

using namespace std;


KinGraph::KinGraph () {
}

KinGraph::~KinGraph () {
	for (auto it=Vertex_map.begin(); it!=Vertex_map.end(); ++it) {
		delete it->second;
	}
}

KinVertex* KinGraph::addVertex(Rigidbody* rb){
	KinVertex* new_vertex = new KinVertex(rb);
  if(rb!=nullptr)
    Vertex_map[rb->id()] = new_vertex;

	return new_vertex;
}

// Add a directed edge from rb_id1 to rb_id2
KinEdge* KinGraph::addEdgeDirected(KinVertex *vertex1, KinVertex *vertex2, Bond * bond)
{
  if(bond!=nullptr) {
    //log("debugRas")<<"KinGraph::addEdgeDirected("<<vertex1->m_rigidbody<<", "<<vertex2->m_rigidbody<<", "<<bond<<"..)"<<endl;
    Atom *atom2, *atom3, *atom4;
    Bond *bond_copy = new Bond(*bond);
    atom2 = bond_copy->Atom1;
    atom3 = bond_copy->Atom2;
    atom4 = nullptr;

    // Find out the atom that covalently bonded to atom3 with smallest Id. It participates in the definition of the torsional angle.
    for (vector<Atom *>::iterator aitr = atom3->Cov_neighbor_list.begin();
         aitr != atom3->Cov_neighbor_list.end(); ++aitr) {
      if ((*aitr) == atom2) continue;
      if (atom4 == nullptr || (*aitr)->getId() < atom4->getId()) {
        atom4 = *aitr;
      }
    }
    for (vector<Atom *>::iterator aitr = atom3->Hbond_neighbor_list.begin();
         aitr != atom3->Hbond_neighbor_list.end(); ++aitr) {
      if ((*aitr) == atom2) continue;
      if (atom4 == nullptr || (*aitr)->getId() < atom4->getId()) {
        atom4 = *aitr;
      }
    }

    // If atom4 is in vertex1, then should flip atom2 and atom3 so that the bond is pointing from vertex1 to vertex2
    Atom *tmp_atom;
    for (vector<Atom *>::iterator svIT = vertex1->m_rigidbody->Atoms.begin();
         svIT != vertex1->m_rigidbody->Atoms.end(); ++svIT) {
      if ((*svIT) == atom4) {
        tmp_atom = bond_copy->Atom1;
        bond_copy->Atom1 = bond_copy->Atom2;
        bond_copy->Atom2 = tmp_atom;
        break;
      }
    }
  }

	KinEdge *edge1 = new KinEdge(vertex1,vertex2,bond);
	vertex1->addEdge(edge1);
	vertex2->setParent(vertex1);
	Edges.push_back(edge1);
	return edge1;
}

bool KinGraph::hasVertex (int rb_id) {
	return (Vertex_map.find(rb_id)!=Vertex_map.end());
}

KinVertex* KinGraph::getVertex (int rb_id) {
	return Vertex_map.find(rb_id)->second;
}


void KinGraph::print () {
	for (auto it=Vertex_map.begin(); it!=Vertex_map.end(); ++it) {
		KinVertex *vertex = it->second;
//		log() << "KinVertex " << vertex->id << ": " << vertex->m_rigidbody->Atoms.size() << " atoms and " << vertex->m_edges.size() << " m_edges" << endl;
		vertex->print();
		for (auto eit=vertex->m_edges.begin(); eit!=vertex->m_edges.end(); ++eit)
			(*eit)->print();
		log() << endl;
	}

    log() << "Total number of edges = " << Edges.size()/2 << endl;
    log() << "Total number of rigid bodies = " << Vertex_map.size() << endl;
}


