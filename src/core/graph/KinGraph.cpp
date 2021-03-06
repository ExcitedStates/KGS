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
	for (auto it=m_vertices.begin(); it!=m_vertices.end(); ++it) {
		delete it->second;
	}
	for (auto &edge: m_edges){
		delete edge;
	}
}

KinVertex* KinGraph::addVertex(Rigidbody* rb){
	KinVertex* new_vertex = new KinVertex(rb);
  if(rb!=nullptr)
    m_vertices[rb->id()] = new_vertex;

	return new_vertex;
}


KinVertex* KinGraph::getVertex (int rb_id) {
  auto vmapIt = m_vertices.find(rb_id);
  if( vmapIt==m_vertices.end() ) return nullptr;
  return vmapIt->second;
}


void KinGraph::print () {
	for (auto it=m_vertices.begin(); it!=m_vertices.end(); ++it) {
		KinVertex *vertex = it->second;
//		log() << "KinVertex " << vertex->id << ": " << vertex->m_rigidbody->Atoms.size() << " atoms and " << vertex->m_edges.size() << " m_edges" << endl;
		vertex->print();
		for (auto eit=vertex->m_edges.begin(); eit!=vertex->m_edges.end(); ++eit)
			(*eit)->print();
		log() << endl;
	}

    log() << "Total number of edges = " << m_edges.size()/2 << endl;
    log() << "Total number of rigid bodies = " << m_vertices.size() << endl;
}


