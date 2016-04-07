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
#include <algorithm>
#include <list>
#include <queue>

#include "RigidbodyGraph.h"
#include "SugarVertex.h"
#include "Logger.h"

using namespace std;

Edge::Edge (RigidbodyGraphVertex *startv, RigidbodyGraphVertex *endv, Bond *bond) :
    m_bond(bond)
{
	StartVertex = startv;
	EndVertex = endv;
	Shared = false;
	RigidbodyGraphCycleCluster_id = 0;
	DOF_id = -1;
	Cycle_DOF_id = -1;
}

void Edge::print () {
	log() << "(" << StartVertex->id << "," << EndVertex->id;
    log() << ",[" << m_bond->Atom1->getResidue()->getId()<<m_bond->Atom1->getName() << ",";
    log() << m_bond->Atom2->getResidue()->getId() << m_bond->Atom2->getName() << "]";
    log() << ",[" << m_bond->Atom1->getId() << "," << m_bond->Atom2->getId() << ",";
    log() << (m_bond->BondType=="HB"?"HB":"CV") << "]," << (Shared?"Shared":"Indepedent");
    log() << ",ClusterId=" << RigidbodyGraphCycleCluster_id << ",DOF_id=" << DOF_id << ",Cycle_DOF_id=" << Cycle_DOF_id << ")" << endl;
}

void Edge::printVerbose () {
  log() << "\t Edge: " << endl;
  log() << "\t\t StartVertex: " << StartVertex->id << "  EndVertex: " << EndVertex->id << endl;
  log() << "\t\t StartVertex: " << (StartVertex->Rb_ptr->isMainchainRb()?"isOnMainchain":"isNOTonMainchain") << endl;
  //log() << "\t\t Direction: Res1_ID Atom1_Name Atom1_ID ---------> Res2_ID Atom2_Name Atom2_ID " << endl;
  log() << "\t\t Direction: " << m_bond->Atom1->getResidue()->getId() << " " << m_bond->Atom1->getId() << " " << m_bond->Atom1->getName() <<
  "\t ---------> \t"
  << m_bond->Atom2->getResidue()->getId() << " " << m_bond->Atom2->getId() << " " << m_bond->Atom2->getName() << endl;
  log() << "\t\t Bond type: " << (m_bond->BondType=="HB"?"HB":"CV") << ", " << (Shared?"Shared":"Indepedent") << endl;
  log() << "\t\t DOF_id= " << DOF_id << "  Cycle_DOF_id= " << Cycle_DOF_id << endl;
  log() << "---------------------------" << endl;
}

void Edge::printShort () {
  log() << "\t Edge: " << endl;
  log() << "\t\t StartVertex: " << StartVertex->id << "  EndVertex: " << EndVertex->id << endl;
  log() << "\t\t Direction: " << m_bond->Atom1->getResidue()->getId() << " " << m_bond->Atom1->getId() << " " << m_bond->Atom1->getName() <<
  "\t ---------> \t"
  << m_bond->Atom2->getResidue()->getId() << " " << m_bond->Atom2->getId() << " " << m_bond->Atom2->getName() << endl;
  log() << "---------------------------" << endl;
}

void Edge::printHTML () {
  int label = EndVertex->Rb_ptr->size();
  //int label = EndVertex->id;
  log() << "var vertex_" << EndVertex->id << " = graph.newNode({label: \"" << label << "\", fill: \"#000000\"});" << endl;
  log() << "graph.newEdge(vertex_" << StartVertex->id << ", vertex_" << EndVertex->id << ");" << endl;
}

void Edge::printHTMLRoot () {
	log() << "var vertex_" << StartVertex->id << " = graph.newNode({label: \"*"
	<< StartVertex->id << "*\", fill: \"#EE2222\"});" << endl;
}

Bond *Edge::getBond() const {
  return m_bond;
}



//bool Edge::containsAtomAtPosition ( const clipper::Coord_orth& pos ) const {
//	if( ( pos[0] >= m_bond->Atom1->m_Position[0] - m_bond->Atom1->Vdw_radius / 100 &&
//	      pos[0] <= m_bond->Atom1->m_Position[0] + m_bond->Atom1->Vdw_radius / 100 &&
//	      pos[1] >= m_bond->Atom1->m_Position[1] - m_bond->Atom1->Vdw_radius / 100 &&
//	      pos[1] <= m_bond->Atom1->m_Position[1] + m_bond->Atom1->Vdw_radius / 100 &&
//	      pos[2] >= m_bond->Atom1->m_Position[2] - m_bond->Atom1->Vdw_radius / 100 &&
//	      pos[2] <= m_bond->Atom1->m_Position[2] + m_bond->Atom1->Vdw_radius / 100 ) ||
//	    ( pos[0] >= m_bond->Atom2->m_Position[0] - m_bond->Atom2->Vdw_radius / 100 &&
//	      pos[0] <= m_bond->Atom2->m_Position[0] + m_bond->Atom2->Vdw_radius / 100 &&
//	      pos[1] >= m_bond->Atom2->m_Position[1] - m_bond->Atom2->Vdw_radius / 100 &&
//	      pos[1] <= m_bond->Atom2->m_Position[1] + m_bond->Atom2->Vdw_radius / 100 &&
//	      pos[2] >= m_bond->Atom2->m_Position[2] - m_bond->Atom2->Vdw_radius / 100 &&
//	      pos[2] <= m_bond->Atom2->m_Position[2] + m_bond->Atom2->Vdw_radius / 100 ) )
//		return true;
//	return false;
//}

//bool Edge::isDOF () {
//	return !Hbond;
//}

RigidbodyGraphVertex::RigidbodyGraphVertex ():
    id(0),
    Rb_ptr(NULL)
{
	Parent = NULL;
	isRibose = false;
}

RigidbodyGraphVertex::RigidbodyGraphVertex (int id_, Rigidbody* rb_ptr):
    id(id_),
    Rb_ptr(rb_ptr)
{
	Parent = NULL;
	//IsSpanned = false;
	Visited = false;
	isRibose = false;
	rb_ptr->setVertex(this);
}

RigidbodyGraphVertex::~RigidbodyGraphVertex () {
	//map<unsigned int,Edge*>::iterator eit;
	for (auto eit=edges.begin(); eit!=edges.end(); ++eit) {
		delete *eit;
	}
	Rb_ptr->setVertex(NULL);
}

void RigidbodyGraphVertex::setParent(RigidbodyGraphVertex* v) {
	Parent = v;
}

void RigidbodyGraphVertex::addEdge (unsigned int neighbor_vertex_id, Edge *edge) {
	//Edges.insert( make_pair(neighbor_vertex_id,edge) );
	edges.push_back( edge );
}

void RigidbodyGraphVertex::print () {
	log() << "Rigidbody_" << id << ", id ";
        for (vector<Atom*>::iterator it=Rb_ptr->Atoms.begin(); it!=Rb_ptr->Atoms.end(); ++it)
                log() << (*it)->getId() << "+";
        log() << endl;
}

void RigidbodyGraphVertex::TransformAtomPosition(RigidTransform *trsfm){
	Atom *pAtom; 
	Vector3 newPos; 
  for (vector<Atom*>::iterator it=Rb_ptr->Atoms.begin(); it!=Rb_ptr->Atoms.end(); ++it){
		pAtom = *it;


		// since an atom can belong more than one rigid body group,
		// only one rigid body transformation can be applied
		newPos = trsfm->R * pAtom->m_Position;

		newPos.x += trsfm->t.x;
		newPos.y += trsfm->t.y;
		newPos.z += trsfm->t.z;

		pAtom->m_Position.x = newPos.x;
		pAtom->m_Position.y = newPos.y;
		pAtom->m_Position.z = newPos.z;
	}
}

RigidbodyGraphCycle::RigidbodyGraphCycle () {
	IndependentLength = 0;

}



void RigidbodyGraphCycleCluster::print () {
	log() << "Number of cycles: " << Cycles.size() << endl;
	for (vector< RigidbodyGraphCycle >::iterator cycle_it=Cycles.begin(); cycle_it!=Cycles.end(); ++cycle_it) {
		log() << cycle_it->CycleEdges.size() << " rigidbody with " << cycle_it->IndependentLength << " independent edges: ";
		for ( vector<Edge*>::iterator it=cycle_it->CycleEdges.begin(); it!=cycle_it->CycleEdges.end(); ++it) {
			log() << (*it)->StartVertex->id << " ";
			//if ( (*it)->Hbond )
			//  log() << "HB" << " ";
		}
		log() << endl;
	}
}

RigidbodyGraph::RigidbodyGraph () {
	MaxCycleClusterId = 0;
}

RigidbodyGraph::~RigidbodyGraph () {
	map<unsigned int, RigidbodyGraphVertex*>::iterator it;
	m_sortedVertices.clear();
	for (it=Vertex_map.begin(); it!=Vertex_map.end(); ++it) {
		delete it->second;
	}

}

RigidbodyGraphVertex* RigidbodyGraph::addVertex (unsigned int rb_id, Rigidbody* rb, bool flexibleRibose){
    //log("debugRas")<<"addVertex("<<rb_id<<" , "<<rb<<" )"<<endl;
	RigidbodyGraphVertex* new_vertex;
	//if( (rb->getAtom("O4'") && rb->getAtom("C3'")) || (rb->getAtom("CB") && rb->getAtom("CB")->getResidue()->getName()=="PRO") ){
  if( flexibleRibose && (rb->getAtom("O4'") && rb->getAtom("C3'")) ){
		new_vertex = new SugarVertex(rb_id,rb);
	}else{
		new_vertex = new RigidbodyGraphVertex(rb_id,rb);
	}
	Vertex_map.insert(make_pair(rb_id,new_vertex));
	return new_vertex;
}

Edge* RigidbodyGraphVertex::findEdge(RigidbodyGraphVertex* v) const
{
  for(auto const& edge: edges){
    if( edge->EndVertex==v )
      return edge;
  }
  return NULL;
}

// Add a directed edge from rb_id1 to rb_id2
Edge* RigidbodyGraph::addEdgeDirected (RigidbodyGraphVertex *vertex1, RigidbodyGraphVertex *vertex2, Bond * bond, int DOF_id)
{
    //log("debugRas")<<"RigidbodyGraph::addEdgeDirected("<<vertex1->Rb_ptr<<", "<<vertex2->Rb_ptr<<", "<<bond<<"..)"<<endl;
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

bool RigidbodyGraph::hasVertex (int rb_id) {
	return (Vertex_map.find(rb_id)!=Vertex_map.end());
}

void RigidbodyGraph::print () {
	map<unsigned int, RigidbodyGraphVertex*>::iterator it;
	for (it=Vertex_map.begin(); it!=Vertex_map.end(); ++it) {
		RigidbodyGraphVertex *vertex = it->second;
		log() << "Rigidbody " << vertex->id << ": " << vertex->Rb_ptr->Atoms.size() << " atoms and " << vertex->edges.size() << " edges" << endl;
		vertex->print();
		//for (map<unsigned int,Edge*>::iterator eit=vertex->Edges.begin(); eit!=vertex->Edges.end(); ++eit) {
		for (auto eit=vertex->edges.begin(); eit!=vertex->edges.end(); ++eit) {
			(*eit)->print();
		}
		log() << endl;
	}

/*	log() << "There are " << CycleClusters.size() << " cycle clusters." << endl;
//	log() << "There are " << (Edges.size()/2-Vertex_map.size()+1) << " cycle clusters." << endl; // Euler rule to determine the number of cycles
	for (map<unsigned int,RigidbodyGraphCycleCluster>::iterator cluster_it=CycleClusters.begin(); cluster_it!=CycleClusters.end(); ++cluster_it) {
		cluster_it->second.print();
	}
    log() << endl;
*/
    log() << "Total number of edges = " << Edges.size()/2 << endl;
    log() << "Total number of rigid bodies = " << Vertex_map.size() << endl;
}

bool compareCycles (RigidbodyGraphCycle c1, RigidbodyGraphCycle c2) {
	if ( c1.IndependentLength < c2.IndependentLength )
		return true;
	else if ( c1.IndependentLength == c2.IndependentLength )
		if ( c1.CycleEdges.size() < c2.CycleEdges.size() )
			return true;
	return false;
}

void RigidbodyTree::printForSpringy () {
	queue<RigidbodyGraphVertex*> node_queue;
	node_queue.push(root);
    //log()<<"var nd_"<<root->id<<" = graph.newNode({label: \"*"<<root->Rb_ptr<<"*\"});"<<endl;
    log() << "var nd_" << root->id << " = graph.newNode({label: \"*" << root->id << "*\"});" << endl;
	while ( node_queue.size()>0 ) {
		// get the first element in the queue
		RigidbodyGraphVertex *cur_node = node_queue.front();
		// for each edge, print it and insert the child into the queue
		//for (map<unsigned int,Edge*>::iterator eit=cur_node->Edges.begin(); eit!=cur_node->Edges.end(); ++eit) {
    for (auto eit=cur_node->edges.begin(); eit!=cur_node->edges.end(); ++eit) {
			RigidbodyGraphVertex* end_node = (*eit)->EndVertex;
      //RigidbodyGraphVertex* end_node = eit->second->EndVertex;
			node_queue.push(end_node);
			log() << "var nd_" << end_node->id << " = graph.newNode({label: \"" << end_node->Rb_ptr << "\"});" << endl;
			log() << "graph.newEdge(nd_" << cur_node->id << ", nd_" << end_node->id << "); //DOF:" << (*eit)->DOF_id << endl;
		}
		node_queue.pop();
	}

    log() << "Go to http://getspringy.com/, download the demo and paste the above into the javascript" << endl;
}

void RigidbodyTree::print() {
	//return;
	// breadth-first-traverse
    log() << "Breadth-first-traversal of the tree:" << endl;
	queue<RigidbodyGraphVertex*> node_queue;
	node_queue.push(root);
	while ( node_queue.size()>0 ) {
		// get the first element in the queue
		RigidbodyGraphVertex *cur_node = node_queue.front();
		// for each edge, print it and insert the child into the queue
    //for (map<unsigned int,Edge*>::iterator eit=cur_node->Edges.begin(); eit!=cur_node->Edges.end(); ++eit) {
		for (auto eit=cur_node->edges.begin(); eit!=cur_node->edges.end(); ++eit) {
			//eit->second->print();
      //node_queue.push(eit->second->EndVertex);
      (*eit)->print();
      node_queue.push((*eit)->EndVertex);
		}
		node_queue.pop();
	}

	// print the edges closing cycles and common ancestors for the anchors in each edge
    log() << "Edges closing cycles:" << endl;
	for (vector< pair<Edge*,RigidbodyGraphVertex*> >::iterator pit=CycleAnchorEdges.begin(); pit!=CycleAnchorEdges.end(); ++pit) {
		pit->first->print();
        log() << "Common ancestor: ";
		pit->second->print();

		// trace the rigid bodies in the cycle
		RigidbodyGraphVertex *start, *end, *cur;
		start = pit->first->StartVertex;
		end = pit->first->EndVertex;
		int length = 1;
        log() << "Left cycle:";
		for (cur=start; cur!=pit->second; cur=cur->Parent) {
            log() << " " << cur->id;
			++length;
		}
        log() << endl;
        log() << "Right cycle:";
		for (cur=end; cur!=pit->second; cur=cur->Parent) {
            log() << " " << cur->id;
			++length;
		}
        log() << endl;
        log() << "Cycle length = " << length << endl;
	}

}

RigidbodyGraphVertex* RigidbodyTree::findCommonAncestor (RigidbodyGraphVertex *v1, RigidbodyGraphVertex *v2) {
	// traverse from v1 to root, and mark every vertex along the way to be Visited
	RigidbodyGraphVertex *cur_node = v1;
	do {
		cur_node->Visited = true;
        //log("debug")<<"Cur node [1] : "<<cur_node->Rb_ptr<<endl;
		if (cur_node == root)
			break;
		else{
			if(cur_node->Parent==NULL){
				cerr<<"RigidbodyTree::findCommonAncestor("<<v1->Rb_ptr<<","<<v2->Rb_ptr<<") node has no m_parent: "<<cur_node->Rb_ptr<<endl;
				cerr<<"You might see this error because of multiple occupancy atoms in the structure"<<endl;
				exit(-1);
			}
			cur_node = cur_node->Parent;
		}
	} while (true);
	// traverse from v2 to root, and stop until meeting a Visited vertex
	cur_node = v2;
    //log("debug")<<"Cur node [2] : "<<cur_node->Rb_ptr<<endl;
	while ( !cur_node->Visited ) {
        //log("debug")<<"Cur node [2] : "<<cur_node->Rb_ptr<<endl;
		if(cur_node->Parent==NULL){
			cerr<<"RigidbodyTree::findCommonAncestor("<<v1->Rb_ptr<<","<<v2->Rb_ptr<<") node has no m_parent: "<<cur_node->Rb_ptr<<endl;
			cerr<<"You might see this error because of multiple occupancy atoms in the structure"<<endl;
			exit(-1);
		}
		cur_node = cur_node->Parent;
	}
	RigidbodyGraphVertex *ancestor = cur_node;
    //log("debug")<<"done [3] "<<endl;
	// unmark all the Visited nodes
	cur_node = v1;
	do {
		cur_node->Visited = false;
		if (cur_node == root)
			break;
		else
			cur_node = cur_node->Parent;
	} while (true);
	return ancestor;
}

RigidbodyGraphVertex* RigidbodyGraph::getVertex (int rb_id) {
	return Vertex_map.find(rb_id)->second;
}

RigidbodyTree::RigidbodyTree(): RigidbodyGraph(){
	num_DOFs = 0;
	Cycle_DOF_num = 0;
}
RigidbodyTree::~RigidbodyTree () {
	for (vector< pair<Edge*,RigidbodyGraphVertex*> >::iterator it=CycleAnchorEdges.begin(); it!=CycleAnchorEdges.end(); ++it) {
		delete it->first;
	}
}


ostream& operator<<(ostream& os, const Edge& e){
	os<<"Edge["<<e.getBond()->Atom1->getName()<<", "<<e.getBond()->Atom2->getName()<<"]";
	return os;
}

ostream& operator<<(ostream& os, const Edge* e){
	os<<*e;
	return os;
}
