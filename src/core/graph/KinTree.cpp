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



#include "KinTree.h"

#include <queue>
#include <cassert>
#include <core/dofs/GlobalTranslateDOF.h>
#include <core/dofs/GlobalRotateDOF.h>
#include <core/dofs/TorsionDOF.h>

#include "Logger.h"

using std::endl;
using std::cout;
using std::cerr;


KinTree::KinTree( const std::vector<Rigidbody*>& rigidbodies, const std::vector<Atom*>& roots ):
    KinGraph()
{
  // Add the super-root
  m_root = addVertex(nullptr);

  // add all rigid bodies as vertices into Rigidbody_graph
  for( auto const &rb: rigidbodies ) {
    addVertex(rb);
  }

  std::list<KinEdge *> cycleEdges;

  //Initialize chain-roots by adding them to the queue and setting up edges from the super-root
  auto my_comp = [](const KinVertex* v1, const KinVertex* v2){
    int id1 = v1->m_rigidbody == nullptr ? 0: v1->m_rigidbody->id();
    int id2 = v2->m_rigidbody == nullptr ? 0: v2->m_rigidbody->id();
    return id1>id2;
  };
  std::priority_queue<KinVertex*, std::vector<KinVertex*>,decltype(my_comp)> queue(my_comp);


  //Build new vector with chain roots given by the user at the front
  //Iterate through this vector instead of vertex_map to use correct chain roots
  std::map<unsigned int,KinVertex*> secondPriorityMap;
  std::vector<KinVertex*> vertexPriorityList;

  //Push the chain roots into Vertex Map, others to second prio map
  for( auto const &id_vertex_pair: Vertex_map ) {
    KinVertex *vertex=id_vertex_pair.second;

    bool foundChainRoot = false;
    for(auto const &rootAtom: roots) {
      if (vertex->m_rigidbody->containsAtom(rootAtom)) {
        foundChainRoot=true;
        break;
      }
    }
    if(foundChainRoot) { vertexPriorityList.push_back(vertex); }
    else { secondPriorityMap[id_vertex_pair.first] = vertex; }
  }

  // add all rigid bodies as vertices into Rigidbody_graph
  for( auto const &id_vertex_pair: secondPriorityMap ) {
    vertexPriorityList.push_back(id_vertex_pair.second);
  }

  //Perform breadth-first-search from all queue vertices and construct KinEdges
  //A variant of Prims algorithm is used to get the orientation of the tree correct in the first go
  std::set<Bond *> visitedBonds;
  std::set<KinVertex *> visitedVertices;

  for(auto const& chainRoot: vertexPriorityList) {
    //If the vertex has been visited before, ignore it, otherwise make it the root of a new tree
    if( visitedVertices.count(chainRoot)>0) continue;

    //Connect to super-root with six global chain dofs
    KinVertex *v2 = addVertex(nullptr);
    KinVertex *v3 = addVertex(nullptr);
    KinVertex *v4 = addVertex(nullptr);
    KinVertex *v5 = addVertex(nullptr);
    KinVertex *v6 = addVertex(nullptr);

    KinEdge *e1 = addEdgeDirected(m_root, v2, nullptr);
    KinEdge *e2 = addEdgeDirected(v2, v3, nullptr);
    KinEdge *e3 = addEdgeDirected(v3, v4, nullptr);
    KinEdge *e4 = addEdgeDirected(v4, v5, nullptr);
    KinEdge *e5 = addEdgeDirected(v5, v6, nullptr);
    KinEdge *e6 = addEdgeDirected(v6, chainRoot, nullptr);

    e1->setDOF(new GlobalTranslateDOF(e1, 0));
    e2->setDOF(new GlobalTranslateDOF(e2, 1));
    e3->setDOF(new GlobalTranslateDOF(e3, 2));
    e4->setDOF(new GlobalRotateDOF(e4, 0));
    e5->setDOF(new GlobalRotateDOF(e5, 1));
    e6->setDOF(new GlobalRotateDOF(e6, 2));
    log("debug") << "KinTree::KinTree(..) - Connecting " << chainRoot->m_rigidbody->Atoms[0]
                 << " to super-root using 6 dofs" << endl;

    queue.push(chainRoot);

    while (!queue.empty()) {
      KinVertex *current_vertex = queue.top();
      queue.pop();
      visitedVertices.insert(current_vertex);
      log("debug") << "KinTree::KinTree(..) - Visiting vertex of size " <<
                   current_vertex->m_rigidbody->Atoms.size() << ", rbID: " << current_vertex->m_rigidbody->id() << ", "
                   <<
                   current_vertex->m_rigidbody->m_bonds.size() << " bonds" << endl;

      for (Bond *const &bond: current_vertex->m_rigidbody->m_bonds) {
        // Determine which other rigid body bond is connected to
        KinVertex *bonded_vertex = bond->m_atom2->getRigidbody()->getVertex();
        if (bonded_vertex == current_vertex)
          bonded_vertex = bond->m_atom1->getRigidbody()->getVertex();

        if (current_vertex == bonded_vertex) {
          log("debug") << "KinTree::KinTree(..) - Bond connecting same rigid body " << bond << endl;
          continue;
        }
        if (visitedBonds.count(bond) > 0) {
          log("debug") << "KinTree::KinTree(..) - Already visited bond " << bond << endl;
          continue;
        }

        visitedBonds.insert(bond);

        if (bond->isHBond()) {
          // If it's an H-bond, it closes a cycle. Add it in m_cycleAnchorEdges.
          KinEdge *edge = new KinEdge(current_vertex, bonded_vertex, bond);
          cycleEdges.push_back(edge);
          log("debug") << "KinTree::KinTree(..) - Adding cycle-edge from h-bond " << edge << endl;
//          cout<< "Molecule::buildSpanningTree() - Adding cycle-edge from h-bond " << edge->getBond()->Atom1->getId() <<", "<<edge->getBond()->m_atom2->getId() << endl;
        } else if (bond->isDBond()) {
          KinEdge *edge = new KinEdge(current_vertex, bonded_vertex, bond);
          cycleEdges.push_back(edge);
          log("debug") << "KinTree::KinTree(..) - Adding cycle-edge from d-bond " << edge << endl;
        } else {
          if (visitedVertices.count(bonded_vertex) > 0) {
            KinEdge *edge = new KinEdge(current_vertex, bonded_vertex, bond);
            cycleEdges.push_back(edge);
            log("debug") << "KinTree::KinTree(..) - Adding cycle-edge from covalent bond " << edge << endl;

          } else {
            // If it's a covalent bond, add it into the tree m_edges
            visitedVertices.insert(bonded_vertex);
            queue.push(bonded_vertex);
            KinEdge *edge = addEdgeDirected(current_vertex, bonded_vertex, bond);
            edge->setDOF(new TorsionDOF(edge));
            log("debug") << "KinTree::KinTree(..) - Adding torsion-edge " << edge << endl;
          }
        }
      }
    } // end while
  }

  collectDOFs();

  //Sort cycle anchor edges to maintain constant row order for different roots
  //This is important for the hydrogen-bond hierarchy analysis!
  cycleEdges.sort(KinEdge::compareIDs);

  // For each hbond KinEdge, find the lowest common ancestor (LCA) of its end-vertices and put all DOFs from the
  // end-points to the LCA into m_spanningTree->m_cycleDOFs.
  for (auto const &h_edge : cycleEdges) {
    KinVertex *lca = findCommonAncestor(h_edge->StartVertex, h_edge->EndVertex);
    m_cycleAnchorEdges.push_back(std::make_pair(h_edge, lca));

    for( KinVertex *v = h_edge->StartVertex; v != lca; v = v->m_parent ) {
      KinEdge *edge = v->m_parent->findEdge(v);
      addCycleDOF(edge->getDOF());
    }

    for( KinVertex *v = h_edge->EndVertex; v != lca; v = v->m_parent ) {
      KinEdge *edge = v->m_parent->findEdge(v);
      addCycleDOF(edge->getDOF());
    }
  }
}

KinTree::~KinTree () {
//  for (std::vector< std::pair<KinEdge*,KinVertex*> >::iterator it=m_cycleAnchorEdges.begin(); it!=m_cycleAnchorEdges.end(); ++it) {
//  delete it->first;
//}
  for (auto edge_nca_pair: m_cycleAnchorEdges) {
    delete edge_nca_pair.first;
  }
}

size_t KinTree::getNumDOFs() const      { return m_dofs.size(); }

size_t KinTree::getNumCycleDOFs() const { return m_cycleDOFs.size(); }

void KinTree::print() const {
  //return;
  // breadth-first-traverse
  log() << "Breadth-first-traversal of the tree:" << endl;
  std::queue<KinVertex*> node_queue;
  node_queue.push(m_root);
  while ( node_queue.size()>0 ) {
    // get the first element in the queue
    KinVertex *cur_node = node_queue.front();
    // for each edge, print it and insert the child into the queue
    //for (map<unsigned int,KinEdge*>::iterator eit=cur_node->Edges.begin(); eit!=cur_node->Edges.end(); ++eit) {
    for (auto eit=cur_node->m_edges.begin(); eit!=cur_node->m_edges.end(); ++eit) {
      //eit->second->print();
      //node_queue.push(eit->second->EndVertex);
      (*eit)->print();
      node_queue.push((*eit)->EndVertex);
    }
    node_queue.pop();
  }

  // print the m_edges closing cycles and common ancestors for the anchors in each edge
  log() << "Edges closing cycles:" << endl;
  for (auto pit=m_cycleAnchorEdges.begin(); pit!=m_cycleAnchorEdges.end(); ++pit) {
    pit->first->print();
    log() << "Common ancestor: ";
    pit->second->print();

    // trace the rigid bodies in the cycle
    KinVertex *start, *end, *cur;
    start = pit->first->StartVertex;
    end = pit->first->EndVertex;
    int length = 1;
    log() << "Left cycle:";
    for (cur=start; cur!=pit->second; cur=cur->m_parent) {
      if(cur->m_rigidbody==nullptr)
        log() << " global";
      else
        log() << " " << cur->m_rigidbody->id();
      ++length;
    }
    log() << endl;
    log() << "Right cycle:";
    for (cur=end; cur!=pit->second; cur=cur->m_parent) {
      if(cur->m_rigidbody==nullptr)
        log() << " global";
      else
        log() << " " << cur->m_rigidbody->id();
      ++length;
    }
    log() << endl;
    log() << "Cycle length = " << length << endl;
  }

}

KinEdge* KinTree::addEdgeDirected(KinVertex *vertex1, KinVertex *vertex2, Bond * bond)
{
  if(bond!=nullptr) {
    //log("debugRas")<<"KinGraph::addEdgeDirected("<<vertex1->m_rigidbody<<", "<<vertex2->m_rigidbody<<", "<<bond<<"..)"<<endl;
    Atom *atom2, *atom3, *atom4;
    Bond *bond_copy = new Bond(*bond);
    atom2 = bond_copy->m_atom1;
    atom3 = bond_copy->m_atom2;
    atom4 = nullptr;

    // Find out the atom that covalently bonded to atom3 with smallest m_id. It participates in the definition of the torsional angle.
    for (std::vector<Atom *>::iterator aitr = atom3->Cov_neighbor_list.begin();
         aitr != atom3->Cov_neighbor_list.end(); ++aitr) {
      if ((*aitr) == atom2) continue;
      if (atom4 == nullptr || (*aitr)->getId() < atom4->getId()) {
        atom4 = *aitr;
      }
    }
    for (std::vector<Atom *>::iterator aitr = atom3->Hbond_neighbor_list.begin();
         aitr != atom3->Hbond_neighbor_list.end(); ++aitr) {
      if ((*aitr) == atom2) continue;
      if (atom4 == nullptr || (*aitr)->getId() < atom4->getId()) {
        atom4 = *aitr;
      }
    }

    // If atom4 is in vertex1, then should flip atom2 and atom3 so that the bond is pointing from vertex1 to vertex2
    Atom *tmp_atom;
    for (std::vector<Atom *>::iterator svIT = vertex1->m_rigidbody->Atoms.begin();
         svIT != vertex1->m_rigidbody->Atoms.end(); ++svIT) {
      if ((*svIT) == atom4) {
        tmp_atom = bond_copy->m_atom1;
        bond_copy->m_atom1 = bond_copy->m_atom2;
        bond_copy->m_atom2 = tmp_atom;
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

KinVertex* KinTree::findCommonAncestor (KinVertex *v1, KinVertex *v2) {
  // traverse from v1 to m_root, and mark every vertex along the way to be Visited
  KinVertex *cur_node = v1;
  do {
    cur_node->Visited = true;
//    log("debug")<<"Cur node [1] : "<<cur_node->m_rigidbody<<endl;
    if (cur_node == m_root)
      break;
    else{
      if(cur_node->m_parent==nullptr){
        cerr<<"KinTree::findCommonAncestor("<<v1->m_rigidbody<<","<<v2->m_rigidbody<<") node has no m_parent: "<<cur_node->m_rigidbody<<endl;
        cerr<<"You might see this error because of multiple occupancy atoms in the structure"<<endl;
        exit(-1);
      }
      cur_node = cur_node->m_parent;
    }
  } while (true);
  // traverse from v2 to m_root, and stop until meeting a Visited vertex
  cur_node = v2;
//  log("debug")<<"Cur node [2] : "<<cur_node->m_rigidbody<<endl;
  while ( !cur_node->Visited ) {
//    log("debug")<<"Cur node [2] : "<<cur_node->m_rigidbody<<endl;
    if(cur_node->m_parent==nullptr){
      cerr<<"KinTree::findCommonAncestor("<<v1->m_rigidbody<<","<<v2->m_rigidbody<<") node has no m_parent: "<<cur_node->m_rigidbody<<endl;
      cerr<<"You might see this error because of multiple occupancy atoms in the structure"<<endl;
      exit(-1);
    }
    cur_node = cur_node->m_parent;
  }
  KinVertex *ancestor = cur_node;
  // log("debug")<<"done [3] "<<endl;
  // unmark all the Visited nodes
  cur_node = v1;
  do {
    cur_node->Visited = false;
    if (cur_node == m_root)
      break;
    else
      cur_node = cur_node->m_parent;
  } while (true);
  return ancestor;
}

void KinTree::collectDOFs()
{
  collectDOFs(m_root);
}

void KinTree::collectDOFs(KinVertex* v)
{
  for (auto const& edge : v->m_edges) {
    if( std::find(m_dofs.begin(), m_dofs.end(), edge->getDOF()) == m_dofs.end() ) {
      DOF* dof = edge->getDOF();
      dof->setIndex(m_dofs.size());
      m_dofs.push_back(dof);
    }

    collectDOFs(edge->EndVertex);
  }
}


DOF* KinTree::getDOF(unsigned int idx) const{
  assert(idx < m_dofs.size());
  return m_dofs[idx];
}

DOF* KinTree::getCycleDOF(unsigned int idx) const
{
  assert(idx < m_cycleDOFs.size());
  return m_cycleDOFs[idx];
}

void KinTree::addCycleDOF(DOF* dof)
{
  if (std::find(m_cycleDOFs.begin(), m_cycleDOFs.end(), dof) == m_cycleDOFs.end()) {
    // DOF is not already in m_cycleDOFs
    dof->setCycleIndex(m_cycleDOFs.size());
    m_cycleDOFs.push_back(dof);
  }
}

