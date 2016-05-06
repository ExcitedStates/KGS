//
// Created by Rasmus Fonseca on 4/7/16.
//

#include "KinTree.h"

#include <queue>
#include <cassert>

#include "Logger.h"

using namespace std;
size_t KinTree::getNumDOFs() const      { return m_dofs.size(); }
size_t KinTree::getNumCycleDOFs() const { return m_cycleDOFs.size(); }

void KinTree::print() const {
  //return;
  // breadth-first-traverse
  log() << "Breadth-first-traversal of the tree:" << endl;
  queue<KinVertex*> node_queue;
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
  for (auto pit=CycleAnchorEdges.begin(); pit!=CycleAnchorEdges.end(); ++pit) {
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

KinVertex* KinTree::findCommonAncestor (KinVertex *v1, KinVertex *v2) {
  // traverse from v1 to m_root, and mark every vertex along the way to be Visited
  KinVertex *cur_node = v1;
  do {
    cur_node->Visited = true;
    //log("debug")<<"Cur node [1] : "<<cur_node->m_rigidbody<<endl;
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
  //log("debug")<<"Cur node [2] : "<<cur_node->m_rigidbody<<endl;
  while ( !cur_node->Visited ) {
    //log("debug")<<"Cur node [2] : "<<cur_node->m_rigidbody<<endl;
    if(cur_node->m_parent==nullptr){
      cerr<<"KinTree::findCommonAncestor("<<v1->m_rigidbody<<","<<v2->m_rigidbody<<") node has no m_parent: "<<cur_node->m_rigidbody<<endl;
      cerr<<"You might see this error because of multiple occupancy atoms in the structure"<<endl;
      exit(-1);
    }
    cur_node = cur_node->m_parent;
  }
  KinVertex *ancestor = cur_node;
  //log("debug")<<"done [3] "<<endl;
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
  for(auto const& edge: v->m_edges){
    if( std::find(m_dofs.begin(),m_dofs.end(), edge->getDOF())==m_dofs.end() ) {
      DOF* dof = edge->getDOF();
      dof->setIndex(m_dofs.size());
      m_dofs.push_back(dof);
    }

    collectDOFs(edge->EndVertex);
  }
}

KinTree::KinTree():
    KinGraph()
{
}

KinTree::~KinTree () {
  for (vector< pair<KinEdge*,KinVertex*> >::iterator it=CycleAnchorEdges.begin(); it!=CycleAnchorEdges.end(); ++it) {
    delete it->first;
  }
}


DOF* KinTree::getDOF(unsigned int idx) const{
  assert(idx<m_dofs.size());
  return m_dofs[idx];
}

DOF* KinTree::getCycleDOF(unsigned int idx) const
{
  assert(idx<m_cycleDOFs.size());
  return m_cycleDOFs[idx];
}

void KinTree::addCycleDOF(DOF* dof)
{
  if(std::find(m_cycleDOFs.begin(), m_cycleDOFs.end(), dof)==m_cycleDOFs.end()) {
    //DOF is not already in m_cycleDOFs
    dof->setCycleIndex(m_cycleDOFs.size());
    m_cycleDOFs.push_back(dof);
  }
}


