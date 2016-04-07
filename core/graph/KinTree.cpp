//
// Created by Rasmus Fonseca on 4/7/16.
//

#include "KinTree.h"

#include <queue>

#include "Logger.h"

using namespace std;

void KinTree::printForSpringy () {
  queue<KinVertex*> node_queue;
  node_queue.push(root);
  //log()<<"var nd_"<<root->id<<" = graph.newNode({label: \"*"<<root->m_rigidbody<<"*\"});"<<endl;
  log() << "var nd_" << root->id << " = graph.newNode({label: \"*" << root->id << "*\"});" << endl;
  while ( node_queue.size()>0 ) {
    // get the first element in the queue
    KinVertex *cur_node = node_queue.front();
    // for each edge, print it and insert the child into the queue
    //for (map<unsigned int,KinEdge*>::iterator eit=cur_node->Edges.begin(); eit!=cur_node->Edges.end(); ++eit) {
    for (auto eit=cur_node->m_edges.begin(); eit!=cur_node->m_edges.end(); ++eit) {
      KinVertex* end_node = (*eit)->EndVertex;
      //KinVertex* end_node = eit->second->EndVertex;
      node_queue.push(end_node);
      log() << "var nd_" << end_node->id << " = graph.newNode({label: \"" << end_node->m_rigidbody << "\"});" << endl;
      log() << "graph.newEdge(nd_" << cur_node->id << ", nd_" << end_node->id << "); //DOF:" << (*eit)->DOF_id << endl;
    }
    node_queue.pop();
  }

  log() << "Go to http://getspringy.com/, download the demo and paste the above into the javascript" << endl;
}

void KinTree::print() {
  //return;
  // breadth-first-traverse
  log() << "Breadth-first-traversal of the tree:" << endl;
  queue<KinVertex*> node_queue;
  node_queue.push(root);
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
  for (vector< pair<KinEdge*,KinVertex*> >::iterator pit=CycleAnchorEdges.begin(); pit!=CycleAnchorEdges.end(); ++pit) {
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
      log() << " " << cur->id;
      ++length;
    }
    log() << endl;
    log() << "Right cycle:";
    for (cur=end; cur!=pit->second; cur=cur->m_parent) {
      log() << " " << cur->id;
      ++length;
    }
    log() << endl;
    log() << "Cycle length = " << length << endl;
  }

}

KinVertex* KinTree::findCommonAncestor (KinVertex *v1, KinVertex *v2) {
  // traverse from v1 to root, and mark every vertex along the way to be Visited
  KinVertex *cur_node = v1;
  do {
    cur_node->Visited = true;
    //log("debug")<<"Cur node [1] : "<<cur_node->m_rigidbody<<endl;
    if (cur_node == root)
      break;
    else{
      if(cur_node->m_parent==NULL){
        cerr<<"KinTree::findCommonAncestor("<<v1->m_rigidbody<<","<<v2->m_rigidbody<<") node has no m_parent: "<<cur_node->m_rigidbody<<endl;
        cerr<<"You might see this error because of multiple occupancy atoms in the structure"<<endl;
        exit(-1);
      }
      cur_node = cur_node->m_parent;
    }
  } while (true);
  // traverse from v2 to root, and stop until meeting a Visited vertex
  cur_node = v2;
  //log("debug")<<"Cur node [2] : "<<cur_node->m_rigidbody<<endl;
  while ( !cur_node->Visited ) {
    //log("debug")<<"Cur node [2] : "<<cur_node->m_rigidbody<<endl;
    if(cur_node->m_parent==NULL){
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
    if (cur_node == root)
      break;
    else
      cur_node = cur_node->m_parent;
  } while (true);
  return ancestor;
}

void KinTree::collectDOFs()
{
  collectDOFs(root);
}

void KinTree::collectDOFs(KinVertex* v)
{
  for(auto const& edge: v->m_edges){
    if( std::find(m_dofs.begin(),m_dofs.end(), edge->getDOF())==m_dofs.end() )
      m_dofs.push_back(edge->getDOF());

    collectDOFs(edge->EndVertex);
  }
}

KinTree::KinTree(): KinGraph(){
  m_numDOFs = 0;
  m_numCycleDOFs = 0;
}
KinTree::~KinTree () {
  for (vector< pair<KinEdge*,KinVertex*> >::iterator it=CycleAnchorEdges.begin(); it!=CycleAnchorEdges.end(); ++it) {
    delete it->first;
  }
}
