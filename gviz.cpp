#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "gviz.hpp"
#include "shell.hpp"
using namespace std;

gviz::gviz(string file){
  filename=file;
  buffer<<"digraph all_nodes {\n";
  buffer<<"size=\"100,100\"\n";
  //  buffer=buffer+"a->b [label = \"1\"]\n";
  //  buffer=buffer+"}\n";
};

void gviz::add_connection(node *nd, node_connection *cn, int id){
  buffer<<nd->gid()<<" -> "<<cn->l1->gid()<<"[color=blue,label="<<id<<"]\n";
  buffer<<nd->gid()<<" -> "<<cn->l2->gid()<<"[color=red,label="<<id<<"]\n";
}

void gviz::shelltosubgraph(shell *sh){
  buffer<<"subgraph shell_"<<sh->get_id()<<" {\n";
  vector<node*> nodes=sh->get_nodes();
  for (vector<node*>::const_iterator nd = nodes.begin(); nd != nodes.end(); ++nd){
    buffer<<(*nd)->gid()<<" ";
  }
  buffer<<";}\n";
}
