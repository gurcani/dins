#ifndef __gviz_hpp
#define __gviz_hpp
#include <vector>
#include <string>
#include <sstream>
#include "node.hpp"
#include "shell.hpp"

using namespace std;

class gviz{
private:
  string filename;
  ostringstream buffer;
public:
  gviz(string file);
  ~gviz(){};
  void print(){ cout<<buffer.str()<<"}\n";};
  void shelltosubgraph(shell *sh);
  void add_connection(node *nd, node_connection *cn,int id);
};
#endif
