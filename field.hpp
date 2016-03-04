#ifndef __field_hpp
#define __field_hpp
#include <vector>
#include <math.h>
#include <string>
#include <iostream>

#include "node.hpp"
#include "shell.hpp"
#include "shell_network.hpp"

using namespace std;

class vecfield{
  vector <vector<complex<double> > > val;
  shell_network* network;
private:
  vecfield();
  ~vecfield();
  void initialize_from_file(string filename):
  void assign_shell_network(shell_network snw){network=snw;};
public:
};



#endif
