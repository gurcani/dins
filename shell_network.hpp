#ifndef __shell_network_h
#define __shell_network_h

#include <vector>
#include "node.hpp"
#include "shell.hpp"
#include <boost/numeric/odeint.hpp>

using namespace std;
using boost::numeric::ublas::matrix;

class shell_network{
private:
  vector<node*> nodes;
  vector<shell*> shells;
  matrix<double> kni;
public:
  shell_network();
  ~shell_network();
  void add_shell(shell*);
  vector<shell*> get_shells(){return shells;};
  vector<node*> get_nodes(){return nodes;};
  vector<shell*>* get_shellsp(){return &shells;};
  vector<node*>* get_nodesp(){return &nodes;};
  void setup_nodes(int N, double g, double k0, bool icosa_first=true);
  void setup_connections();
  matrix<double>& get_kni(){return kni;};
  void symmetrize_shells(matrix<complex<double > > &mat){for(int l=0;l<shells.size();l++) shells[l]->symmetrize(mat);}
  void force_shell(int l,complex<double>f, matrix<complex<double> > &fm){shells[l]->force(f,fm);};
};
#endif
