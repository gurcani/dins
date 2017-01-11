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
  void remove_divergence(matrix<complex<double > > &mat){
    double kn;
    vector<double> kni;
    for(int n=0;n<nodes.size();n++){
      kn=nodes[n]->k();
      kni=nodes[n]->kni();
      complex<double> k_dot_mat=kni[0]*mat(n,0)+kni[1]*mat(n,1)+kni[2]*mat(n,2);
      for(int i=0;i<3;i++){
	mat(n,i)=mat(n,i)-k_dot_mat*kni[i]/kn/kn;
      }
    }
  }
  void force_shell(int l,complex<double>f, matrix<complex<double> > &fm){shells[l]->force(f,fm);};
};
#endif
