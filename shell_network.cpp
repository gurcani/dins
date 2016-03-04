#include <vector>
#include <math.h>
#include <iostream>
#include "node.hpp"
#include "shell.hpp"
#include "shell_network.hpp"
#include <boost/numeric/odeint.hpp>
using namespace std;
using boost::numeric::ublas::matrix;

shell_network::shell_network(){
}

shell_network::~shell_network(){
}

void shell_network::setup_nodes(int N, double g, double k0, bool icosa_first){
  shell* tsh;
  double kn;
  for(int l=0;l<N;l++){
    kn=k0*pow(g,l);
    if( (((bool)(l%2)) ^ (icosa_first)) )
      tsh=new icosa_shell(l,kn*sqrt(sqrt(5)/3));
    else
      tsh=new dodeca_shell(l,kn);
    add_shell(tsh);
  }
  kni=matrix<double>(nodes.size(),3,0);
  vector<double> k;
  for(int l=0;l<kni.size1();l++){
    k=nodes[l]->kni();
    for(int m=0;m<kni.size2();m++){
      kni(l,m)=k[m];
    }
  }
}

void shell_network::setup_connections(){
  int N=shells.size();
  for(int l=0; l<N;l++){
    if(l==0){
      shells[l]->connect_shells(shells[l+1], shells[l+2],0);
    }
    else if(l==1){
      shells[l]->connect_shells(shells[l-1], shells[l+1],1);
      shells[l]->connect_shells(shells[l+1], shells[l+2],0);
    }
    else if(l==N-1){
      shells[l]->connect_shells(shells[l-2], shells[l-1],2);
    }
    else if(l==N-2){
      shells[l]->connect_shells(shells[l-2], shells[l-1],2);
      shells[l]->connect_shells(shells[l-1], shells[l+1],1);
    }
    else{
      shells[l]->connect_shells(shells[l-2], shells[l-1],2);
      shells[l]->connect_shells(shells[l-1], shells[l+1],1);
      shells[l]->connect_shells(shells[l+1], shells[l+2],0);
    }
  }
}

void shell_network::add_shell(shell *sh){
  int id;
  if(shells.empty())
    id=-1;
  else
    id=shells.back()->get_id();
  sh->set_id(id+1);
  shells.push_back(sh);
  vector<node*> nds=sh->get_nodes();
  for (vector<node*>::const_iterator nd = nds.begin(); nd != nds.end(); ++nd){
    if(nodes.empty())
      id=-1;
    else
      id=nodes.back()->gid();
    (*nd)->set_global_id(id+1);
    nodes.push_back(*nd);
    //    kni.push_back(nodes.kni());
  }
}
