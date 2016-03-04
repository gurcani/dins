#include <vector>
#include <math.h>
#include <iostream>
#include "node.hpp"
#include "shell.hpp"
#include <boost/numeric/odeint.hpp>

using namespace std;
using boost::numeric::ublas::matrix;

shell::shell(int n, double k){
  kn=k;
  global_id=n;
}
shell::~shell(){ 
}

void shell::print_nodes(){
  for (vector<node*>::const_iterator nd = nodes.begin(); nd != nodes.end(); ++nd)
    //  cout<<"";
    cout<<"gid="<<(*nd)->gid()<<",sid="<<(*nd)->id()<<",theta="<<(*nd)->theta()*180/PI<<", phi="<<(*nd)->phi()*180/PI<<", k="<<(*nd)->k()<<"\n";
}

void shell::connect_shells(shell *sh1, shell *sh2, int order_sh0){
  int ***cons;
  int n1,n2;
  bool ico;
  if(order_sh0>2 || order_sh0<0){
    cerr<<"ERROR:order_sh0 can not be "<<order_sh0<<"!, assuming 1\n";
    order_sh0=1;
  }
  ico=is_icosa();
  if(!ico && !is_dodeca()){
    cerr<<"ERROR:shell is neither icosahedral nor dedocahedral, assuming icosahedral!\n";
    cerr<<"number of nodes = "<< nodes.size()<<" \n";
    ico=true;
  }
  if(ico){
    for (int l=0;l<12;l++){
      for(int m=0;m<5;m++){
	if(order_sh0==0){
	  n1=connections_icosa012[l%6][m][0];
	  n2=connections_icosa012[l%6][m][1];
	  if (l>=6){
	    n1=(n1+10)%20;
	    n2=(n2+6)%12;
	  }
	}
	else if(order_sh0==1){
	  n1=connections_icosa102[l%6][m][0];
	  n2=connections_icosa102[l%6][m][1];
	  if (l>=6){
	    n1=(n1+10)%20;
	    n2=(n2+10)%20;
	  }
	}
	else if(order_sh0==2){
	  n1=connections_icosa120[l%6][m][0];
	  n2=connections_icosa120[l%6][m][1];
	  if (l>=6){
	    n1=(n1+6)%12;
	    n2=(n2+10)%20;
	  }
	}
	nodes[l]->add_connection(sh1->nodes[n1],sh2->nodes[n2]);
      }
    }
  }
  else{
    for (int l=0;l<20;l++){
      for(int m=0;m<3;m++){
	if(order_sh0==0){
	  n1=connections_dodeca012[l%10][m][0];
	  n2=connections_dodeca012[l%10][m][1];
	  if (l>=10){
	    n1=(n1+6)%12;
	    n2=(n2+10)%20;
	  }
	}
	else if(order_sh0==1){
	  n1=connections_dodeca102[l%10][m][0];
	  n2=connections_dodeca102[l%10][m][1];
	  if (l>=10){
	    n1=(n1+6)%12;
	    n2=(n2+6)%12;
	  }
	}
	else if(order_sh0==2){
	  n1=connections_dodeca120[l%10][m][0];
	  n2=connections_dodeca120[l%10][m][1];
	  if (l>=10){
	    n1=(n1+10)%20;
	    n2=(n2+6)%12;
	  }
	}
	nodes[l]->add_connection(sh1->nodes[n1],sh2->nodes[n2]);
      }
    }
  }
}

void shell::symmetrize(matrix<complex<double> > & mat){
  int num=nodes.size();
  int i;
  for(int l=num/2;l<num;l++){
    i=nodes[l]->gid();
    for(int j=0;j<mat.size2();j++){
      mat(i,j)=conj(mat(i-num/2,j));
    }
  }
}

icosa_shell::icosa_shell(int n, double k) : shell(n,k){
  node_icosa *nd;
  vector<int> *f;
  for(int l=0;l<12;l++){
    nd=new node_icosa(l,k);
    push_node(nd);
  }
  for(int l=0;l<20;l++){
    f=new vector<int>(faces_icosa[l],faces_icosa[l]+sizeof(faces_icosa[l])/sizeof(faces_icosa[l][0]));
    push_face(*f);
  }
}

dodeca_shell::dodeca_shell(int n, double k) : shell(n,k){
  node_dodeca *nd;
  vector<int> *f;
  for(int l=0;l<20;l++){
    nd=new node_dodeca(l,k);
    push_node(nd);
  }
  for(int l=0;l<12;l++){
    f=new vector<int>(faces_dodeca[l],faces_dodeca[l]+sizeof(faces_dodeca[l])/sizeof(faces_dodeca[l][0]));
    push_face(*f);
  }
}
