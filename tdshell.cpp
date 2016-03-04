#include <vector>
#include <string>
#include <iostream>
#include <complex>
#include <cmath>
#include "tdshell.hpp"
#include "shell_network.hpp"
#include "shell.hpp"
#include "node.hpp"
#include <boost/timer/timer.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::timer;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::identity_matrix;

tdshell::tdshell(){
}

tdshell::tdshell(int N, double k0){
  double g=sqrt((1.0+sqrt(5.0))/2.0);
  snw=new shell_network();
  snw->setup_nodes(N,g,k0);
  snw->setup_connections();
  assign_network(snw);
}
tdshell::tdshell(string filename,int N, double k0){
  double g=sqrt((1.0+sqrt(5.0))/2.0);
  snw=new shell_network();
  snw->setup_nodes(N,g,k0);
  snw->setup_connections();
  assign_network(snw);
  datfile=new hdf5datfile(filename.c_str());
  datfile->init_write("u",ui);
}
void tdshell::set_params(double alphapar,double nupar, double nuLpar){
  double k;
  alpha=alphapar;
  nu=nupar;
  nuL=nuLpar;
  for(int l=0;l<d.size1();l++){
    k=(*(snw->get_nodesp()))[l]->k();
    for(int m=0;m<d.size2();m++){
      d(l,m)=pow(k,4)*nu;
    }
  }
}

void tdshell::assign_network(shell_network *nw, bool init_vectors){
  snw=nw;
  vector<node*>* nodes=snw->get_nodesp();
  int num_nodes=nodes->size();
  //  vector<vector<vector<vector<double> > > > M0(num_nodes,vector<vector<vector<double> > >(3,vector<vector<double> >(3,vector<double>(3,0.0))));
  Mijk=init_Mijk();
  //  vector<vector<complex<double> > > v0(num_nodes,vector<complex <double> >(3,complex<double>(0.0,0.0)));
  matrix<complex<double> > v0(num_nodes,3,complex<double>(0.0,0.0));
  ui=v0;
  f=v0;
  d=v0;
  double k;
  for(int l=0;l<d.size1();l++){
    k=(*nodes)[l]->k();
    for(int m=0;m<d.size2();m++){
      d(l,m)=pow(k,4)*nu;
    }
  }
       
}

vector<vector<vector<vector<double> > > > tdshell::init_Mijk(){
  vector<node*> *nodes=snw->get_nodesp();
  double kn;
  vector<double> kni;
  identity_matrix<double> delt(3,3);
  vector<vector<vector<vector<double> > > > M0(nodes->size(),vector<vector<vector<double> > >(3,vector<vector<double> >(3,vector<double>(3,0.0))));
  for(int n=0;n<nodes->size();n++){
    kn=(*nodes)[n]->k();
    kni=(*nodes)[n]->kni();
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	for(int k=0;k<3;k++){
	  M0[n][i][j][k]=kni[k]*(delt(i,j)-kni[i]*kni[j]/kn/kn);
	}
      }
    }
  }
  return M0;
}


void tdshell::set_integration_params(double t0par,double tmaxpar, double dtpar, double dtoutpar){
  t0=t0par;
  tmax=tmaxpar;
  dt=dtpar;
  dtout=dtoutpar;
}

void tdshell::set_forcing(matrix<complex<double> > const &forcing){
  /*  for(int l; l<forcing.size();l++){
    f[l]=forcing[l];
    }*/
  f=forcing;
}


void tdshell::datwrite(){
  datfile->append(ui);
}
/*
void tdshell::operator()(vector<vector<complex<double> > > const &u, vector<vector<complex<double> > > &dudt, const double t){
  dudt=u;
}
*/
void tdshell::operator()(matrix<complex<double> > const &u, matrix<complex<double> > &dudt, const double t){
  vector<node*> *nodes;
  node_connection *ncon;
  int np,npp;
  complex<double> I(0.0,-1.0);
  nodes=snw->get_nodesp();
  // sum over j,k,l
  for(int n=0;n<nodes->size();n++){
    for(int i=0;i<3;i++){
      dudt(n,i)=-d(n,i)*u(n,i)+f(n,i);
      for(int l=0;l<(*nodes)[n]->num_connections();l++){
	ncon=(*nodes)[n]->connection(l);
	np=ncon->l1->gid();
	npp=ncon->l2->gid();
	for(int j=0;j<3;j++){
	  for(int k=0;k<3;k++){
	    dudt(n,i)+=I*Mijk[n][i][j][k]*(conj(u(np,j))*conj(u(npp,k)) + conj(u(npp,j))*conj(u(np,k)));
	  }
	}
      }
    }
  }
}
void tdshell::run(){
  t=t0;
  int ntmax=(tmax-t0)/dtout;
  //  runge_kutta4<matrix<complex<double> > > stepper;
  adams_bashforth_moulton< 4 , matrix<complex<double> > > stepper;
  //  runge_kutta_fehlberg78<matrix<complex<double> > > stepper;
  datfile->append(ui);
  cout<<"t="<<t<<"\n";
  const complex<double> i=1i;
  double th;
  auto_cpu_timer tim;
  //  vector<complex<double> > vi=flatten(ui);
  for(int l=0;l<ntmax;l++){
    integrate_const(stepper,(*this),ui,t,t+dtout,dt);
    tim.report();
    t+=dtout;
    cout<<"t="<<t<<"\n";
    datfile->append(ui);
  }
}

