#include <vector>
#include <math.h>
#include <iostream>

#include "node.hpp"
#include "shell.hpp"
#include "shell_network.hpp"
#include "gviz.hpp"
#include "obj3d.hpp"
#include "tdshell.hpp"
#include <boost/numeric/odeint.hpp>


using namespace std;
using boost::numeric::ublas::matrix;

int main(void){
  tdshell tds("test.h5",40,1.0);
  tds.set_params(1.0,1.0e-13,1.0e3);
  tds.set_integration_params(0.0,100.0,1.0e-5,1.0e-2);
  matrix<complex<double> > *vec=tds.get_ui();
  tds.hdf5fp()->write("k",tds.network()->get_kni());
    // initialization has to be symmetric...
  /*  for(int l=0;l<vec->size1();l++){
    for(int m=0;m<vec->size2();m++){
      (*vec)(l,m)=complex<double>(0.1,0.1);
    }
    }*/
  
  (*vec)=matrix<complex<double> >(vec->size1(),vec->size2(),complex<double>(0.0001,0.0001));
  cout<<"sz1="<<vec->size1()<<",sz2="<<vec->size2()<<"\n";
  tds.network()->symmetrize_shells((*vec));
  //  tds.force_shell(2,complex<double>(0.1,0.1));
  matrix<complex<double> > f = matrix<complex<double> >(vec->size1(),vec->size2());
  tds.network()->force_shell(2,complex<double>(0.01,0.01),f);
  tds.network()->force_shell(3,complex<double>(0.01,-0.01),f);
  tds.set_forcing(f);
  cout<<"f(22,0)="<<tds.get_forcing()(22,0)<<"\n";
  //  tds.datwrite();
  tds.run();
  tds.close();
  /*  shell_network mesh;
  gviz graph("test.gv");
  obj3d obj("test.obj");
  double g=sqrt((1.0+sqrt(5.0))/2.0);
  mesh.setup_nodes(10,g,1.0);
  mesh.setup_connections();
  
  vector<shell*> shells=mesh.get_shells();
  vector<node*> nodes=shells[0]->get_nodes();
  graph.shelltosubgraph(shells[0]);
  graph.shelltosubgraph(shells[1]);
  graph.shelltosubgraph(shells[2]);
  for(vector<shell*>::const_iterator sh = shells.begin(); sh != shells.end();++sh){
    obj.shell_to_obj3d(*sh);
  }
  for (vector<node*>::const_iterator nd = nodes.begin(); nd != nodes.end(); ++nd){
    for(int l=0;l<(*nd)->num_connections();l++){
      //cout<<"{"<<(*nd)->gid()<<","<<(*nd)->connection(l)->l1->gid()<<","<<(*nd)->connection(l)->l2->gid()<<"}\n";
      graph.add_connection((*nd),(*nd)->connection(l),l);
    }
    //cout<<"bla";
  }
  /*    for (vector<shell*>::const_iterator sh = shells.begin(); sh != shells.end(); ++sh){
    (*sh)->print_nodes();
    cout<<"==============================\n";
    }*/
  //  obj.print();
  //  graph.print();
    //  ish.print_nodes();
}
