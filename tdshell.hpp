#ifndef __tdshell_hpp
#define __tdshell_hpp

#include <vector>
#include <complex>
#include "hdf5wrap.hpp"
#include "node.hpp"
#include "shell.hpp"
#include "shell_network.hpp"
#include <boost/numeric/odeint.hpp>

using namespace std;
using boost::numeric::ublas::matrix;

class tdshell{
private:
  shell_network *snw;
  hdf5datfile *datfile;
  double alpha;
  double nu,nuL;
  vector<vector<vector<vector<double> > > > Mijk;
  matrix<complex<double> > f,d;
  double t;
  double t0;
  double tmax;
  double dt;
  double dtout;
  //  vector <vector<complex<double> > > ui;
  matrix<complex<double> > ui;
public:
  tdshell();
  tdshell(int N, double k0);
  tdshell(string filename, int N, double k0);
  ~tdshell(){};
  void assign_network(shell_network *nw, bool init_vectors=true);
  void assign_hdfdatfile(hdf5datfile *file){datfile=file;};
  void set_params(double alpha,double nu, double nuL);
  vector<vector<vector<vector<double> > > > init_Mijk();
  void operator()(matrix<complex<double> > const &u, matrix<complex<double> > &dudt, const double t);
//void operator()(vector<vector<complex<double> > > const &u, vector<vector<complex<double> > > &dudt, const double t);
  void set_integration_params(double t0par,double tmaxpar, double dtpar, double dtoutpar);
  void set_forcing(matrix<complex<double> > const &);
  matrix<complex<double> > & get_forcing(){return f;};
  matrix<complex<double> >* get_ui(){return &ui;};
  void datwrite();
  void close(){datfile->close();}
  hdf5datfile *hdf5fp(){return datfile;};
  void run();
  shell_network* network(){return snw;};
  void force_shell(int l, complex<double> ff){snw->force_shell(l,ff,f);}
  //complex<double> nlterm(vector<complex<double> > const &,int, int, int);
};

#endif
