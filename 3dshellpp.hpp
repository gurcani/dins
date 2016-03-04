#ifndef __3dshellpp_hpp
#define __3dshellpp_hpp

#include <vector>
#include <complex>
#include "hdf5wrap.hpp"
#include "node.hpp"
#include "shell.hpp"
#include "shell_network.hpp"

using namespace std;
class 3dshellpp {
private:
shell_network *nw;
double alpha;
double nu,nuL;
vector<double> a,b,c;
vector<complex<double> > f;
double t;
double t0;
double tmax;
double dt;
double dtout;
vector <vector<complex<double> > > u;
public:
3dshellhpp();
~3dshellhpp();
void init_wparams(double alpha,double nu, double nuL);
void operator()(vector<vector<complex<double> > > const &u, vector<vector<complex<double> > > &dudt, const double t);
void set_integration_params(double t0par,double tmaxpar, double dtpar, double dtoutpar);
void set_forcing(vector<vector<complex<double> > > const &);
};

#endif
