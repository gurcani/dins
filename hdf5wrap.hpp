#ifndef HDF5WRAP_HPP
#define HDF5WRAP_HPP
#include <hdf5.h>
#include <vector>
#include <complex>
#include <boost/numeric/odeint.hpp>
using boost::numeric::ublas::matrix;
using namespace std;

class hdf5datfile{
private:
  hid_t file;
  hid_t group;
  hid_t dtype;
  hid_t data;
  hid_t par_id;
  int rank;
  int nt;
  int szx,szy;
public:
  hdf5datfile(const char *filename);
  ~hdf5datfile();
  void open_file(const char *filename);
  //  void write(vector<double>)
  void init_write(vector<double> const &, const char *);
  void init_write(const char *, int, int);
  void init_write(const char *, vector<vector<complex<double> > > const &);
  void init_write(const char *, matrix<complex<double> > const &);
  void write(const char *, matrix<double> const &);
  void append(vector<double> const &vec);
  void append(matrix<complex <double> > const &vec);
  void close(){H5Fclose(file);}
  void append(vector<vector<complex <double> > > const &vec);
};

vector <complex<double > > flatten(vector<vector<complex<double> > > vec);
complex<double> * flatten(matrix<complex<double> > const &mat);
double * flatten(matrix<double > const &mat);
#endif
