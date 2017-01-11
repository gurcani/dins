#include <iostream>
#include <vector>
#include <complex>
#include <boost/type_traits.hpp>
#include <hdf5wrap.hpp>
#include <hdf5.h>
#include <boost/numeric/odeint.hpp>


using namespace std;
using boost::numeric::ublas::matrix;

hdf5datfile::hdf5datfile(const char *filename){
  file=H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  group=H5Gcreate(file,"fields",H5P_DEFAULT);  
}
hdf5datfile::~hdf5datfile(){
  herr_t hstat = H5Fclose(file);
}

void hdf5datfile::open_file(const char *filename){
  if((file=H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT))<0){
    cout<<"error: can not open file: " << filename <<"\n";
  }
}
/*
void hdf5datfile::init_write(vector<complex<double> > const &vec,const char*varname,hsize_t szxpar,hsize_t szypar){
  szx=szxpar;
  szy=szypar;
  dtype=H5Tcreate (H5T_COMPOUND, sizeof (complex<double>));
  H5Tinsert (dtype, "real", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert (dtype, "imag", sizeof(complex<double>)/2, H5T_NATIVE_DOUBLE);
  //  int sz=vec.size();
  hsize_t dims[]={1,szy,szx};
  hsize_t maxdims[]={H5S_UNLIMITED,szy,szx};
  hid_t space_id=H5Screate_simple(3,dims,maxdims);
  hid_t par_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(par_id, 3, dims);
  data=H5Dcreate(group,varname,dtype,space_id,par_id);
  nt=0;
  }


void hdf5datfile::init_write(const char *varname, vector<vector<complex<double > > > const &vec){
  szx=vec.size();
  szy=vec[0].size();
  dtype=H5Tcreate (H5T_COMPOUND, sizeof (complex<double>));
  H5Tinsert (dtype, "real", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert (dtype, "imag", sizeof(complex<double>)/2, H5T_NATIVE_DOUBLE);
  hsize_t dims[]={1,szy,szx};
  hsize_t maxdims[]={H5S_UNLIMITED,szy,szx};
  hid_t space_id=H5Screate_simple(3,dims,maxdims);
  hid_t par_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(par_id, 3, dims);
  data=H5Dcreate(group,varname,dtype,space_id,par_id);
  nt=0;
}
*/

void hdf5datfile::init_write(const char*varname,hsize_t szxpar, hsize_t szypar){
  szx=szxpar;
  szy=szypar;
  dtype=H5Tcreate (H5T_COMPOUND, sizeof (complex<double>));
  H5Tinsert (dtype, "real", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert (dtype, "imag", sizeof(complex<double>)/2, H5T_NATIVE_DOUBLE);
  hsize_t dims[]={1,szy,szx};
  hsize_t maxdims[]={H5S_UNLIMITED,szy,szx};
  hid_t space_id=H5Screate_simple(3,dims,maxdims);
  par_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(par_id, 3, dims);
  data=H5Dcreate(group,varname,dtype,space_id,par_id);
  nt=0;
}

void hdf5datfile::init_write(const char *varname, matrix<complex<double > > const &vec){
  szx=vec.size1();
  szy=vec.size2();
  dtype=H5Tcreate (H5T_COMPOUND, sizeof (complex<double>));
  H5Tinsert (dtype, "real", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert (dtype, "imag", sizeof(complex<double>)/2, H5T_NATIVE_DOUBLE);
  hsize_t dims[]={1,szy,szx};
  hsize_t maxdims[]={H5S_UNLIMITED,szy,szx};
  hid_t space_id=H5Screate_simple(3,dims,maxdims);
  par_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(par_id, 3, dims);
  data=H5Dcreate(group,varname,dtype,space_id,par_id);
  nt=0;
}

void hdf5datfile::write(const char *varname, matrix<double> const &vec){
  hsize_t wx=vec.size1();
  hsize_t wy=vec.size2();
  hsize_t dims[]={1,wy,wx};
  hsize_t maxdims[]={H5S_UNLIMITED,wy,wx};
  hid_t space_id=H5Screate_simple(3,dims,maxdims);
  hid_t wpar_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(wpar_id, 3, dims);
  hid_t wdata=H5Dcreate(group,varname,H5T_NATIVE_DOUBLE,space_id,wpar_id);
  hid_t filespace = H5Dget_space(wdata);
  
  //  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, NULL, NULL,dims, NULL);
  
  H5Dwrite(wdata,H5T_NATIVE_DOUBLE,space_id,H5S_ALL,H5P_DEFAULT,flatten(vec));
  herr_t hstat = H5Fflush(file,H5F_SCOPE_GLOBAL);
  nt=0;
}

void hdf5datfile::append(vector<double> const &vec){
  hid_t space_id,filespace;
  hsize_t dims[]={nt+1,szy,szx};
  hsize_t offset[]={nt,0,0};
  hsize_t dims_sel[]={1,szy,szx};
  hsize_t maxdims[]={H5S_UNLIMITED,szy,szx};
  H5Pset_chunk(par_id, 3, dims_sel);
  space_id=H5Screate_simple(3,dims_sel,maxdims);
  H5Dset_extent(data,dims);
  filespace = H5Dget_space(data);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,dims_sel, NULL);
  H5Dwrite(data,dtype,space_id,filespace,H5P_DEFAULT,&vec[0]);
  herr_t hstat = H5Fflush(file,H5F_SCOPE_GLOBAL);
  nt++;
}

void hdf5datfile::append(vector<vector<complex <double> > > const &vec){
  hid_t space_id,filespace;
  hsize_t dims[]={nt+1,szy,szx};
  hsize_t offset[]={nt,0,0};
  hsize_t dims_sel[]={1,szy,szx};
  hsize_t maxdims[]={H5S_UNLIMITED,szy,szx};
  H5Pset_chunk(par_id, 3, dims_sel);
  space_id=H5Screate_simple(3,dims_sel,maxdims);
  H5Dset_extent(data,dims);
  filespace = H5Dget_space(data);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,dims_sel, NULL);
  vector <complex<double> > dvec=flatten(vec);
  H5Dwrite(data,dtype,space_id,filespace,H5P_DEFAULT,&dvec[0]);
  herr_t hstat = H5Fflush(file,H5F_SCOPE_GLOBAL);
  nt++;
}

void hdf5datfile::append(matrix<complex <double> > const &vec){
  hid_t space_id,filespace;
  hsize_t dims[]={nt+1,szy,szx};
  hsize_t offset[]={nt,0,0};
  hsize_t dims_sel[]={1,szy,szx};
  hsize_t maxdims[]={H5S_UNLIMITED,szy,szx};
  dtype=H5Tcreate (H5T_COMPOUND, sizeof (complex<double>));
  H5Tinsert (dtype, "real", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert (dtype, "imag", sizeof(complex<double>)/2, H5T_NATIVE_DOUBLE);
  H5Pset_chunk(par_id, 3, dims);
  space_id=H5Screate_simple(3,dims_sel,maxdims);
  H5Dset_extent(data,dims);
  filespace = H5Dget_space(data);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,dims_sel, NULL);
  //  vector <complex<double> > dvec=flatten(vec);
  //  complex<double> *dvec=&vec.data()[0];
  H5Dwrite(data,dtype,space_id,filespace,H5P_DEFAULT,flatten(vec));
  herr_t hstat = H5Fflush(file,H5F_SCOPE_GLOBAL);
  nt++;
}

vector <complex<double > > flatten(vector<vector<complex<double> > > vec){
  vector<complex<double> > ovec;
  for(vector<vector<complex<double> > >::const_iterator iv=vec.begin();iv!=vec.end();++iv){
    ovec.insert(ovec.end(),(*iv).begin(),(*iv).end());
  }
  return ovec;
}

complex<double> * flatten(matrix<complex<double> > const &mat){
  complex<double> * dat= new complex<double>[mat.size1()*mat.size2()];
  for(int l=0;l<mat.size1();l++){
    for(int m=0;m<mat.size2();m++){
      dat[m*mat.size1()+l]=mat(l,m);
    }
  }
  return dat;
}

double * flatten(matrix<double > const &mat){
  double * dat= new double[mat.size1()*mat.size2()];
  for(int l=0;l<mat.size1();l++){
    for(int m=0;m<mat.size2();m++){
      dat[m*mat.size1()+l]=mat(l,m);
    }
  }
  return dat;
}
