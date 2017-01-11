#ifndef __shell_hpp
#define __shell_hpp

#include <vector>
#include <iostream>
#include <complex>
#include "node.hpp"
#include <boost/numeric/odeint.hpp>

using namespace std;
using boost::numeric::ublas::matrix;

class shell{
private:
  int global_id;
  vector <node *> nodes;
  vector <vector<int> > faces;
  double kn;
public:
  shell(int n, double k);
  ~shell();
  vector<node*> get_nodes(){return nodes;};
  vector<vector<int> > get_faces(){return faces;};
  void set_nodes( vector <node*> nds){nodes=nds;};
  void set_id(int n){global_id=n;};
  int get_id(){return global_id;};
  void push_node( node *nd){nodes.push_back(nd);};
  void push_face( vector<int> f){faces.push_back(f);};
  void print_nodes();
  int size(){return nodes.size();};
  virtual bool is_icosa(){return false;};
  virtual bool is_dodeca(){return false;};
  void symmetrize( matrix<complex<double > > &);
  void force(complex<double> f, matrix<complex<double> > &fm){
    cout<<"nodes.size()="<<nodes.size()<<"\n";
    for(int l=0;l<nodes.size()/2;l++){
      cout<<"gid="<<nodes[l]->gid()<<"\n";
      fm(nodes[l]->gid(),0)=fm(nodes[l]->gid(),1)=fm(nodes[l]->gid(),2)=f;
    }
    for(int l=nodes.size()/2;l<nodes.size();l++){
      cout<<"gid*="<<nodes[l]->gid()<<"\n";
      fm(nodes[l]->gid(),0)=fm(nodes[l]->gid(),1)=fm(nodes[l]->gid(),2)=conj(f);
    }
  };
  void connect_shells(shell *sh1, shell *sh2, int order_sh0=1);
};

class icosa_shell : public shell{
private:
public:
  icosa_shell(int n, double k);
  bool is_icosa(){return true;};
  //  void symmetrize( matrix<complex<double > > &);
};

class dodeca_shell : public shell{
private:
public:
  dodeca_shell(int n, double k);
  bool is_dodeca(){return true;};
  //  void symmetrize( matrix<complex<double > > &);
};

const int connections_icosa102[6][5][2]={{{5,10},{6,11},{7,12},{8,13},{9,14}},{{1,10},{18,15},{12,7},{16,19},{3,14}},{{4,10},{17,15},{13,8},{19,16},{2,11}},{{0,11},{18,16},{14,9},{15,17},{3,12}},{{1,12},{19,17},{10,5},{16,18},{4,13}},{{2,13},{15,18},{11,6},{17,19},{0,14}}};

const int connections_icosa012[6][5][2]={{{15,10},{16,11},{17,7},{18,8},{19,9}},{{11,3},{2,6},{13,4},{6,8},{8,11}},{{12,4},{3,6},{14,5},{7,9},{9,7}},{{10,1},{13,5},{4,6},{5,8},{8,10}},{{0,6},{11,2},{14,1},{6,9},{9,11}},{{10,2},{1,6},{12,3},{5,7},{7,10}}};

const int connections_icosa120[6][5][2]={{{10,10},{11,11},{7,12},{8,13},{9,14}},{{3,10},{4,14},{11,15},{6,7},{8,19}},{{5,10},{4,11},{9,15},{7,16},{6,8}},{{1,11},{5,12},{10,16},{8,17},{6,9}},{{2,12},{1,13},{6,5},{11,17},{9,18}},{{3,13},{2,14},{6,6},{7,18},{10,19}}};

const int connections_dodeca102[10][3][2]={{{4,6},{9,7},{11,8}},{{5,6},{7,9},{10,8}},{{1,6},{11,9},{8,10}},{{2,6},{9,11},{7,10}},{{3,6},{8,11},{10,7}},{{3,8},{5,7},{6,4}},{{1,8},{4,9},{6,5}},{{2,9},{5,10},{6,1}},{{3,10},{1,11},{6,2}},{{2,7},{4,11},{6,3}}};

const int connections_dodeca012[10][3][2]={{{3,11},{10,15},{5,14}},{{1,10},{4,12},{11,16}},{{7,17},{2,11},{5,13}},{{1,14},{8,18},{3,12}},{{2,10},{9,19},{4,13}},{{0,10},{9,7},{11,8}},{{0,11},{7,9},{10,8}},{{0,12},{8,5},{11,9}},{{0,13},{7,5},{9,6}},{{0,14},{8,6},{10,7}}};

const int connections_dodeca120[10][3][2]={{{15,6},{11,7},{14,8}},{{16,6},{12,8},{10,9}},{{17,6},{13,9},{11,10}},{{18,6},{14,10},{12,11}},{{19,6},{13,7},{10,11}},{{8,7},{7,8},{10,4}},{{9,8},{8,9},{11,5}},{{12,1},{5,9},{9,10}},{{13,2},{6,10},{5,11}},{{6,7},{14,3},{7,11}}};

const int faces_icosa[20][3]={{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},{1,10,2},{2,11,3},{3,7,4},{4,8,5},{5,9,1},{10,11,2},{11,7,3},{7,8,4},{8,9,5},{9,10,1},{11,10,6},{7,11,6},{8,7,6},{9,8,6},{10,9,6}};

const int faces_dodeca[12][5]={{0,1,2,3,4},{0,5,18,6,1},{1,6,19,7,2},{2,7,15,8,3},{3,8,16,9,4},{4,9,17,5,0},{5,17,12,13,18},{6,18,13,14,19},{7,19,14,10,15},{8,15,10,11,16},{9,16,11,12,17},{12,11,10,14,13}};


//const int connections_icosa_dodeca120[6][5][2]={{{10,10},{11,11},{12,7},{13,8},{14,9}},{{10,3},{14,4},{15,11},{7,6},{19,8}},{{10,5},{11,4},{15,9},{16,7},{8,6}},{{11,1},{12,5},{16,10},{17,8},{9,6}},{{12,2},{13,1},{5,6},{17,11},{18,9}},{{13,3},{14,2},{6,6},{18,7},{19,10}}};

#endif
