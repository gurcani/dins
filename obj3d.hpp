#ifndef __obj3d_hpp
#define __obj3d_hpp
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "shell.hpp"
#include "node.hpp"
using namespace std;

class obj3d{
private:
  string filename;
  ostringstream buffer;
  vector<vector<int> > faces;
  vector<cartesian_pos> verts;
  vector<cartesian_pos> normals;
public:
  obj3d(string file);
  ~obj3d(){};
  //  void shells_to_obj3d(vector<shell*> shells);
  void shell_to_obj3d(shell* sh);
  void print();
};

/*
faces_icosa[20][3]={{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},{1,10,2},{2,11,3},{3,7,4},{4,8,5},{5,9,1},{10,11,2},{11,7,3},{7,8,4},{8,9,5},{9,10,1},{11,10,6},{7,11,6},{8,7,6},{9,8,6},{10,9,6}};
faces_dodeca[12][5]={{0,1,2,3,4},{0,5,18,6,1},{1,6,19,7,2},{2,7,15,8,3},{3,8,16,9,4},{4,9,17,5,0},{5,17,12,13,18},{6,18,13,14,19},{7,19,14,10,15},{8,15,10,11,16},{9,16,11,12,17},{12,11,10,14,13}};
*/
cartesian_pos cross_prod(cartesian_pos *a, cartesian_pos *b);
vector<cartesian_pos> compute_normals(vector<cartesian_pos> verts, vector<vector<int> > faces);

#endif
