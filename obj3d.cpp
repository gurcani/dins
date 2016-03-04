#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include "shell.hpp"
#include "node.hpp"
#include "obj3d.hpp"
using namespace std;

obj3d::obj3d(string file){
  filename=file;
}

void obj3d::shell_to_obj3d(shell *sh){
  vector<node*> nodes=sh->get_nodes();
  vector<vector<int> > sfaces=sh->get_faces();
  vector<int> sf;
  int if0;
  //vector<cartesian_pos > shvs;
  //vector<cartesian_pos > shnms;
  if0=verts.size();
  for (vector<node*>::const_iterator nd = nodes.begin(); nd != nodes.end(); ++nd){
    verts.push_back((*nd)->node_position());
  }
  for (int l=0;l<sfaces.size();l++){
    sf=sfaces[l];
    for(vector<int>::iterator fi = sf.begin(); fi!=sf.end(); ++fi)
      (*fi) += if0;
    faces.push_back(sf);
  }
  normals=compute_normals(verts,faces);
}

void obj3d::print(){
  buffer<<"o composite\n";
  for (vector<cartesian_pos>::const_iterator v = verts.begin(); v != verts.end(); ++v){
    buffer<<"v "<<(*v).x<<" "<<(*v).y<<" "<<(*v).z<<"\n";
  }
  for (vector<cartesian_pos>::const_iterator vn = normals.begin(); vn != normals.end(); ++vn){
    buffer<<"vn "<<(*vn).x<<" "<<(*vn).y<<" "<<(*vn).z<<"\n";
  }
  for (int l=0;l<faces.size();l++){
    buffer<<"f ";
    for(int m=0;m<faces[l].size();m++){
      buffer<<faces[l][m]+1<<"//"<<l+1<<" ";
    }
    buffer<<"\n";
  }
  cout<<buffer.str();
};
/*
void obj3d::shell_to_obj3d(shell *sh){
  vector<node*> nodes=sh->get_nodes();
  vector<vector<int> > faces=sh->get_faces();
  vector<cartesian_pos > verts;
  vector<cartesian_pos > normals;
  buffer<<"o shell"<<sh->get_id()<<"\n";
  for (vector<node*>::const_iterator nd = nodes.begin(); nd != nodes.end(); ++nd){
    buffer<<"v "<<(*nd)->posx()<<" "<<(*nd)->posy()<<" "<<(*nd)->posz()<<"\n";
    verts.push_back((*nd)->node_position());
  }
  normals=compute_normals(verts,faces);
  for (vector<cartesian_pos>::const_iterator vn = normals.begin(); vn != normals.end(); ++vn){
    buffer<<"vn "<<(*vn).x<<" "<<(*vn).y<<" "<<(*vn).z<<"\n";
  }
  for (int l=0;l<faces.size();l++){
    buffer<<"f ";
    for(int m=0;m<faces[l].size();m++){
      buffer<<faces[l][m]+1<<"//"<<l+1<<" ";
    }
    buffer<<"\n";
  }
  /*
  if(sh->is_icosa()){
    for(l=0;l<20;l++){
      buffer<<"f "faces_icosa[l][0]-1<<"//"<<l<<" ";
      buffer<<faces_icosa[l][1]-1<<"//"<<l<<" ";
      buffer<<faces_icosa[l][2]-1<<"//"<<l<<"\n";
    }
  }
  else{
    for(l=0;l<12;l++){
      buffer<<"f "faces_dodeca[l][0]-1<<"//"<<l<<" ";
      buffer<<faces_dodeca[l][1]-1<<"//"<<l<<" ";
      buffer<<faces_dodeca[l][2]-1<<"//"<<l<<" ";
      buffer<<faces_dodeca[l][3]-1<<"//"<<l<<" ";
      buffer<<faces_dodeca[l][4]-1<<"//"<<l<<"\n";
    }    
  }
}  */

vector<cartesian_pos> compute_normals(vector<cartesian_pos> verts, vector<vector<int> > faces){
  cartesian_pos v1,v2,vn,p0,p1,p2;
  double norm;
  vector<cartesian_pos> vns;
  for (vector<vector<int> >::const_iterator fc = faces.begin(); fc != faces.end(); ++fc){
    p0=verts[(*fc)[0]];
    p1=verts[(*fc)[1]];
    p2=verts[(*fc)[2]];
    v1=(cartesian_pos){p1.x-p0.x,p1.y-p0.y,p1.z-p0.z};
    v2=(cartesian_pos){p2.x-p0.x,p2.y-p0.y,p2.z-p0.z};
    vn=cross_prod(&v1,&v2);
    norm=sqrt(vn.x*vn.x+vn.y*vn.y+vn.z*vn.z);
    vn=(cartesian_pos){vn.x/norm,vn.y/norm,vn.z/norm};
    vns.push_back(vn);
  }
  return vns;
};

cartesian_pos cross_prod(cartesian_pos *a, cartesian_pos *b){
  cartesian_pos c;
  c.x=a->y*b->z-a->z*b->y;
  c.y=a->z*b->x-a->x*b->z;
  c.z=a->x*b->y-a->y*b->x;
  return c;
}
