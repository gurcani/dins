#include <vector>
#include <math.h>
#include <iostream>

#include "node.hpp"
using namespace std;

node::node(int n,spherical_angle thphi,double r){
  shell_id=n;
  set_position_spherical(thphi,r);
}

node::node(int n,cartesian_pos pos){
  shell_id=n;
  set_position(pos);
}

node::node(int n){
  shell_id=n;
}

node::~node(){
}

void node::add_connection(node *nd1,node *nd2){
  node_connection *con=new node_connection();
  con->l1=nd1;
  con->l2=nd2;
  connected_nodes.push_back(con);
}
void node::set_position(cartesian_pos pos){
  position=pos;
  radius=sqrt(pos.x*pos.x+pos.y*pos.y+pos.z*pos.z);
  angles.theta=acos(pos.z/radius);
  angles.phi=atan(pos.y/pos.z);
}
void node::set_position(double px, double py, double pz){
  set_position((cartesian_pos){px,py,pz});
}
void node::set_position_spherical(spherical_angle thphi, double r){
  angles=thphi;
  radius=r;
  position.x=r*sin(thphi.theta)*cos(thphi.phi);
  position.y=r*sin(thphi.theta)*sin(thphi.phi);
  position.z=r*cos(thphi.theta);  
}
void node::set_position_spherical(double theta, double phi,double r){
  set_position_spherical((spherical_angle){theta,phi},r);
}
node_icosa::node_icosa(int n, double r) : node(n){
  int nbase=n%6;
  int nhalf=(n>=6);
  double th,ph;
  if(nbase==0){
    th=0.0+PI*nhalf;
    ph=0.0;
  }
  else{
    th=0.5*PI-atan(0.5)*(0.5-nhalf)*2.0;
    ph=(2.0*PI/5.0*(nbase-1)+PI*nhalf);
    ph=ph-(2.0*PI)*(ph>=2.0*PI);
  }
  set_position_spherical(th,ph,r);
  if(n>11 || n<0)
    cerr<<"ERROR: We return something, but the node number is out of bounds!\n";
}

node_dodeca::node_dodeca(int n, double r) : node(n){
  int nbase=n%10;
  int nhalf=(n>=10);
  double gr=(1.0+sqrt(5.0))/2.0;
  double alp=asin(gr/sqrt(3.0))-acos(gr/sqrt(gr+2.0));
  double beta=atan(2.0*gr*gr);
  double th,ph;
  th=(nbase<5)*alp+(nbase>=5)*beta;
  th=((nhalf) ? (PI-th) : (th) );
  ph=PI/5.0+2.0*PI/5.0*(nbase%5)+PI*nhalf;  
  ph=ph-(2.0*PI)*(ph>=2.0*PI);
  set_position_spherical(th,ph,r);
  if(n>19 || n<0)
    cerr<<"ERROR: We return something, but the node number is out of bounds!\n";
}
