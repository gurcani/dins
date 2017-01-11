#ifndef __node_hpp
#define __node_hpp
#include <vector>
using namespace std;

class node;

typedef struct node_connection_{
  node *l1,*l2;
}node_connection;

typedef struct spherical_angle_{
  double theta,phi;
}spherical_angle;

typedef struct cartesian_pos_{
  double x,y,z;
} cartesian_pos;

class node{
  int global_id;
  int shell_id;
  int anti_node_gid;
  spherical_angle angles;
  double radius;
  cartesian_pos position;
  vector<node_connection*> connected_nodes;
  //  double theta,phi;
  //  double posx,posy,posz;
private:
public:
  node(int n);
  node(int n, cartesian_pos pos);
  node(int n, spherical_angle thphi, double r);
  ~node();
  void set_position(cartesian_pos pos);
  void set_position(double px,double py,double pz);
  void set_position_spherical(spherical_angle thphi,double r=1.0);
  void set_position_spherical(double theta,double phi,double r=1.0);
  spherical_angle node_angles(){ return angles; }
  cartesian_pos node_position(){ return position; }
  void set_global_id(int n){global_id=n;};
  void set_anti_gid(int n){anti_node_gid=n;};
  int gid(){return global_id;};
  int id(){return shell_id;};
  double theta(){ return angles.theta;}
  double phi(){return angles.phi;}
  double k(){return radius;}
  vector<double> kni(){double ar[3]={position.x,position.y,position.z};
    return vector<double>(ar,ar+sizeof(ar)/sizeof(ar[0]));}; //probably this is messed up!!!
  double posx(){ return position.x;}
  double posy(){return position.y;}
  double posz(){return position.z;}
  int num_connections(){return connected_nodes.size();}
  void add_connection(node *, node *);
  node_connection* connection(int n){return connected_nodes[n];};
  int anti_gid(){return anti_node_gid;};
};

#define PI 3.141592653589793238463

class node_icosa : public node {
private:
public:
  node_icosa(int n,double r=1.0);
};
class node_dodeca : public node {
private:
public:
  node_dodeca(int n, double r=1.0);
};
#endif
