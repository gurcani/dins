#include <vector>
#include <node.hpp>
using namespace std;

class shell_network{
private:
  vector<node*> nodes;
  vector<shell*> shells;
public:
  shell_network();
  ~shell_network();
  icosa_dodeca_setup_nodes(int N, double g, double k0);
  dodeca_icosa_setup_nodes(int N, double g, double k0);
}
