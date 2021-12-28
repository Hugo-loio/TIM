// basic file operations
#include <iostream>
#include "odata.h"
#include "TBmod.h"
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {

  cout << "SSH" << endl;

  int L[1] = {2};
  TBmod SSH(1, 2, L);

  int n1[1] = {0};
  int n2[1] = {1};
  SSH.set_hop(0,1, n1, -1);
  SSH.set_hop(1,0, n2, -2);

  cout << "PBCs\n" <<  SSH.get_rH() << endl;

  int bc[1] = {0};
  SSH.set_bc(bc);

  cout << "No PBCs\n" << SSH.get_rH() << endl;

  bc[0] = 2;
  SSH.set_bc(bc);

  double twisted[1] = {M_PI/2};
  SSH.set_twists(twisted);

  cout << "Pi twisted BCs\n" << SSH.get_rH() << endl;

  cout << "SSH2D" << endl;


  int L2[2] = {2,2};
  TBmod SSH2D(2,4,L2);

  int n21[2] = {0,0};
  int n22[2] = {1,0};
  int n23[2] = {0,1};
  //Intracell hoppings
  SSH2D.set_hop(0,1, n21, -1);
  SSH2D.set_hop(0,2, n21, -1);
  SSH2D.set_hop(1,3, n21, -1);
  SSH2D.set_hop(2,3, n21, -1);
  //Intercell hoppings
  SSH2D.set_hop(1,0, n22, -2);
  SSH2D.set_hop(3,2, n22, -2);
  SSH2D.set_hop(2,0, n23, -2);
  SSH2D.set_hop(3,1, n23, -2);

  ofstream f;
  f.open("test.txt");

  //f << "PBCs\n" <<  SSH2D.get_rH() << endl;
  SSH2D.get_rH().print(f,"PBCs\n");

  int bc2[2] = {0,0};
  SSH2D.set_bc(bc2);

  SSH2D.get_rH().print(f,"No PBCs\n");

  f.close();

  /*
     bc[0] = 2;
     SSH.set_bc(bc);

     double twisted[1] = {M_PI/2};
     SSH.set_twists(twisted);

     cout << "Pi twisted BCs\n" << SSH.get_rH() << endl;
     */

  return 0;
}
