// basic file operations
#include <iostream>
#include "odata.h"
#include "TBmod.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  /*
     cx_mat A(4,4, fill::randu);

     cout << A << endl;

     A.resize(4,2);

     cout << A << endl;

     int a = A.size();

     int b = size(A)[0];

     cout << A.size() << " " << size(A) << " " << a << " " << b << " " << size(A)[1] <<  endl;
     */

  int L[1] = {2};
  TBmod SSH(1, 2, L);

  int n1[1] = {0};
  int n2[1] = {1};
  SSH.set_hop(0,1, n1, 1);
  SSH.set_hop(1,0, n2, 2);

  cout << "ola" << endl;

  cout << "PBCs\n" <<  SSH.get_rH() << endl;

  int bc[1] = {0};
  SSH.set_bc(bc);

  cout << "No PBCs\n" << SSH.get_rH() << endl;

  bc[0] = 2;
  SSH.set_bc(bc);

  double twisted[1] = {M_PI/2};
  SSH.set_twists(twisted);

  cout << "Pi twisted BCs\n" << SSH.get_rH() << endl;

  return 0;
}
