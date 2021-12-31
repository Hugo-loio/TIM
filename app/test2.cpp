#include "TBmodOp.h"

using namespace std;
using namespace arma;

int main(){
  TBmodOp SSH(1, 2);

  int n1[1] = {0};
  int n2[1] = {1};
  SSH.set_hop(0,1,n1, -1);
  SSH.set_hop(1,0,n2, -2);

  SSH.set_sparse(false);
  cout << SSH.berry_phase_k(10,1) << endl;

  int L2[2] = {2,2};
  TBmodOp SSH2D(2,4,L2);

  int n21[2] = {0,0};
  int n22[2] = {1,0};
  int n23[2] = {0,1};

  double lambda1 = -1;
  double lambda2 = -1.5;

  //Intracell hoppings
  SSH2D.set_hop(0,1, n21, lambda1);
  SSH2D.set_hop(0,2, n21, lambda1);
  SSH2D.set_hop(1,3, n21, lambda1);
  SSH2D.set_hop(2,3, n21, lambda1);
  //Intercell hoppings
  SSH2D.set_hop(1,0, n22, lambda2);
  SSH2D.set_hop(3,2, n22, lambda2);
  SSH2D.set_hop(2,0, n23, lambda2);
  SSH2D.set_hop(3,1, n23, lambda2);

  int bc2[2] = {0,1};
  SSH2D.set_bc(bc2);

  SSH2D.set_dir(1);
  SSH2D.set_sparse(true);
  cout << SSH2D.berry_phase_k(10,2) << " << Berry phase" << endl;

  SSH2D.set_sparse(false);
  cout <<  SSH2D.berry_phase_k(10,2) << " << Berry phase" << endl;
}

