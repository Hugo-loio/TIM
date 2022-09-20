#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main (int argc, char ** argv) {

  sp_mat mat = sprandu<sp_mat>(1000, 1000, 0.1);
  cx_vec eigVal;
  cx_mat eigVec;

  //Shift-invert
  eigs_gen(eigVal, eigVec, mat, 10, 0.0);

  cout << eigVal << endl;

  //Test with MPI
}
