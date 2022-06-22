#include <iostream>
#include <armadillo>
#include "OData.h"
#include "BBH3D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  cout << "3D BBH" << endl;
  BBH3D bbh3D(2,1);

  int l[3] = {1,1,1};
  int n[2] = {10,10};
  bbh3D.getSupercellNestedWannierBands(argv[0], "SupercellNestedWannierBandsBBH3D_intra2_1x1x1_z.dat", l, n);
  l[2] = {2};
  bbh3D.getSupercellNestedWannierBands(argv[0], "SupercellNestedWannierBandsBBH3D_intra2_1x1x2_z.dat", l, n);


  return 0;
}
