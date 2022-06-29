#include <iostream>
#include <armadillo>
#include "OData.h"
#include "BBH2D.h"
#include "BBH3D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  cout << "BBH 2D" << endl;
  BBH2D bbh2(0.5,1);
  int l2[2] = {10,100};
  bbh2.getBoundaryHam(l2, 0, argv[0], "BoundaryHamxBBH2D_10x10_intra0.5.dat");
  l2[0] = 100;
  l2[1] = 10;
  bbh2.getBoundaryHam(l2, 1, argv[0], "BoundaryHamyBBH2D_10x10_intra0.5.dat");
  bbh2.setIntraHop(2);
  l2[0] = 10;
  l2[1] = 100;
  bbh2.getBoundaryHam(l2, 0, argv[0], "BoundaryHamxBBH2D_10x10_intra2.dat");
  l2[0] = 100;
  l2[1] = 10;
  bbh2.getBoundaryHam(l2, 1, argv[0], "BoundaryHamyBBH2D_10x10_intra2.dat");
  bbh2.setIntraHop(0.5);
  l2[0] = 10;
  l2[1] = 5;
  bbh2.getBoundaryHam(l2, 0, argv[0], "BoundaryHamxBBH2D_10x5_intra0.5.dat");
  l2[1] = 1;
  bbh2.getBoundaryHam(l2, 0, argv[0], "BoundaryHamxBBH2D_10x1_intra0.5.dat");
  l2[1] = 3;
  bbh2.getBoundaryHam(l2, 0, argv[0], "BoundaryHamxBBH2D_10x3_intra0.5.dat");
  bbh2.setIntraHop(2);
  l2[1] = 5;
  bbh2.getBoundaryHam(l2, 0, argv[0], "BoundaryHamxBBH2D_10x5_intra2.dat");
  l2[1] = 1;
  bbh2.getBoundaryHam(l2, 0, argv[0], "BoundaryHamxBBH2D_10x1_intra2.dat");
  l2[1] = 3;
  bbh2.getBoundaryHam(l2, 0, argv[0], "BoundaryHamxBBH2D_10x3_intra2.dat");

  cout << "BBH 3D" << endl;
  BBH3D bbh3(0.5,1);
  int l3[3] = {4,4,4};
  bbh3.getBoundaryHam(l3, 0, argv[0], "BoundaryHamxBBH3D_4x4x4_intra0.5.dat");
  bbh3.getBoundaryHam(l3, 1, argv[0], "BoundaryHamyBBH3D_4x4x4_intra0.5.dat");
  bbh3.getBoundaryHam(l3, 2, argv[0], "BoundaryHamzBBH3D_4x4x4_intra0.5.dat");
  bbh3.setIntraHop(2);
  bbh3.getBoundaryHam(l3, 0, argv[0], "BoundaryHamxBBH3D_4x4x4_intra2.dat");
  bbh3.getBoundaryHam(l3, 1, argv[0], "BoundaryHamyBBH3D_4x4x4_intra2.dat");
  bbh3.getBoundaryHam(l3, 2, argv[0], "BoundaryHamzBBH3D_4x4x4_intra2.dat");

  return 0;
}
