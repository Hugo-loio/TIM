#include <iostream>
#include <armadillo>
#include "OData.h"
#include "BBH3D.h"
#include "BBH2D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  cout << "3D BBH" << endl;
  BBH3D bbh3D(1,2);
  bbh3D.getNestedWannierBands(argv[0], "NestedWannierBandsBBH3D_inter2_z.dat"); 
  bbh3D.setInterHop(1);
  bbh3D.getNestedWannierBands(argv[0], "NestedWannierBandsBBH3D_inter1_z.dat"); 
  bbh3D.setInterHop(0.5);
  bbh3D.getNestedWannierBands(argv[0], "NestedWannierBandsBBH3D_inter0.5_z.dat"); 

  BBH2D bbh2D(1,2);
  bbh2D.test();
  return 0;
}
