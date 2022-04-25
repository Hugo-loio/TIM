#include <iostream>
#include <armadillo>
#include "OData.h"
#include "BBH2D.h"
#include "BBH3D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  /*
     cout << "2D BBH" << endl;
     BBH2D bbh2D(1,2);
     bbh2D.setInterHop(2);
     bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter2_y.dat", 5, 5, 1);

     bbh2D.setInterHop(0.99);
     bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter2_x.dat", 5, 5, 0);

     bbh2D.setInterHop(1.01);
     bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter2_x.dat", 5, 5, 0);

     bbh2D.setInterHop(1);
     bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter1_x.dat", 5, 5, 0);

     bbh2D.setInterHop(0.5);
     bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter0.5_x.dat", 4, 4, 0);
     */

  cout << "3D BBH" << endl;
  BBH3D bbh3D(2,1);

  int l[3] = {1,1,1};
  bbh3D.getSupercellNestedWannierBands(argv[0], "SupercellNestedWannierBandsBBH3D_intra2_1x1x1_z.dat", l);
  l[2] = {2};
  bbh3D.getSupercellNestedWannierBands(argv[0], "SupercellNestedWannierBandsBBH3D_intra2_1x1x2_z.dat", l);


  return 0;
}
