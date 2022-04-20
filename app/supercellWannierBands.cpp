#include <iostream>
#include <armadillo>
#include "OData.h"
#include "BBH2D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  cout << "2D BBH" << endl;
  BBH2D bbh2D(1,2);
  bbh2D.setInterHop(2);
  bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter2_x.dat", 1, 2, 0);
  //  bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter2_y.dat", 5, 5, 1);

  /*j
    bbh2D.setInterHop(0.99);
    bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter2_x.dat", 5, 5, 0);

    bbh2D.setInterHop(1.01);
    bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter2_x.dat", 5, 5, 0);
    */

  bbh2D.setInterHop(1);
  // bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter1_x.dat", 5, 5, 0);

  bbh2D.setInterHop(0.5);
  //bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsBBH2D_inter0.5_x.dat", 5, 5, 0);

  return 0;
}
