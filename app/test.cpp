#include <iostream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  cout << "2D BBH" << endl;
  BBH2D bbh2D(1,2);
  bbh2D.getChargeDensity(argv[0], "ChargeDensityBBH2D_inter2.dat", 20, 20, 2);
  bbh2D.setInterHop(1);
  bbh2D.getChargeDensity(argv[0], "ChargeDensityBBH2D_inter1.dat", 20, 20, 2);
  bbh2D.setInterHop(0.5);
  bbh2D.getChargeDensity(argv[0], "ChargeDensityBBH2D_inter0.5.dat", 20, 20, 2);

  bbh2D.setOnSite(1e-3);
  bbh2D.setInterHop(2);
  bbh2D.getChargeDensity(argv[0], "ChargeDensityBBH2D_inter2_delta.dat", 20, 20, 2);
  bbh2D.setInterHop(1);
  bbh2D.getChargeDensity(argv[0], "ChargeDensityBBH2D_inter1_delta.dat", 20, 20, 2);
  bbh2D.setInterHop(0.5);
  bbh2D.getChargeDensity(argv[0], "ChargeDensityBBH2D_inter0.5_delta.dat", 20, 20, 2);

  return 0;
}
