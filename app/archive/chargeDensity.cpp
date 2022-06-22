#include <iostream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include "SOTAI.h"

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
  bbh2D.setInterHop(0.99);
  bbh2D.getChargeDensity(argv[0], "ChargeDensityBBH2D_inter0.99_delta.dat", 20, 20, 2);
  bbh2D.setInterHop(1.01);
  bbh2D.getChargeDensity(argv[0], "ChargeDensityBBH2D_inter1.01_delta.dat", 20, 20, 2);

  cout << "3D BBH" << endl;
  BBH3D bbh3D(1,2);
  int l[3] = {8,8,8};
  bbh3D.getChargeDensity(argv[0], "ChargeDensityBBH3D_inter2.dat", l, 4);
  bbh3D.setInterHop(1);
  bbh3D.getChargeDensity(argv[0], "ChargeDensityBBH3D_inter1.dat", l, 4);
  bbh3D.setInterHop(0.5);
  bbh3D.getChargeDensity(argv[0], "ChargeDensityBBH3D_inter0.5.dat", l, 4);

  bbh3D.setOnSite(1e-3);
  bbh3D.setInterHop(2);
  bbh3D.getChargeDensity(argv[0], "ChargeDensityBBH3D_inter2_delta.dat", l, 4);
  bbh3D.setInterHop(1);
  bbh3D.getChargeDensity(argv[0], "ChargeDensityBBH3D_inter1_delta.dat", l, 4);
  bbh3D.setInterHop(0.5);
  bbh3D.getChargeDensity(argv[0], "ChargeDensityBBH3D_inter0.5_delta.dat", l, 4);

  cout << "SOTAI" << endl;
  SOTAI sotai(2);
  sotai.setOnSite(1e-3);
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m2.dat", 20, 20, 2);
  sotai.setM(0.5);
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m0.5.dat", 20, 20, 2);
  sotai.setM(1);
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m1.dat", 20, 20, 2);
  sotai.setM(0.99);
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m0.99.dat", 20, 20, 2);
  sotai.setM(1.01);
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m1.01.dat", 20, 20, 2);

  return 0;
}
