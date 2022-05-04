#include <iostream>
#include "DisorderedBBH2D.h"
#include "DisorderedBBH3D.h"
#include "DisorderedSOTAI.h"

using namespace std;

int main (int arc, char ** argv) {
  int l[2] = {4,4};
  int n[2] = {10,10};

  cout << "SOTAI" << endl;
  DisorderedSOTAI sotai(1.1);
  sotai.setOnSite(1e-2);
  cout << "W = 0" << endl;
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m1.1_w0.dat", 10, 10, 2);
  sotai.setW(1);
  cout << "W = 1" << endl;
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m1.1_w1.dat", 10, 10, 2);
  sotai.setW(2.5);
  cout << "W = 2.5" << endl;
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m1.1_w2.5.dat", 10, 10, 2);
  sotai.setW(3);
  cout << "W = 3" << endl;
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m1.1_w3.dat", 10, 10, 2);

  sotai.setM(0.5);
  sotai.setW(0);
  cout << sotai.getQuadrupoleNestedSupercell(l,n) << endl;
  cout << "W = 0" << endl;
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m0.5_w0.dat", 10, 10, 2);
  sotai.setW(1);
  cout << "W = 1" << endl;
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m0.5_w1.dat", 10, 10, 2);
  sotai.setW(2);
  cout << "W = 2" << endl;
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m0.5_w2.dat", 10, 10, 2);
  sotai.setW(3);
  cout << "W = 3" << endl;
  sotai.getChargeDensity(argv[0], "ChargeDensitySOTAI_m0.5_w3.dat", 10, 10, 2);
}
