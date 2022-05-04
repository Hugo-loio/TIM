#include <iostream>
#include "DisorderedBBH2D.h"
#include "DisorderedBBH3D.h"
#include "DisorderedSOTAI.h"

using namespace std;

void printQuadrupole(DisorderedBBH2D & mod, int * l, int *n){
  for(int i = 0; i < 5; i++){
    cout << i << ": " << mod.getQuadrupoleNestedSupercell(l,n) << endl;
  }
}

int main (int arc, char ** argv) {
  DisorderedBBH2D bbh2D(2,1);
  int l[2] = {4,4};
  int n[2] = {10,10};

  cout << "p = 0" << endl;
  printQuadrupole(bbh2D, l, n);
  cout << "p = 0.1" << endl;
  bbh2D.setProbDisorder(0.1);
  printQuadrupole(bbh2D, l, n);
  cout << "p = 0.2" << endl;
  bbh2D.setProbDisorder(0.2);
  printQuadrupole(bbh2D, l, n);
  cout << "p = 0.5" << endl;
  bbh2D.setProbDisorder(0.5);
  printQuadrupole(bbh2D, l, n);
  cout << "p = 0.8" << endl;
  bbh2D.setProbDisorder(0.8);
  printQuadrupole(bbh2D, l, n);
  cout << "p = 0.9" << endl;
  bbh2D.setProbDisorder(0.9);
  printQuadrupole(bbh2D, l, n);
  cout << "p = 1" << endl;
  bbh2D.setProbDisorder(1);
  printQuadrupole(bbh2D, l, n);

  bbh2D.setProbDisorder(0.8);
  bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsDisorderedBBH2D_p0.8_4x4_x.dat", 4, 4, 0);
  bbh2D.setProbDisorder(0.9);
  bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsDisorderedBBH2D_p0.9_4x4_x.dat", 4, 4, 0);
  bbh2D.setProbDisorder(1);
  bbh2D.getSupercellWannierBands(argv[0], "SupercellWannierBandsDisorderedBBH2D_p1_4x4_x.dat", 4, 4, 0);
}
