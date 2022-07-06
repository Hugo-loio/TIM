#include <iostream>
#include "DisorderedBBH3D.h"
#include "DisorderedSOTAI.h"

using namespace std;

void printQuadrupole(DisorderedSOTAI & mod, int * l, int *n){
  for(int i = 0; i < 5; i++){
    cout << i << ": " << mod.getQuadrupoleNestedSupercell(l,n) << endl;
  }
}

int main (int arc, char ** argv) {
  int l[2] = {4,4};
  int n[2] = {10,10};

  cout << "SOTAI" << endl;
  DisorderedSOTAI sotai(1.1);
  cout << "W = 0" << endl;
  printQuadrupole(sotai, l, n);
  cout << "W = 1" << endl;
  sotai.setW(1);
  printQuadrupole(sotai, l, n);
  cout << "W = 2.5" << endl;
  sotai.setW(2.5);
  printQuadrupole(sotai, l, n);
  cout << "W = 3" << endl;
  sotai.setW(3);
  printQuadrupole(sotai, l, n);
}
