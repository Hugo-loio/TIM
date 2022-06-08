#include <iostream>
#include <armadillo>
#include "OData.h"
#include "DisorderFunctions.h"
#include "DisorderedSOTAI.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  /*
     for(int i = 0; i < 100; i++){
     cout << uniform(0,1) << endl;
     }
     */

  DisorderedSOTAI sotai(2);
  sotai.setW(1);

  int l[2] = {5,5};
  cx_mat h = sotai.getHam(l);

  OData o(argv[0], "testH.dat");
  o.matrixWeights(h);

  return 0;
}
