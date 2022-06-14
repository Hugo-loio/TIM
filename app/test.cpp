#include <iostream>
#include <armadillo>
#include "OData.h"
#include "DisorderFunctions.h"
#include "DisorderedSOTAI.h"
#include "BBH2D.h"
#include <thread>
#include <cstdlib>

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  /*
     for(int i = 0; i < 100; i++){
     cout << uniform(0,1) << endl;
     }
     */

  DisorderedSOTAI sotai(1.1);
  sotai.setW(1);

  int l[2] = {5,5};
  cx_mat h = sotai.getHam(l);

  OData o(argv[0], "testH.dat");
  o.matrixWeights(h);

  //int l[2] = {2,10};
  //vec eigVal;
  //cx_mat eigVec;
  //eig_sym(eigVal, eigVec, sotai.getHam(l));
  //sotai.getQuadrupoleManyBody(l);
  //cout << sotai.getHam(l) << endl;
  //sotai.getTopInv(l);

  //unsigned int n = thread::hardware_concurrency();
  //cout << n << " concurrent threads are supported.\n";

  BBH2D bbh;
  bbh.test();

  return 0;
}
