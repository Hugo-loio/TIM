#include <iostream>
#include <armadillo>
#include "OData.h"
#include "DisorderFunctions.h"
#include "DisorderedSOTAI.h"
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

  DisorderedSOTAI sotai(2);
  sotai.setW(1);

  int l[2] = {40,40};
  vec eigVal;
  cx_mat eigVec;
  //eig_sym(eigVal, eigVec, sotai.getHam(l));
  sotai.getQuadrupoleManyBody(l);

  unsigned int n = thread::hardware_concurrency();
  cout << n << " concurrent threads are supported.\n";

  return 0;
}
