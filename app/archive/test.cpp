#include <iostream>
#include <armadillo>
#include "OData.h"
#include "DisorderFunctions.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  for(int i = 0; i < 100; i++){
    cout << uniform(0,1) << endl;
  }

  int * a = new int[0];
  delete[] a;
  return 0;
}
