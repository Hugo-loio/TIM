#include <iostream>
#include <time.h>
#include <armadillo>
#include "OData.h"
#include "DisorderedSOTAI.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include <thread>
#include <cstdlib>

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  DisorderedSOTAI sotai(1.1);
  int l[2] = {20,20};
  sotai.setSize(l);
  sotai.setW(0);
  for(int i = 0; i < 10; i++){
    sotai.generateDisorder();
    sotai.getIPR(10);
  }
  //sotai.getTMM(100, 0, 10);

  return 0;
}
