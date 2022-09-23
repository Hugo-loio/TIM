#include <iostream>
#include <time.h>
#include <armadillo>
#include "OData.h"
#include "DisorderedSOTAI.h"
#include "DisorderedSSH.h"
#include "DisorderedBBH3D.h"
#include "Anderson1D.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include <thread>
#include <cstdlib>
#include "ParallelMPI.h"

using namespace std;
using namespace arma;

int main (int argc, char ** argv) {

  DisorderedSOTAI sotai(1.1);
  int l[2] = {10,10};
  sotai.setSize(l);

  sotai.setW(0);
  sotai.generateDisorder();
  cout << "w = 0" << endl;
  for(int i = 0; i < 10 ; i++){
    sotai.generateDisorder();
    sotai.test(argv[0]);
  }

  cout << "w = 1" << endl;
  sotai.setW(1);
  for(int i = 0; i < 10 ; i++){
    sotai.generateDisorder();
    sotai.test(argv[0]);
  }

  cout << "w = 2" << endl;
  sotai.setW(2);
  for(int i = 0; i < 10 ; i++){
    sotai.generateDisorder();
    sotai.test(argv[0]);
  }

  cout << "w = 3" << endl;
  sotai.setW(3);
  for(int i = 0; i < 10 ; i++){
    sotai.generateDisorder();
    sotai.test(argv[0]);
  }

  cout << "w = 4" << endl;
  sotai.setW(4);
  for(int i = 0; i < 10 ; i++){
    sotai.generateDisorder();
    sotai.test(argv[0]);
  }

  cout << "w = 5" << endl;
  sotai.setW(5);
  for(int i = 0; i < 10 ; i++){
    sotai.generateDisorder();
    sotai.test(argv[0]);
  }

  cout << "w = 9" << endl;
  sotai.setW(9);
  for(int i = 0; i < 10 ; i++){
    sotai.generateDisorder();
    sotai.test(argv[0]);
  }

  return 0;
}
