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

  Anderson1D a(-1);
  a.setW(3);
  int l[1] = {10};
  a.setSize(l);
  a.generateDisorder();
  //a.test();
  //cout << "E = 0: " << a.getTMM(1, 0) << endl;
  //cout << "E = 1: " << a.getTMM(10, 1) << endl;
  //cout << "E = -1: " << a.getTMM(10, -1) << endl;

  a.setSize(l);
  a.generateDisorder();
  //cout << a.getHam() << endl;
  //

  DisorderedSOTAI sotai(1.1);
  sotai.setW(5);
  int l2[2] = {100, 100};
  sotai.setSize(l);
  sotai.generateDisorder();
  cout << sotai.getLSR(10) << endl;

  return 0;
}
