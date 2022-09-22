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

  Anderson1D a(1);
  a.setW(3);
  int l[1] = {2};
  a.setSize(l);
  //a.test();
  cout << "E = 1:\n" << 1/a.getTMM(1, 1) << endl;
  cout << "E = 1:\n" << 1/a.getTMM(10, 1) << endl;
  l[0] = {4};
  a.setSize(l);
  cout << "4 layers" << endl;
  cout << "E = 1:\n" << 1/a.getTMM(1, 1) << endl;
  //cout << "E = 1:\n" << 1/a.getTMM(5, 1) << endl;
  //cout << "E = 1:\n" << 1/a.getTMM(10, 1) << endl;

  return 0;
}
