#include <iostream>
#include <time.h>
#include <armadillo>
#include "OData.h"
#include "DisorderedSOTAI.h"
#include "DisorderedSSH.h"
#include "DisorderedBBH3D.h"
#include "Anderson1D.h"
#include "Anderson3D.h"
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

  sotai.setW(3);
  sotai.generateDisorder();
  //cout << "w = 3" << endl;
  //cout << sotai.getTMM(10, 0, 8)[0] << endl;

  DisorderedBBH3D bbh3d(1.1);
  int l2[3] = {4,4,4};
  //sotai.setSize(l);

  sotai.setW(3);
  sotai.generateDisorder();
  //cout << sotai.getTMM(10, 0, 8, 0)[0] << endl;
  //cout << "w = 3" << endl;
  //cout <<  bbh3d.getTMM(10,0,6)[0] << endl;
  
  Anderson3D and3d(1);
  and3d.test();
  return 0;
}
