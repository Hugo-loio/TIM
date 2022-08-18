#include <iostream>
#include <time.h>
#include <armadillo>
#include "OData.h"
#include "DisorderedSOTAI.h"
#include "DisorderedSSH.h"
#include "DisorderedBBH3D.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include <thread>
#include <cstdlib>
#include "ParallelMPI.h"

using namespace std;
using namespace arma;

int sampPerJob = 4;
double intra = 1.1;
double weight = 3;

int main (int argc, char ** argv) {

  DisorderedBBH3D bbh3d(intra, 1);
  int l[3] = {5,5,5};
  bbh3d.setSize(l);
  bbh3d.setW(weight);

  double quadyz;
  double quadxz;
  double quadxy;
  for(int i = 0; i < sampPerJob; i++){
    try{
      bbh3d.generateDisorder();
      quadyz = bbh3d.getBoundQuadrupole(0);
      quadxz = bbh3d.getBoundQuadrupole(1);
      quadxy = bbh3d.getBoundQuadrupole(2);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at disorder weight = "  << weight << endl;
      i--;
    }
  }

  return 0;
}
