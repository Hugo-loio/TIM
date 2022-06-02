#include "BoundaryGreenH.h"
#include <iostream>

using namespace arma;
using namespace std;

BoundaryGreenH::BoundaryGreenH(Hamiltonian * ham, int nOrb, int * l) : ham(ham), Hamiltonian(ham->getNDim()-1){
  blockSize = nOrb;
  for(int i = 0; i < nDim; i++){
    blockSize *= l[i];
  }
  nLayers = l[nDim];
  isSparse = false;

  //cout << blockSize << " " << nLayers << " " << nDim << endl;
}

BoundaryGreenH::~BoundaryGreenH(){
}

cx_mat BoundaryGreenH::boundaryGreenFunc(double * k){
}

cx_mat BoundaryGreenH::H(double * k){
  return boundaryGreenFunc(k).i();
}
