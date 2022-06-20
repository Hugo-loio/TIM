#include "BoundaryGreenH.h"
#include <iostream>

using namespace arma;
using namespace std;

BoundaryGreenH::BoundaryGreenH(Hamiltonian * ham, int blockSize, int nLayers) : ham(ham), Hamiltonian(ham->getNDim()-1){
  this->blockSize = blockSize;
  this->nLayers = nLayers;
  isSparse = false;
}

BoundaryGreenH::~BoundaryGreenH(){
}

cx_mat BoundaryGreenH::boundaryGreenFunc(double * k){
  cx_mat res = -(ham->blockH(0,0,k)).i();
  cx_mat v;
  for(int i = 1; i < nLayers; i++){
    v = ham->blockH(i-1,i,k);
    res = (-ham->blockH(i,i,k) - v*res*(v.t())).i();
  }
  return res;
}

cx_mat BoundaryGreenH::H(double * k){
  return boundaryGreenFunc(k).i();
}
