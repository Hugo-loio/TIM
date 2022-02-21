#include "Hop.h"
#include <iostream>

using namespace std;

Hop::Hop(int nOrb1, int nOrb2, int * n, complex<double> hop, int nDim){
  this->nOrb1 = nOrb1;
  this->nOrb2 = nOrb2;
  this->n = new int[nDim];
  for(int i = 0; i < nDim; i++){
    this->n[i] = n[i];
  }
  this->hop = hop;
  this->nDim = nDim;
}

Hop::Hop(const Hop & copy){
  nOrb1 = copy.nOrb1;
  nOrb2 = copy.nOrb2;
  nDim = copy.nDim;
  hop = copy.hop;
  n = new int[nDim];
  for(int i = 0; i < nDim; i++){
    n[i] = copy.n[i];
  }
}

Hop::~Hop(){
  delete[] n;
}

int Hop::getMaxN(){
  int max = 0;
  for(int i = 0; i < nDim; i++){
    if(n[i] > max){
      max = n[i];
    }
  }
  return max;
}
