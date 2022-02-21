#include "Hop_old.h"
#include <iostream>

using namespace std;

Hop::Hop(int norb1, int norb2, int * n, complex<double> hop, int ndim){
  this->norb1 = norb1;
  this->norb2 = norb2;
  this->n = new int[ndim];
  for(int i = 0; i < ndim; i++){
    this->n[i] = n[i];
  }
  this->hop = hop;
  this->ndim = ndim;
}

Hop::Hop(const Hop & copy){
  norb1 = copy.norb1;
  norb2 = copy.norb2;
  ndim = copy.ndim;
  hop = copy.hop;
  n = new int[ndim];
  for(int i = 0; i < ndim; i++){
    n[i] = copy.n[i];
  }
}

Hop::~Hop(){
  delete[] n;
}

int Hop::get_maxn(){
  int max = 0;
  for(int i = 0; i < ndim; i++){
    if(n[i] > max){
      max = n[i];
    }
  }
  return max;
}
