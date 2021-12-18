#include "hop.h"

hop::hop(int norb1, int norb2, int * n, double hop, int ndim){
  this->norb1 = norb1;
  this->norb2 = norb2;
  this->n = new int[ndim];
  for(int i = 0; i < ndim; i++){
    this->n[i] = n[i];
  }
  this->hop = hop;
  this->ndim = ndim;
}

hop::~hop(){
  delete[] n;
}

int hop::get_maxn(){
  int max = 0;
  for(int i = 0; i < ndim; i++){
    if(n[i] > max){
      max = n[i];
    }
  }
  return max;
}
