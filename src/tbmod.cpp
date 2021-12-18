#include "tbmod.h"
#include <iostream>

using namespace arma;
using namespace std;

tbmod::tbmod(int ndim, int norb, int * L = NULL){
  this->ndim = ndim;
  this->norb = norb;
  if(check_dim(L)){
    cout << "Array of system size doesn't correspond to the system dimensionality " << endl;
    this->L = NULL;
  }
  this->pcb = NULL;
}

tbmod::~tbmod(){
  if(L != NULL){
    delete[] L;
  }
  if(pcb != NULL){
    delete[] pcb;
  }
}

template<typename T>
bool tbmod::check_dim(T * vec){
  bool res = false;
  if(sizeof(vec)/sizeof(vec[0]) == ndim){
    res = true;
  }
  return res;
}
