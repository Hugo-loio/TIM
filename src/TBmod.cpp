#include "TBmod.h"
#include <iostream>

using namespace arma;
using namespace std;

TBmod::TBmod(int ndim, int norb, int * L){
  this->ndim = ndim;
  this->norb = norb;
  if(!check_dim(L)){
    cout << "Array of system size doesn't correspond to the system dimensionality " << endl;
    this->L = NULL;
  }
  this->pcb = new bool[ndim];
  this->theta = new double[ndim];
  for(int i = 0; i < ndim; i++){
    this->pcb[i] = true; 
    this->theta[i] = 0;
  }
}

TBmod::~TBmod(){
  if(L != NULL){
    delete[] L;
  }
  delete[] pcb;
  delete[] theta;
}

template<typename T>
bool TBmod::check_dim(T * vec){
  bool res = false;
  if(sizeof(vec)/sizeof(vec[0]) == ndim){
    res = true;
  }
  return res;
}

void TBmod::set_hop(int norb1, int norb2, int * n, double hop){
  if(check_dim(n)){
    this->hop.push_back(Hop(norb1, norb2, n, hop, ndim));
  }
  else{
    cout << "Neighbour cell array doesn't correspond to system dimension. Hopping not added." << endl;
  }
}

void TBmod::set_onsite(int norb, double en){
  int * n = new int[ndim];
  for(int i = 0; i < ndim; i++){
    n[i] = 0;
  }
  hop.push_back(Hop(norb, norb, n, en, ndim));
  delete[] n;
}

void TBmod::set_periodic(bool * pbc){
  if(check_dim(pcb)){
    for(int i = 0; i < ndim; i++){
      this->pcb[i] = pbc[i];
    }
  }
  else{
    cout << "Array for periodic boundary conditions doesn't correspond to system dimensionality. Boundary conditions not changed" << endl;
  }
}

void TBmod::set_twisted(double * theta){
  if(check_dim(theta)){
    for(int i = 0; i < ndim; i++){
      this->theta[i] = theta[i];
    }
  }
}
