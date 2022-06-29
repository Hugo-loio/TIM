#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(int nDim){
  this->nDim = nDim;
  theta = new double[nDim];
  for(int i = 0; i < nDim; i++){
    theta[i] = 0;
  }
}

Hamiltonian::~Hamiltonian(){
  delete[] theta;
}

Hamiltonian::Hamiltonian(const Hamiltonian & copy){
  cout << __PRETTY_FUNCTION__ << endl;
  nDim = copy.nDim;
  theta = new double[nDim];
  for(int i = 0; i < nDim; i++){
    theta[i] = copy.theta[i];
  }
  isSparse = copy.isSparse;
}

void Hamiltonian::setTwists(double * theta){
  for(int i = 0; i < nDim; i++){
    this->theta[i] = theta[i];
  }
}
