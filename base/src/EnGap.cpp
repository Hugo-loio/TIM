#include "EnGap.h"

EnGap::EnGap(Hamiltonian * ham){
  this->ham = ham;
}

EnGap::~EnGap(){
}

double EnGap::getGap(double en, double * k){
  if(ham->getIsSparse()){
    cx_vec eigVal;
    vec realEig;
    int nStates = 10;
    double res;

    while(1){
      eigs_gen(eigVal, ham->spH(k), nStates, en);
      if(size(eigVal)[0] != nStates){
	cout << "Found " << size(eigVal)[0] << " states" << endl;
	throw runtime_error("Diagonalization failed.");
      }
      realEig = vec(nStates, fill::zeros);
      for(int i = 0; i < nStates; i++){
	realEig[i] = eigVal[i].real();
	if(realEig[i] == en){
	  return 0;
	}
      }
      realEig = sort(realEig);

      for(int i = 0; i < nStates - 1; i++){
	if(realEig[i+1] > en && realEig[i] < en){
	  if(nStates != 10){
	    cout << __PRETTY_FUNCTION__ << " needed " << nStates << " states" << endl;
	  }
	  return realEig[i+1] - realEig[i];
	}
      }
      nStates*= 2;
    }
    return -1;
  }
  else{
    cout << __PRETTY_FUNCTION__ << " hasn't been implemented for dense Hamiltonians" << endl;
    return -1;
  }
}
