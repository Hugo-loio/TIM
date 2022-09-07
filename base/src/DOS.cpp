#include "DOS.h"

DOS::DOS(Hamiltonian * ham){
  this->ham = ham;
}

DOS::~DOS(){
  if(momentsFound){
    delete[] mu;
  }
}

DOS::kpm(double en, int nMoments, int nRandVecs){
  if(ham->getIsSparse()){
    findRescaling(k);
    if(momentsFound){
      if(this->nMoments != nMoments){
	momentsFound = false;
      }
      if(this->nRandVecs != nRandVecs){
	momentsFound = false;
      }
      if(!momentsFound){
	delete[] mu;
      }
    }
    if(!momentsFound){
      mu = new double[nMoments];
      this->nMoments = nMoments;
      this->nRandVecs = nRandVecs;
      calculateMoments()
    }
  } 
  else{
    cout << __PRETTY_FUNCTION__ << " hasn't been implemented for dense Hamiltonians" << endl;
    return 0;
  }
}

DOS::jacksonKernel(int n){
}

DOS::findRescaling(double * k){
  cx_vec eigVal;
  eigs_gen(eigVal, ham->spH(k), 1, "lr");
  double eMax = eigVal[0].real();
  eigs_gen(eigVal, ham->spH(k), 1, "sr");
  double eMin = eigVal[0].real();

  a = (eMax - eMin)/(2 - err);
  b = (eMax + eMin)/2;
}

