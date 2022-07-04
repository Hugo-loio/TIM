#include "BoundaryGreenH.h"
#include <iostream>
#include <chrono>

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
    /*
       if(i == 1){
       auto t1 = chrono::high_resolution_clock::now();
       v = ham->blockH(i-1,i,k);
       auto t2 = chrono::high_resolution_clock::now();
       res = (-ham->blockH(i,i,k) - v*res*(v.t())).i();
       auto t3 = chrono::high_resolution_clock::now();
       auto d1 = chrono::duration_cast<chrono::microseconds>(t2 - t1);
       auto d2 = chrono::duration_cast<chrono::microseconds>(t3 - t2);
       cout << "Block ham: " << d1.count() << endl; 
       cout << "iteration: " << d2.count() << endl; 
       }
       else{
       */
    v = ham->blockH(i-1,i,k);
    res = (-ham->blockH(i,i,k) - v*res*(v.t())).i();
    //}
  }
  if(lowerBound){
    int begin = blockSize - lowerBlockSize -1;
    int end = blockSize -1;
    return res.submat(begin, begin, end, end);
  }
  else{
    return res;
  }
}

cx_mat BoundaryGreenH::H(double * k){
  return boundaryGreenFunc(k).i();
}

void BoundaryGreenH::setLowerBound(bool lowerBound, int blockSize){
  this->lowerBound = lowerBound;
  lowerBlockSize = blockSize;
}
