#include "BoundaryGreenH.h"
#include <iostream>

using namespace arma;
using namespace std;

BoundaryGreenH::BoundaryGreenH(Hamiltonian * ham, int nOrb, int * l) : ham(ham), Hamiltonian(ham->getNDim()-1){
  blockSize = nOrb/2;
  for(int i = 0; i < nDim; i++){
    blockSize *= l[i];
  }
  nLayers = 2*l[nDim];
  isSparse = false;

  cout << blockSize << " " << nLayers << " " << nDim << endl;
}

BoundaryGreenH::~BoundaryGreenH(){
}

cx_mat BoundaryGreenH::boundaryGreenFunc(double * k){
  if(ham->getIsSparse() == true){
    sp_cx_mat H = ham->spH(k);
    cx_mat m(2*blockSize, 2*blockSize, fill::zeros);
    m.submat(blockSize, 0, 2*blockSize -1, blockSize -1) = eye<cx_mat>(blockSize, blockSize);
    cx_mat mN = eye<cx_mat>(2*blockSize, 2*blockSize);
    cx_mat v1;
    cx_mat v2 = cx_mat(H.submat((nLayers-2)*blockSize, (nLayers-1)*blockSize, (nLayers-1)*blockSize - 1, nLayers*blockSize - 1));
    cx_mat vN = v2;
    cx_mat g;
    int n = blockSize*nLayers;
    for(int i = nLayers; i > 0; i--){
      cout << "i " << i << "\n " << v2 << endl;
      v1 = v2.t().i();
      v2 = cx_mat(H.submat((i-3)*blockSize, (i-2)*blockSize, (i-2)*blockSize - 1, (i-1)*blockSize - 1)).t().i();
      g = -cx_mat(H.submat((i-1)*blockSize, (i-1)*blockSize, i*blockSize - 1, i*blockSize - 1)).i();
      m.submat(0, 0, blockSize - 1, blockSize -1) = v1*g;
      m.submat(0, blockSize, blockSize - 1, 2*blockSize -1) = -v1*v2;
      mN *= m;
    }
    return mN.submat(blockSize, 0, 2*blockSize -1, blockSize -1)*(mN.submat(0,0,blockSize - 1, blockSize - 1).i())*(vN.t().i());
  }
}

cx_mat BoundaryGreenH::H(double * k){
  return boundaryGreenFunc(k).i();
}
