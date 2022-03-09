#include "BoundaryWilsonH.h"

BoundaryWilsonH::BoundaryWilsonH(Hamiltonian * h, int dir, int n, int m) : wilson(h), Hamiltonian(h->getNDim() -1){
  ham = h;
  this->dir = dir;
  wilson.setLoopDir(dir);
  nOcc = m;
  nK = n;
  isSparse = false;
}

BoundaryWilsonH::~BoundaryWilsonH(){
}

cx_mat BoundaryWilsonH::H(double * k){
  wilson.setLoopStart(k);
  complex<double> ii(0,1);
  return (-ii)*logmat(wilson.wilsonLoop(nK, nOcc));
}

sp_cx_mat BoundaryWilsonH::spH(double * k){
  cout << "Hamiltonian resulting of Wilson loop is not a sparse matrix" << endl;
  return sp_cx_mat();
}
