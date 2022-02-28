#include "Wilson.h"
#include <iostream>

using namespace std;
using namespace arma;

Wilson::Wilson(int dim){
  k0 = new double[dim];
  for (int i = 0; i < dim; i++){
    k0[i] = 0;
  }
  this->dim = dim;
  dir = 0;
}

Wilson::~Wilson(){
  delete[] k0;
}

void Wilson::setLoopDir(int dir){
  if(dir < dim){
    this->dir = dir;
  }
  else{
    cout << "Error: Direction is larger than the number of spatial dimensions" << endl;
  }
}

void Wilson::setLoopStart(double * k0){
  for(int i = 0; i < dim; i++){
    this->k0[i] = k0[i];
  }
}

cx_mat Wilson::evOcc(int m, cx_mat h){
  vec eigVal;
  cx_mat eigVec;
  int n = size(h)[0];
  eig_sym(eigVal, eigVec, h);
  eigVec.resize(n,m);
  eigVal.resize(m);
  //cout << eigVal << eigVec << endl;
  return eigVec;
}

cx_mat Wilson::wilsonLoop(Hamiltonian & ham, int n, int m){
  double deltaK = 2*M_PI/(double)n;
  double * k = new double[dim];
  for(int i = 0; i < dim; i++){
    k[i] = k0[i];
  }
  k[dir] += deltaK;

  if(ham.getIsSparse()){
    cx_vec eigVal;
    cx_mat eigVec0, eigVec1;
    eigs_gen(eigVal, eigVec0, ham.spH(k0), m, "sr");
    eigs_gen(eigVal, eigVec1, ham.spH(k), m, "sr");
    cx_mat eigVec2 = eigVec1;

    cx_mat u = eigVec0.t() * eigVec1;

    for(int i = 2; i < n; i++){
      eigVec1 = eigVec2;
      k[dir] = k0[dir] + i*deltaK;
      eigs_gen(eigVal, eigVec2, ham.spH(k), m , "sr");
      //cout << i << "\n" << "eigVal\n" << eigVal << "eigVec\n" << eigVec2 << endl;
      u = u* (eigVec1.t() * eigVec2);
    }
    u = u * (eigVec2.t() * eigVec0);

    delete[] k;
    return u;
  }
  else{
    cx_mat eigVec0 = evOcc(m, ham.H(k0));
    cx_mat eigVec1 = evOcc(m, ham.H(k));
    cx_mat eigVec2 = eigVec1;

    cx_mat u = eigVec0.t() * eigVec1;

    for(int i = 2; i < n; i++){
      eigVec1 = eigVec2;
      k[dir] = k0[dir] + i*deltaK;
      eigVec2 = evOcc(m, ham.H(k));
      u = u* (eigVec1.t() * eigVec2);
    }
    u = u * (eigVec2.t() * eigVec0);

    delete[] k;
    return u;
  }
}

double Wilson::berryPhase(Hamiltonian & ham, int n, int m){
  return - log_det(this->wilsonLoop(ham,n,m)).imag();
}
