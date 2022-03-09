#include "Wilson.h"
#include <iostream>

using namespace std;
using namespace arma;

Wilson::Wilson(Hamiltonian * ham){
  this->dim = ham->getNDim();
  k0 = new double[dim];
  for (int i = 0; i < dim; i++){
    k0[i] = 0;
  }
  dir = 0;
  this->ham = ham;
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
  return eigVec;
}

cx_mat Wilson::wilsonLoop(int n, int m){
  double deltaK = 2*M_PI/(double)n;
  double * k = new double[dim];
  for(int i = 0; i < dim; i++){
    k[i] = k0[i];
  }
  k[dir] += deltaK;
  cx_mat u;

  if(ham->getIsSparse()){
    cx_vec eigVal;
    cx_mat eigVec0, eigVec1;
    eigs_gen(eigVal, eigVec0, ham->spH(k0), m, "sr");
    eigs_gen(eigVal, eigVec1, ham->spH(k), m, "sr");
    cx_mat eigVec2 = eigVec1;

    u = eigVec0.t() * eigVec1;

    for(int i = 2; i < n; i++){
      eigVec1 = eigVec2;
      k[dir] = k0[dir] + i*deltaK;
      eigs_gen(eigVal, eigVec2, ham->spH(k), m , "sr");
      u = u* (eigVec1.t() * eigVec2);
    }
    u = u * (eigVec2.t() * eigVec0);
  }
  else{
    cx_mat eigVec0 = evOcc(m, ham->H(k0));
    cx_mat eigVec1 = evOcc(m, ham->H(k));
    cx_mat eigVec2 = eigVec1;

    u = eigVec0.t() * eigVec1;

    for(int i = 2; i < n; i++){
      eigVec1 = eigVec2;
      k[dir] = k0[dir] + i*deltaK;
      eigVec2 = evOcc(m, ham->H(k));
      u = u* (eigVec1.t() * eigVec2);
    }
    u = u * (eigVec2.t() * eigVec0);
  }
  delete[] k;
  return u;
}

double Wilson::berryPhase(int n, int m){
  return - log_det(this->wilsonLoop(n,m)).imag();
}

cx_mat Wilson::wilsonLoopSupercell(int n, int m, double * k){
  double deltaTheta = 2*M_PI/(double)n;
  double * theta = new double[dim];
  double * theta0 = new double[dim];
  for(int i = 0; i < dim; i++){
    theta[i] = ham->getTwists()[i];
    theta0[i] = ham->getTwists()[i];
  }
  theta[dir] += deltaTheta;
  cx_mat u;

  if(ham->getIsSparse()){
    cx_vec eigVal;
    cx_mat eigVec0, eigVec1;
    eigs_gen(eigVal, eigVec0, ham->spH(k), m, "sr");
    ham->setTwists(theta);
    eigs_gen(eigVal, eigVec1, ham->spH(k), m, "sr");
    cx_mat eigVec2 = eigVec1;

    u = eigVec0.t() * eigVec1;

    for(int i = 2; i < n; i++){
      eigVec1 = eigVec2;
      theta[dir] = theta0[dir] + i*deltaTheta;
      ham->setTwists(theta);
      eigs_gen(eigVal, eigVec2, ham->spH(k), m , "sr");
      u = u* (eigVec1.t() * eigVec2);
    }
    u = u * (eigVec2.t() * eigVec0);
  }
  else{
    cx_mat eigVec0 = evOcc(m, ham->H(k));
    ham->setTwists(theta);
    cx_mat eigVec1 = evOcc(m, ham->H(k));
    cx_mat eigVec2 = eigVec1;

    u = eigVec0.t() * eigVec1;

    for(int i = 2; i < n; i++){
      eigVec1 = eigVec2;
      theta[dir] = theta0[dir] + i*deltaTheta;
      ham->setTwists(theta);
      eigVec2 = evOcc(m, ham->H(k));
      u = u* (eigVec1.t() * eigVec2);
    }
    u = u * (eigVec2.t() * eigVec0);
  }

  //cout << "ola2" << endl;
  delete[] theta;
  return u;
}

double Wilson::berryPhaseSupercell(int n, int m, double * k){
  return - log_det(this->wilsonLoopSupercell(n,m,k)).imag();
}
