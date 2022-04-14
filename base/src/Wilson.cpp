#include "Wilson.h"
#include <iostream>
#include <complex>

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
    vec eigVal;
    cx_mat eigVec0, eigVec1;
    eig_sym(eigVal, eigVec0, ham->H(k0));
    eigVec0 = eigVec0.cols(0,m-1);
    eig_sym(eigVal, eigVec1, ham->H(k));
    eigVec1 = eigVec1.cols(0,m-1);
    cx_mat eigVec2 = eigVec1;

    u = eigVec0.t() * eigVec1;

    for(int i = 2; i < n; i++){
      eigVec1 = eigVec2;
      k[dir] = k0[dir] + i*deltaK;
      eig_sym(eigVal, eigVec2, ham->H(k));
      eigVec2 = eigVec2.cols(0,m-1);
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

vec Wilson::wilsonPhases(int n, int m){
  cx_vec eigVal;
  vec phase(m);
  eig_gen(eigVal, this->wilsonLoop(n,m));
  complex<double> ii(0,1);  
  for(int i = 0; i < m; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  return sort(phase);
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
    vec eigVal;
    cx_mat eigVec0, eigVec1;
    eig_sym(eigVal, eigVec0, ham->H(k));
    eigVec0 = eigVec0.cols(0,m-1);
    ham->setTwists(theta);
    eig_sym(eigVal, eigVec1, ham->H(k));
    eigVec1 = eigVec1.cols(0,m-1);
    cx_mat eigVec2 = eigVec1;

    u = eigVec0.t() * eigVec1;

    for(int i = 2; i < n; i++){
      eigVec1 = eigVec2;
      theta[dir] = theta0[dir] + i*deltaTheta;
      ham->setTwists(theta);
      eig_sym(eigVal, eigVec2, ham->H(k));
      eigVec2 = eigVec2.cols(0,m-1);
      u = u* (eigVec1.t() * eigVec2);
    }
    u = u * (eigVec2.t() * eigVec0);
  }

  delete[] theta;
  return u;
}

double Wilson::berryPhaseSupercell(int n, int m, double * k){
  return - log_det(this->wilsonLoopSupercell(n,m,k)).imag();
}

vec Wilson::supercellWilsonPhases(int n, int m, double * k){
  cx_vec eigVal;
  vec phase(m);
  eig_gen(eigVal, this->wilsonLoopSupercell(n,m,k));
  complex<double> ii(0,1);  
  for(int i = 0; i < m; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  return sort(phase);
}

cx_mat Wilson::wilsonEigVec(int n, int m){
  cx_vec eigVal;
  cx_mat eigVec;
  vec phase(m);
  eig_gen(eigVal, eigVec, this->wilsonLoop(n,m));
  complex<double> ii(0,1);  
  for(int i = 0; i < m; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  uvec sort = sort_index(phase);
  cx_mat sortedEigVec(m,m);
  for(int i = 0; i < m; i++){
    sortedEigVec.col(i) = eigVec.col(sort[i]);
  }
  //cout << "k: " << k0[0] << " " << k0[1] << "\n " << sortedEigVec << "\n" << eigVec << "\n " << sort << "\n" << phase << endl;
  return sortedEigVec;
}

cx_mat Wilson::nestedWilsonLoop(int * n, int * dir, int m){
  double deltaK = 2*M_PI/(double)n[dir[1]];
  double * k = new double[dim];
  for(int i = 0; i < dim; i++){
    k[i] = k0[i];
  }
  cx_mat u;

  setLoopDir(dir[0]);

  if(ham->getIsSparse()){
  }
  else{
    vec eigVal;
    cx_mat eigVec;
    cx_mat wannier0, wannier1;
    eig_sym(eigVal, eigVec, ham->H(k0));
    wannier0 = (eigVec.cols(0,m-1))*(wilsonEigVec(n[dir[0]],m).cols(0,m/2-1));
    k0[dir[1]] += deltaK;
    eig_sym(eigVal, eigVec, ham->H(k0));
    wannier1 = (eigVec.cols(0,m-1))*(wilsonEigVec(n[dir[0]],m).cols(0,m/2-1));
    cx_mat wannier2 = wannier1;

    u = wannier0.t() * wannier1;

    for(int i = 2; i < n[dir[1]]; i++){
      wannier1 = wannier2;
      k0[dir[1]] = k[dir[1]] + i*deltaK;
      eig_sym(eigVal, eigVec, ham->H(k0));
      wannier2 = (eigVec.cols(0,m-1))*(wilsonEigVec(n[dir[0]],m).cols(0,m/2-1));
      u = u* (wannier1.t() * wannier2);
    }
    u = u * (wannier2.t() * wannier0);
  }

  k0[dir[1]] = k[dir[1]];
  delete[] k;
  return u;
}

vec Wilson::nestedWilsonPhases(int * n, int * dir, int m){
  cx_vec eigVal;
  vec phase(m/2);
  eig_gen(eigVal, this->nestedWilsonLoop(n,dir,m));
  complex<double> ii(0,1);  
  for(int i = 0; i < m/2; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  return sort(phase);
}

