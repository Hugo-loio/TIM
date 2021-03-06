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
    u0 = eigVec0;

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
    u0 = eigVec0;

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
  return sortedEigVec;
}

cx_mat Wilson::nestedWilsonLoop(int * n, int * dir, int m){
  double deltaK = 2*M_PI/(double)n[dir[1]];
  double * k = new double[dim];
  for(int i = 0; i < dim; i++){
    k[i] = k0[i];
  }
  cx_mat w;

  setLoopDir(dir[0]);

  if(ham->getIsSparse()){
  }
  else{
    cx_mat wannier0, wannier1;
    cx_mat nu = wilsonEigVec(n[dir[0]],m).cols(0,m/2-1);

    wannier0 = u0*nu;
    w0 = wannier0;
    k0[dir[1]] += deltaK;
    nu = wilsonEigVec(n[dir[0]],m).cols(0,m/2-1);
    wannier1 = u0*nu;
    cx_mat wannier2 = wannier1;

    w = wannier0.t() * wannier1;

    for(int i = 2; i < n[dir[1]]; i++){
      wannier1 = wannier2;
      k0[dir[1]] = k[dir[1]] + i*deltaK;
      nu = wilsonEigVec(n[dir[0]],m).cols(0,m/2-1);
      wannier2 = u0*nu;
      w = w* (wannier1.t() * wannier2);
    }
    w = w * (wannier2.t() * wannier0);
  }

  k0[dir[1]] = k[dir[1]];
  delete[] k;
  return w;
}

vec Wilson::nestedWilsonPhases(int * n, int * dir, int m){
  cx_vec eigVal;
  vec phase(m/2);
  eig_gen(eigVal, this->nestedWilsonLoop(n,dir,m));
  complex<double> ii(0,1);  
  for(int i = 0; i < m/2; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  //cout << k0[0] << " " << k0[1] << " " << k0[2] << " " << eigVal[0] << " " << eigVal[1] << endl;
  return sort(phase);
}

cx_mat Wilson::nestedWilsonEigVec(int * n, int * dir, int m){
  cx_vec eigVal;
  cx_mat eigVec;
  vec phase(m/2);
  eig_gen(eigVal, eigVec, this->nestedWilsonLoop(n, dir, m));
  complex<double> ii(0,1);
  for(int i = 0; i < m/2; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  uvec sort = sort_index(phase);
  cx_mat sortedEigVec(m/2,m/2);
  for(int i = 0; i < m/2; i++){
    sortedEigVec.col(i) = eigVec.col(sort[i]);
  }
  return sortedEigVec;
}

cx_mat Wilson::nestedNestedWilsonLoop(int * n, int * dir, int m){
  double deltaK = 2*M_PI/(double)n[dir[2]];
  double * k = new double[dim];
  for(int i = 0; i < dim; i++){
    k[i] = k0[i];
  }
  cx_mat u;

  setLoopDir(dir[0]);

  if(ham->getIsSparse()){
  }
  else{
    cx_mat wannier;
    cx_mat nestedW0, nestedW1;
    cx_mat eta = nestedWilsonEigVec(n,dir,m).cols(0,m/4-1);
    nestedW0 = w0*eta;
    k0[dir[2]] += deltaK;
    eta = nestedWilsonEigVec(n,dir,m).cols(0,m/4-1);
    nestedW1 = w0*eta;
    cx_mat nestedW2 = nestedW1;

    u = nestedW0.t() * nestedW1;

    for(int i = 2; i < n[dir[2]]; i++){
      nestedW1 = nestedW2;
      k0[dir[2]] = k[dir[2]] + i*deltaK;
      eta = nestedWilsonEigVec(n,dir,m).cols(0,m/4-1);
      nestedW2 = w0*eta;
      u = u* (nestedW1.t() * nestedW2);
    }
    u = u * (nestedW2.t() * nestedW0);

    //cout << "w: " << size(w0) << " nw: " << size(nestedW0) << " loop: " << size(wilsonEigVec(n[dir[0]], m)) << " nloop: " << size(nestedWilsonEigVec(n,dir,m)) << " nnloop: " << size(u) << endl;
  }

  k0[dir[2]] = k[dir[2]];
  delete[] k;
  return u;
}

//Supercell

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
    u0 = eigVec0;

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

  ham->setTwists(theta0);
  delete[] theta;
  delete[] theta0;
  return u;
}

double Wilson::berryPhaseSupercell(int n, int m, double * k){
  return - log_det(this->wilsonLoopSupercell(n,m,k)).imag();
}

vec Wilson::wilsonPhasesSupercell(int n, int m, double * k){
  cx_vec eigVal;
  vec phase(m);
  eig_gen(eigVal, this->wilsonLoopSupercell(n,m,k));
  complex<double> ii(0,1);  
  for(int i = 0; i < m; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  return sort(phase);
}

cx_mat Wilson::wilsonEigVecSupercell(int n, int m, double * k){
  cx_vec eigVal;
  cx_mat eigVec;
  vec phase(m);
  eig_gen(eigVal, eigVec, this->wilsonLoopSupercell(n,m,k));
  complex<double> ii(0,1);  
  for(int i = 0; i < m; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  uvec sort = sort_index(phase);
  cx_mat sortedEigVec(m,m);
  for(int i = 0; i < m; i++){
    sortedEigVec.col(i) = eigVec.col(sort[i]);
  }
  return sortedEigVec;
}


cx_mat Wilson::nestedWilsonLoopSupercell(int * n, int * dir, int * m, double * k){
  double deltaTheta = 2*M_PI/(double)n[dir[1]];
  double * theta = new double[dim];
  double * theta0 = new double[dim];
  for(int i = 0; i < dim; i++){
    theta[i] = ham->getTwists()[i];
    theta0[i] = ham->getTwists()[i];
  }
  theta[dir[1]] += deltaTheta;
  cx_mat u;

  setLoopDir(dir[0]);

  if(ham->getIsSparse()){
  }
  else{
    vec eigVal;
    cx_mat eigVec;
    cx_mat wannier0, wannier1;
    eig_sym(eigVal, eigVec, ham->H(k));
    wannier0 = (eigVec.cols(0,m[0]-1))*(wilsonEigVecSupercell(n[dir[0]],m[0],k).cols(0,m[1]-1));
    ham->setTwists(theta);
    eig_sym(eigVal, eigVec, ham->H(k));
    wannier1 = (eigVec.cols(0,m[0]-1))*(wilsonEigVecSupercell(n[dir[0]],m[0],k).cols(0,m[1]-1));
    cx_mat wannier2 = wannier1;

    u = wannier0.t() * wannier1;

    for(int i = 2; i < n[dir[1]]; i++){
      wannier1 = wannier2;
      theta[dir[1]] = theta0[dir[1]] + i*deltaTheta;
      ham->setTwists(theta);
      eig_sym(eigVal, eigVec, ham->H(k));
      wannier2 = (eigVec.cols(0,m[0]-1))*(wilsonEigVecSupercell(n[dir[0]],m[0],k).cols(0,m[1]-1));
      u = u* (wannier1.t() * wannier2);
    }
    u = u * (wannier2.t() * wannier0);
  }

  ham->setTwists(theta0);

  delete[] theta;
  delete[] theta0;
  return u;
}

cx_mat Wilson::nestedWilsonEigVecSupercell(int * n, int * dir, int * m, double * k){
  cx_vec eigVal;
  cx_mat eigVec;
  vec phase(m[1]);
  eig_gen(eigVal, eigVec, this->nestedWilsonLoopSupercell(n,dir,m,k));
  complex<double> ii(0,1);  
  for(int i = 0; i < m[1]; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  uvec sort = sort_index(phase);
  cx_mat sortedEigVec(m[1],m[1]);
  for(int i = 0; i < m[1]; i++){
    sortedEigVec.col(i) = eigVec.col(sort[i]);
  }
  return sortedEigVec;
}

vec Wilson::nestedWilsonPhasesSupercell(int * n, int * dir, int * m, double * k){
  cx_vec eigVal;
  vec phase(m[1]);
  eig_gen(eigVal, this->nestedWilsonLoopSupercell(n,dir,m,k));
  complex<double> ii(0,1);  
  for(int i = 0; i < m[1]; i++){
    phase[i] = (-ii*log(eigVal[i])).real()/(2*M_PI);
  }
  return sort(phase);
}

cx_mat Wilson::nestedNestedWilsonLoopSupercell(int * n, int * dir, int * m, double * k){
  double deltaTheta = 2*M_PI/(double)n[dir[2]];
  double * theta = new double[dim];
  double * theta0 = new double[dim];
  for(int i = 0; i < dim; i++){
    theta[i] = ham->getTwists()[i];
    theta0[i] = ham->getTwists()[i];
  }
  theta[dir[1]] += deltaTheta;
  cx_mat u;

  setLoopDir(dir[0]);

  if(ham->getIsSparse()){
  }
  else{
    vec eigVal;
    cx_mat eigVec;
    cx_mat wannier;
    cx_mat nestedW0, nestedW1;
    eig_sym(eigVal, eigVec, ham->H(k));
    wannier = (eigVec.cols(0,m[0]-1))*(wilsonEigVecSupercell(n[dir[0]],m[0],k).cols(0,m[1]-1));
    nestedW0 = (wannier.cols(0,m[1]-1))*(nestedWilsonEigVecSupercell(n,dir,m,k).cols(0,m[2]-1));
    ham->setTwists(theta);
    eig_sym(eigVal, eigVec, ham->H(k));
    wannier = (eigVec.cols(0,m[0]-1))*(wilsonEigVecSupercell(n[dir[0]],m[0],k).cols(0,m[1]-1));
    nestedW1 = (wannier.cols(0,m[1]-1))*(nestedWilsonEigVecSupercell(n,dir,m,k).cols(0,m[2]-1));
    cx_mat nestedW2 = nestedW1;

    u = nestedW0.t() * nestedW1;

    for(int i = 2; i < n[dir[2]]; i++){
      nestedW1 = nestedW2;
      theta[dir[2]] = theta0[dir[2]] + i*deltaTheta;
      ham->setTwists(theta);
      eig_sym(eigVal, eigVec, ham->H(k));
      wannier = (eigVec.cols(0,m[0]-1))*(wilsonEigVecSupercell(n[dir[0]],m[0],k).cols(0,m[1]-1));
      nestedW2 = (wannier.cols(0,m[1]-1))*(nestedWilsonEigVecSupercell(n,dir,m,k).cols(0,m[2]-1));
      u = u* (nestedW1.t() * nestedW2);
    }
    u = u * (nestedW2.t() * nestedW0);

    //cout << "w: " << size(wannier) << " nw: " << size(nestedW0) << " loop: " << size(wilsonEigVecSupercell(n[dir[0]], m[0],k)) << " nloop: " << size(nestedWilsonEigVecSupercell(n,dir,m,k)) << " nnloop: " << size(u) << endl;
  }

  ham->setTwists(theta0);

  delete[] theta;
  delete[] theta0;
  return u;
}
