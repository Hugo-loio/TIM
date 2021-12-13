#include "SSH_1D.h"
#include "tbhop.h"
#include <complex>

SSH_1D::SSH_1D(double t1, double t2){
  this->t1 = t1;
  this->t2 = t2;
  L = 0;
}

SSH_1D::~SSH_1D(){
}

void SSH_1D::set_intrahop(double t1){
  this->t1 = t1;
}

void SSH_1D::set_interhop(double t2){
  this->t2 = t2;
}

double SSH_1D::get_intrahop(){
  return t1;
}

double SSH_1D::get_interhop(){
  return t2;
}

cx_mat SSH_1D::kH(double k){
  cx_mat H(2, 2, fill::zeros);
  complex<double> i(0,1);
  H(0,1) = -t1-t2*exp(-i*k);
  H(1,0) = -t1-t2*exp(i*k);

  return H;
}

vec SSH_1D::eigenvalk(double k){
  cx_mat H = this->kH(k);
  vec eigval;
  eig_sym(eigval, H);
  return eigval;
}

cx_mat SSH_1D::eigenveck(double k){
  cx_mat H = this->kH(k);
  vec eigval;
  cx_mat eigvec;
  eig_sym(eigval,eigvec, H);
  return eigvec;
}

double SSH_1D::berry_kspace(int N){
  return abs(this->berry_phase(N, 1));
}

void SSH_1D::set_rH(int L, bool closed){
  rH = sp_mat(2*L,2*L);
  rH(1,0) = -t1;
  rH(0,1) = -t1;
  rH(2*L-1,2*L-2) = -t1;
  rH(2*L-2,2*L-1) = -t1;
  for(int i = 1; i < 2*L-1; i++){
    rH(i,i+1) = (i%2)*(-t2)+((i+1)%2)*(-t1);
    rH(i,i-1) = (i%2)*(-t1)+((i+1)%2)*(-t2);
  }
  if(closed){
    rH(0,2*L-1) = -t1;
    rH(2*L-1,0) = -t1;
  }
}


