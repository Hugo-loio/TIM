#include "SSH_1D.h"
#include <complex>

SSH_1D::SSH_1D(double t1, double t2){
  this->t1 = t1;
  this->t2 = t2;
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
  cx_vec u0 = eigenveck(0).col(0);
  cx_vec u1 = u0;
  cx_vec u2 = eigenveck(2*M_PI/N).col(0);
  complex<double> res = cdot(u1,u2);
  for(int i = 2; i < N; i++){
    u1 = u2;
    u2 = eigenveck(2*M_PI*(double)i/N).col(0);
    res *= cdot(u1,u2);
  }
  res *= cdot(u2,u0);
  return abs(log(res).imag());
}
