#include "tbhop.h"

using namespace arma;

cx_mat ev_occ(int M, cx_mat H){
  vec eigval;
  cx_mat eigvec;
  int N = size(H)[0];
  eig_sym(eigval, eigvec, H);
  eigvec.resize(N,M);
  return eigvec;
}

cx_mat Kloop::wilson_loop(int N, int M){
  double deltak = 2*M_PI/(double)N;

  cx_mat eigvec0 = ev_occ(M, this->hloop(0));
  cx_mat eigvec1 = ev_occ(M, this->hloop(deltak));
  cx_mat eigvec2 = eigvec1;

  cx_mat U = eigvec0.t() * eigvec1;

  for(int i = 2; i < N; i++){
    eigvec1 = eigvec2;
    eigvec2 = ev_occ(M, this->hloop(i*deltak));
    U = U* (eigvec1.t() * eigvec2);
  }
  U = U * (eigvec2.t() * eigvec0);
  return U;
}

double Kloop::berry_phase(int N, int M){
  return - log_det(this->wilson_loop(N,M)).imag();
}
