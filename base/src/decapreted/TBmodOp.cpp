#include "TBmodOp.h"

using namespace std;
using namespace arma;

TBmodOp::TBmodOp(int ndim, int norb, int * L) : TBmod(ndim, norb, L){
  sparse = true;
  dir = 0;
  k0 = new double[ndim];
  for(int i = 0; i < ndim; i++){
    k0[i] = 0;
  }
}

TBmodOp::~TBmodOp(){
  delete[] k0;
}

void TBmodOp::set_sparse(bool val){
  sparse = val;
}

void TBmodOp::set_k0(double * k0){
  for(int i = 0; i < ndim; i++){
    this->k0[i] = k0[i];
  }
}

void TBmodOp::set_dir(int dir){
  this->dir = dir;
}

cx_mat TBmodOp::ev_occ(int M, cx_mat H){
  vec eigval;
  cx_mat eigvec;
  int N = size(H)[0];
  eig_sym(eigval, eigvec, H);
  eigvec.resize(N,M);
  eigval.resize(M);
  cout << eigval << eigvec << endl;
  return eigvec;
}

cx_mat TBmodOp::hloop(double k){
  double * kvec = new double[ndim];
  for(int i = 0; i < ndim; i++){
    kvec[i] = k0[i];
  }
  kvec[dir] += k;
  cx_mat res = get_H(kvec);
  delete[] kvec;
  return res;
}

sp_cx_mat TBmodOp::sp_hloop(double k){
  double * kvec = new double[ndim];
  for(int i = 0; i < ndim; i++){
    kvec[i] = k0[i];
  }
  kvec[dir] += k;
  sp_cx_mat res = get_spH(kvec);
  delete[] kvec;
  return res;
}

cx_mat TBmodOp::wilson_loop_k(int N, int M){
  double deltak = 2*M_PI/(double)N;

  if(sparse){
    cx_vec eigval;
    cx_mat eigvec0, eigvec1;
    cout << sp_hloop(0) << hloop(0) << "\n" << M << endl;
    eigs_gen(eigval, eigvec0, sp_hloop(0), M, "sr");
    cout << "eigval\n" << eigval << "eigvec\n" << eigvec0 << endl;
    eigs_gen(eigval, eigvec1, sp_hloop(deltak), M, "sr");
    cx_mat eigvec2 = eigvec1;

    cx_mat U = eigvec0.t() * eigvec1;

    for(int i = 2; i < N; i++){
      eigvec1 = eigvec2;
      eigs_gen(eigval, eigvec2, sp_hloop(i*deltak), M , "sr");
      cout << i << "\n" << "eigval\n" << eigval << "eigvec\n" << eigvec2 << endl;
      U = U* (eigvec1.t() * eigvec2);
    }
    U = U * (eigvec2.t() * eigvec0);
    return U;
  }
  else{
    cx_mat eigvec0 = ev_occ(M, this->hloop(0));
    cx_mat eigvec1 = ev_occ(M, this->hloop(deltak));
    cx_mat eigvec2 = eigvec1;

    cx_mat U = eigvec0.t() * eigvec1;

    for(int i = 2; i < N; i++){
      eigvec1 = eigvec2;
      cout << i << endl;
      eigvec2 = ev_occ(M, this->hloop(i*deltak));
      U = U* (eigvec1.t() * eigvec2);
    }
    U = U * (eigvec2.t() * eigvec0);
    return U;
  }
}

double TBmodOp::berry_phase_k(int N, int M){
  return - log_det(this->wilson_loop_k(N,M)).imag();
}
