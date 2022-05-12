#include "MultipoleOp.h"

using namespace std;
using namespace arma;

MultipoleOp::MultipoleOp(Hamiltonian * ham, int * l,int dim, int nOrb) : ham(ham), dim(dim), nOrb(nOrb){
  lAccum = new int[dim];
  this->l = new int[dim];
  for(int i = 0; i < dim; i++){
    this->l[i] = l[i];
  }
  lAccum[0] = l[0];
  for(int i = 1; i < dim; i++){
    lAccum[i] = lAccum[i-1]*l[i];
  }
}

MultipoleOp::~MultipoleOp(){
  delete[] lAccum;
}

int MultipoleOp::getN(int * n){
  int res = 1;
  for(int i = 0; i < dim; i++){
    if(i == 0){
      res = n[0];
    }
    else{
      res += lAccum[i-1]*n[i];
    }
  }
  return res;
}

void MultipoleOp::nextPoint(int depth, int * point, bool up){
  if(depth == dim-1){
    //increase coordinate
    if(point[dim-1] == l[dim -1] - 1){
      point[dim -1] = 0; 
      nextPoint(--depth, point, true);
    }
    else{
      point[dim-1]++;
    }
  }
  else if(depth != -1){
    if(up){
      //increase coordinate
      if(point[depth] == l[depth] -1 ){
	point[depth] = 0;
	nextPoint(--depth, point, true);
      }
      else{
	point[depth]++;
      }
    }
    else{
      //Go to next direction
      nextPoint(++depth, point, false);
    }
  }
}

double MultipoleOp::polarization(int a, double * k){
  int * point = new int[dim];
  for(int i = 0; i < dim; i++){
    point[i] = 0;
  }
  cx_mat psi;
  cx_mat psi_tilde;
  if(!ham->getIsSparse()){
    vec eigVal;
    eig_sym(eigVal, psi, ham->H(k));
    psi = psi.cols(0,nOcc-1);
    psi_tilde = psi;
    int m;
    for(int i = 0; i < lAccum[dim-1]; i++){
      complex<double> ii(0,1);
      complex<double> u = exp(ii*(double)2*M_PI*(double)point[a]/(double)l[a]);
      m = nOrb*getN(point);
      psi_tilde.rows(m, m + nOrb -1) *= u;
      nextPoint(0, point, false);
    }
  }
  return fmod((l[a]/(2*M_PI*lAccum[dim-1]))*log_det(psi.t()*psi_tilde).imag(),1);
}
