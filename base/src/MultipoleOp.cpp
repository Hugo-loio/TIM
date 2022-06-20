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
  nOcc = lAccum[dim-1]*nOrb/2;
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

double MultipoleOp::chop(double d){
  double res = fmod(abs(d),1);
  double err = 1e-9;
  if(res > 1-err){
    res = 0;
  }
  return res;
}

double MultipoleOp::polarization(int a, double * k){
  int * point = new int[dim];
  for(int i = 0; i < dim; i++){
    point[i] = 0;
  }
  cx_mat psi;
  cx_mat psi_tilde;
  double p0 = ((double)1/2)*(nOcc/l[a])*(1+l[a]);
  if(!ham->getIsSparse()){
    vec eigVal;
    eig_sym(eigVal, psi, ham->H(k));
    psi = psi.cols(0,nOcc-1);
    psi_tilde = psi;
    int m;
    complex<double> ii(0,1);
    complex<double> u;
    for(int i = 0; i < lAccum[dim-1]; i++){
      u = exp(ii*(double)2*M_PI*(double)point[a]/(double)l[a]);
      m = nOrb*getN(point);
      psi_tilde.rows(m, m + nOrb -1) *= u;
      nextPoint(0, point, false);
    }
  }
  return chop((l[a]/(2*M_PI*lAccum[dim-1]))*log_det(psi.t()*psi_tilde).imag() - p0);
}

double MultipoleOp::quadrupole(int a, int b, double * k){
  int * point = new int[dim];
  for(int i = 0; i < dim; i++){
    point[i] = 0;
  }
  cx_mat psi;
  cx_mat psi_tilde;
  double p0 = ((double)1/4)*(nOcc/(l[a]*l[b]))*(1+l[a])*(1+l[b]);
  if(!ham->getIsSparse()){
    vec eigVal;
    eig_sym(eigVal, psi, ham->H(k));
    psi = psi.cols(0,nOcc-1);
    psi_tilde = psi;
    int m;
    complex<double> ii(0,1);
    double phase;
    complex<double> u;
    for(int i = 0; i < lAccum[dim-1]; i++){
      phase = (double)(point[a]+1)*(point[b]+1)/(double)(l[a]*l[b]);
      u = exp(ii*(double)2*M_PI*phase);
      m = nOrb*getN(point);
      //cout << m << " " << phase << endl;
      psi_tilde.rows(m, m + nOrb -1) *= u;
      nextPoint(0, point, false);
    }
  }
  return chop((l[a]*l[b]/(2*M_PI*lAccum[dim-1]))*log_det(psi.t()*psi_tilde).imag() - p0);
}

double MultipoleOp::octupole(int a, int b, int c, double * k){
  int * point = new int[dim];
  for(int i = 0; i < dim; i++){
    point[i] = 0;
  }
  cx_mat psi;
  cx_mat psi_tilde;
  double p0 = ((double)1/8)*(nOcc/(l[a]*l[b]*l[c]))*(1+l[a])*(1+l[b])*(1+l[c]);
  if(!ham->getIsSparse()){
    vec eigVal;
    eig_sym(eigVal, psi, ham->H(k));
    psi = psi.cols(0,nOcc-1);
    psi_tilde = psi;
    int m;
    complex<double> ii(0,1);
    double phase;
    complex<double> u;
    for(int i = 0; i < lAccum[dim-1]; i++){
      phase = (double)(point[a]*point[b]*point[c])/(double)(l[a]*l[b]*l[c]);
      u = exp(ii*(double)2*M_PI*phase);
      m = nOrb*getN(point);
      psi_tilde.rows(m, m + nOrb -1) *= u;
      nextPoint(0, point, false);
    }
  }
  //cout << "p0 : " << p0 << endl;
  return chop((l[a]*l[b]*l[c]/(2*M_PI*lAccum[dim-1]))*log_det(psi.t()*psi_tilde).imag() - p0);
}
