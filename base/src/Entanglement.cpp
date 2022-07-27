#include "Entanglement.h"

Entanglement::Entanglement(Hamiltonian * ham, int nOcc, double * k) : ham(ham){
  this->nOcc = nOcc;
  diagonalizeH(k);
}

Entanglement::~Entanglement(){
}

void Entanglement::diagonalizeH(double * k){
  vec eigVal;
  eig_sym(eigVal, eigStates, ham->H(k)); 
  eigStates = eigStates.cols(0, nOcc - 1);
  for(int i = 0; i < nOcc; i++){
    eigStates.col(i) = normalise(eigStates.col(i));
  }
}

void Entanglement::setOcc(int nOcc){
  this->nOcc = nOcc;
  diagonalizeH();
}

double Entanglement::bipEntropy(uvec cut){
  cx_mat cutEigStates = eigStates.rows(cut);
  cx_mat corr = conj(cutEigStates)*(cutEigStates.st());
  vec zeta;
  eig_sym(zeta, corr);
  double z;
  double s = 0;
  double err = 1e-30;
  for(int i = 0; i < size(zeta)[0]; i++){
    z = chop(zeta[i]);
    s+= -z*log(z + err) - (1-z)*log(1-z + err);
    //cout << "z: " << z << " s: " << s << endl;
  }
  //cout << "s final: " << s << endl;
  return s;
}

double Entanglement::chop(double d){
  if(d < 0){
    return 0;
  }
  else if(d > 1){
    return 1;
  }
  return d;
}

