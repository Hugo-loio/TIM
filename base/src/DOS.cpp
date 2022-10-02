#include "DOS.h"
#include "AuxFunctions.h"

DOS::DOS(Hamiltonian * ham){
  this->ham = ham;
}

DOS::~DOS(){
  if(momentsFound){
    delete[] mu;
  }
}

double DOS::kpm(double en, int nMoments, int nRandVecs, double * k){
  if(ham->getIsSparse()){
    if((k != NULL || !rescalingFound) && !customERange){
      findRescaling(k);
      rescalingFound = true;
    }
    if(momentsFound){
      if(this->nMoments != nMoments || this->nRandVecs != nRandVecs || k != NULL){
	momentsFound = false;
      }
      if(!momentsFound){
	delete[] mu;
      }
    }
    if(!momentsFound){
      mu = new double[nMoments];
      for(int i = 0; i < nMoments; i++){
	mu[i] = 0;
      }
      this->nMoments = nMoments;
      this->nRandVecs = nRandVecs;
      if(ham->getIsReal()){
	calculateMoments<sp_mat>(real(ham->spH(k)));
      }
      else{
	calculateMoments<sp_cx_mat>(ham->spH(k));
      }
      momentsFound = true;
    }
    double res = jacksonKernel(0)*mu[0];
    //cout << res << endl;
    en = (en - b)/a;
    for(int i = 1; i < nMoments; i++){
      res += 2*jacksonKernel(i)*mu[i]*chebyshev(i, en);
      //cout << i << " " << res << " " << jacksonKernel(i) << " " << mu[i] << " " << chebyshev(i, en) << endl;
    }
    return (1/a)*res/(M_PI*sqrt(1 - en*en));
  } 
  else{
    cout << __PRETTY_FUNCTION__ << " hasn't been implemented for dense Hamiltonians" << endl;
    return 0;
  }
}

void DOS::findRescaling(double * k){
  cx_vec eigVal;
  eigs_gen(eigVal, ham->spH(k), 1, "lr");
  double eMax = eigVal[0].real();
  eigs_gen(eigVal, ham->spH(k), 1, "sr");
  double eMin = eigVal[0].real();

  a = (eMax - eMin)/(2 - err);
  b = (eMax + eMin)/2;
}

double DOS::jacksonKernel(int n){
  double aux = nMoments + 1;
  return ((aux - n)*cos((M_PI*n)/aux) + sin((M_PI*n)/aux)*cot(M_PI/aux))/aux;
}

double DOS::chebyshev(int n, double x){
  return cos(n * acos(x));
}

template <class mat> void DOS::calculateMoments(mat h){
  int d = size(h)[0];
  h = (h - b*speye<mat>(size(h)))/a;
  cx_vec rand(d);
  int i,e;
  int n = nMoments/2;
  mu[0] = 1;
  for(i = 0; i < nRandVecs; i++){
    randomize(rand, d);
    cx_vec a1 = h*rand;
    cx_vec a2 = rand;
    cx_vec a3;
    double mu1 = cdot(rand, a1).real();
    mu[1] += mu1;
    for(e = 1; e < n; e++){
      a3 = 2*h*a1-a2;
      a2 = a1;
      a1 = a3;
      mu[2*e] += 2*cdot(a2,a2).real() - 1;
      mu[2*e + 1] += 2*cdot(a1,a2).real() - mu1;
    }
    if(nMoments % 2 == 1){
      mu[nMoments - 1] += 2*cdot(a1,a1).real() - 1;
    } 
  }
  for(i = 1; i < nMoments; i++){
    mu[i] *= (1/(double)nRandVecs);
  }
}

void DOS::randomize(cx_vec & rand, int d){
  random_device dev;
  mt19937 generator(dev());
  uniform_real_distribution<double> uni(0, 2*M_PI);
  double norm = 1/sqrt((double)d);
  complex<double> ii(0,1);
  for(int i = 0; i < d; i++){
    rand[i] = norm*exp(ii*uni(generator));
  }
}

void DOS::setKpmERange(double eMin, double eMax){
  customERange = true;
  a = (eMax - eMin)/(2 - err);
  b = (eMax + eMin)/2;
}
