#include "LocalizationStats.h"

LocalizationStats::LocalizationStats(Hamiltonian * ham){
  this->ham = ham;
}

LocalizationStats::~LocalizationStats(){
}

double LocalizationStats::ipr(int vol, int nOrb, int nStates, double * k){
  if(ham->getIsSparse()){
    cx_vec eigVal;
    cx_mat eigVec;
    eigs_gen(eigVal, eigVec, ham->spH(k), nStates, 0.0);
    double res = 0;
    double amp = 0;
    if(size(eigVal)[0] != nStates){
      throw runtime_error("Diagonalization failed.");
    }
    for(int i = 0; i < nStates; i++){
      //cout << eigVal[i] << endl;
      eigVec.col(i) = normalise(eigVec.col(i));
      for(int e = 0; e < vol; e++){
	amp = 0;
	for(int j = 0; j < nOrb; j++){
	  amp += (eigVec(e*nOrb + j,i)*conj(eigVec(e*nOrb + j,i))).real();
	}
	res += amp*amp;
      }
    }
    return res/(double)nStates;
  }
  else{
    cout << __PRETTY_FUNCTION__ << " hasn't been implemented for dense Hamiltonians" << endl;
    return -1;
  }
}

double LocalizationStats::tmm(int nLayers, int qrIt, double en, double * k){
  int size = ham->blockH(0,0,k).n_cols;
  int tSize = 2*size;
  double * c = new double[tSize];
  double * d = new double[tSize];
  for(int i = 0; i < tSize; i++){
    c[i] = 0;
    d[i] = 0;
  }
  cx_mat tAux(tSize, tSize, fill::zeros);
  tAux.submat(size, 0, tSize -1, size -1) = eye<cx_mat>(size, size);
  cx_mat enMat = en*eye<cx_mat>(size,size);
  cx_mat v = ham->blockH(0,1,k).i();
  tAux.submat(0, 0, size -1, size -1) = v*(enMat - ham->blockH(0,0,k));
  tAux.submat(0, size, size -1, tSize -1) = -v;
  cx_mat t = tAux;
  cx_mat q,r;

  int i,e,minIndex;
  double rTemp;

  for(i = 1; i < maxItTMM; i++){
    if(i % qrIt == 0 && i/qrIt > 1){
      qr(q, r, tAux);
      tAux = q;
      for(e = 0; e < tSize; e++){
	//TODO: do this if r is complex?
	rTemp = log(abs(r(e,e)));
	c[e] += rTemp;
	d[e] += rTemp*rTemp;
      }
      if(testTmmConv(c, d, tSize, i/qrIt, minIndex)){
	break;
      }
    }
    e = i % nLayers;
    if(e != 0){
      ham->generateDisorder();
      //cout << "aqui 3" << endl;
      //cout << ham->blockH(e, e+1, k) << endl;
      v = ham->blockH(e,e+1,k).i();
      //cout << "aqui 4" << endl;
      t.submat(0, 0, size -1, size -1) = v*(enMat - ham->blockH(e,e,k));
      t.submat(0, size, size -1, tSize -1) = -v*ham->blockH(e,e-1,k);
      tAux = t*tAux;
    }
    else{
      t.submat(0, 0, size -1, size -1) = enMat - ham->blockH(nLayers-1,nLayers-1,k);
      t.submat(0, size, size -1, tSize -1) = -ham->blockH(nLayers-1,nLayers-2,k);
      ham->generateDisorder();
      //cout << "aqui 1" << endl;
      v = ham->blockH(0,1,k).i();
      //cout << "aqui 2" << endl;
      t.submat(0, 0, size -1, size -1) = v*t.submat(0, 0, size -1, size -1);
      t.submat(0, size, size -1, tSize -1) = v*t.submat(0, size, size -1, tSize -1);
      tAux = t*tAux;
    }
  }
  //cout << "Iterations: " << i << endl;

  cout << "iterations: " << i << endl;
  if(i == maxItTMM){
    cout << "Warning: transfer matrix method stopped at iteration limit" << endl;
  }

  for(e = 0; e < tSize; e++){
    cout << abs((double)i/c[e]) << endl;
  }

  double res = abs((double)i/c[minIndex]);
  delete[] c;
  delete[] d;
  return res;
}

bool LocalizationStats::testTmmConv(double * c, double * d, int size, int nQR, int & minIndex){
  minIndex = 0; 
  double minAbs = abs(c[0]);
  double thisAbs = 0;
  for(int e = 1; e < size; e++){
    thisAbs = abs(c[e]);
    if(thisAbs < minAbs){
      minIndex = e;
      minAbs = thisAbs;
    }
  }
  double err = abs(sqrt(d[minIndex]/(double)nQR - (c[minIndex]/(double)nQR)*(c[minIndex]/(double)nQR))*sqrt((double)nQR)/c[minIndex]);
  /*
     if(nQR % 1000 == 0){
     cout << "nQR: " << nQR << " err " << err << " d " << d[minIndex] << " c " << c[minIndex] <<  endl;
     cout << d[minIndex]/(double)nQR << " " << (c[minIndex]/(double)nQR)*(c[minIndex]/(double)nQR) << endl;
     }
     */
  if(err < tmmErr){
    return true;
  }
  return false;
}

double LocalizationStats::lsr(int nStates, double * k){
  if(ham->getIsSparse()){
    cx_vec eigVal;
    eigs_gen(eigVal, ham->spH(k), nStates, 0.0);
    vec realEig = vec(nStates, fill::zeros);
    for(int i = 0; i < nStates; i++){
      realEig[i] = eigVal[i].real();
    }
    realEig = sort(realEig);
    double res = 0;
    double max, min;
    if(size(eigVal)[0] != nStates){
      throw runtime_error("Diagonalization failed.");
    }
    double * s = new double[nStates-2];
    for(int i = 0; i < nStates - 2; i++){
      if(realEig[i+1] <= 0){
	s[i] = realEig[i+1] - realEig[i];
      }
      else{
	s[i] = realEig[i+2] - realEig[i+1];
      }
    }
    for(int i = 0; i < nStates - 3; i++){
      if(s[i] < s[i+1]){
	min = s[i];
	max = s[i+1];
      }
      else{
	min = s[i+1];
	max = s[i];
      }
      res += min/max;
    }
    delete s;
    return res/((double)nStates-3);
  }
  else{
    cout << __PRETTY_FUNCTION__ << " hasn't been implemented for dense Hamiltonians" << endl;
    return -1;
  }
}
