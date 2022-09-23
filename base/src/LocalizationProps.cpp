#include "LocalizationProps.h"

LocalizationProps::LocalizationProps(Hamiltonian * ham){
  this->ham = ham;
}

LocalizationProps::~LocalizationProps(){
}

double LocalizationProps::ipr(int vol, int nOrb, int nStates, double en, double * k){
  if(ham->getIsSparse()){
    sparseDiag(nStates, en, k);
    double res = 0;
    double amp = 0;
    for(int i = 0; i < nStates; i++){
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

double LocalizationProps::lsr(int nStates, double en, double * k){
  if(ham->getIsSparse()){
    sparseDiag(nStates, en, k);
    vec realEig = sort(eigVal.subvec(0,nStates-1));
    double res = 0;
    double max, min;
    double * s = new double[nStates-2];
    for(int i = 0; i < nStates - 2; i++){
      if(realEig[i+1] <= en){
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
      if(max == 0){
	res += 1;
      }
      else{
	res += min/max;
      } 
    }
    delete s;
    return res/((double)nStates-3);
  }
  else{
    cout << __PRETTY_FUNCTION__ << " hasn't been implemented for dense Hamiltonians" << endl;
    return -1;
  }
}

double LocalizationProps::gap(double en, double * k){
  if(ham->getIsSparse()){
    vec realEig;
    int nStates = 10;
    if(nStatesFound > nStates){
      nStates = nStatesFound;
    }
    double res;
    bool extraStates = false;

    while(1){
      sparseDiag(nStates, en, k);
      realEig = sort(eigVal);
      for(int i = 0; i < nStates; i++){
	if(realEig[i] == en){
	  return 0;
	}
      }
      for(int i = 0; i < nStates - 1; i++){
	if(realEig[i+1] > en && realEig[i] < en){
	  if(extraStates){
	    cout << __PRETTY_FUNCTION__ << " needed " << nStates << " states" << endl;
	  }
	  return realEig[i+1] - realEig[i];
	}
      }
      nStates*= 2;
      extraStates = true;
    }
    return -1;
  }
  else{
    cout << __PRETTY_FUNCTION__ << " hasn't been implemented for dense Hamiltonians" << endl;
    return -1;
  }
}

void LocalizationProps::sparseDiag(int nStates, double en, double * k){
  if(forceDiag == false){
    if(k == NULL && nStates <= nStatesFound && lastEn == en){
      return;
    }
  }
  else{
    forceDiag = false;
  }
  cx_vec eigValTemp;
  eigs_gen(eigValTemp, eigVec, ham->spH(k), nStates, en);
  if(size(eigValTemp)[0] != nStates){
    cout << "Found " << size(eigVal)[0] << " states out of " << nStates << endl;
    throw runtime_error("Diagonalization failed.");
  }
  eigVal = vec(nStates, fill::zeros);
  for(int i = 0; i < nStates; i++){
    eigVal[i] = eigValTemp[i].real();
  }
  nStatesFound = nStates;
  lastEn = en;
}

double LocalizationProps::tmm(int nLayers, int qrIt, double en, double * k){
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
  double rTemp, err;

  for(i = 1; i < maxItTMM; i++){
    if(i % 100000 == 0){
      cout << __PRETTY_FUNCTION__ << " in iteration " << i << endl;
    }
    if(i % qrIt == 0){
      qr(q, r, tAux);
      tAux = q;
      for(e = 0; e < tSize; e++){
	//TODO: do this if r is complex?
	rTemp = log(abs(r(e,e)));
	c[e] += rTemp;
	d[e] += rTemp*rTemp;
      }
      if(i/qrIt > 1){
	if(testTmmConv(c, d, tSize, i/qrIt, minIndex, err)){
	  break;
	}
      }
    }
    e = i % (nLayers - 1);
    if(e != 0){
      ham->generateDisorder();
      v = ham->blockH(e,e+1,k).i();
      t.submat(0, 0, size -1, size -1) = v*(enMat - ham->blockH(e,e,k));
      t.submat(0, size, size -1, tSize -1) = -v*ham->blockH(e,e-1,k);
      tAux = t*tAux;
    }
    else{
      t.submat(0, 0, size -1, size -1) = enMat - ham->blockH(nLayers-1,nLayers-1,k);
      t.submat(0, size, size -1, tSize -1) = -ham->blockH(nLayers-1,nLayers-2,k);
      ham->generateDisorder();
      v = ham->blockH(0,1,k).i();
      t.submat(0, 0, size -1, size -1) = v*t.submat(0, 0, size -1, size -1);
      t.submat(0, size, size -1, tSize -1) = v*t.submat(0, size, size -1, tSize -1);
      tAux = t*tAux;
    }
  }


  if(i == maxItTMM){
    cout << "Warning: transfer matrix method stopped at iteration limit " << maxItTMM << " with error " << err << endl;
  }
  else{
    cout << "Iterations: " << i << endl;
  }

  for(e = 0; e < tSize; e++){
    //cout << abs((double)i/c[e]) << endl;
  }

  double res = abs((double)i/c[minIndex]);
  delete[] c;
  delete[] d;
  return res;
}

bool LocalizationProps::testTmmConv(double * c, double * d, int size, int nQR, int & minIndex, double & err){
  minIndex = 0; 
  double minAbs = abs(c[0]);
  double thisAbs;
  for(int e = 1; e < size; e++){
    thisAbs = abs(c[e]);
    if(thisAbs < minAbs){
      minIndex = e;
      minAbs = thisAbs;
    }
  }
  err = abs(sqrt(d[minIndex]/(double)nQR - (c[minIndex]/(double)nQR)*(c[minIndex]/(double)nQR))*sqrt((double)nQR)/c[minIndex]);
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
