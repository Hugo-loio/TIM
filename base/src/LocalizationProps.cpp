#include "LocalizationProps.h"
#include "AuxFunctions.h"

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
    //cout << nStates << realEig << endl;
    double res = 0;
    double max, min;
    double * s = new double[nStates-2];
    int countS = 0;
    for(int i = 0; i < nStates - 2; i++){
      if(realEig[i+1] <= en){
	s[i] = realEig[i+1] - realEig[i];
	countS++;
      }
    }
    for(int i = 0; i < countS - 1; i++){
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
    return res/((double)countS-1);
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
  if(ham->getIsReal()){
    mat eigVecTemp;
    eigs_sym(eigVal, eigVecTemp, real(ham->spH(k)), nStates, en);
    if(size(eigVal)[0] != nStates){
      cout << "Found " << size(eigVal)[0] << " states out of " << nStates << endl;
      throw runtime_error("Diagonalization failed.");
    }
    eigVec = cx_mat(eigVecTemp, mat(size(eigVecTemp)[0], size(eigVecTemp)[1], fill::zeros));
  }
  else{
    cx_vec eigValTemp;
    eigs_gen(eigValTemp, eigVec, ham->spH(k), nStates, en);
    if(size(eigValTemp)[0] != nStates){
      cout << "Found " << size(eigValTemp)[0] << " states out of " << nStates << endl;
      throw runtime_error("Diagonalization failed.");
    }
    eigVal = vec(nStates, fill::zeros);
    for(int i = 0; i < nStates; i++){
      eigVal[i] = eigValTemp[i].real();
    }
  }
  //Sort eigVals and eigVecs
  uvec sortIndex(nStates);
  vector<int> unSorted;
  for(int i = 0; i < nStates; i++){
    unSorted.push_back(i);
  }
  double absMin;
  int index;
  for(int i = 0; i < nStates; i++){
    absMin = abs(eigVal[unSorted[0]]);
    index = 0;
    for(int e = 1; e < unSorted.size(); e++){
      if(abs(eigVal[unSorted[e]]) < absMin){
	absMin = abs(eigVal[unSorted[e]]);
	index = e;
      }
    }
    sortIndex(i) = unSorted[index];
    unSorted.erase(unSorted.begin() + index);
  }
  cx_mat eigVecTemp = eigVec;
  vec eigValTemp = eigVal;
  for(int i = 0; i < nStates; i++){
    eigVal(i) = eigValTemp(sortIndex(i));
    eigVec.col(i) = eigVecTemp.col(sortIndex(i));
  }
  nStatesFound = nStates;
  lastEn = en;
}

vector<double> LocalizationProps::tmm(int nLayers, int qrIt, double en, double * k){
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
  int printIt = maxItTMM/10;

  for(i = 1; i < maxItTMM; i++){
    if(i % printIt == 0){
      cout << __PRETTY_FUNCTION__ << " in iteration " << i << endl;
    }
    if(i % qrIt == 0){
      qr(q, r, tAux);
      tAux = q;
      for(e = 0; e < tSize; e++){
	rTemp = log(abs(r(e,e)));
	c[e] += rTemp;
	d[e] += rTemp*rTemp;
      }
      if(i/qrIt > 10){
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

  vector<double> res;
  res.push_back(abs((double)i/c[minIndex]));
  res.push_back(err);
  delete[] c;
  delete[] d;
  return res;
}

vector<double> LocalizationProps::tmmSpecial(int nLayers, int qrIt, double en, double * k){
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
  cx_mat v = ham->blockH(0,1,k); 
  diagInverse<cx_mat>(v, size);
  tAux.submat(0, 0, size -1, size -1) = diagLeftMult<cx_mat>(v, enMat - ham->blockH(0,0,k), size);
  tAux.submat(0, size, size -1, tSize -1) = -v;
  cx_mat t = tAux;
  cx_mat q,r;

  int i,e,minIndex;
  double rTemp, err;
  double deltaErr = tmmErr;
  int printIt = maxItTMM/10;

  for(i = 1; i < maxItTMM; i++){
    if(i % printIt == 0){
      tmmErr += deltaErr;
      cout << __PRETTY_FUNCTION__ << " in iteration " << i << endl;
    }
    if(i % qrIt == 0){
      qr(q, r, tAux);
      tAux = q;
      for(e = 0; e < tSize; e++){
	rTemp = log(abs(r(e,e)));
	c[e] += rTemp;
	d[e] += rTemp*rTemp;
      }
      if(i/qrIt > 10){
	if(testTmmConv(c, d, tSize, i/qrIt, minIndex, err)){
	  break;
	}
      }
    }
    e = i % (nLayers - 1);
    if(e != 0){
      if(e % 2 == 1){
	ham->generateDisorder();
	t.submat(0, 0, size -1, size -1) = enMat - ham->blockH(e,e,k);
	t.submat(0, size, size -1, tSize -1) = -ham->blockH(e,e-1,k);
	tAux = updateT<cx_mat>(t, tAux, size);
      }
      else{
	ham->generateDisorder();
	cx_mat v = ham->blockH(e,e+1,k);
	diagInverse(v, size);
	t.submat(0, 0, size -1, size -1) = diagLeftMult<cx_mat>(v, enMat - ham->blockH(e,e,k), size);
	t.submat(0, size, size -1, tSize -1) = -v;
	tAux = updateT<cx_mat>(t, tAux, size);
      }
    }
    else{
      t.submat(0, 0, size -1, size -1) = enMat - ham->blockH(nLayers-1,nLayers-1,k);
      ham->generateDisorder();
      v = ham->blockH(0,1,k); 
      diagInverse<cx_mat>(v, size);
      t.submat(0, 0, size -1, size -1) = diagLeftMult<cx_mat>(v, t.submat(0, 0, size-1, size-1), size);
      t.submat(0, size, size -1, tSize -1) = -v;
      tAux = updateT<cx_mat>(t, tAux, size);
    }
  }


  if(i == maxItTMM){
    cout << "Warning: transfer matrix method stopped at iteration limit " << maxItTMM << " with error " << err << endl;
  }
  else{
    cout << "Iterations: " << i << endl;
  }

  vector<double> res;
  res.push_back(abs((double)i/c[minIndex]));
  res.push_back(err);
  delete[] c;
  delete[] d;
  return res;
}

vector<double> LocalizationProps::tmmSpecial2(int nLayers, int qrIt, double en, double * k){
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
  cx_mat v = ham->blockH(0,1,k); 
  diagInverse<cx_mat>(v, size);
  tAux.submat(0, 0, size -1, size -1) = diagLeftMult<cx_mat>(v, enMat - ham->blockH(0,0,k), size);
  tAux.submat(0, size, size -1, tSize -1) = -v;
  cx_mat t = tAux;
  cx_mat q,r;

  cx_mat interV = ham->blockH(1,2,k);

  int i,e,minIndex;
  double rTemp, err;
  double deltaErr = tmmErr;
  int printIt = maxItTMM/10;

  for(i = 1; i < maxItTMM; i++){
    if(i % printIt == 0){
      tmmErr += deltaErr;
      cout << __PRETTY_FUNCTION__ << " in iteration " << i << endl;
    }
    if(i % qrIt == 0){
      qr(q, r, tAux);
      tAux = q;
      for(e = 0; e < tSize; e++){
	rTemp = log(abs(r(e,e)));
	c[e] += rTemp;
	d[e] += rTemp*rTemp;
      }
      if(i/qrIt > 10){
	if(testTmmConv(c, d, tSize, i/qrIt, minIndex, err)){
	  break;
	}
      }
    }
    e = i % (nLayers - 1);
    if(e != 0){
      if(e % 2 == 1){
	ham->generateDisorder();
	t.submat(0, 0, size -1, size -1) = diagLeftMult<cx_mat>(interV, enMat - ham->blockH(e,e,k), size);
	t.submat(0, size, size -1, tSize -1) = -diagDoubleMult<cx_mat>(interV, ham->blockH(e,e-1,k), size);
	tAux = updateT<cx_mat>(t, tAux, size);
      }
      else{
	ham->generateDisorder();
	cx_mat v = ham->blockH(e,e+1,k);
	diagInverse(v, size);
	t.submat(0, 0, size -1, size -1) = diagLeftMult<cx_mat>(v, enMat - ham->blockH(e,e,k), size);
	t.submat(0, size, size -1, tSize -1) = -diagDoubleMult<cx_mat>(v, interV, size);
	tAux = updateT<cx_mat>(t, tAux, size);
      }
    }
    else{
      t.submat(0, 0, size -1, size -1) = enMat - ham->blockH(nLayers-1,nLayers-1,k);
      t.submat(0, size, size -1, tSize -1) = -ham->blockH(nLayers-1,nLayers-2,k);
      ham->generateDisorder();
      v = ham->blockH(0,1,k); 
      diagInverse<cx_mat>(v, size);
      t.submat(0, 0, size -1, size -1) = diagLeftMult<cx_mat>(v, t.submat(0, 0, size-1, size-1), size);
      t.submat(0, size, size -1, tSize -1) = diagDoubleMult<cx_mat>(v, t.submat(0, size, size -1, tSize -1), size);
      tAux = updateT<cx_mat>(t, tAux, size);
    }
  }


  if(i == maxItTMM){
    cout << "Warning: transfer matrix method stopped at iteration limit " << maxItTMM << " with error " << err << endl;
  }
  else{
    cout << "Iterations: " << i << endl;
  }

  vector<double> res;
  res.push_back(abs((double)i/c[minIndex]));
  res.push_back(err);
  delete[] c;
  delete[] d;
  return res;
}

vector<double> LocalizationProps::tmmSpecialReal(int nLayers, int qrIt, double en, double * k){
  int size = ham->blockH(0,0,k).n_cols;
  int tSize = 2*size;
  double * c = new double[tSize];
  double * d = new double[tSize];
  for(int i = 0; i < tSize; i++){
    c[i] = 0;
    d[i] = 0;
  }
  mat tAux(tSize, tSize, fill::zeros);
  tAux.submat(size, 0, tSize -1, size -1) = eye<mat>(size, size);
  mat enMat = en*eye<mat>(size,size);
  mat v = real(ham->blockH(0,1,k)); 
  diagInverse<mat>(v, size);
  tAux.submat(0, 0, size -1, size -1) = diagLeftMult<mat>(v, enMat - real(ham->blockH(0,0,k)), size);
  tAux.submat(0, size, size -1, tSize -1) = -v;
  mat t = tAux;
  mat q,r;

  int i,e,minIndex;
  double rTemp, err;
  double deltaErr = tmmErr;
  int printIt = maxItTMM/10;

  for(i = 1; i < maxItTMM; i++){
    if(i % printIt == 0){
      tmmErr += deltaErr;
      cout << __PRETTY_FUNCTION__ << " in iteration " << i << endl;
    }
    if(i % qrIt == 0){
      qr(q, r, tAux);
      tAux = q;
      for(e = 0; e < tSize; e++){
	rTemp = log(abs(r(e,e)));
	c[e] += rTemp;
	d[e] += rTemp*rTemp;
      }
      if(i/qrIt > 10){
	if(testTmmConv(c, d, tSize, i/qrIt, minIndex, err)){
	  break;
	}
      }
    }
    e = i % (nLayers - 1);
    if(e != 0){
      if(e % 2 == 1){
	ham->generateDisorder();
	t.submat(0, 0, size -1, size -1) = enMat - real(ham->blockH(e,e,k));
	t.submat(0, size, size -1, tSize -1) = -real(ham->blockH(e,e-1,k));
	tAux = updateT<mat>(t, tAux, size);
      }
      else{
	ham->generateDisorder();
	mat v = real(ham->blockH(e,e+1,k));
	diagInverse(v, size);
	t.submat(0, 0, size -1, size -1) = diagLeftMult<mat>(v, enMat - real(ham->blockH(e,e,k)), size);
	t.submat(0, size, size -1, tSize -1) = -v;
	tAux = updateT<mat>(t, tAux, size);
      }
    }
    else{
      t.submat(0, 0, size -1, size -1) = enMat - real(ham->blockH(nLayers-1,nLayers-1,k));
      ham->generateDisorder();
      v = real(ham->blockH(0,1,k)); 
      diagInverse<mat>(v, size);
      t.submat(0, 0, size -1, size -1) = diagLeftMult<mat>(v, t.submat(0, 0, size-1, size-1), size);
      t.submat(0, size, size -1, tSize -1) = -v;
      tAux = updateT<mat>(t, tAux, size);
    }
  }


  if(i == maxItTMM){
    cout << "Warning: transfer matrix method stopped at iteration limit " << maxItTMM << " with error " << err << endl;
  }
  else{
    cout << "Iterations: " << i << endl;
  }

  vector<double> res;
  res.push_back(abs((double)i/c[minIndex]));
  res.push_back(err);
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
  if(err < tmmErr){
    return true;
  }
  return false;
}

template <class mat> void LocalizationProps::diagInverse(mat & v, int size){
  for(int i = 0; i < size; i++){
    v(i,i) = (double)1/v(i,i);
  }
}

template <class mat> mat LocalizationProps::diagLeftMult(mat & d, mat m, int size){
  mat res(size, size);
  for(int i = 0; i < size; i++){
    res.row(i) = d(i,i)*m.row(i);
  }
  return res;
}

template <class mat> mat LocalizationProps::diagDoubleMult(mat & d1, mat d2, int size){
  mat res(size, size, fill::zeros);
  for(int i = 0; i < size; i++){
    res(i,i) = d1(i,i)*d2(i,i);
  }
  return res;
}

template <class mat> mat LocalizationProps::updateT(mat & t, mat & tAux, int size){
  mat res = mat(2*size, 2*size);
  res.submat(size, 0, 2*size-1, size-1) = tAux.submat(0, 0, size-1, size-1);
  res.submat(size, size, 2*size-1, 2*size-1) = tAux.submat(0, size, size-1, 2*size-1);
  mat d = t.submat(0,size,size-1,2*size-1);
  res.submat(0, 0, size-1, size-1) = t.submat(0,0,size-1,size-1)*tAux.submat(0,0,size-1,size-1) + diagLeftMult<mat>(d, tAux.submat(size,0,2*size-1,size-1), size);
  res.submat(0, size, size-1, 2*size-1) = t.submat(0,0,size-1,size-1)*tAux.submat(0,size,size-1,2*size-1) + diagLeftMult<mat>(d, tAux.submat(size,size,2*size-1,2*size-1), size);
  return res;
}
