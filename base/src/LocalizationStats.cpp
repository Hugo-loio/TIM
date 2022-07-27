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

double LocalizationStats::tmm(int nLayers, int nIt, double en, double * k){
  int size = ham->blockH(0,0,k).n_cols;
  cout << size << endl;
  cx_mat tN(2*size, 2*size, fill::zeros);
  tN.submat(size, 0, 2*size -1, size -1) = eye<cx_mat>(size, size);
  cx_mat enMat = en*eye<cx_mat>(size,size);
  cx_mat v = ham->blockH(0,1,k).i();
  tN.submat(0, 0, size -1, size -1) = v*(enMat - ham->blockH(0,0,k));
  tN.submat(0, size, size -1, 2*size -1) = -v;
  cx_mat t = tN;

  ham->generateDisorder();
  for(int i = 0; i < nIt; i++){
    for(int e = 1; e < nLayers; e++){
      v = ham->blockH(e,e+1,k).i();
      t.submat(0, 0, size -1, size -1) = v*(enMat - ham->blockH(e,e,k));
      t.submat(0, size, size -1, 2*size -1) = -v*ham->blockH(e,e-1,k);
      tN = t*tN;
    }
    t.submat(0, 0, size -1, size -1) = enMat - ham->blockH(nLayers-1,nLayers-1,k);
    t.submat(0, size, size -1, 2*size -1) = -ham->blockH(nLayers-1,nLayers-2,k);
    ham->generateDisorder();
    v = ham->blockH(0,1,k).i();
    t.submat(0, 0, size -1, size -1) = v*t.submat(0, 0, size -1, size -1);
    t.submat(0, size, size -1, 2*size -1) = v*t.submat(0, size, size -1, 2*size -1);
    tN = t*tN;
  }
  cout << tN << endl;

  cx_mat gamma = tN*tN.t();
  gamma = powmat(gamma, 1);
  cout << gamma << endl;
  vec eigVal = eig_sym(gamma);
  cout << nLayers*nIt << endl;
  for(int i = 0; i < eigVal.size(); i++){
    cout << pow(eigVal(i), 1/(double)(nLayers*nIt)) << endl;
  }
  return 0;
}
