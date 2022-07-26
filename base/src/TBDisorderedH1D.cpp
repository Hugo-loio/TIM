#include "TBDisorderedH1D.h"

TBDisorderedH1D::TBDisorderedH1D(TBModel model) : TBCleanH(model){
  nDisHop = 0;
  nDisOnSite = 0;
  isDisHop = new bool[nHop];
  indexDisHop = new int[nHop];
  for(int i = 0; i < nHop; i++){
    isDisHop[i] = false;
    indexDisHop[i] = 0;
  }
  isDisOnSite = new bool[nOnSite];
  indexDisOnSite = new int[nOnSite];
  for(int i = 0; i < nOnSite; i++){
    isDisOnSite[i] = false;
    indexDisOnSite[i] = 0;
  }
}

TBDisorderedH1D::~TBDisorderedH1D(){
  deleteDisArrays();
  delete[] isDisHop;
  delete[] isDisOnSite;
  delete[] indexDisHop;
  delete[] indexDisOnSite;
}

void TBDisorderedH1D::setDisHop(int nDisHop, int * hop){
  deleteDisArrays();
  for(int i = 0; i < nHop; i++){
    isDisHop[i] = false;
    indexDisHop[i] = 0;
  }
  for(int i = 0; i < nDisHop; i++){
    isDisHop[hop[i]] = true;
    indexDisHop[hop[i]] = i;
  }
  this->nDisHop = nDisHop;
}

void TBDisorderedH1D::setDisOnSite(int nDisOnSite, int * onSite){
  deleteDisArrays();
  for(int i = 0; i < nOnSite; i++){
    isDisOnSite[i] = false;
    indexDisOnSite[i] = 0;
  }
  for(int i = 0; i < nDisOnSite; i++){
    isDisOnSite[onSite[i]] = true;
    indexDisOnSite[onSite[i]] = i;
  }
  this->nDisOnSite = nDisOnSite;
}

void TBDisorderedH1D::setSize(int * l){
  deleteDisArrays();
  TBCleanH::setSize(l);
}

void TBDisorderedH1D::deleteDisArrays(){
  if(disHop != NULL && disOnSite != NULL){
    for(int i = 0; i < nDisHop; i++){
      delete[] disHop[i];
    }
    for(int i = 0; i < nDisOnSite; i++){
      delete[] disOnSite[i];
    }
    delete[] disHop;
    delete[] disOnSite;
    disHop = NULL;
    disOnSite = NULL;
  }
}

void TBDisorderedH1D::createDisArrays(){
  disHop = new complex<double> * [nDisHop];
  for(int i = 0; i < nDisHop; i++){
    disHop[i] = new complex<double> [l[0]];
  }

  disOnSite = new complex<double> * [nDisOnSite];
  for(int i = 0; i < nDisOnSite; i++){
    disOnSite[i] = new complex<double> [l[0]];
  }
}

cx_mat TBDisorderedH1D::H(double* k){
  if(nRDim != 2 || nDim != 2){
    cout << "1D disordered class doesn't have a 1D TB model in real space. Expect things to go wrong, not my problem." << endl;
  }
  return TBCleanH::H(k);
}

sp_cx_mat TBDisorderedH1D::spH(double* k){
  if(nRDim != 2 || nDim != 2){
    cout << "1D disordered class doesn't have a 1D TB model in real space. Expect things to go wrong, not my problem." << endl;
  }
  return TBCleanH::spH(k);
}

cx_mat TBDisorderedH1D::blockH(int line, int col, double* k){
  if(nRDim != 2 || nDim != 2){
    cout << "1D disordered class doesn't have a 1D TB model in real space. Expect things to go wrong, not my problem." << endl;
  }
  return TBCleanH::blockH(line, col, k);
}

void TBDisorderedH1D::fillHopWrapper(cx_mat & res, int hopIndex, double *k, int dim, int addI, int addJ, int addN, int * startI){
  fillHop<cx_mat>(res, hopIndex, k, dim, addI, addJ, addN, startI);
}

void TBDisorderedH1D::fillHopWrapper(sp_cx_mat & res, int hopIndex, double *k, int dim, int addI, int addJ, int addN, int * startI){
  fillHop<sp_cx_mat>(res, hopIndex, k, dim, addI, addJ, addN, startI);
}

void TBDisorderedH1D::fillOnSiteWrapper(cx_mat & res, int onSiteIndex, double * k, int dim, int addN, int * startI){
  fillOnSite<cx_mat>(res, onSiteIndex, k, dim, addN, startI);
}

void TBDisorderedH1D::fillOnSiteWrapper(sp_cx_mat & res, int onSiteIndex, double * k, int dim, int addN, int * startI){
  fillOnSite<sp_cx_mat>(res, onSiteIndex, k, dim, addN, startI);
}

template <class mat> void TBDisorderedH1D::fillOnSite(mat & res, int onSiteIndex, double * k, int dim, int addN, int * startI){
  int e = onSiteIndex;
  int j,n;

  int * start = new int[nRDim + 1];
  int * end = new int[dim + 1];
  end[dim] = 1;
  if(startI != NULL){
    end[dim] = startI[0] + 1;
  }
  for(j = 0; j < nRDim + 1; j++){
    start[j] = 0;
  }

  for(j = 0; j < dim; j++){
    end[j] = l[rIndex[j]] - 1;
  }

  n = flatten2(model.getOnSite(e).getNOrb(), start) + addN;
  if(isDisOnSite[e]){
    fillDis<mat>(res, disOnSite[indexDisOnSite[e]], 1, n, incNOnSite[e], start, end, dim, 0, 0, startI);
  }
  else{
    fill<mat>(res, model.getOnSite(e).getEn(), n, incNOnSite[e], start, end, dim, 0, 0, startI);
  }

  delete[] start;
  delete[] end;
}

template <class mat> void TBDisorderedH1D::fillHop(mat & res, int hopIndex, double * k, int dim, int addI, int addJ, int addN, int * startI){
  complex<double> t;
  int j,n,q,r;
  complex<double> phase;
  complex<double> ii(0,1);

  int * b = new int[dim + 1];
  int * start = new int[nRDim + 1];
  int * end = new int[dim + 1];
  if(startI != NULL){
    end[dim] = startI[0] + 1;
  }
  else{
    end[dim] = 1;
  }
  for(j = dim; j < nRDim + 1; j++){
    start[j] = 0;
  }
  int e = hopIndex;

  t = model.getHop(e).getHop();

  //Bulk lattice
  if(startHopUCBulk[e][nRDim] != 1){
    for(j = 0; j < dim; j++){
      end[j] = endHopUCBulk[e][j];
    }

    n = flatten2(model.getHop(e).getNOrb1(), startHopUCBulk[e]) + addN;

    if(isDisHop[e]){
      fillDis<mat>(res, disHop[indexDisHop[e]], 1, n, incNHopBulk[e], startHopUCBulk[e], end, dim, addI, nHopBulk[e] + addJ, startI);
    }
    else{
      fill<mat>(res, t, n, incNHopBulk[e], startHopUCBulk[e], end, dim, addI, nHopBulk[e] + addJ, startI);
    }
  }

  //Boundary lattice
  for(j = 0; j < dim + 1; j++){
    b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
  }

  b[0]++;
  if(dim != 0){
    r = 0;
    while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
      b[r] = 0;
      b[++r]++;
      if(r == dim){
	break;
      }
    }
  }

  while(b[dim] == 0){

    q = 0;
    phase = 1;
    for(j = 0; j < dim; j++){
      if(b[j] == 0){
	start[j] = startHopUCBulk[e][j];
	end[j] = endHopUCBulk[e][j];
      }
      else{
	start[j] = startHopUCBound[e][j];
	end[j] = endHopUCBound[e][j];
	q += pow(2,rIndex[j]);
	phase *= exp(-ii*theta[rIndex[j]]);
      }
    }

    t = model.getHop(e).getHop()*phase;

    n = flatten2(model.getHop(e).getNOrb1(), start) + addN;

    if(isDisHop[e]){
      fillDis<mat>(res, disHop[indexDisHop[e]], phase, n, incNHopBound[e][q], start, end, dim, addI, nHopBound[e][q] + addJ, startI);
    }
    else{
      fill<mat>(res, t, n, incNHopBound[e][q], start, end, dim, addI, nHopBound[e][q] + addJ, startI);
    }

    b[0]++;
    r = 0;
    while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
      b[r] = 0;
      b[++r]++;
      if(r == dim){
	break;
      }
    }
  }

  delete[] b;
  delete[] start;
  delete[] end;
} 

template <class mat> void TBDisorderedH1D::fillDis(mat & res, complex<double> * w, complex<double> phase, int n, int * incN, int * start, int * end, int dim, int addI, int addJ, int * startI){
  int * i = new int[nRDim + 1];
  for(int j = 0; j < nRDim + 1; j++){
    i[j] = start[j];
  }
  if(startI != NULL){
    for(int j = dim; j < nRDim; j++){
      i[j] = startI[j - dim];
    }
  }
  int val = i[dim];
  int p;
  while(i[dim] == val){
    //cout << n << " " << addI << " " << addJ << " " << end[0] << " " << end[1] << " " << end[2] << " " << dim << " " << i[dim] << " " << start[dim] << endl;
    res(n + addI, n + addJ) += w[i[0]]*phase;

    i[0] += 1;
    n += incN[0];
    p = 0;
    while(i[p] > end[p]){
      i[p] = start[p];
      i[++p] += 1;
      n += incN[p];
    }
  }
  delete[] i;
}
