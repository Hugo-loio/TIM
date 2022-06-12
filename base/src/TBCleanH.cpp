#include "TBCleanH.h"
#include <iostream>

using namespace arma;
using namespace std;

TBCleanH::TBCleanH(TBModel model) : model(model), Hamiltonian(model.getNDim()){
  nOrb = model.getNOrb();
  nHop = model.getNHop();
  nOnSite = model.getNOnSite();

  nRDim = 0;
  bC = new int[nDim];
  l = new int[nDim];
  lAccum = new int[nDim];
  lBound = new int * [nDim];
  rIndex = new int[nRDim];
  nBulk = NULL;
  nBound = NULL;

  for(int i = 0; i < nDim; i++){
    bC[i] = 1;
    l[i] = 1;
    lBound[i] = new int[2];
    lBound[i][0] = 0;
    lBound[i][1] = 0;
    int n = 0;
    for(int e = 0; e < nHop; e++){
      n = model.getHop(e).getN(i);
      if(abs(n) > l[i]){
	l[i] = abs(n);
      }
      if(n > 0 && n > lBound[i][1]){
	lBound[i][1] = n;
      }
      else if(n < 0 && n < -lBound[i][0]){
	lBound[i][0] = - n;
      }
    }
  } 
  calcVol();
}

TBCleanH::~TBCleanH(){
  delete[] bC;
  delete[] l;
  delete[] lAccum;
  for(int i = 0; i < nDim; i++){
    delete[] lBound[i];
  }
  delete[] lBound;
  delete[] rIndex;
  delete_meshes();
}

void TBCleanH::delete_meshes(){
  if(nBulk != NULL){
    for(int i = 0; i < vBulk; i++){
      delete nBulk[i];
    }
    delete[] nBulk;
  }
  if(nBound != NULL){
    for(int i = 0; i < vBound; i++){
      delete nBound[i];
    }
    delete[] nBound;
  }
}

void TBCleanH::setSparse(bool val){
  isSparse = val;
}

void TBCleanH::setSize(int * l){
  if(l != NULL){
    for(int i = 0; i < nDim; i++){
      this->l[i] = l[i];
    }
  }
  checkSize();
  calcVol();
}

void TBCleanH::setBC(int * bC){
  nRDim = 0;
  for(int i = 0; i < nDim; i++){
    this->bC[i] = bC[i];
    if(bC[i] != 1){
      nRDim++;
    }
  }

  delete[] rIndex;
  rIndex = new int[nRDim];
  int e = 0;
  for(int i = 0; i < nDim; i++){
    if(bC[i] != 1){
      rIndex[e] = i;
      e++;
    }
  }

  delete[] lAccum;
  lAccum = new int[nRDim];

  calcVol();
}

void TBCleanH::checkSize(){
  bool changed = false;
  for(int e = 0; e < nHop; e++){
    for(int i = 0; i < nDim; i++){
      if(abs(model.getHop(e).getN(i)) > l[i]){
	l[i] = abs(model.getHop(e).getN(i));
	changed = true;
      }
    }
  }
  if(changed){
    cout << "System size wasn't large enough for the hopping terms\nNew system size:" << endl;
    for(int i = 0; i < nDim; i++){
      cout << "l[" << i << "] = " << l[i] << endl;
    }
  }
}

void TBCleanH::calcVol(){
  //Delete meshes before calculating volumes
  delete_meshes();

  //Calculate lAccum
  for(int i = 0; i < nRDim; i++){
    lAccum[i] = 1;
    for(int e = 0; e <= i; e++){
      lAccum[i] *= l[rIndex[e]];
    }
  }

  vBulk = 1;
  vBound = 0;
  int e;
  for(int i = 0; i < nRDim; i++){
    e = rIndex[i];
    vBulk *= l[e] -(lBound[e][0] + lBound[e][1]);
  }
  if(vBulk < 0){
    vBulk = 0;
  }
  if(nRDim != 0){
    vBound = lAccum[nRDim - 1] - vBulk;
  }

  //Calculate meshes
  calcN();
}

int TBCleanH::getN(int * n){
  int res = 1;
  for(int i = 0; i < nRDim; i++){
    if(i == 0){
      res = n[0];
    }
    else{
      res += lAccum[i-1]*n[i];
    }
  }
  return res;
}

void TBCleanH::calcN(){
  nBulk = new int * [vBulk];
  nBound = new int * [vBound];

  int * point = new int[nRDim];
  for(int i = 0; i < nRDim; i++){
    point[i] = 0;
  }


  int e,j;
  int countBulk = 0;
  int countBound = 0;
  bool bound = false;
  if(nRDim != 0){
    for(int i = 0; i < lAccum[nRDim -1] ;i++){
      for(e = 0; e < nRDim; e++){
	j = rIndex[e];
	if(point[e] > l[j] - lBound[j][1] - 1 || point[e] < lBound[j][0]){
	  bound = true;
	}
      }
      if(bound){
	nBound[countBound] = new int[nRDim + 1];
	for(e = 0; e < nRDim; e++){
	  nBound[countBound][e] = point[e];
	}
	nBound[countBound][nRDim] = getN(point);
	countBound++;
      }
      else{
	nBulk[countBulk] = new int[nRDim + 1];
	for(e = 0; e < nRDim; e++){
	  nBulk[countBulk][e] = point[e];
	}
	nBulk[countBulk][nRDim] = getN(point);
	countBulk++;
      }
      bound = false;
      nextPoint(0, point, false);
    }
  }
  else{
    nBulk[0] = new int[2];
    nBulk[0][0] = 0;
    nBulk[0][1] = 0;
  }
  delete[] point;
}

void TBCleanH::nextPoint(int depth, int * point, bool up){
  if(depth == nRDim-1){
    //increase coordinate
    if(point[nRDim-1] == l[rIndex[nRDim -1]] - 1){
      point[nRDim -1] = 0; 
      nextPoint(--depth, point, true);
    }
    else{
      point[nRDim-1]++;
    }
  }
  else if(depth != -1){
    if(up){
      //increase coordinate
      if(point[depth] == l[rIndex[depth]] -1 ){
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

cx_mat TBCleanH::H(double * k){
  int size; 

  if(nRDim != 0){
    size = nOrb*lAccum[nRDim-1];
  }
  else{
    size = nOrb;
  }
  cx_mat res(size, size, fill::zeros);

  complex<double> t;
  int i,j,h;
  complex<double> phase;
  complex<double> kPhase;
  complex<double> ii(0,1);

  int * newNBound = new int[nRDim];
  int incN;

  //Hopping terms
  for(int e = 0; e < nHop; e++){
    //Increase in n to get neighbour cell in hop (in case of bulk)
    incN = 0;
    j = 0;
    kPhase = 1;
    for(i = 0; i < nDim; i++){
      if(bC[i] != 1){
	if(i == rIndex[0]){
	  incN = model.getHop(e).getN(i);
	}
	else{
	  incN += model.getHop(e).getN(i)*lAccum[j-1];
	}
	j++;
      }
      else{
	//Extra phase due to k space hamiltonian
	kPhase *= exp(ii*k[i]*(double)model.getHop(e).getN(i));
      }
    }

    //Go through bulk mesh
    for(i = 0; i < vBulk; i++){
      res(model.getHop(e).getNOrb1() + nOrb*nBulk[i][nRDim], model.getHop(e).getNOrb2() + nOrb*(nBulk[i][nRDim] + incN)) += model.getHop(e).getHop()*kPhase;
    }
    //Go through boundary mesh and check BCs
    for(i = 0; i < vBound; i++){
      phase = 1;
      for(j = 0; j < nRDim; j++){
	h = rIndex[j];
	newNBound[j] = (nBound[i][j] + model.getHop(e).getN(h) + l[h]) % l[h];
	if(nBound[i][j] + model.getHop(e).getN(h) > l[h] - 1){
	  switch(bC[h]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(-ii*theta[j]);
	  }
	}
	else if(nBound[i][j] + model.getHop(e).getN(h) < 0){
	  switch(bC[h]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(ii*theta[j]);
	  }
	}
      }
      res(model.getHop(e).getNOrb1() + nOrb*nBound[i][nRDim], model.getHop(e).getNOrb2() + nOrb*getN(newNBound)) += kPhase*phase*model.getHop(e).getHop();
    }
  }

  delete[] newNBound;

  res = res + res.t();

  //On-site terms
  int m;
  for(int e  = 0; e < nOnSite; e++){
    //Go through bulk mesh
    for(i = 0; i < vBulk; i++){
      m = model.getOnSite(e).getNOrb() + nOrb*nBulk[i][nRDim];
      res(m,m) += model.getOnSite(e).getEn();
    }
    //Go through boundary mesh
    for(i = 0; i < vBound; i++){
      m = model.getOnSite(e).getNOrb() + nOrb*nBound[i][nRDim];
      res(m,m) += model.getOnSite(e).getEn();
    }
  }

  return res;
}

sp_cx_mat TBCleanH::spH(double * k){
  int size; 

  if(nRDim != 0){
    size = nOrb*lAccum[nRDim-1];
  }
  else{
    size = nOrb;
  }
  sp_cx_mat res(size, size);

  complex<double> t;
  int i,j,h;
  complex<double> phase;
  complex<double> kPhase;
  complex<double> ii(0,1);

  int * newNBound = new int[nRDim];
  int incN;

  //Hopping terms
  for(int e = 0; e < nHop; e++){
    //Increase in n to get neighbour cell in hop (in case of bulk)
    incN = 0;
    j = 0;
    kPhase = 1;
    for(i = 0; i < nDim; i++){
      if(bC[i] != 1){
	if(i == rIndex[0]){
	  incN = model.getHop(e).getN(i);
	}
	else{
	  incN += model.getHop(e).getN(i)*lAccum[j-1];
	}
	j++;
      }
      else{
	//Extra phase due to k space hamiltonian
	kPhase *= exp(ii*k[i]*(double)model.getHop(e).getN(i));
      }
    }

    //Go through bulk mesh
    for(i = 0; i < vBulk; i++){
      res(model.getHop(e).getNOrb1() + nOrb*nBulk[i][nRDim], model.getHop(e).getNOrb2() + nOrb*(nBulk[i][nRDim] + incN)) += model.getHop(e).getHop()*kPhase;
    }
    //Go through boundary mesh and check BCs
    for(i = 0; i < vBound; i++){
      phase = 1;
      for(j = 0; j < nRDim; j++){
	h = rIndex[j];
	newNBound[j] = (nBound[i][j] + model.getHop(e).getN(h) + l[h]) % l[h];
	if(nBound[i][j] + model.getHop(e).getN(h) > l[h] - 1){
	  switch(bC[h]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(-ii*theta[j]);
	  }
	}
	else if(nBound[i][j] + model.getHop(e).getN(h) < 0){
	  switch(bC[h]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(ii*theta[j]);
	  }
	}
      }
      res(model.getHop(e).getNOrb1() + nOrb*nBound[i][nRDim], model.getHop(e).getNOrb2() + nOrb*getN(newNBound)) += kPhase*phase*model.getHop(e).getHop();
    }
  }

  delete[] newNBound;

  res = res + res.t();

  //On-site terms
  int m;
  for(int e  = 0; e < nOnSite; e++){
    //Go through bulk mesh
    for(i = 0; i < vBulk; i++){
      m = model.getOnSite(e).getNOrb() + nOrb*nBulk[i][nRDim];
      res(m,m) += model.getOnSite(e).getEn();
    }
    //Go through boundary mesh
    for(i = 0; i < vBound; i++){
      m = model.getOnSite(e).getNOrb() + nOrb*nBound[i][nRDim];
      res(m,m) += model.getOnSite(e).getEn();
    }
  }

  return res;
}
