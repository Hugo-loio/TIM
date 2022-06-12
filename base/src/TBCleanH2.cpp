#include "TBCleanH2.h"
#include <iostream>

using namespace arma;
using namespace std;

TBCleanH2::TBCleanH2(TBModel model) : model(model), Hamiltonian(model.getNDim()){
  nOrb = model.getNOrb();
  nHop = model.getNHop();
  nOnSite = model.getNOnSite();

  nRDim = 0;
  bC = new int[nDim];
  l = new int[nDim];
  lAccum = new int[nDim];
  rIndex = new int[nRDim];
  rOrder = new int[nDim];
  lOrb = new int * [nOrb];
  mOrb = new int [nOrb];
  nu = new int[nDim];
  nuAccum = new int[nRDim];
  nHopBulk = new int [nHop];
  nHopBound = new int * [nHop];
  incHop = new int * [nHop];
  incNHop = new int * [nHop];
  startHopBulk = new int * [nHop];
  startHopBound = new int * [nHop];
  endHopBulk = new int * [nHop];
  endHopBound = new int * [nHop];
  incOnSite = new int * [nOnSite];
  incNOnSite = new int * [nOnSite];
  startOnSite = new int * [nOnSite];
  endOnSite = new int * [nOnSite];

  for(int i = 0; i < nDim; i++){
    bC[i] = 1;
    l[i] = 1;
    rOrder[i] = i;
    nu[i] = 1;
    int n = 0;
    for(int e = 0; e < nHop; e++){
      n = model.getHop(e).getN(i);
      if(abs(n) > l[i]){
	l[i] = abs(n);
      }
    }
  } 

  for(int i = 0; i < nOrb; i++){
    lOrb[i] = new int[nDim];
    for(int e = 0; e < nDim; e++){
      lOrb[i][e] = 0;
    }
    mOrb[i] = i;
  }

  for(int i = 0; i < nHop; i++){
    nHopBound[i] = new int[(int)pow(2,nDim)];
    incHop[i] = new int[nRDim + 1];
    incNHop[i] = new int[nRDim + 1];
    startHopBulk[i] = new int[nRDim];
    startHopBound[i] = new int[nRDim];
    endHopBulk[i] = new int[nRDim];
    endHopBound[i] = new int[nRDim];
  }
  for(int i = 0; i < nOnSite; i++){
    incOnSite[i] = new int[nRDim + 1];
    incNOnSite[i] = new int[nRDim + 1];
    startOnSite[i] = new int[nRDim];
    endOnSite[i] = new int[nRDim];
  }
}

TBCleanH2::~TBCleanH2(){
  delete[] bC;
  delete[] l;
  for(int i = 0; i < nOrb; i++){
    delete[] lOrb;
  }
  delete[] lOrb;
  delete[] mOrb;
  delete[] rIndex;
  delete[] rOrder;
  delete[] nuAccum;
  delete[] nHopBulk;
  for(int i = 0; i < nHop; i++){
    delete[] nHopBound[i];
    delete[] incHop[i];
    delete[] incNHop[i];
    delete[] startHopBulk[i];
    delete[] startHopBound[i];
    delete[] endHopBulk[i];
    delete[] endHopBound[i];
  }
  delete[] nHopBound;
  delete[] incHop;
  delete[] incNHop;
  delete[] startHopBulk;
  delete[] startHopBound;
  delete[] endHopBulk;
  delete[] endHopBound;
  for(int i = 0; i < nOnSite; i++){
    delete[] incOnSite[i];
    delete[] incNOnSite[i];
    delete[] startOnSite[i];
    delete[] endOnSite[i];
  }
  delete[] incOnSite;
  delete[] incNOnSite;
  delete[] startOnSite;
  delete[] endOnSite;
}

void TBCleanH2::setSparse(bool val){
  isSparse = val;
}

void TBCleanH2::setSize(int * l){
  if(l != NULL){
    for(int i = 0; i < nDim; i++){
      this->l[i] = l[i];
    }
  }

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

void TBCleanH2::setBC(int * bC){
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
  delete[] nuAccum;
  lAccum = new int[nRDim];
  nuAccum = new int[nRDim];

  for(int i = 0; i < nHop; i++){
    delete[] incHop[i];
    delete[] incNHop[i];
    delete[] startHopBulk[i];
    delete[] startHopBound[i];
    delete[] endHopBulk[i];
    delete[] endHopBound[i];
  }
  for(int i = 0; i < nHop; i++){
    incHop[i] = new int[nRDim + 1];
    incNHop[i] = new int[nRDim + 1];
    startHopBulk[i] = new int[nRDim];
    startHopBound[i] = new int[nRDim];
    endHopBulk[i] = new int[nRDim];
    endHopBound[i] = new int[nRDim];
  }
  for(int i = 0; i < nOnSite; i++){
    delete[] incOnSite[i];
    delete[] incNOnSite[i];
    delete[] startOnSite[i];
    delete[] endOnSite[i];
  }
  for(int i = 0; i < nOnSite; i++){
    incOnSite[i] = new int[nRDim + 1];
    incNOnSite[i] = new int[nRDim + 1];
    startOnSite[i] = new int[nRDim];
    endOnSite[i] = new int[nRDim];
  }

  sortRIndex();
}

void TBCleanH2::setOrder(int * o){
  for(int i = 0; i < nDim; i++){
    rOrder[i] = o[i];
  }
  sortRIndex();
}

void TBCleanH2::setOrbLayer(int ** l){
  //Check if coordinates are valid
  int * lDist = new int[nOrb];
  int * nuTemp = new int[nDim];
  for(int i = 0; i < nDim; i++){
    nuTemp[i] = 1;
  }
  int lCoord;
  for(int e = 0; e < nDim; e++){
    for(int i = 0; i < nOrb; i++){
      lDist[i] = 0;
    }
    for(int i = 0; i < nOrb; i++){
      if(l[i][e] < nOrb && l[i][e] >= 0){
	lDist[l[i][e]]++;
      }
      else{
	cout << "Bad layer coordinates, layers not changed" << endl;
	delete[] lDist;
	delete[] nuTemp;
	return;
      }
    }
    for(int i = 1; i < nOrb; i++){
      if(lDist[i] != 0){
	nuTemp[e]++;
      }
    }
    for(int i = 0; i < nuTemp[e]; i++){
      if(lDist[i] != nOrb/nuTemp[e]){
	cout << "Orbitals not equally divides between layers, layers not changed" << endl;
	delete[] lDist;
	delete[] nuTemp;
	return;
      }
    }
  }

  for(int i = 0; i < nOrb; i++){
    for(int e = 0; e < nDim; e++){
      lOrb[i][e] = l[i][e];
    }
  }
  for(int i = 0; i < nDim; i++){
    nu[i] = nuTemp[i];
  }

  map<int *, vector<int>> map; //Maps layer coord to orbitals in layer

  for(int i = 0; i < nOrb; i++){
    if(map.count(l[i]) == 0){
      vector <int> o;
      o.push_back(i);
      map[l[i]] = o;
    }
    else{
      map[l[i]].push_back(i);
    }
  }

  for(std::map<int *, vector<int>>::iterator i = map.begin(); i != map.end(); i++){
    for(int e = 0; e < i->second.size(); e++){
      mOrb[i->second[e]] = e;
    }
  }

  delete[] lDist;
  delete[] nuTemp;
}

void TBCleanH2::sortRIndex(){
  int tempDir;
  for(int i = 0; i < nRDim; i++){
    for(int e = i + 1; e < nRDim; e++){
      if(rOrder[rIndex[i]] > rOrder[rIndex[e]]){
	tempDir = rIndex[e];
	rIndex[e] = rIndex[i];
	rIndex[i] = tempDir;
      }
    }
  }
}

int TBCleanH2::flatten(int alpha, int * n){
  int res = mOrb[alpha];
  for(int i = 0; i < nRDim; i++){
    res += ((nOrb*n[rIndex[i]] + (nOrb*lOrb[alpha][rIndex[i]])/nu[rIndex[i]])*lAccum[i])/nuAccum[i];
  }
  return res;
}

void TBCleanH2::calcAux(){
  //Calculate lAccum
  for(int i = 0; i < nRDim; i++){
    lAccum[i] = 1;
    for(int e = 0; e <= i; e++){
      lAccum[i] *= l[rIndex[e]];
    }
  }
  //Calculate nuAccum
  for(int i = 0; i < nRDim; i++){
    nuAccum[i] = 1;
    for(int e = i + 1; e < nRDim; e++){
      nuAccum[i] *= nu[rIndex[e]];
    }
  }

  //Calculate increments in flattened indices
  int * nAux = new int[nDim];
  int * nAux2 = new int[nDim];
  int * boundDir  = new int[nDim];
  int n;
  for(int i = 0; i < nHop; i++){

    for(int e = 0; e < nDim; e++){
      nAux[e] = 0;
    }

    bool isBulk = true;
    int nDir = 0;
    for(int e = 0; e < nDim; e++){
      n = model.getHop(i).getN()[e];
      if(bC[e] != 1){
	if(abs(n) >= l[e]){
	  isBulk = false;
	}
	else if(n < 0){
	  nAux[e] -= n;
	}
      }
      if(bC[e] == 2 && n != 0){
	nDir++;
	boundDir[e] = 1;
      }
      else{
	boundDir[e] = 0;
      }
      nAux2[e] = nAux[e] + n;
    }


    if(isBulk){
      nHopBulk[i] = flatten(model.getHop(i).getNOrb2(), nAux2) - flatten(model.getHop(i).getNOrb1(), nAux);
    }
    else{
      nHopBulk[i] = 0;
    }

    int * boundOpt = new int[nDim + 1];
    for(int e = 0; e < nDim + 1; e++){
      boundOpt[e] = 0;
    }
    int p;
    bool isBound = false;
    int index;
    while(boundOpt[nDim] != 1){

      index = 0;
      for(int e = 0; e < nDim; e++){
	nAux[e] = 0;
	n = model.getHop(i).getN()[e];
	index += e*pow(2,e);
	if(boundOpt[e] == 1){
	  if(n != 0){
	    isBound = true;
	    if(n > 0){
	      nAux[e] = l[e] - 1;
	    }
	    nAux2[e] = (nAux[e] + n + l[e]) % l[e];
	  }
	}
      } 

      if(isBound){
	nHopBound[i][index] = flatten(model.getHop(i).getNOrb2(), nAux2) - flatten(model.getHop(i).getNOrb1(), nAux);
      }
      else{
	nHopBound[i][index] = nHopBulk[i];
      }

      boundOpt[0]++;
      p = 0;
      while(boundOpt[p] > 1){
	boundOpt[p] = 0;
	boundOpt[++p]++;
      }

    }

  }

  delete[] boundDir;
  delete[] nAux;
  delete[] nAux2;

  //Calculate limits and increments to loop over the lattice
  //Hopping
  for(int i = 0; i < nHop; i++){
    int nOrb1 = model.getHop(i).getNOrb1();
    int nH;
    for(int e = 0; e < nRDim; e++){
      nH = model.getHop(i).getN()[e];
      incHop[i][e] = mOrb[nOrb1] + lOrb[nOrb1][rIndex[e]]*nOrb/nuAccum[nRDim - 1];
      if(n < 0){
	startHopBulk[i][e] = (1-nH)*incHop[i][e];
	startHopBound[i][e] = incHop[i][e];
	endHopBound[i][e] = -nH*incHop[i][e];
	endHopBulk[i][e] = (l[rIndex[e]]*nOrb*nu[rIndex[e]])/nuAccum[nRDim -1] - incHop[i][e];
      }
      else if(n > 0){
	startHopBulk[i][e] = incHop[i][e];
	startHopBound[i][e] = (l[rIndex[e]]*nOrb*nu[rIndex[e]])/nuAccum[nRDim -1] - (nH - 1)*incHop[i][e];
	endHopBound[i][e] = (l[rIndex[e]]*nOrb*nu[rIndex[e]])/nuAccum[nRDim -1] - incHop[i][e];
	endHopBulk[i][e] = (l[rIndex[e]]*nOrb*nu[rIndex[e]])/nuAccum[nRDim -1] - nH*incHop[i][e];
      }
      else{
	startHopBulk[i][e] = incHop[i][e];
	endHopBulk[i][e] = (l[rIndex[e]]*nOrb*nu[rIndex[e]])/nuAccum[nRDim -1] - incHop[i][e];
	startHopBound[i][e] = 0;
	endHopBound[i][e] = 0;
      }
      if(e == 0){
	incNHop[i][e] = incHop[i][e];
      }
      else{
	incNHop[i][e] = incHop[i][e]*lAccum[e]*nOrb/nu[e];
      }
    }
    incHop[i][nRDim] = 1;
    incNHop[i][nRDim] = 1;
  }

  //On-site
  for(int i = 0; i < nOnSite; i++){
    int orb = model.getOnSite(i).getNOrb();
    for(int e = 0; e < nRDim; e++){
      incOnSite[i][e] = mOrb[orb] + lOrb[orb][rIndex[e]]*nOrb/nuAccum[nRDim - 1];
      startOnSite[i][e] = incOnSite[i][e];
      endHopBulk[i][e] = (l[rIndex[e]]*nOrb*nu[rIndex[e]])/nuAccum[nRDim -1] - incOnSite[i][e];
      if(e == 0){
	incNOnSite[i][e] = incOnSite[i][e];
      }
      else{
	incNOnSite[i][e] = incOnSite[i][e]*lAccum[e]*nOrb/nu[e];
      }
    }
    incOnSite[i][nRDim] = 1;
    incNOnSite[i][nRDim] = 1;
  }
}

cx_mat TBCleanH2::H(double * k){
  int size; 

  if(nRDim != 0){
    size = nOrb*lAccum[nRDim-1];
  }
  else{
    size = nOrb;
  }
  cx_mat res(size, size, fill::zeros);

  complex<double> t;
  int j,n,p;
  complex<double> phase;
  complex<double> kPhase;
  complex<double> ii(0,1);

  int * i = new int[nRDim + 1];

  //Hopping terms
  for(int e = 0; e < nHop; e++){
    kPhase = 1;
    for(j = 0; j < nDim; j++){
      if(bC[j] == 1){
	//Extra phase due to k space hamiltonian
	kPhase *= exp(ii*k[j]*(double)model.getHop(e).getN(j));
      }
    }

    t = model.getHop(e).getHop()*kPhase;

    //Bulk lattice
    for(j = 0; j < nRDim; j++){
      i[j] = startHopBulk[e][j];
    }
    n = i[0];
    for(j = 1; j < nRDim; j++){
      n += (i[j]*nOrb*lAccum[j])/nu[j];
    }
    i[nRDim] = 0;

    while(i[nRDim] == 0){
      res(n, n + nHopBulk[e]) += t;

      p = 0;
      while(i[p] > endHopBulk[e][p]){
	i[p] = startHopBulk[e][p];
	i[++p] += incHop[e][p];
	n += incNHop[e][p];
      }
    }

    //Boundary lattice

    /*
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
    */
  }
  res = res + res.t();

  //On-site terms
  for(int e = 0; e < nOnSite; e++){

    for(j = 0; j < nRDim; j++){
      i[j] = startOnSite[e][j];
    }
    n = i[0];
    for(j = 1; j < nRDim; j++){
      n += (i[j]*nOrb*lAccum[j])/nu[j];
    }
    i[nRDim] = 0;

    while(i[nRDim] == 0){
      res(n, n) += model.getOnSite(e).getEn();

      p = 0;
      while(i[p] > endHopBulk[e][p]){
	i[p] = startHopBulk[e][p];
	i[++p] += incOnSite[e][p];
	n += incNOnSite[e][p];
      }
    }
  }


  delete[] i;

  return res;
}

/*
   sp_cx_mat TBCleanH2::spH(double * k){
   int size; 

   if(nRDim != 0){
   size = nOrb*lAccum[nRDim-1];
   }
   else{
   size = nOrb;
   }
   sp_cx_mat res(size, size);


   return res;
   }
   */
