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
  rIndex = new int[nRDim];
  rOrder = new int[nDim];
  lOrb = new int * [nOrb];
  mOrb = new int [nOrb];
  nu = new int[nDim];
  nuAccum = new int[nRDim];
  nHopBulk = new int [nHop];
  nHopBound = new int * [nHop];
  inc = new int [nRDim + 1];
  incNHopBulk = new int * [nHop];
  incNHopBound = new int ** [nHop];
  startHopBulk = new int * [nHop];
  startHopBound = new int * [nHop];
  endHopBulk = new int * [nHop];
  endHopBound = new int * [nHop];
  startHopUCBulk = new int * [nHop];
  startHopUCBound = new int * [nHop];
  endHopUCBulk = new int * [nHop];
  endHopUCBound = new int * [nHop];
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
    incNHopBulk[i] = new int[nRDim + 1];
    incNHopBound[i] = new int * [(int)pow(2,nDim)];
    startHopBulk[i] = new int[nRDim];
    startHopBound[i] = new int[nRDim];
    endHopBulk[i] = new int[nRDim + 1];
    endHopBound[i] = new int[nRDim + 1];
    startHopUCBulk[i] = new int [nRDim + 1];
    startHopUCBound[i] = new int [nRDim];
    endHopUCBulk[i] = new int [nRDim + 1];
    endHopUCBound[i] = new int [nRDim + 1];
    for(int e = 0; e < pow(2,nDim); e++){
      incNHopBound[i][e] = new int[nRDim + 1];
    }
  }
  for(int i = 0; i < nOnSite; i++){
    incNOnSite[i] = new int[nRDim + 1];
    startOnSite[i] = new int[nRDim];
    endOnSite[i] = new int[nRDim + 1];
  }

}

TBCleanH::~TBCleanH(){
  delete[] bC;
  delete[] l;
  for(int i = 0; i < nOrb; i++){
    delete[] lOrb[i];
  }
  delete[] lOrb;
  delete[] mOrb;
  delete[] rIndex;
  delete[] rOrder;
  delete[] nuAccum;
  delete[] nHopBulk;
  for(int i = 0; i < nHop; i++){
    for(int e = 0; e < pow(2,nDim); e++){
      delete[] incNHopBound[i][e];
    }
    delete[] nHopBound[i];
    delete[] incNHopBulk[i];
    delete[] incNHopBound[i];
    delete[] startHopBulk[i];
    delete[] startHopBound[i];
    delete[] endHopBulk[i];
    delete[] endHopBound[i];
    delete[] startHopUCBulk[i];
    delete[] startHopUCBound[i];
  }
  delete[] nHopBound;
  delete[] inc;
  delete[] incNHopBulk;
  delete[] incNHopBound;
  delete[] startHopBulk;
  delete[] startHopBound;
  delete[] endHopBulk;
  delete[] endHopBound;
  delete[] startHopUCBulk;
  delete[] startHopUCBound;
  delete[] endHopUCBulk;
  delete[] endHopUCBound;
  for(int i = 0; i < nOnSite; i++){
    delete[] incNOnSite[i];
    delete[] startOnSite[i];
    delete[] endOnSite[i];
  }
  delete[] incNOnSite;
  delete[] startOnSite;
  delete[] endOnSite;
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
  isUpdated = false;
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

  bDim = nRDim - 1;

  delete[] lAccum;
  delete[] nuAccum;
  delete[] inc;
  lAccum = new int[nRDim];
  nuAccum = new int[nRDim];
  inc = new int[nRDim + 1];

  for(int i = 0; i < nHop; i++){
    for(int e = 0; e < pow(2,nDim); e++){
      delete[] incNHopBound[i][e];
    }
    delete[] incNHopBulk[i];
    delete[] startHopBulk[i];
    delete[] startHopBound[i];
    delete[] endHopBulk[i];
    delete[] endHopBound[i];
    delete[] startHopUCBulk[i];
    delete[] startHopUCBound[i];
    delete[] endHopUCBulk[i];
    delete[] endHopUCBound[i];
  }
  for(int i = 0; i < nHop; i++){
    incNHopBulk[i] = new int[nRDim + 1];
    startHopBulk[i] = new int[nRDim];
    startHopBound[i] = new int[nRDim];
    endHopBulk[i] = new int[nRDim + 1];
    endHopBound[i] = new int[nRDim + 1];
    startHopUCBulk[i] = new int [nRDim + 1];
    startHopUCBound[i] = new int [nRDim];
    endHopUCBulk[i] = new int [nRDim + 1];
    endHopUCBound[i] = new int [nRDim + 1];
    for(int e = 0; e < pow(2,nDim); e++){
      incNHopBound[i][e] = new int[nRDim + 1];
    }
  }
  for(int i = 0; i < nOnSite; i++){
    delete[] incNOnSite[i];
    delete[] startOnSite[i];
    delete[] endOnSite[i];
  }
  for(int i = 0; i < nOnSite; i++){
    incNOnSite[i] = new int[nRDim + 1];
    startOnSite[i] = new int[nRDim];
    endOnSite[i] = new int[nRDim + 1];
  }
  isUpdated = false;
}

void TBCleanH::setOrder(int * o){
  for(int i = 0; i < nDim; i++){
    rOrder[i] = o[i];
  }
  isUpdated = false;
}

void TBCleanH::setOrbLayer(int ** l){
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


  delete[] lDist;
  delete[] nuTemp;
  isUpdated = false;
}

int TBCleanH::flatten(int alpha, int * n){
  int res = mOrb[alpha];
  for(int i = 0; i < nRDim; i++){
    if(i == 0){
      res += (nOrb*n[rIndex[i]] + (nOrb*lOrb[alpha][rIndex[i]])/nu[rIndex[i]])/nuAccum[i];
    }
    else{
      res += ((nOrb*n[rIndex[i]] + (nOrb*lOrb[alpha][rIndex[i]])/nu[rIndex[i]])*lAccum[i - 1])/nuAccum[i];
    }
  }
  return res;
}

int TBCleanH::flatten2(int alpha, int * n){
  int res = mOrb[alpha];
  for(int i = 0; i < nRDim; i++){
    if(i == 0){
      res += (nOrb*n[i] + (nOrb*lOrb[alpha][rIndex[i]])/nu[rIndex[i]])/nuAccum[i];
    }
    else{
      res += ((nOrb*n[i] + (nOrb*lOrb[alpha][rIndex[i]])/nu[rIndex[i]])*lAccum[i - 1])/nuAccum[i];
    }
  }
  return res;
}

void TBCleanH::calcAux(){
  //Sort rIndex
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

  //Calculate mOrb
  map<int *, vector<int>> map; //Maps layer coord to orbitals in layer
  for(int i = 0; i < nOrb; i++){
    int j = i;
    for(int e = 0; e < i; e++){
      bool diff = false;
      for(int k = 0; k < nRDim; k++){
	if(lOrb[e][rIndex[k]] != lOrb[i][rIndex[k]]){
	  diff = true;
	}
      }
      if(!diff){
	j = e;
	break;
      }
    }

    if(map.count(lOrb[j]) == 0){
      vector <int> o;
      o.push_back(i);
      map[lOrb[j]] = o;
    }
    else{
      map[lOrb[j]].push_back(i);
    }
  }
  for(std::map<int *, vector<int>>::iterator i = map.begin(); i != map.end(); i++){
    for(int e = 0; e < i->second.size(); e++){
      mOrb[i->second[e]] = e;
    }
  }

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
  int * boundOpt = new int[nDim + 1];
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

    for(int e = 0; e < nDim + 1; e++){
      boundOpt[e] = 0;
    }
    boundOpt[0]++;
    nHopBound[i][0] = nHopBulk[i];
    int p;
    bool isBound;
    int index;
    while(boundOpt[nDim] != 1){

      isBound = true;
      index = 0;
      for(int e = 0; e < nDim; e++){
	nAux[e] = 0;
	n = model.getHop(i).getN()[e];
	index += boundOpt[e]*pow(2,e);
	if(boundOpt[e] == 1){
	  if(n != 0){
	    if(n > 0){
	      nAux[e] = l[e] - 1;
	    }
	    nAux2[e] = (nAux[e] + n + l[e]) % l[e];
	  }
	  else{
	    isBound = false;
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

  int nL;
  if(nRDim != 0){
    nL = nuAccum[0]*nu[rIndex[0]];
  }

  //Calculate limits and increments to loop over the lattice
  //Hopping
  for(int i = 0; i < nRDim; i++){
    if(i == 0){
      inc[i] = (nOrb*nu[rIndex[i]])/nL;
    }
    else{
      inc[i] = nu[rIndex[i]];
    }
  }
  inc[nRDim] = 1;
  int * start = new int[nRDim];
  int * end = new int[nRDim];
  for(int i = 0; i < nHop; i++){
    int nOrb1 = model.getHop(i).getNOrb1();
    int nH;
    if(nRDim != 0){
      nH = model.getHop(i).getN()[rIndex[0]];
      if(nH < 0){
	startHopBulk[i][0] = mOrb[nOrb1] + lOrb[nOrb1][rIndex[0]]*nOrb/nL - nH*inc[0];
	startHopBound[i][0] = mOrb[nOrb1] + lOrb[nOrb1][rIndex[0]]*nOrb/nL;
	endHopBound[i][0] = startHopBound[i][0] - (nH+1)*inc[0];
	endHopBulk[i][0] = (l[rIndex[0]]*nOrb*nu[rIndex[0]])/nL - (inc[0] - startHopBound[i][0]);
      }
      else if(nH > 0){
	startHopBulk[i][0] = mOrb[nOrb1] + lOrb[nOrb1][rIndex[0]]*nOrb/nL;
	startHopBound[i][0] = (l[rIndex[0]]*nOrb*nu[rIndex[0]])/nL - nH*inc[0] + startHopBulk[i][0];
	endHopBound[i][0] = (l[rIndex[0]]*nOrb*nu[rIndex[0]])/nL - (inc[0] - startHopBulk[i][0]);
	endHopBulk[i][0] = (l[rIndex[0]]*nOrb*nu[rIndex[0]])/nL - (nH + 1)*inc[0] + startHopBulk[i][0];
      }
      else{
	startHopBulk[i][0] = mOrb[nOrb1] + lOrb[nOrb1][rIndex[0]]*nOrb/nL;
	endHopBulk[i][0] = (l[rIndex[0]]*nOrb*nu[rIndex[0]])/nL - (inc[0] - startHopBulk[i][0]);
	startHopBound[i][0] = 0;
	endHopBound[i][0] = 0;
      }
      incNHopBulk[i][0] = inc[0];
    }
    for(int e = 1; e < nRDim; e++){
      nH = model.getHop(i).getN()[rIndex[e]];
      if(nH < 0){
	startHopBulk[i][e] = lOrb[nOrb1][rIndex[e]] - nH*inc[e];
	startHopBound[i][e] =  lOrb[nOrb1][rIndex[e]];
	endHopBound[i][e] = startHopBound[i][e] - (nH+1)*inc[e];
	endHopBulk[i][e] = l[rIndex[e]]*nu[rIndex[e]] - (inc[e] - startHopBound[i][e]);
      }
      else if(nH > 0){
	startHopBulk[i][e] =  lOrb[nOrb1][rIndex[e]];
	startHopBound[i][e] = l[rIndex[e]]*nu[rIndex[e]] - nH*inc[e] + startHopBulk[i][e];
	endHopBound[i][e] = l[rIndex[e]]*nu[rIndex[e]] - (inc[e] - startHopBulk[i][e]);
	endHopBulk[i][e] = l[rIndex[e]]*nu[rIndex[e]] - (nH + 1)*inc[e] + startHopBulk[i][e];
      }
      else{
	startHopBulk[i][e] =  lOrb[nOrb1][rIndex[e]];
	endHopBulk[i][e] = l[rIndex[e]]*nu[rIndex[e]] - (inc[e] - startHopBulk[i][e]);
	startHopBound[i][e] = 0;
	endHopBound[i][e] = 0;
      }

      //Bulk incN
      incNHopBulk[i][e] = - endHopBulk[i][0] + startHopBulk[i][0] + inc[e]*lAccum[e-1]*nOrb/nuAccum[e-1];
      for(int j = 1; j < e; j++){
	incNHopBulk[i][e] += (startHopBulk[i][j] - endHopBulk[i][j])*lAccum[j-1]*nOrb/nuAccum[j-1];
      }
      for(int j = e-1; j >= 0;j--){
	incNHopBulk[i][e] -= incNHopBulk[i][j];
      }
    }

    //Bound incN
    for(int e = 0; e < nDim + 1; e++){
      boundOpt[e] = 0;
    }
    boundOpt[0]++;
    for(int e = 0; e < nRDim; e++){
      incNHopBound[i][0][e] = incNHopBulk[i][e];
    }
    incNHopBound[i][0][nRDim] = 1;
    int p;
    bool isBound;
    int index;
    while(boundOpt[nDim] != 1){

      isBound = true;
      index = 0;
      for(int e = 0; e < nDim; e++){
	n = model.getHop(i).getN()[e];
	index += boundOpt[e]*pow(2,e);
	if(boundOpt[e] == 1 && n == 0){
	  isBound = false;
	}
      } 

      if(nRDim != 0){
	incNHopBound[i][index][0] = inc[0];
      }
      incNHopBound[i][index][nRDim] = 1;

      if(isBound){
	for(int e = 0; e < nRDim; e++){
	  if(boundOpt[rIndex[e]] == 0){
	    start[e] = startHopBulk[i][e];
	    end[e] = endHopBulk[i][e];
	  }
	  else{
	    start[e] = startHopBound[i][e];
	    end[e] = endHopBound[i][e];
	  }
	}
	for(int e = 1; e < nRDim; e++){
	  incNHopBound[i][index][e] = - end[0] + start[0] + inc[e]*lAccum[e-1]*nOrb/nuAccum[e-1];
	  for(int j = 1; j < e; j++){
	    incNHopBound[i][index][e] += (start[j] - end[j])*lAccum[j-1]*nOrb/nuAccum[j-1];
	  }
	  for(int j = e-1; j >= 0;j--){
	    incNHopBound[i][index][e] -= incNHopBound[i][index][j];
	  }
	}
      }
      else{
	for(int e = 0; e < nRDim; e++){
	  incNHopBound[i][index][e] = incNHopBulk[i][e];
	}
      }

      boundOpt[0]++;
      p = 0;
      while(boundOpt[p] > 1){
	boundOpt[p] = 0;
	boundOpt[++p]++;
      }
    }

    incNHopBulk[i][nRDim] = 1;
    endHopBound[i][nRDim] = 1;
    endHopBulk[i][nRDim] = 1;

    //Limits in unit cells
    startHopUCBulk[i][nRDim] = 0;
    for(int e = 0; e < nRDim; e++){
      nH = model.getHop(i).getN()[rIndex[e]];
      if(abs(nH) >= l[rIndex[e]]){
	startHopUCBulk[i][nRDim] = 1;
      }
      if(nH < 0){
	startHopUCBulk[i][e] = nH;
	startHopUCBound[i][e] = 0;
	endHopUCBulk[i][e] = l[rIndex[e]] - 1;
	endHopUCBound[i][e] = nH-1;
      }
      else if(nH > 0){
	startHopUCBulk[i][e] = 0;
	startHopUCBound[i][e] = l[rIndex[e]] - nH;
	endHopUCBulk[i][e] = l[rIndex[e]] - 1 - nH;
	endHopUCBound[i][e] = l[rIndex[e]] - 1;
      }
      else{
	startHopUCBulk[i][e] = 0;
	startHopUCBound[i][e] = 0;
	endHopUCBulk[i][e] = l[rIndex[e]] - 1;
	endHopUCBound[i][e] = 0;
      }
    }

    endHopUCBulk[i][nRDim] = 1;
    endHopUCBound[i][nRDim] = 1;

  }

  delete[] boundOpt;
  delete[] start;
  delete[] end;

  //On-site
  for(int i = 0; i < nOnSite; i++){
    int orb = model.getOnSite(i).getNOrb();
    if(nRDim != 0){
      startOnSite[i][0] = mOrb[orb] + lOrb[orb][rIndex[0]]*nOrb/nL;
      endOnSite[i][0] = (l[rIndex[0]]*nOrb*nu[rIndex[0]])/nL - inc[0] + startOnSite[i][0];
      incNOnSite[i][0] = inc[0];
    }
    for(int e = 1; e < nRDim; e++){
      startOnSite[i][e] =  lOrb[orb][rIndex[e]];
      endOnSite[i][e] = l[rIndex[e]]*nu[rIndex[e]] - (inc[e] - startOnSite[i][e]);
      incNOnSite[i][e] = - endOnSite[i][0] + startOnSite[i][0] + inc[e]*lAccum[e-1]*nOrb/nuAccum[e-1];
      for(int j = 1; j < e; j++){
	incNOnSite[i][e] += (startOnSite[i][j] - endOnSite[i][j])*lAccum[j-1]*nOrb/nuAccum[j-1];
      }
      for(int j = e-1; j >= 0;j--){
	incNOnSite[i][e] -= incNOnSite[i][j];
      }
    }
    incNOnSite[i][nRDim] = 1;
    endOnSite[i][nRDim] = 1;
  }
}

cx_mat TBCleanH::H(double * k){
  if(!isUpdated){
    calcAux();
    isUpdated = true;
  }

  int size; 
  if(nRDim != 0){
    size = nOrb*lAccum[nRDim-1];
  }
  else{
    size = nOrb;
  }
  cx_mat res(size, size, fill::zeros);

  complex<double> t;
  int j,n,p,q,r;
  complex<double> phase;
  complex<double> kPhase;
  complex<double> ii(0,1);

  int * i = new int[nRDim + 1];
  int * b = new int[nRDim + 1];
  int * start = new int[nRDim];
  int * end = new int[nRDim + 1];
  end[nRDim] = 1;

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
    for(j = 0; j <= nRDim; j++){
      i[j] = startHopUCBulk[e][j];
    }

    n = flatten2(model.getHop(e).getNOrb1(), i);

    while(i[nRDim] == 0){
      res(n, n + nHopBulk[e]) += t;

      i[0] += 1;
      n += incNHopBulk[e][0];
      p = 0;
      while(i[p] > endHopUCBulk[e][p]){
	i[p] = startHopUCBulk[e][p];
	i[++p] += 1;
	n += incNHopBulk[e][p];
      }
    }

    //Boundary lattice
    for(j = 0; j < nRDim + 1; j++){
      b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
    }

    b[0]++;
    if(nRDim != 0){
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == nRDim){
	  break;
	}
      }
    }

    while(b[nRDim] == 0){

      q = 0;
      phase = 1;
      for(j = 0; j < nRDim; j++){
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
	i[j] = start[j];
      }
      i[nRDim] = 0;

      t = model.getHop(e).getHop()*kPhase*phase;

      n = flatten2(model.getHop(e).getNOrb1(), i);

      while(i[nRDim] == 0){
	res(n, n + nHopBound[e][q]) += t;

	i[0] += 1;
	n += incNHopBound[e][q][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = start[p];
	  i[++p] += 1;
	  n += incNHopBound[e][q][p];
	}
      }


      b[0]++;
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == nRDim){
	  break;
	}
      }
    }

  }
  res = res + res.t();

  //On-site terms

  for(int e = 0; e < nRDim; e++){
    end[e] = l[rIndex[e]] - 1;
  }

  for(int e = 0; e < nOnSite; e++){

    for(j = 0; j < nRDim; j++){
      i[j] = 0;
    }
    n = flatten2(model.getOnSite(e).getNOrb(), i);
    i[nRDim] = 0;

    while(i[nRDim] == 0){
      res(n, n) += model.getOnSite(e).getEn();

      i[0] += 1;
      n += incNOnSite[e][0];
      p = 0;
      while(i[p] > end[p]){
	i[p] = 0;
	i[++p] += 1;
	n += incNOnSite[e][p];
      }
    }
  }

  delete[] i;
  delete[] b;
  delete[] start;
  delete[] end;

  return res;
}

sp_cx_mat TBCleanH::spH(double * k){
  if(!isUpdated){
    calcAux();
    isUpdated = true;
  }

  int size; 

  if(nRDim != 0){
    size = nOrb*lAccum[nRDim-1];
  }
  else{
    size = nOrb;
  }
  sp_cx_mat res(size, size);

  complex<double> t;
  int j,n,p,q,r;
  complex<double> phase;
  complex<double> kPhase;
  complex<double> ii(0,1);

  int * i = new int[nRDim + 1];
  int * b = new int[nRDim + 1];
  int * start = new int[nRDim];
  int * end = new int[nRDim + 1];
  end[nRDim] = 1;

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
    for(j = 0; j <= nRDim; j++){
      i[j] = startHopUCBulk[e][j];
    }

    n = flatten2(model.getHop(e).getNOrb1(), i);

    while(i[nRDim] == 0){
      res(n, n + nHopBulk[e]) += t;

      i[0] += 1;
      n += incNHopBulk[e][0];
      p = 0;
      while(i[p] > endHopUCBulk[e][p]){
	i[p] = startHopUCBulk[e][p];
	i[++p] += 1;
	n += incNHopBulk[e][p];
      }
    }

    //Boundary lattice
    for(j = 0; j < nRDim + 1; j++){
      b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
    }

    b[0]++;
    if(nRDim != 0){
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == nRDim){
	  break;
	}
      }
    }

    while(b[nRDim] == 0){

      q = 0;
      phase = 1;
      for(j = 0; j < nRDim; j++){
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
	i[j] = start[j];
      }
      i[nRDim] = 0;

      t = model.getHop(e).getHop()*kPhase*phase;

      n = flatten2(model.getHop(e).getNOrb1(), i);

      while(i[nRDim] == 0){

	res(n, n + nHopBound[e][q]) += t;

	i[0] += 1;
	n += incNHopBound[e][q][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = start[p];
	  i[++p] += 1;
	  n += incNHopBound[e][q][p];
	}
      }


      b[0]++;
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == nRDim){
	  break;
	}
      }
    }

  }
  res = res + res.t();

  //On-site terms

  for(int e = 0; e < nRDim; e++){
    end[e] = l[rIndex[e]] - 1;
  }

  for(int e = 0; e < nOnSite; e++){

    for(j = 0; j < nRDim; j++){
      i[j] = 0;
    }
    n = flatten2(model.getOnSite(e).getNOrb(), i);
    i[nRDim] = 0;

    while(i[nRDim] == 0){
      res(n, n) += model.getOnSite(e).getEn();

      i[0] += 1;
      n += incNOnSite[e][0];
      p = 0;
      while(i[p] > end[p]){
	i[p] = 0;
	i[++p] += 1;
	n += incNOnSite[e][p];
      }
    }
  }

  delete[] i;
  delete[] b;
  delete[] start;
  delete[] end;

  return res;
}

void TBCleanH::setBlockDim(int bDim){
  if(bDim >= 0 && bDim < nRDim){
    this->bDim = bDim;
  }
}

cx_mat TBCleanH::blockH(int line, int col, double * k){
  if(!isUpdated){
    calcAux();
    isUpdated = true;
  }

  int size; 
  if(nRDim == 0){
    cout << "__PRETTY_FUNCTION__" << " can't divide systems in spatial layers, returning empty matrix" << endl;
    return cx_mat(0,0);
  }
  else{ 
    int nL = nuAccum[0]*nu[rIndex[0]];
    if(bDim == 0){
      size = nOrb/nL;
    }
    else{
      size = l[rIndex[0]]*(nOrb/nuAccum[0]);
      for(int i = 1; i < bDim; i++){
	size *= l[rIndex[i]]*nu[rIndex[i]];
      }
    }
  }
  cx_mat res(size, size, fill::zeros);

  int * startNL = new int[nRDim - bDim];
  int * startL = new int[nRDim - bDim];
  int * diffNL = new int[nRDim - bDim];
  int * diffN = new int[nRDim - bDim];
  int * diffL = new int[nRDim - bDim];

  //This might not be necessary
  bool transpose = false;
  if(col < line){
    int temp = col;
    col = line;
    line = temp;
    transpose = true;
  }

  diffNL[0] = (col-line) % (l[rIndex[bDim]]*nu[rIndex[bDim]]);
  startNL[0] = line % (l[rIndex[nDim]]*nu[rIndex[bDim]]);
  for(int i = 1; i < nRDim-bDim; i++){
    diffNL[i] = (col - line - diffNL[i-1])/((l[rIndex[bDim + i -1]]*nu[rIndex[bDim + i - 1]]));
    startNL[i] = line/((l[rIndex[bDim + i -1]]*nu[rIndex[bDim + i - 1]]));
    if(bDim + i < nRDim){
      diffNL[i] = diffNL[i] % (l[rIndex[bDim + i]]*nu[rIndex[bDim + i]]);
      startNL[i] = startNL[i] % (l[rIndex[bDim + i]]*nu[rIndex[bDim + i]]);
    }
  }

  for(int i = 0; i < nRDim - bDim; i++){
    startL[i] = startNL[i] % nu[rIndex[bDim + i]];
    diffN[i] = (startL[i] + diffNL[i])/nu[rIndex[bDim + i]];
    diffL[i] = diffNL[i] % nu[rIndex[bDim + i]];
    if((startL[i] + diffL[i]) % nu[rIndex[bDim + i]] < startL[i]){
      diffL[i] = -diffL[i];
    }
    cout << "n " << diffN[i] << " l " << diffL[i] << " nl " << diffNL[i] << " " << l[rIndex[bDim + i]] << endl;
  }

  cout << "size " << size << endl;

  bool isHop;
  bool isHopHerm;
  vector<int> h;
  vector<int> hHerm; //Conjugate transpose the matrix comming from these
  int n, l1, l2;
  for(int e = 0; e < nHop; e++){
    isHop = true;
    isHopHerm = true;
    for(int i = 0; i < nRDim - bDim; i++){
      n = model.getHop(e).getN(rIndex[bDim + i]);
      l2 = lOrb[model.getHop(e).getNOrb2()][rIndex[bDim + i]];
      l1 = lOrb[model.getHop(e).getNOrb1()][rIndex[bDim + i]];
      cout << e << " " << i << " " << n << " " << l2 - l1 << endl;
      if(n != diffN[i] || (l2 - l1) != diffL[i] || l1 != startL[i]){
	isHop = false;
      }
      if(n != - diffN[i] || (l2 - l1) != - diffL[i] || l2 != startL[i]){
	isHopHerm = false;
      }
    }
    if(isHop == true){
      h.push_back(e);
    }
    else if(isHopHerm == true){
      hHerm.push_back(e);
    }
  }

  for(int i = 0; i < h.size(); i++){
    cout << "h: " << h[i] << endl;
  }

  for(int i = 0; i < hHerm.size(); i++){
    cout << "hHerm: " << hHerm[i] << endl;
  }

  delete diffL;
  delete diffN;
  delete diffNL;
  delete startNL;
  delete startL;

  //TODO: until here some things can be precalculated to increase efficiency

  complex<double> t;
  int j,p,q,r;
  complex<double> phase;
  complex<double> kPhase;
  complex<double> ii(0,1);

  int * i = new int[nRDim + 1];
  int * b = new int[bDim + 1];
  int * start = new int[bDim];
  int * end = new int[bDim + 1];
  end[bDim] = 1;

  int sub = (col-line)*size;


  //Hopping terms (hermitian conjugate)
  for(int m = 0; m < hHerm.size(); m++){
    int e = hHerm[m];
    kPhase = 1;
    for(j = 0; j < nDim; j++){
      if(bC[j] == 1){
	//Extra phase due to k space hamiltonian
	kPhase *= exp(ii*k[j]*(double)model.getHop(e).getN(j));
      }
    }

    t = model.getHop(e).getHop()*kPhase;

    //Bulk lattice
    for(j = 0; j < bDim; j++){
      i[j] = startHopUCBulk[e][j];
    }
    for(j = bDim; j <= nRDim; j++){
      i[j] = 0;
    }

    n = flatten2(model.getHop(e).getNOrb1(), i);

    while(i[bDim] == 0){
      cout << n << " "  << n-sub << " " << n + nHopBulk[e] << endl;
      res(n - sub, n + nHopBulk[e]) += t;

      i[0] += 1;
      n += incNHopBulk[e][0];
      p = 0;
      while(i[p] > endHopUCBulk[e][p]){
	i[p] = startHopUCBulk[e][p];
	i[++p] += 1;
	n += incNHopBulk[e][p];
      }
    }

    //Boundary lattice
    for(j = 0; j < bDim + 1; j++){
      b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
    }

    b[0]++;
    if(bDim != 0){
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == bDim){
	  break;
	}
      }
    }

    while(b[bDim] == 0){

      q = 0;
      phase = 1;
      for(j = 0; j < bDim; j++){
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
	i[j] = start[j];
      }
      for(j = bDim; j <= nRDim; j++){
	i[j] = 0;
      }

      t = model.getHop(e).getHop()*kPhase*phase;

      n = flatten2(model.getHop(e).getNOrb1(), i);

      while(i[bDim] == 0){

	res(n, n + nHopBound[e][q]) += t;

	i[0] += 1;
	n += incNHopBound[e][q][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = start[p];
	  i[++p] += 1;
	  n += incNHopBound[e][q][p];
	}
      }


      b[0]++;
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == bDim){
	  break;
	}
      }
    }

  }

  res = res.t();

  //Hopping terms 
  for(int m = 0; m < h.size(); m++){
    int e = h[m];
    kPhase = 1;
    for(j = 0; j < nDim; j++){
      if(bC[j] == 1){
	//Extra phase due to k space hamiltonian
	kPhase *= exp(ii*k[j]*(double)model.getHop(e).getN(j));
      }
    }

    t = model.getHop(e).getHop()*kPhase;

    //Bulk lattice
    for(j = 0; j < bDim; j++){
      i[j] = startHopUCBulk[e][j];
    }
    for(j = bDim; j <= nRDim; j++){
      i[j] = 0;
    }

    n = flatten2(model.getHop(e).getNOrb1(), i);

    while(i[bDim] == 0){
      cout << e << " " << n << " " << sub << " " << nHopBulk[e] << endl;
      res(n, n + nHopBulk[e] - sub) += t;

      i[0] += 1;
      n += incNHopBulk[e][0];
      p = 0;
      while(i[p] > endHopUCBulk[e][p]){
	i[p] = startHopUCBulk[e][p];
	i[++p] += 1;
	n += incNHopBulk[e][p];
      }
    }

    //Boundary lattice
    for(j = 0; j < bDim + 1; j++){
      b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
    }

    b[0]++;
    if(bDim != 0){
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == bDim){
	  break;
	}
      }
    }

    while(b[bDim] == 0){

      q = 0;
      phase = 1;
      for(j = 0; j < bDim; j++){
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
	i[j] = start[j];
      }
      for(j = bDim; j <= nRDim; j++){
	i[j] = 0;
      }

      t = model.getHop(e).getHop()*kPhase*phase;

      n = flatten2(model.getHop(e).getNOrb1(), i);

      while(i[bDim] == 0){

	res(n, n + nHopBound[e][q]) += t;

	i[0] += 1;
	n += incNHopBound[e][q][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = start[p];
	  i[++p] += 1;
	  n += incNHopBound[e][q][p];
	}
      }


      b[0]++;
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == bDim){
	  break;
	}
      }
    }

  }

  if(transpose){
    res = res + res.t();
  }

  //On-site terms

  if(sub == 0){
    for(int e = 0; e < nRDim; e++){
      end[e] = l[rIndex[e]] - 1;
    }

    for(int e = 0; e < nOnSite; e++){

      for(j = 0; j < nRDim; j++){
	i[j] = 0;
      }
      n = flatten2(model.getOnSite(e).getNOrb(), i);
      i[nRDim] = 0;

      while(i[bDim] == 0){
	res(n, n) += model.getOnSite(e).getEn();

	i[0] += 1;
	n += incNOnSite[e][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = 0;
	  i[++p] += 1;
	  n += incNOnSite[e][p];
	}
      }
    }
  }

  delete[] i;
  delete[] b;
  delete[] start;
  delete[] end;

  return res;
}

