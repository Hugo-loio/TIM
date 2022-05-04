#include "TBDisorderedH.h"
#include <iostream>

using namespace arma;
using namespace std;

TBDisorderedH::TBDisorderedH(TBModel model) : TBCleanH(model){
  hopBulk = new complex<double> * [vBulk];
  hopBound = new complex<double> * [vBound];
  onSiteBulk = new complex<double> * [vBulk];
  onSiteBound = new complex<double> * [vBound];

  for(int i = 0; i < vBulk; i++){
    hopBulk[i] = new complex<double>[nHop];
    onSiteBulk[i] = new complex<double>[nHop];
  }

  for(int i = 0; i < vBound; i++){
    hopBound[i] = new complex<double>[nOnSite];
    onSiteBound[i] = new complex<double>[nOnSite];
  }

  generateDisorder();
}

TBDisorderedH::~TBDisorderedH(){
  delete_arrays();
}

void TBDisorderedH::delete_arrays(){
  for(int i = 0; i < vBulk; i++){
    delete[] hopBulk[i];
    delete[] onSiteBulk[i];
  }
  for(int i = 0; i < vBound; i++){
    delete[] hopBound[i];
    delete[] onSiteBound[i];
  }
  delete[] hopBulk;
  delete[] hopBound;
  delete[] onSiteBulk;
  delete[] onSiteBound;
}

void TBDisorderedH::resize(){
  hopBulk = new complex<double> * [vBulk];
  hopBound = new complex<double> * [vBound];
  onSiteBulk = new complex<double> * [vBulk];
  onSiteBound = new complex<double> * [vBound];

  for(int i = 0; i < vBulk; i++){
    hopBulk[i] = new complex<double>[nHop];
    onSiteBulk[i] = new complex<double>[nOnSite];
  }

  for(int i = 0; i < vBound; i++){
    hopBound[i] = new complex<double>[nHop];
    onSiteBound[i] = new complex<double>[nOnSite];
  }
}

void TBDisorderedH::setBC(int * bC){
  delete_arrays();
  TBCleanH::setBC(bC);
  resize();
  generateDisorder();
}

void TBDisorderedH::setSize(int * l){
  delete_arrays();
  TBCleanH::setSize(l);
  resize();
  generateDisorder();
}

void TBDisorderedH::generateDisorder(){
  int e;
  if(hopFunc != NULL){
    for(int i = 0; i < vBulk; i++){
      for(e = 0; e < nHop; e++){
	hopBulk[i][e] = hopFunc(model.getHop(e), wHop);
      }
    }
    for(int i = 0; i < vBound; i++){
      for(e = 0; e < nHop; e++){
	hopBound[i][e] = hopFunc(model.getHop(e), wHop);
      }
    }
  }
  else{
    for(int i = 0; i < vBulk; i++){
      for(e = 0; e < nHop; e++){
	hopBulk[i][e] = model.getHop(e).getHop();
      }
    }
    for(int i = 0; i < vBound; i++){
      for(e = 0; e < nHop; e++){
	hopBound[i][e] = model.getHop(e).getHop();
      }
    }
  }

  if(onSiteFunc != NULL){
    for(int i = 0; i < vBulk; i++){
      for(e = 0; e < nOnSite; e++){
	onSiteBulk[i][e] = onSiteFunc(model.getOnSite(e), wOnSite);
      }
    }
    for(int i = 0; i < vBound; i++){
      for(e = 0; e < nOnSite; e++){
	onSiteBound[i][e] = onSiteFunc(model.getOnSite(e), wOnSite);
      }
    }
  }
  else{
    for(int i = 0; i < vBulk; i++){
      for(e = 0; e < nOnSite; e++){
	onSiteBulk[i][e] = model.getOnSite(e).getEn();
      }
    }
    for(int i = 0; i < vBound; i++){
      for(e = 0; e < nOnSite; e++){
	onSiteBound[i][e] = model.getOnSite(e).getEn();
      }
    }
  }
}

cx_mat TBDisorderedH::H(double *k){
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
      res(model.getHop(e).getNOrb1() + nOrb*nBulk[i][nRDim], model.getHop(e).getNOrb2() + nOrb*(nBulk[i][nRDim] + incN)) += hopBulk[i][e]*kPhase;
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
      res(model.getHop(e).getNOrb1() + nOrb*nBound[i][nRDim], model.getHop(e).getNOrb2() + nOrb*getN(newNBound)) += kPhase*phase*hopBound[i][e];
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
      res(m,m) += onSiteBulk[i][e];
    }
    //Go through boundary mesh
    for(i = 0; i < vBound; i++){
      m = model.getOnSite(e).getNOrb() + nOrb*nBound[i][nRDim];
      res(m,m) += onSiteBound[i][e];
    }
  }

  return res;
}

sp_cx_mat TBDisorderedH::spH(double * k){
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
      res(model.getHop(e).getNOrb1() + nOrb*nBulk[i][nRDim], model.getHop(e).getNOrb2() + nOrb*(nBulk[i][nRDim] + incN)) += hopBulk[i][e]*kPhase;
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
      res(model.getHop(e).getNOrb1() + nOrb*nBound[i][nRDim], model.getHop(e).getNOrb2() + nOrb*getN(newNBound)) += kPhase*phase*hopBound[i][e];
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
      res(m,m) += onSiteBulk[i][e];
    }
    //Go through boundary mesh
    for(i = 0; i < vBound; i++){
      m = model.getOnSite(e).getNOrb() + nOrb*nBound[i][nRDim];
      res(m,m) += onSiteBound[i][e];
    }
  }

  return res;
}