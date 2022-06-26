#include "TBDisorderedH2D.h"

TBDisorderedH2D::TBDisorderedH2D(TBModel model) : TBCleanH(model){
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

TBDisorderedH2D::~TBDisorderedH2D(){
  deleteDisArrays();
  delete[] isDisHop;
  delete[] isDisOnSite;
  delete[] indexDisHop;
  delete[] indexDisOnSite;
}

void TBDisorderedH2D::setDisHop(int nDisHop, int * hop){
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

void TBDisorderedH2D::setDisOnSite(int nDisOnSite, int * onSite){
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

void TBDisorderedH2D::setSize(int * l){
  deleteDisArrays();
  TBCleanH::setSize(l);
}

void TBDisorderedH2D::deleteDisArrays(){
  if(disHop != NULL && disOnSite != NULL){
    int e;
    for(int i = 0; i < nDisHop; i++){
      for(e = 0; e < l[0]; e++){
	delete[] disHop[i][e];
      }
      delete[] disHop[i];
    }
    for(int i = 0; i < nDisOnSite; i++){
      for(e = 0; e < l[0]; e++){
	delete[] disOnSite[i][e];
      }
      delete[] disOnSite[i];
    }
    delete[] disHop;
    delete[] disOnSite;
    disHop = NULL;
    disOnSite = NULL;
  }
}

void TBDisorderedH2D::createDisArrays(){
  disHop = new complex<double> ** [nDisHop];
  int e;
  for(int i = 0; i < nDisHop; i++){
    disHop[i] = new complex<double> * [l[0]];
    for(e = 0; e < l[1]; e++){
      disHop[i][e] = new complex<double> [l[1]];
    }
  }

  disOnSite = new complex<double> ** [nDisOnSite];
  for(int i = 0; i < nDisOnSite; i++){
    disOnSite[i] = new complex<double> * [l[0]];
    for(e = 0; e < l[1]; e++){
      disOnSite[i][e] = new complex<double> [l[1]];
    }
  }
}

/*
  if(nRDim != 2 || nDim != 2){
    cout << "2D disordered class doesn't have a 2D TB model in real space. Expect things to go wrong, not my problem." << endl;
  }
  */
