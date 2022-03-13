#include "TBModel.h"
#include <iostream>

using namespace arma;
using namespace std;

TBModel::TBModel(int nDim, int nOrb){
  this->nDim = nDim;
  this->nOrb = nOrb;
}

TBModel::TBModel(const TBModel & model){
  nDim = model.nDim;
  nOrb = model.nOrb;
  hop = model.hop;
  onSite = model.onSite;
}

TBModel::~TBModel(){
}

void TBModel::setHop(int nOrb1, int nOrb2, int * n, complex<double> hop){
  if(nOrb1 >= nOrb || nOrb2 >= nOrb){
    cout << __PRETTY_FUNCTION__ << "\nError: orbital index not compatible with number of orbitals" << endl;
  }
  else{
    this->hop.push_back(Hop(nOrb1, nOrb2, n, hop, nDim));
  }
}

void TBModel::setHop(Hop hop){
  if(hop.getNOrb1() >= nOrb || hop.getNOrb2() >= nOrb){
    cout << __PRETTY_FUNCTION__ << "\nError: orbital index not compatible with number of orbitals" << endl;
  }
  else{
    this->hop.push_back(hop);
  }
}

void TBModel::setOnSite(int nOrb, complex<double> en){
  if(nOrb >= this->nOrb){
    cout << __PRETTY_FUNCTION__ << "\nError: orbital index not compatible with number of orbitals" << endl;
  }
  else{
    this->onSite.push_back(OnSite(nOrb, en));
  }
}

void TBModel::setOnSite(OnSite onSite){
  if(onSite.getNOrb() >= nOrb){
    cout << __PRETTY_FUNCTION__ << "\nError: orbital index not compatible with number of orbitals" << endl;
  }
  else{
    this->onSite.push_back(onSite);
  }
}

Hop & TBModel::getHop(int nHop){
  return hop[nHop];
}

OnSite & TBModel::getOnSite(int nOnSite){
  return onSite[nOnSite];
}
