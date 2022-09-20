#include "Anderson1D.h"
#include "OData.h"
#include "LocalizationProps.h"

Anderson1D::Anderson1D(double t){
  model = new TBModel(1, 1);
  int n[1] = {1};
  model->setHop(0, 0, n, t);
  model->setOnSite(0,0);
  ham = new DisorderedAnderson1D(*model);
}

Anderson1D::~Anderson1D(){
  delete ham;
  delete model;
}

void Anderson1D::setHop(double t){
  model->getHop(0).setHop(t);
  delete ham;
  ham = new DisorderedAnderson1D(*model);
}

void Anderson1D::setSize(int * l){
  ham->setSize(l);
}


void Anderson1D::setW(double w){
  ham->setWeight(w);
}

void Anderson1D::generateDisorder(){
  ham->generateDisorder();
}

double Anderson1D::ipr(int nStates, double en){
  int bC[2] = {2};
  ham->setBC(bC);
  ham->setSparse(true);
  int vol = ham->getSize()[0];

  LocalizationProps loc(ham);
  return loc.ipr(vol, 1, nStates, en);
}

cx_mat Anderson1D::getHam(){
  int bC[1] = {0};
  bool layers[1] = {true};
  int order[1] = {0};
  ham->setBC(bC);
  ham->setSparse(false);
  //setLayers(layers);
  //ham->setOrder(order);
  return ham->H(NULL);
}

double Anderson1D::getTMM(int qrIt, double en){
  int bC[2] = {0};
  int lVec[2] = {3};
  setSize(lVec);
  ham->setBC(bC);
  generateDisorder();

  LocalizationProps loc(ham);
  return loc.tmm(2, qrIt, en);
}

void Anderson1D::test(char * argv0){
  int bC[2] = {0};
  int lVec[2] = {3};
  setSize(lVec);
  ham->setBC(bC);
  generateDisorder();

  cout << ham->blockH(0,1) << endl;
  cout << ham->blockH(1,2) << endl;
  //cout << ham->blockH(2,3) << endl;
  cout << ham->H() << endl;
}
