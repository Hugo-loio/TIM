#include "Anderson3D.h"
#include "OData.h"
#include "LocalizationProps.h"

Anderson3D::Anderson3D(double t){
  model = new TBModel(3, 1);
  int nX[3] = {1,0,0};
  int nY[3] = {0,1,0};
  int nZ[3] = {0,0,1};
  model->setHop(0, 0, nX, t);
  model->setHop(0, 0, nY, t);
  model->setHop(0, 0, nZ, t);
  model->setOnSite(0,0);
  ham = new DisorderedAnderson3D(*model);
}

Anderson3D::~Anderson3D(){
  delete ham;
  delete model;
}

void Anderson3D::setHop(double t){
  model->getHop(0).setHop(t);
  model->getHop(1).setHop(t);
  model->getHop(2).setHop(t);
  delete ham;
  ham = new DisorderedAnderson3D(*model);
}

void Anderson3D::setSize(int * l){
  ham->setSize(l);
}


void Anderson3D::setW(double w){
  ham->setWeight(w);
}

void Anderson3D::generateDisorder(){
  ham->generateDisorder();
}

double Anderson3D::ipr(int nStates, double en){
  int bC[3] = {2,2,2};
  ham->setBC(bC);
  ham->setSparse(true);
  int vol = ham->getSize()[0]*ham->getSize()[1]*ham->getSize()[2];

  LocalizationProps loc(ham);
  return loc.ipr(vol, 1, nStates, en);
}

cx_mat Anderson3D::getHam(){
  int bC[3] = {0,0,0};
  int order[3] = {0, 1, 2};
  ham->setBC(bC);
  ham->setSparse(false);
  ham->setOrder(order);
  return ham->H(NULL);
}

vector<double> Anderson3D::getTMM(int qrIt, double en){
  int bC[1] = {0};
  ham->setBC(bC);
  generateDisorder();

  LocalizationProps loc(ham);
  return loc.tmm(ham->getSize()[0], qrIt, en);
}

void Anderson3D::test(char * argv0){
  int bC[3] = {0,0,0};
  int order[3] = {0, 1, 2};
  ham->setBC(bC);
  ham->setSparse(false);
  ham->setOrder(order);
  int lVec[3] = {2,2,3};
  setSize(lVec);
  generateDisorder();

  cout << ham->blockH(0,1) << endl;
  cout << ham->blockH(1,2) << endl;
  cout << ham->H() << endl;
}
