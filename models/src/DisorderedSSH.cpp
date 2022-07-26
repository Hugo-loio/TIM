#include "DisorderedSSH.h"
#include "MultipoleOp.h"
#include "OData.h"
#include "LocalizationStats.h"

DisorderedSSH::DisorderedSSH(double t1, double t2){
  model = new TBModel(1, 2);
  int n1[1] = {0};
  int n2[1] = {1};
  model->setHop(0,1,n1, t1);
  model->setHop(0,1,n2, t2);
  ham = new DisorderedHopH1D(*model);
}

DisorderedSSH::~DisorderedSSH(){
  delete ham;
  delete model;
}

void DisorderedSSH::setIntraHop(double t1){
  model->getHop(0).setHop(t1);
  delete ham;
  ham = new DisorderedHopH1D(*model);
}

void DisorderedSSH::setInterHop(double t2){
  model->getHop(1).setHop(t2);
  delete ham;
  ham = new DisorderedHopH1D(*model);
}

void DisorderedSSH::setW(double w){
  ham->setWeight(w);
}

void DisorderedSSH::generateDisorder(){
  ham->generateDisorder();
}

double DisorderedSSH::ipr(int nStates){
  int bC[2] = {2};
  ham->setBC(bC);
  ham->setSparse(true);
  int vol = ham->getSize()[0];

  LocalizationStats loc(ham);
  return loc.ipr(vol, 2, nStates);
}

double DisorderedSSH::polarization(){
  int bC[1] = {2};
  int l[1] = {ham->getSize()[0]};
  ham->setBC(bC);
  ham->setSparse(false);
  MultipoleOp p(ham, l, 1, 2);
  p.setOcc(l[0]);
  return p.polarization(0);
}

