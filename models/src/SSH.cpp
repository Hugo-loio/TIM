#include "SSH.h"
#include "OData.h"

SSH::SSH(double t1, double t2){
  model = new TBModel(1, 2);
  int n1[1] = {0};
  int n2[1] = {1};
  model->setHop(0,1,n1, t1);
  model->setHop(0,1,n2, t2);
  ham = new TBCleanH(*model);
};

SSH::~SSH(){
  delete model;
  delete ham;
};

void SSH::setIntraHop(double t1){
  model->getHop(0).setHop(t1);
  delete ham;
  ham = new TBCleanH(*model);
}

void SSH::setInterHop(double t2){
  model->getHop(1).setHop(t2);
  delete ham;
  ham = new TBCleanH(*model);
}

void SSH::getBands(char * argv0, string fileName, int n){
  OData o(argv0, fileName);
  o.eBands2D(*ham, n);
}

double SSH::berryPhase(int n){
  int bC[1] = {1};
  ham->setBC(bC);
  ham->setSparse(false);
  Wilson wilson(ham);
  return wilson.berryPhase(n, 1);
}

double SSH::berryPhaseSupercell(int n, int l){
  int bC[1] = {2};
  int lVec[1] = {l};
  ham->setBC(bC);
  ham->setSize(lVec);
  ham->setSparse(true);
  Wilson wilson(ham);
  return wilson.berryPhaseSupercell(n,1*l);
}
