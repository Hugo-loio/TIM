#include "SSH2D.h"
#include "OData.h"

SSH2D::SSH2D(double t1, double t2){
  model = new TBModel(2, 4);
  int n1[2] = {0,0};
  int n2[2] = {1,0};
  int n3[2] = {0,1};
  //Intracell hoppings
  model->setHop(0,1, n1, t1);
  model->setHop(0,2, n1, t1);
  model->setHop(1,3, n1, t1);
  model->setHop(2,3, n1, t1);
  //Intercell hoppings
  model->setHop(1,0, n2, t2);
  model->setHop(3,2, n2, t2);
  model->setHop(2,0, n3, t2);
  model->setHop(3,1, n3, t2);
  ham = new TBCleanH(*model);
};

SSH2D::~SSH2D(){
  delete model;
  delete ham;
};

void SSH2D::setIntraHop(double t1){
  for(int i = 0; i < 4; i++){
    model->getHop(i).setHop(t1);
  }
  delete ham;
  ham = new TBCleanH(*model);
}

void SSH2D::setInterHop(double t2){
  for(int i = 4; i < 8; i++){
    model->getHop(i).setHop(t2);
  }
  delete ham;
  ham = new TBCleanH(*model);
}

void SSH2D::getBands(char * argv0, string fileName, int nx, int ny){
  OData o(argv0, fileName);
  o.eBands3D(*ham, nx, ny);
}

double SSH2D::berryPhase(int n, int dir, double * k0){
  int bC[2] = {1,1};
  ham->setBC(bC);
  ham->setSparse(true);
  Wilson wilson(ham);
  wilson.setLoopDir(dir);
  if(k0 != NULL){
    wilson.setLoopStart(k0);
  }
  return wilson.berryPhase(n, 1);
}

double SSH2D::berryPhaseSupercell(int n, int dir, int * bC, int * lVec, double * k){
  ham->setBC(bC);
  ham->setSize(lVec);
  ham->setSparse(false);
  double theta[2] = {0,0};
  ham->setTwists(theta);
  Wilson wilson(ham);
  wilson.setLoopDir(dir);
  int m = lVec[dir];
  for(int i = 0; i < ham->getNDim(); i++){
    if(i != dir && bC[i] != 1){
      m*=lVec[i];
    }
  }
  return wilson.berryPhaseSupercell(n,m,k);
}

cx_mat SSH2D::getH(int * l, double * k){
  double theta[2] = {0,0};
  ham->setTwists(theta);
  ham->setSize(l);
  return ham->H(k);
}
