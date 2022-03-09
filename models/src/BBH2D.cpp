#include "BBH2D.h"
#include "OData.h"

BBH2D::BBH2D(double t1, double t2){
  model = new TBModel(2, 4);
  int n1[2] = {0,0};
  int n2[2] = {1,0}; int n3[2] = {0,1};
  //Intracell hoppings
  model->setHop(0,1, n1, t1);
  model->setHop(0,2, n1, -t1);
  model->setHop(1,3, n1, t1);
  model->setHop(2,3, n1, t1);
  //Intercell hoppings
  model->setHop(1,0, n2, t2);
  model->setHop(3,2, n2, t2);
  model->setHop(2,0, n3, -t2);
  model->setHop(3,1, n3, t2);
  ham = new TBCleanH(*model);
};

BBH2D::~BBH2D(){
  delete model;
  delete ham;
};

void BBH2D::setIntraHop(double t1){
  model->getHop(0).setHop(t1);
  model->getHop(1).setHop(-t1);
  model->getHop(2).setHop(t1);
  model->getHop(3).setHop(t1);
  delete ham;
  ham = new TBCleanH(*model);
}

void BBH2D::setInterHop(double t2){
  model->getHop(4).setHop(t2);
  model->getHop(5).setHop(t2);
  model->getHop(6).setHop(-t2);
  model->getHop(7).setHop(t2);
  delete ham;
  ham = new TBCleanH(*model);
}

void BBH2D::getBands(char * argv0, string fileName, int nx, int ny){
  OData o(argv0, fileName);
  o.eBands3D(*ham, nx, ny);
}

double BBH2D::berryPhase(int n, int dir, double * k0){
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

void BBH2D::getWannierBands(char * argv0, string fileName, int dir, int n){
  if(boundH != NULL){
    delete boundH;
  }
  int bC[2] = {1,1};
  ham->setBC(bC);
  ham->setSparse(true);
  boundH = new BoundaryWilsonH(ham, dir, n, 2);
  OData o(argv0, fileName);
  double k[2] = {0,0};
  if(dir == 0){
    o.eBands2D(*boundH, 100, 1, k);
  }
  else{
    o.eBands2D(*boundH, 100, 0, k);
  }
}

void BBH2D::test(){
}
