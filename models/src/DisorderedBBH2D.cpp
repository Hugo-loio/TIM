#include "DisorderedBBH2D.h"
#include "OData.h"

DisorderedBBH2D::DisorderedBBH2D(double t1, double t2, double delta){
  model = new TBModel(2, 4);
  int n1[2] = {0,0};
  int n2[2] = {1,0}; 
  int n3[2] = {0,1};
  //Intracell hoppings
  model->setHop(0,3, n1, t1);
  model->setHop(0,2, n1, t1);
  model->setHop(1,2, n1, -t1);
  model->setHop(1,3, n1, t1);
  //Intercell hoppings
  model->setHop(0,2, n2, t2);
  model->setHop(3,1, n2, t2);
  model->setHop(0,3, n3, t2);
  model->setHop(2,1, n3, -t2);
  //OnSite
  model->setOnSite(0,delta);
  model->setOnSite(1,delta);
  model->setOnSite(2,-delta);
  model->setOnSite(3,-delta);
  ham = new TBDisorderedH(*model);
};

DisorderedBBH2D::~DisorderedBBH2D(){
  delete model;
  delete ham;
};

void DisorderedBBH2D::setIntraHop(double t1){
  model->getHop(0).setHop(t1);
  model->getHop(1).setHop(t1);
  model->getHop(2).setHop(-t1);
  model->getHop(3).setHop(t1);
  delete ham;
  ham = new TBDisorderedH(*model);
}

void DisorderedBBH2D::setInterHop(double t2){
  model->getHop(4).setHop(t2);
  model->getHop(5).setHop(t2);
  model->getHop(6).setHop(t2);
  model->getHop(7).setHop(-t2);
  delete ham;
  ham = new TBDisorderedH(*model);
}

void DisorderedBBH2D::setOnSite(double delta){
  model->getOnSite(0).setEn(delta);
  model->getOnSite(1).setEn(delta);
  model->getOnSite(2).setEn(-delta);
  model->getOnSite(3).setEn(-delta);
}

void DisorderedBBH2D::setProbDisorder(double p){
  ham->setHopDisorderFunction(&probDisorder2D);
  ham->setHopDisorderWeight(p);
}

void DisorderedBBH2D::generateDisorder(){
  ham->generateDisorder();
}

void DisorderedBBH2D::getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled){
  int bC[2] = {0,0};
  ham->setBC(bC);
  ham->setSparse(false);
  int l[2] = {nx, ny};
  ham->setSize(l);
  OData o(argv0, fileName);
  o.chargeDensity(*ham, 4, nOrbFilled, l);
}

double DisorderedBBH2D::getQuadrupoleNestedSupercell(int * l, int * n){
  int bc[2] = {2,2};
  ham->setBC(bc);
  ham->setSize(l);
  ham->setSparse(false);
  Wilson wilson(ham);
  int dir[2] = {0,1};
  int m[2] = {l[0]*l[1]*2, l[1]};


  return fmod(abs(log_det(wilson.nestedWilsonLoopSupercell(n, dir, m)).imag())/(2*M_PI), 1);
}

void DisorderedBBH2D::getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson){
  int nPoints[1] = {100};
  int bC[2] = {2,2};
  ham->setBC(bC);
  ham->setSparse(false);
  int size[2] = {nx,ny};
  ham->setSize(size);
  OData o(argv0, fileName);
  o.supercellWannierBands(*ham, nPoints, 10, dirWilson, nx*ny*2);
}
