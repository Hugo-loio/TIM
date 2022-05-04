#include "DisorderedBBH3D.h"
#include "OData.h"

DisorderedBBH3D::DisorderedBBH3D(double t1, double t2, double delta){
  model = new TBModel(3, 8);
  int n0[3] = {0,0,0};
  int nX[3] = {1,0,0}; 
  int nY[3] = {0,1,0};
  int nZ[3] = {0,0,1};
  //Intracell hoppings
  model->setHop(0,2, n0, t1);
  model->setHop(0,3, n0, t1);
  model->setHop(0,4, n0, t1);
  model->setHop(1,2, n0, -t1);
  model->setHop(1,3, n0, t1);
  model->setHop(1,5, n0, t1);
  model->setHop(2,6, n0, t1);
  model->setHop(3,7, n0, t1);
  model->setHop(4,6, n0, -t1);
  model->setHop(4,7, n0, -t1);
  model->setHop(5,7, n0, -t1);
  model->setHop(5,6, n0, t1);
  //Intercell hoppings
  model->setHop(0,2, nX, t2);
  model->setHop(3,1, nX, t2);
  model->setHop(4,6, nX, -t2);
  model->setHop(7,5, nX, -t2);
  model->setHop(0,3, nY, t2);
  model->setHop(2,1, nY, -t2);
  model->setHop(4,7, nY, -t2);
  model->setHop(6,5, nY, t2);
  model->setHop(0,4, nZ, t2);
  model->setHop(1,5, nZ, t2);
  model->setHop(2,6, nZ, t2);
  model->setHop(3,7, nZ, t2);
  //OnSite
  model->setOnSite(0,delta);
  model->setOnSite(1,delta);
  model->setOnSite(2,-delta);
  model->setOnSite(3,-delta);
  model->setOnSite(4,-delta);
  model->setOnSite(5,-delta);
  model->setOnSite(6,delta);
  model->setOnSite(7,delta);
  ham = new TBCleanH(*model);
};

DisorderedBBH3D::~DisorderedBBH3D(){
  delete model;
  delete ham;
};

void DisorderedBBH3D::setIntraHop(double t1){
  for(int i = 0; i < 12; i++){
    if(i == 3 || i == 8 || i == 9 || i == 10){
      model->getHop(i).setHop(-t1);
    }
    else{
      model->getHop(i).setHop(t1);
    }
  }
  delete ham;
  ham = new TBCleanH(*model);
}

void DisorderedBBH3D::setInterHop(double t2){
  for(int i = 12; i < 24; i++){
    if(i == 14 || i == 15 || i == 17 || i == 18){
      model->getHop(i).setHop(-t2);
    }
    else{
      model->getHop(i).setHop(t2);
    }
  }
  delete ham;
  ham = new TBCleanH(*model);
}

void DisorderedBBH3D::setOnSite(double delta){
  model->getOnSite(0).setEn(delta);
  model->getOnSite(1).setEn(delta);
  model->getOnSite(2).setEn(-delta);
  model->getOnSite(3).setEn(-delta);
  model->getOnSite(4).setEn(-delta);
  model->getOnSite(5).setEn(-delta);
  model->getOnSite(6).setEn(delta);
  model->getOnSite(7).setEn(delta);
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