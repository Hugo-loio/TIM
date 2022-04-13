#include "BBH2D.h"
#include "OData.h"

BBH2D::BBH2D(double t1, double t2, double delta){
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
  ham = new TBCleanH(*model);
};

BBH2D::~BBH2D(){
  delete model;
  delete ham;
};

void BBH2D::setIntraHop(double t1){
  model->getHop(0).setHop(t1);
  model->getHop(1).setHop(t1);
  model->getHop(2).setHop(-t1);
  model->getHop(3).setHop(t1);
  delete ham;
  ham = new TBCleanH(*model);
}

void BBH2D::setInterHop(double t2){
  model->getHop(4).setHop(t2);
  model->getHop(5).setHop(t2);
  model->getHop(6).setHop(t2);
  model->getHop(7).setHop(-t2);
  delete ham;
  ham = new TBCleanH(*model);
}

void BBH2D::setOnSite(double delta){
  model->getOnSite(0).setEn(delta);
  model->getOnSite(1).setEn(delta);
  model->getOnSite(2).setEn(-delta);
  model->getOnSite(3).setEn(-delta);
}

void BBH2D::getBands(char * argv0, string fileName, int nx, int ny){
  OData o(argv0, fileName);
  o.eBands3D(*ham, nx, ny);
}

double BBH2D::berryPhase(int n, int dir, double * k0){
  int bC[2] = {1,1};
  ham->setBC(bC);
  ham->setSparse(false);
  Wilson wilson(ham);
  wilson.setLoopDir(dir);
  if(k0 != NULL){
    wilson.setLoopStart(k0);
  }
  return wilson.berryPhase(n, 2);
}

void BBH2D::getWannierBands(char * argv0, string fileName, int dir){
  int bC[2] = {1,1};
  ham->setBC(bC);
  ham->setSparse(false);
  OData o(argv0, fileName);
  int nVec[1] = {500};
  o.wannierBands(*ham, nVec ,50, dir, 2);
}

void BBH2D::getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled){
  int bC[2] = {0,0};
  ham->setBC(bC);
  ham->setSparse(false);
  int l[2] = {nx, ny};
  ham->setSize(l);
  OData o(argv0, fileName);
  o.chargeDensity(*ham, 4, nOrbFilled, l);
}

double BBH2D::getQuadrupoleNested(int nx, int ny, double * k0){
  int bC[2] = {1,1};
  ham->setBC(bC);
  ham->setSparse(false);
  int n[2] = {nx,ny};
  int dir[2] = {0,1};
  Wilson wilson(ham);
  double k[2] = {0,0};
  if(k0 != NULL){
    wilson.setLoopStart(k0);
    k[0] = k0[0];
    k[1] = k0[1];
  }
  double quad = 0;
  double nw;
  for(int i = 0; i < n[dir[0]]; i++){
    wilson.setLoopStart(k);
    nw = abs(log_det(wilson.nestedWilsonLoop(n, dir, 2)).imag())/(2*M_PI);
    quad += nw;
    //    cout << "i: " << i << " k[dir[0]]: " << k[dir[0]] << " nested:" << nw << endl;
    k[dir[0]] += 2*M_PI/(double)n[dir[0]];
  }
  return quad/(double)n[dir[0]];
}
