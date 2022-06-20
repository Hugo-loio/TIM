#include "BBH2D.h"
#include "OData.h"
#include "MultipoleOp.h"
#include "BoundaryGreenH.h"

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

double BBH2D::getQuadrupoleNestedSupercell(int * l, int * n){
  int bc[2] = {2,2};
  ham->setBC(bc);
  ham->setSize(l);
  ham->setSparse(false);
  Wilson wilson(ham);
  int dir[2] = {0,1};
  int m[2] = {l[0]*l[1]*2, l[1]};

  return fmod(abs(log_det(wilson.nestedWilsonLoopSupercell(n, dir, m)).imag())/(2*M_PI), 1);
}

void BBH2D::getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson){
  int nPoints[1] = {100};
  int bC[2] = {2,2};
  ham->setBC(bC);
  ham->setSparse(false);
  int size[2] = {nx,ny};
  ham->setSize(size);
  OData o(argv0, fileName);
  o.supercellWannierBands(*ham, nPoints, 10, dirWilson, nx*ny*2);
}

void BBH2D::test(char * argv0){
  int bC[2] = {0,0};
  //int layers[4][2] = {{1,1},{0,0},{0,1},{1,0}};
  int ** layers = new int * [4];
  for(int i = 0; i < 4; i++){
    layers[i] = new int[2];
  }
  layers[0][0] = 1;
  layers[0][1] = 1;
  layers[1][0] = 0;
  layers[1][1] = 0;
  layers[2][0] = 0;
  layers[2][1] = 1;
  layers[3][0] = 1;
  layers[3][1] = 0;
  double k[2] = {M_PI/2,-M_PI/2};
  int l[2] = {2,2};
  int order[2] = {1,0};

  ham->setSize(l);
  ham->setBC(bC);
  ham->setSparse(false);
  ham->setOrbLayer(layers);
  ham->setOrder(order);

  cx_mat h = ham->H(k);
  cout << h << endl;

  OData o(argv0, "testH.dat");
  o.matrixWeights(h);

  cx_mat h2 = ham->blockH(3,1);
  cout << h2 << endl;
  OData o2(argv0, "testH2.dat");
  o2.matrixWeights(h2);

  for(int i = 0; i < 4; i++){
    delete[] layers[i];
  }
  delete[] layers;
}

double BBH2D::getQuadrupoleManyBody(int * l){
  int bC[2] = {2,2};
  ham->setBC(bC);
  ham->setSize(l);
  ham->setSparse(false);
  MultipoleOp q(ham, l, 2, 4);
  q.setOcc(2*l[0]*l[1]);
  return q.quadrupole(0,1);
}

double BBH2D::getBoundPolarization(int * l, int dir){
  int ** layers = new int * [4];
  for(int i = 0; i < 4; i++){
    layers[i] = new int[2];
  }
  layers[0][0] = 1;
  layers[0][1] = 1;
  layers[1][0] = 0;
  layers[1][1] = 0;
  layers[2][0] = 0;
  layers[2][1] = 1;
  layers[3][0] = 1;
  layers[3][1] = 0;
  int order[2] = {0,1};
  int bC[2] = {2,0};
  if(dir == 1){
    bC[0] = 0;
    bC[1] = 2;
    order[0] = 1;
    order[1] = 0;
  }
  ham->setOrbLayer(layers);
  ham->setOrder(order);
  ham->setBC(bC);
  ham->setSize(l);
  ham->setSparse(false);
  BoundaryGreenH green(ham, l[order[0]]*2, l[order[1]]*2);
  
  int lBound[1] = {l[order[0]]};
  MultipoleOp o(&green, lBound, 1, 2);
  return o.polarization(0);
}
