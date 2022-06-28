#include "DisorderedSOTAI.h"
#include "MultipoleOp.h"
#include "BoundaryGreenH.h"
#include "Wilson.h"
#include "OData.h"

DisorderedSOTAI::DisorderedSOTAI(double m,double delta){
  this->m = m;
  this->delta = delta;
  model = new TBModel(2,4);
  int n1[2] = {0,0};
  int n2[2] = {1,0}; 
  int n3[2] = {0,1};
  complex<double> ii(0,1);
  //Intracell hoppings
  model->setHop(0,1, n1, -ii*m);
  model->setHop(0,2, n1, -ii*m);
  model->setHop(1,3, n1, ii*m);
  model->setHop(2,3, n1, -ii*m);
  //Intercell hoppings
  model->setHop(2,0, n2, 1);
  model->setHop(3,1, n2, -1);
  model->setHop(1,0, n3, 1);
  model->setHop(3,2, n3, 1);
  //OnSite
  if(delta != 0){
    model->setOnSite(0,-delta);
    model->setOnSite(1,delta);
    model->setOnSite(2,delta);
    model->setOnSite(3,-delta);
  }
  ham = new DisorderedHopH2D(*model);
}

DisorderedSOTAI::~DisorderedSOTAI(){
  delete model;
  delete ham;
}

void DisorderedSOTAI::setM(double m){
  this->m = m;
  this->~DisorderedSOTAI();
  DisorderedSOTAI(m,delta);
}

void DisorderedSOTAI::setOnSite(double delta){
  this->delta = delta;
  this->~DisorderedSOTAI();
  DisorderedSOTAI(m,delta);
}


void DisorderedSOTAI::setW(double w){
  ham->setWeight(w);
}

void DisorderedSOTAI::generateDisorder(){
  ham->generateDisorder();
}

void DisorderedSOTAI::getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled){
  int bC[2] = {0,0};
  ham->setBC(bC);
  ham->setSparse(false);
  int l[2] = {nx, ny};
  ham->setSize(l);
  OData o(argv0, fileName);
  o.chargeDensity(*ham, 4, nOrbFilled, l);
}

double DisorderedSOTAI::getQuadrupoleNestedSupercell(int * l, int * n){
  int bc[2] = {2,2};
  ham->setBC(bc);
  ham->setSize(l);
  ham->setSparse(false);
  Wilson wilson(ham);
  int dir[2] = {0,1};
  int m[2] = {l[0]*l[1]*2, l[1]};


  return fmod(abs(log_det(wilson.nestedWilsonLoopSupercell(n, dir, m)).imag())/(2*M_PI), 1);
}

void DisorderedSOTAI::getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson){
  int nPoints[1] = {100};
  int bC[2] = {2,2};
  ham->setBC(bC);
  ham->setSparse(false);
  int size[2] = {nx,ny};
  ham->setSize(size);
  OData o(argv0, fileName);
  o.supercellWannierBands(*ham, nPoints, 10, dirWilson, nx*ny*2);
}

double DisorderedSOTAI::getQuadrupoleManyBody(){
  int bC[2] = {2,2};
  int l[2] = {ham->getSize()[0], ham->getSize()[1]};
  ham->setBC(bC);
  ham->setSparse(false);
  MultipoleOp q(ham, l, 2, 4);
  q.setOcc(2*l[0]*l[1]);
  return q.quadrupole(0,1);
}

cx_mat DisorderedSOTAI::getHam(int * l){
  int bC[2] = {0,0};
  ham->setBC(bC);
  ham->setSize(l);
  ham->setSparse(false);
  return ham->H(NULL);
}

double DisorderedSOTAI::getBoundPolarization(int dir){
  int order[2] = {0,1};
  int bC[2] = {2,0};
  if(dir == 1){
    bC[0] = 0;
    bC[1] = 2;
    order[0] = 1;
    order[1] = 0;
  }
  bool layerDir[2] = {true, true};
  setLayers(layerDir);
  ham->setOrder(order);
  ham->setBC(bC);
  ham->setSparse(false);
  int l[2] = {ham->getSize()[0], ham->getSize()[1]};
  BoundaryGreenH green(ham, l[order[0]]*2, l[order[1]]*2);

  int lBound[1] = {l[order[0]]};
  MultipoleOp o(&green, lBound, 1, 2);

  return o.polarization(0);
}


void DisorderedSOTAI::test(char * argv0){
  int bC[2] = {2,2};
  int l[2] = {20,20};
  int order[2] = {0,1};

  bool layerDir[2] = {true, true};
  setLayers(layerDir);
  ham->setSize(l);
  ham->setBC(bC);
  ham->setSparse(false);
  ham->setOrder(order);

  cx_mat h = ham->H();

  cx_mat h2 = ham->blockH(0,1);
  OData o2(argv0, "testH2.dat");
  o2.matrixWeights(h2);

  setW(8);
  generateDisorder();

  h = ham->H();

  OData o(argv0, "testH.dat");
  o.matrixWeights(h);
}

void DisorderedSOTAI::setLayers(bool * layerDir){
  int ** layers = new int * [4];
  for(int i = 0; i < 4; i++){
    layers[i] = new int[2];
    for(int e = 0; e < 2; e++){
      layers[i][e] = 0;
    }
  }
  if(layerDir[0]){
    layers[2][0] = 1;
    layers[3][0] = 1;
  }
  if(layerDir[1]){
    layers[1][1] = 1;
    layers[3][1] = 1;
  }

  ham->setOrbLayer(layers);

  for(int i = 0; i < 4; i++){
    delete[] layers[i];
  }
  delete[] layers;
}
