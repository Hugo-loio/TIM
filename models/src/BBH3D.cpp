#include "BBH3D.h"
#include "OData.h"
#include "MultipoleOp.h"
#include "BoundaryGreenH.h"

BBH3D::BBH3D(double t1, double t2, double delta){
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

BBH3D::~BBH3D(){
  delete model;
  delete ham;
};

void BBH3D::setIntraHop(double t1){
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

void BBH3D::setInterHop(double t2){
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

void BBH3D::setOnSite(double delta){
  model->getOnSite(0).setEn(delta);
  model->getOnSite(1).setEn(delta);
  model->getOnSite(2).setEn(-delta);
  model->getOnSite(3).setEn(-delta);
  model->getOnSite(4).setEn(-delta);
  model->getOnSite(5).setEn(-delta);
  model->getOnSite(6).setEn(delta);
  model->getOnSite(7).setEn(delta);
}

void BBH3D::getBands(char * argv0, string fileName){
  OData o(argv0, fileName);
  int nSegs = 6;
  double ** kI = new double * [nSegs];
  double ** kF = new double * [nSegs];
  //Gamma
  double * Gamma = new double[3];
  Gamma[0] = 0;
  Gamma[1] = 0;
  Gamma[2] = 0;
  //R
  double * R = new double[3];
  R[0] = M_PI;
  R[1] = M_PI;
  R[2] = M_PI;
  //M
  double * M = new double[3];
  M[0] = M_PI;
  M[1] = M_PI;
  M[2] = 0;
  //X
  double * X = new double[3];
  X[0] = 0;
  X[1] = M_PI;
  X[2] = 0;
  //kI
  kI[0] = Gamma;
  kI[1] = X;
  kI[2] = M;
  kI[3] = Gamma;
  kI[4] = R;
  kI[5] = M;
  //kF
  kF[0] = X;
  kF[1] = M;
  kF[2] = Gamma;
  kF[3] = R;
  kF[4] = X;
  kF[5] = R;
  int * segPoints = new int[nSegs];
  for(int i = 0; i < nSegs; i++){
    segPoints[i] = 10;
  }
  int bC[3] = {1,1,1};
  ham->setBC(bC);
  o.eBandsPath(*ham,nSegs,kI,kF,segPoints);
  delete[] kI;
  delete[] kF;
  delete[] segPoints;
}

cx_mat BBH3D::getH(double * k){
  return ham->H(k);
}

void BBH3D::getChargeDensity(char * argv0, string fileName, int * l, int nOrbFilled){
  int bC[3] = {0,0,0};
  ham->setBC(bC);
  ham->setSparse(false);
  ham->setSize(l);
  OData o(argv0, fileName);
  o.chargeDensity(*ham, 8, nOrbFilled, l);
}

void BBH3D::getWannierBands(char * argv0, string fileName, int dir){
  int bC[3] = {1,1,1};
  ham->setBC(bC);
  ham->setSparse(false);
  OData o(argv0, fileName);
  int nVec[2] = {20,20};
  o.wannierBands(*ham, nVec ,10, dir, 4);
}

void BBH3D::getNestedWannierBands(char * argv0, string fileName){
  int bC[3] = {1,1,1};
  ham->setBC(bC);
  ham->setSparse(false);
  OData o(argv0, fileName);
  int dirWilson[2] = {0,1};
  int nWilson[2] = {10,10};

  o.nestedWannierBands(*ham, 100, 2, nWilson, dirWilson, 4);
}

double BBH3D::getOctupoleNested(int nx, int ny, int nz, double * k0){
  int bC[3] = {1,1,1};
  ham->setBC(bC);
  ham->setSparse(false);
  int n[3] = {nx,ny,nz};
  int dir[3] = {0,1,2};
  Wilson wilson(ham);
  double k[3] = {0,0,0};
  if(k0 != NULL){
    wilson.setLoopStart(k0);
    k[0] = k0[0];
    k[1] = k0[1];
    k[2] = k0[2];
  }

  return fmod(abs(log_det(wilson.nestedNestedWilsonLoop(n, dir, 4)).imag())/(2*M_PI), 1);
}

void BBH3D::getSupercellNestedWannierBands(char * argv0, string fileName, int * l, int * n){
  int bC[3] = {2,2,2};
  ham->setBC(bC);
  ham->setSparse(false);
  ham->setSize(l);
  OData o(argv0, fileName);
  int m[2] = {l[0]*l[1]*l[2]*4, l[1]*l[2]*2};
  int dir[2] = {0,1};
  o.supercellNestedWannierBands(*ham, 100, 2, n, dir, m);
}

double BBH3D::getOctupoleNestedSupercell(int * l, int * n){
  int bc[3] = {2,2,2};
  ham->setBC(bc);
  ham->setSize(l);
  ham->setSparse(false);
  Wilson wilson(ham);
  int dir[3] = {0,1,2};
  int m[3] = {l[0]*l[1]*l[2]*4, l[1]*l[2]*2, l[2]};

  return fmod(abs(log_det(wilson.nestedNestedWilsonLoopSupercell(n, dir, m)).imag())/(2*M_PI), 1);
}

double BBH3D::getOctupoleManyBody(int * l){
  int bC[3] = {2,2,2};
  ham->setBC(bC);
  ham->setSize(l);
  ham->setSparse(false);
  MultipoleOp o(ham, l, 3, 8);
  o.setOcc(4*l[0]*l[1]*l[2]);
  return o.octupole(0,1,2);
}

void BBH3D::test(char * argv0){
  int bC[3] = {0,0,0};
  bool layerDir[3] = {true, true, true};

  double k[3] = {M_PI/2,M_PI/2,M_PI/2};
  int l[3] = {2,2,2};
  int order[3] = {2,0,1};
  //int order[3] = {0,2,1};

  ham->setSize(l);
  ham->setBC(bC);
  ham->setSparse(false);
  setLayers(layerDir);
  ham->setOrder(order);

  cx_mat h = ham->H(k);
  //cout << h << endl;

  OData o(argv0, "testH.dat");
  o.matrixWeights(h);

  ham->setBlockDim(2);
  //ham->blockH(0,7);

  cx_mat h2 = ham->blockH(0,0);
  //cout << h2 << endl;
  OData o2(argv0, "testH2.dat");
  o2.matrixWeights(h2);
}

double BBH3D::getBoundQuadrupole(int * l, int dir){
  int bC[3] = {0,0,0};

  int order[3] = {0,1,2};
  bool layerDir[3] = {false,false,false};
  layerDir[dir] = true;
  int a,b;
  if(dir == 0){
    order[0] = 2;
    order[1] = 0;
    order[2] = 1;
    bC[1] = 2;
    bC[2] = 2;
  }
  else if(dir == 1){
    order[0] = 0;
    order[1] = 2;
    order[2] = 1;
    bC[0] = 2;
    bC[2] = 2;
  }
  else{
    bC[0] = 2;
    bC[1] = 2;
  }

  ham->setSize(l);
  ham->setBC(bC);
  ham->setSparse(false);
  setLayers(layerDir);
  ham->setOrder(order);

  BoundaryGreenH green(ham, l[order[0]]*l[order[1]]*4, l[order[2]]*2);
  int lBound[2] = {l[order[0]], l[order[1]]};
  MultipoleOp o(&green, lBound, 2, 4);

  return o.quadrupole(0,1);
}

void BBH3D::setLayers(bool * layerDir){
  int ** layers = new int * [8];
  for(int i = 0; i < 8; i++){
    layers[i] = new int[3];
    for(int e = 0; e < 3; e++){
      layers[i][e] = 0;
    }
  }
  
  if(layerDir[0]){
    layers[0][0] = 1;
    layers[3][0] = 1;
    layers[4][0] = 1;
    layers[7][0] = 1;
  }
  if(layerDir[1]){
    layers[0][1] = 1;
    layers[2][1] = 1;
    layers[4][1] = 1;
    layers[6][1] = 1;
  }
  if(layerDir[2]){
    layers[0][2] = 1;
    layers[1][2] = 1;
    layers[2][2] = 1;
    layers[3][2] = 1;
  }

  ham->setOrbLayer(layers);

  for(int i = 0; i < 4; i++){
    delete[] layers[i];
  }
  delete[] layers;
}

