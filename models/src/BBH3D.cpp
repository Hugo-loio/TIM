#include "BBH3D.h"
#include "OData.h"

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
  int dirWilson[2] = {1,0};
  int nWilson[2] = {10,10};

  o.nestedWannierBands(*ham, 100, 2, nWilson, dirWilson, 4);
}

/*
   double BBH3D::berryPhase(int n, int dir, double * k0){
   int bC[2] = {1,1};
   ham->setBC(bC);
   ham->setSparse(false);
   Wilson wilson(ham);
   wilson.setLoopDir(dir);
   if(k0 != NULL){
   wilson.setLoopStart(k0);
   }
   return wilson.berryPhase(n, 1);
   }
   */
