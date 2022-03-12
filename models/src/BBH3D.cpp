#include "BBH3D.h"
#include "OData.h"

BBH3D::BBH3D(double t1, double t2){
  model = new TBModel(3, 8);
  int n0[3] = {0,0,0};
  int nX[3] = {1,0,0}; 
  int nY[3] = {0,1,0};
  int nZ[3] = {0,0,1};
  //Intracell hoppings
  model->setHop(0,1, n0, -t1);
  model->setHop(0,2, n0, t1);
  model->setHop(0,4, n0, t1);
  model->setHop(1,3, n0, -t1);
  model->setHop(1,5, n0, t1);
  model->setHop(2,3, n0, -t1);
  model->setHop(2,6, n0, t1);
  model->setHop(3,7, n0, t1);
  model->setHop(4,5, n0, t1);
  model->setHop(4,6, n0, -t1);
  model->setHop(5,7, n0, t1);
  model->setHop(6,7, n0, t1);
  //Intercell hoppings
  model->setHop(2,1, nX, -t2);
  model->setHop(4,3, nX, -t2);
  model->setHop(6,5, nX, t2);
  model->setHop(8,7, nX, t2);
  model->setHop(3,1, nY, t2);
  model->setHop(4,2, nY, -t2);
  model->setHop(7,5, nY, -t2);
  model->setHop(8,6, nY, t2);
  model->setHop(5,1, nZ, t2);
  model->setHop(6,2, nZ, t2);
  model->setHop(7,3, nZ, t2);
  model->setHop(8,4, nZ, t2);
  ham = new TBCleanH(*model);
};

BBH3D::~BBH3D(){
  delete model;
  delete ham;
};

void BBH3D::setIntraHop(double t1){
  for(int i = 0; i < 12; i++){
    if(i == 0 || i == 3 || i == 5 || i == 9){
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
    if(i == 12 || i == 13 || i == 17 || i == 18){
      model->getHop(i).setHop(-t2);
    }
    else{
      model->getHop(i).setHop(t2);
    }
  }
  delete ham;
  ham = new TBCleanH(*model);
}


void BBH3D::getBands(char * argv0, string fileName){
  OData o(argv0, fileName);
  int nPoints = 8;
  double ** k = new double * [nPoints];
  for(int i = 0, i < nPoints; i++){
    k[i] = new double[3];
  }
  //Gamma
  k[0][1] = 0;
  k[0][2] = 0;
  k[0][3] = 0;
  k[3][1] = 0;
  k[3][2] = 0;
  k[3][3] = 0;
  //R
  //M
  k[2][1] = M_PI;
  k[2][2] = M_PI;
  k[2][3] = 0;
  k[6][1] = M_PI;
  k[6][2] = M_PI;
  k[6][3] = 0;
  //X
  int segPoints = new int[nPoints-1];
  for(int i = 0; i < nPoints -1; i++){
    segPoints[i] = 100;
  }
  segPoints[5] = 0;
  o.eBandsPath(*ham,nPoints,k,segPoints);
  delete[] k;
  delete[] segPoints;
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

   void BBH3D::getWannierBands(char * argv0, string fileName, int dir, int n){
   if(boundH != NULL){
   delete boundH;
   }
   int bC[2] = {1,1};
   ham->setBC(bC);
   ham->setSparse(false);
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

   void BBH3D::test(){
   }

*/
