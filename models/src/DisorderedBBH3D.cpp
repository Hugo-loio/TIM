#include "DisorderedBBH3D.h"
#include "OData.h"
#include "MultipoleOp.h"
#include "BoundaryGreenH.h"

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
  if(delta != 0){
    model->setOnSite(0,delta);
    model->setOnSite(1,delta);
    model->setOnSite(2,-delta);
    model->setOnSite(3,-delta);
    model->setOnSite(4,-delta);
    model->setOnSite(5,-delta);
    model->setOnSite(6,delta);
    model->setOnSite(7,delta);
  }
  ham = new DisorderedHopH3D(*model);
  ham->setDisType(1);
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
  ham = new DisorderedHopH3D(*model);
  ham->setDisType(1);
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
  ham = new DisorderedHopH3D(*model);
  ham->setDisType(1);
}

void DisorderedBBH3D::setOnSite(double delta){
  if(model->getNOnSite() == 0){
    model->setOnSite(0,delta);
    model->setOnSite(1,delta);
    model->setOnSite(2,-delta);
    model->setOnSite(3,-delta);
    model->setOnSite(4,-delta);
    model->setOnSite(5,-delta);
    model->setOnSite(6,delta);
    model->setOnSite(7,delta);
  }
  else{
    model->getOnSite(0).setEn(delta);
    model->getOnSite(1).setEn(delta);
    model->getOnSite(2).setEn(-delta);
    model->getOnSite(3).setEn(-delta);
    model->getOnSite(4).setEn(-delta);
    model->getOnSite(5).setEn(-delta);
    model->getOnSite(6).setEn(delta);
    model->getOnSite(7).setEn(delta);
  }
  delete ham;
  ham = new DisorderedHopH3D(*model);
  ham->setDisType(1);
}


void DisorderedBBH3D::setW(double w){
  ham->setWeight(w);
}

void DisorderedBBH3D::generateDisorder(){
  ham->generateDisorder();
}

void DisorderedBBH3D::setLayers(bool * layerDir){
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

  for(int i = 0; i < 8; i++){
    delete[] layers[i];
  }
  delete[] layers;
}

void DisorderedBBH3D::getChargeDensity(char * argv0, string fileName, int nOrbFilled){
  int bC[3] = {0,0,0};
  int l[3] = {ham->getSize()[0], ham->getSize()[1], ham->getSize()[2]};
  ham->setBC(bC);
  ham->setSparse(false);
  OData o(argv0, fileName);
  o.chargeDensity(*ham, 8, nOrbFilled, l);
}

double DisorderedBBH3D::getBoundQuadrupole(int dir){
  int bC[3] = {0,0,0};
  int l[3] = {ham->getSize()[0], ham->getSize()[1], ham->getSize()[2]};

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

  ham->setBC(bC);
  ham->setSparse(false);
  setLayers(layerDir);
  ham->setOrder(order);

  BoundaryGreenH green(ham, l[order[0]]*l[order[1]]*4, l[order[2]]*2);
  int lBound[2] = {l[order[0]], l[order[1]]};
  MultipoleOp o(&green, lBound, 2, 4);

  return o.quadrupole(0,1);
}

double DisorderedBBH3D::getOctupoleManyBody(){
  int bC[3] = {2,2,2};
  int l[3] = {ham->getSize()[0], ham->getSize()[1], ham->getSize()[2]};
  ham->setBC(bC);
  ham->setSparse(false);
  MultipoleOp o(ham, l, 3, 8);
  o.setOcc(4*l[0]*l[1]*l[2]);
  return o.octupole(0,1,2);
}

