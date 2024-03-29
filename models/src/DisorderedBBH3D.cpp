#include "DisorderedBBH3D.h"
#include "OData.h"
#include "MultipoleOp.h"
#include "BoundaryGreenH.h"
#include "LDOS.h"

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
  loc = new LocalizationProps(ham);
};

DisorderedBBH3D::~DisorderedBBH3D(){
  delete model;
  delete ham;
  delete loc;
  if(dos != NULL){
    delete dos;
  }
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
  delete loc;
  loc = new LocalizationProps(ham);
  updateDOS = true;
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
  updateDOS = true;
  delete loc;
  loc = new LocalizationProps(ham);
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
  updateDOS = true;
  delete loc;
  loc = new LocalizationProps(ham);
}


void DisorderedBBH3D::setW(double w){
  ham->setWeight(w);
}

void DisorderedBBH3D::generateDisorder(){
  ham->generateDisorder();
  loc->setForceDiag();
  updateDOS = true;
}

void DisorderedBBH3D::setSize(int * l){
  ham->setSize(l);
  loc->setForceDiag();
  updateDOS = true;
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

double DisorderedBBH3D::getIPR(int nStates, double en){
  int order[3] = {0,1,2};
  int bC[3] = {2,2,2};
  bool layerDir[3] = {false,false,false};
  setLayers(layerDir);
  ham->setOrder(order);
  ham->setBC(bC);
  ham->setSparse(true);
  ham->setIsReal(true);
  int vol = ham->getSize()[0]*ham->getSize()[1]*ham->getSize()[2];

  return loc->ipr(vol, 8, nStates, en);
}

vector<double> DisorderedBBH3D::getTMM(int qrIt, double en, int m, int dir){
  int order[3] = {0,1,2};
  int bC[3] = {2,2,0};
  int lVec[3] = {m, m, 2};
  if(dir == 0){
    bC[0] = 0;
    bC[2] = 2;
    lVec[0] = 2;
    lVec[2] = m;
    order[0] = 2;
    order[1] = 0;
    order[2] = 1;
  }
  else if(dir == 1){
    bC[1] = 0;
    bC[2] = 2;
    lVec[1] = 2;
    lVec[2] = m;
    order[1] = 2;
    order[0] = 0;
    order[2] = 1;
  }
  setSize(lVec);
  bool layerDir[3] = {true,true,true};
  setLayers(layerDir);
  ham->setOrder(order);
  ham->setBC(bC);
  generateDisorder();

  vector<double> res;
  if(bC[2] == 0){
    res = loc->tmmSpecialReal(3, qrIt, en);
  }
  else{
    res = loc->tmmSpecialReal2(3, qrIt, en);
  }
  res[0] /= (double)m;

  return res;
}

double DisorderedBBH3D::getLSR(int nStates, double en){
  int order[3] = {0,1,2};
  int bC[3] = {2,2,2};
  bool layerDir[3] = {false,false,false};
  setLayers(layerDir);
  ham->setOrder(order);
  ham->setBC(bC);
  ham->setSparse(true);
  ham->setIsReal(true);

  return loc->lsr(nStates, en);
}

double DisorderedBBH3D::getDOS(double en, int nMoments, int nRandVecs, double eMax){
  if(updateDOS){
    int order[3] = {0,1,2};
    int bC[3] = {2,2,2};
    bool layerDir[3] = {false,false,false};
    setLayers(layerDir);
    ham->setOrder(order);
    ham->setBC(bC);
    ham->setSparse(true);
    ham->setIsReal(true);

    if(dos != NULL){
      delete dos;
    }
    dos = new DOS(ham);
    if(eMax != 0){
      dos->setKpmERange(-eMax, eMax);
    }
    updateDOS = false;
  }
  return dos->kpm(en, nMoments, nRandVecs);
}

double DisorderedBBH3D::getLDOS(int * n, double en, int nMoments, double eMax){
  int order[3] = {0,1,2};
  int bC[3] = {0,0,0};
  bool layerDir[3] = {false,false,false};
  //setLayers(layerDir);
  ham->setOrder(order);
  ham->setBC(bC);
  ham->setSparse(true);
  //int startIndex = 4*(ham->getSize()[0]*n[1] + n[0]);
  int startIndex = ham->getIndex(0, n);
  //cout << startIndex << " " << ham->getIndex(7,n) << endl;

  LDOS ldos(ham, startIndex, startIndex + 8);
  if(eMax != 0){
    ldos.setKpmERange(-eMax, eMax);
  }
  return ldos.kpm(en, nMoments);
}

double DisorderedBBH3D::getEnGap(double en){
  int order[3] = {0,1,2};
  int bC[3] = {2,2,2};
  bool layerDir[3] = {false, false, false};
  setLayers(layerDir);
  ham->setOrder(order);
  ham->setBC(bC);
  ham->setSparse(true);
  ham->setIsReal(false);

  return loc->gap(en);
}

double DisorderedBBH3D::getMaxE(){
  int bC[3] = {2,2,2};
  int order[3] = {0,1,2};

  bool layerDir[3] = {true, true, true};
  setLayers(layerDir);
  ham->setBC(bC);
  ham->setSparse(false);
  ham->setOrder(order);

  vec eigVal;

  //eig_sym(eigVal, ham->H(NULL));
  //eigs_sym(eigVal, real(ham->spH(NULL)), 10);
  eigs_sym(eigVal, real(ham->spH(NULL)), 10);
  eigVal = sort(eigVal);
  return eigVal(size(eigVal)[0] - 1);
}

double DisorderedBBH3D::probDensE0(int * n){
  if(size(eigVecE0)[0] == 0){
    int nStates = 10;

    int order[3] = {0,1,2};
    int bC[3] = {0,0,0};
    bool layerDir[3] = {false,false,false};
    //setLayers(layerDir);
    ham->setOrder(order);
    ham->setBC(bC);
    ham->setSparse(true);

    vec eigVal;
    mat eigVec;
    eigs_sym(eigVal, eigVec, real(ham->spH(NULL)), nStates, 0.0);
    if(size(eigVal)[0] != nStates){
      cout << "Found " << size(eigVal)[0] << " states out of " << nStates << endl;
      throw runtime_error("Diagonalization failed.");
    }

    int absMin = abs(eigVal[0]);
    int indexMin = 0;
    for(int i = 1; i < nStates; i++){
      if(abs(eigVal[i]) < absMin){
	absMin = abs(eigVal[i]);
	indexMin = i;
      }
    }
    eigVecE0 = normalise(eigVec.col(indexMin));
  }

  int startN = ham->getIndex(0, n);
  double res = 0;
  for(int i = startN; i < startN + 8; i++){
    res += eigVecE0[i]*eigVecE0[i];
  }

  return res;
}

vec DisorderedBBH3D::spectrumE0(int nStates){
  int order[3] = {0,1,2};
  int bC[3] = {2,2,2};
  ham->setOrder(order);
  ham->setBC(bC);
  ham->setSparse(true);

  //vec eigVal;
  //eigs_sym(eigVal, real(ham->spH(NULL)), nStates, 0.0);
  cx_vec eigVal;
  eigs_gen(eigVal, ham->spH(NULL), nStates, 0.0);
  if(size(eigVal)[0] != nStates){
    cout << "Found " << size(eigVal)[0] << " states out of " << nStates << endl;
    throw runtime_error("Diagonalization failed.");
  }
  return sort(real(eigVal));
}

void DisorderedBBH3D::test(){
  int order[3] = {0,1,2};
  int bC[3] = {2,2,2};
  ham->setOrder(order);
  ham->setBC(bC);
  ham->setSparse(true);
  cout << ham->spH(NULL) << endl;
  cout << real(ham->spH(NULL)) << endl;
}
