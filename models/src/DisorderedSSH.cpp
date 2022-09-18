#include "DisorderedSSH.h"
#include "MultipoleOp.h"
#include "OData.h"
#include "LocalizationStats.h"

DisorderedSSH::DisorderedSSH(double t1, double t2, double delta){
  model = new TBModel(1, 2);
  int n1[1] = {0};
  int n2[1] = {1};
  model->setHop(0,1,n1, t1);
  model->setHop(1,0,n2, t2);
  if(delta != 0){
    model->setOnSite(0,-delta);
    model->setOnSite(0,delta);
  }
  ham = new DisorderedHopH1D(*model);
}

DisorderedSSH::~DisorderedSSH(){
  delete ham;
  delete model;
  if(ent != NULL){
    delete ent;
  }
}

void DisorderedSSH::setIntraHop(double t1){
  model->getHop(0).setHop(t1);
  delete ham;
  ham = new DisorderedHopH1D(*model);
}

void DisorderedSSH::setInterHop(double t2){
  model->getHop(1).setHop(t2);
  delete ham;
  ham = new DisorderedHopH1D(*model);
}

void DisorderedSSH::setSize(int * l){
  ham->setSize(l);
  diagEnt = true;
}


void DisorderedSSH::setW(double w){
  ham->setDisType(2);
  ham->setWeight(w);
}

void DisorderedSSH::generateDisorder(){
  ham->generateDisorder();
  diagEnt = true;
}

double DisorderedSSH::ipr(int nStates, double en){
  int bC[2] = {2};
  ham->setBC(bC);
  ham->setSparse(true);
  int vol = ham->getSize()[0];

  LocalizationStats loc(ham);
  return loc.ipr(vol, 2, nStates, en);
}

double DisorderedSSH::polarization(){
  int bC[1] = {2};
  int l[1] = {ham->getSize()[0]};
  ham->setBC(bC);
  ham->setSparse(false);
  MultipoleOp p(ham, l, 1, 2);
  p.setOcc(l[0]);
  return p.polarization(0);
}

double DisorderedSSH::entanglement(int type){
  double s = -1;
  int bC[1] = {0};
  bool layers[1] = {true};
  setLayers(layers);
  ham->setBC(bC);
  ham->setSparse(false);
  int l = ham->getSize()[0];
  if(diagEnt){
    if(ent != NULL){
      delete ent;
    }
    ent = new Entanglement(ham, l);
    diagEnt = false;
  }
  switch(type){
    case 0: //real cut
      {
	uvec cut(l, fill::zeros);
	for(int i = 0; i < l; i++){
	  cut[i] = i;
	}
	s = ent->bipEntropy(cut);
	s /= log(2);
      }
      break;
    case 1: //orbital cut
      {
	uvec cut(l, fill::zeros);
	for(int i = 0; i < l; i++){
	  cut[i] = 2*i;
	}
	s = ent->bipEntropy(cut);
	s /= log(2);
	break;
      }
    case 2: //disconnected
      {
	//A
	uvec cut(l, fill::zeros);
	for(int i = 0; i < l; i++){
	  cut[i] = i;
	}
	s = ent->bipEntropy(cut);
	//B
	for(int i = 0; i < l/2; i++){
	  cut[i] = l/2 + i;
	}
	for(int i = 0; i < l/2; i++){
	  cut[i+l/2] = 3*l/2 + i;
	}
	s += ent->bipEntropy(cut);
	//Union
	cut = uvec(3*l/2, fill::zeros);
	for(int i = 0; i < l; i++){
	  cut[i] = i;
	}
	for(int i = 0; i < l/2; i++){
	  cut[i + l] = 3*l/2 + i;
	}
	s -= ent->bipEntropy(cut);
	//Intersection
	cut = uvec(l/2, fill::zeros);
	for(int i = 0; i < l/2; i++){
	  cut[i] = l/2 + i;
	}
	s -= ent->bipEntropy(cut);
	s /= log(2);
	break;
      }
  }
  return s;
}

cx_mat DisorderedSSH::getHam(){
  int bC[1] = {0};
  bool layers[1] = {true};
  int order[1] = {0};
  ham->setBC(bC);
  ham->setSparse(false);
  //setLayers(layers);
  //ham->setOrder(order);
  return ham->H(NULL);
}

void DisorderedSSH::setLayers(bool * layerDir){
  int ** layers = new int * [2];
  for(int i = 0; i < 2; i++){
    layers[i] = new int[1];
    layers[i][0] = 0;
  }
  if(layerDir[0]){
    layers[1][0] = 1;
  }
  ham->setOrbLayer(layers);

  for(int i = 0; i < 2; i++){
    delete[] layers[i];
  }
  delete[] layers;
}
