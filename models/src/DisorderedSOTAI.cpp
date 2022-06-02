#include "DisorderedSOTAI.h"
#include "MultipoleOp.h"
#include "BoundaryGreenH.h"
#include "OData.h"

DisorderedSOTAI::DisorderedSOTAI(double m,double delta){
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
  ham = new TBDisorderedH(*model);
  ham->setHopDisorderFunction(&sotaiDisorder);
  ham->setHopDisorderWeight(w);
  //ham->setMatchingHopDisorder(1,2,-1);
  //ham->setMatchingHopDisorder(0,3);
}

DisorderedSOTAI::~DisorderedSOTAI(){
  delete model;
  delete ham;
}

void DisorderedSOTAI::setM(double m){
  complex<double> ii(0,1);
  model->getHop(0).setHop(-ii*m);
  model->getHop(1).setHop(-ii*m);
  model->getHop(2).setHop(ii*m);
  model->getHop(3).setHop(-ii*m);
  delete ham;
  ham = new TBDisorderedH(*model);
  ham->setHopDisorderFunction(&sotaiDisorder);
  ham->setHopDisorderWeight(w);
  //ham->setMatchingHopDisorder(1,2,-1);
  //ham->setMatchingHopDisorder(0,3);
}

void DisorderedSOTAI::setOnSite(double delta){
  if(model->getNOnSite() == 0){
    model->setOnSite(0,-delta);
    model->setOnSite(1,delta);
    model->setOnSite(2,delta);
    model->setOnSite(3,-delta);
  }
  else{
    model->getOnSite(0).setEn(-delta);
    model->getOnSite(1).setEn(delta);
    model->getOnSite(2).setEn(delta);
    model->getOnSite(3).setEn(-delta);
  }
  delete ham;
  ham = new TBDisorderedH(*model);
  ham->setHopDisorderFunction(&sotaiDisorder);
  ham->setHopDisorderWeight(w);
  //ham->setMatchingHopDisorder(1,2,-1);
  //ham->setMatchingHopDisorder(0,3);
}


void DisorderedSOTAI::setW(double w){
  this->w = w;
  ham->setHopDisorderWeight(w);
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

double DisorderedSOTAI::getQuadrupoleManyBody(int * l){
  int bC[2] = {2,2};
  ham->setBC(bC);
  ham->setSize(l);
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

double DisorderedSOTAI::getTopInv(int * l){
  int bC[2] = {0,0};
  ham->setBC(bC);
  ham->setSize(l);
  ham->setSparse(true);
  BoundaryGreenH green(ham, 4, l);
  int lGreen[1] = {2*l[0]};
  MultipoleOp p(&green, lGreen, 1, 2); 
  return p.polarization(0);
}
