#include "DisorderedHopH2D.h"

DisorderedHopH2D::DisorderedHopH2D(TBModel model) : TBDisorderedH2D(model){
  w = 0;
}

DisorderedHopH2D::~DisorderedHopH2D(){
}

void DisorderedHopH2D::generateDisorder(){
  random_device dev;
  mt19937 generator(dev());
  uniform_real_distribution<double> uni(-0.5, 0.5);

  setIntraDisHop();
  createDisArrays();

  complex<double> ii(0,1);

  int e,j;
  complex<double> hop;
  for(int i = 0; i < nDisHop; i++){
    for(e = 0; e < l[0]; e++){
      for(j = 0; j < l[0]; j++){
	hop = model.getHop(indexDisHop[i]).getHop();
	disHop[i][e][j] = hop + (hop/abs(hop))*w*uni(generator); 
      }
    }
  }
}

void DisorderedHopH2D::setIntraDisHop(){
  bool isIntraHop;
  vector<int> intraHopVec;
  for(int e = 0; e < nHop; e++){
    isIntraHop = true;
    for(int i = 0; i < nDim; i++){
      if(model.getHop(e).getN()[i] != 0){
	isIntraHop = false;
      }
    }
    if(isIntraHop){
      intraHopVec.push_back(e);
    }
  }
  int * intraHop = new int[intraHopVec.size()];
  for(int i = 0; i < intraHopVec.size(); i++){
    intraHop[i] = intraHopVec[i];
  }
  setDisHop(intraHopVec.size(), intraHop);
  delete[] intraHop;
}
