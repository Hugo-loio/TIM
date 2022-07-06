#include "DisorderedHopH3D.h"

DisorderedHopH3D::DisorderedHopH3D(TBModel model) : TBDisorderedH3D(model){
  w = 0;
  disType = 0;
}

DisorderedHopH3D::~DisorderedHopH3D(){
}

void DisorderedHopH3D::setDisType(int type){
  this->disType = type;
  deleteDisArrays();
}

void DisorderedHopH3D::generateDisorder(){
  random_device dev;
  mt19937 generator(dev());
  uniform_real_distribution<double> uni(-0.5, 0.5);

  switch(disType){
    case 0:
      {
	//imaginary intracell hop disorder
	setIntraDisHop();
	if(disHop == NULL){
	  createDisArrays();
	}
	complex<double> ii(0,1);

	int e,j,k;
	complex<double> hop;
	for(int i = 0; i < nDisHop; i++){
	  hop = model.getHop(indexDisHop[i]).getHop();
	  for(e = 0; e < l[0]; e++){
	    for(j = 0; j < l[1]; j++){
	      for(k = 0; k < l[2]; k++){
		disHop[i][e][j][k] = hop + ii*w*uni(generator); 
	      }
	    }
	  }
	}
	break;
      }
    case 1:
      {
	//real intracell hop disorder
	setIntraDisHop();
	if(disHop == NULL){
	  createDisArrays();
	}
	int e,j,k;
	complex<double> hop;
	for(int i = 0; i < nDisHop; i++){
	  hop = model.getHop(indexDisHop[i]).getHop();
	  for(e = 0; e < l[0]; e++){
	    for(j = 0; j < l[1]; j++){
	      for(k = 0; k < l[2]; k++){
		disHop[i][e][j][k] = hop + w*uni(generator); 
	      }
	    }
	  }
	}
	break;
      }
  }
}

void DisorderedHopH3D::setIntraDisHop(){
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
