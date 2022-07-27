#include "DisorderedHopH1D.h"

DisorderedHopH1D::DisorderedHopH1D(TBModel model) : TBDisorderedH1D(model){
  w = 0;
  disType = 0;
}

DisorderedHopH1D::~DisorderedHopH1D(){
}

void DisorderedHopH1D::setDisType(int type){
  this->disType = type;
  deleteDisArrays();
}

void DisorderedHopH1D::generateDisorder(){
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

	int e;
	complex<double> hop;
	for(int i = 0; i < nDisHop; i++){
	  hop = model.getHop(indexDisHop[i]).getHop();
	  for(e = 0; e < l[0]; e++){
	    disHop[i][e] = hop + ii*w*uni(generator); 
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
	int e;
	complex<double> hop;
	for(int i = 0; i < nDisHop; i++){
	  hop = model.getHop(indexDisHop[i]).getHop();
	  for(e = 0; e < l[0]; e++){
	    disHop[i][e] = hop + w*uni(generator); 
	  }
	}
	break;
      }
    case 2:
      {
	//real hop disorder
	int * hopVec = new int[nHop];
	for(int i = 0; i < nHop; i++){
	  hopVec[i] = i;
	}
	setDisHop(nHop, hopVec);
	delete[] hopVec;

	if(disHop == NULL){
	  createDisArrays();
	}
	int e;
	complex<double> hop;
	for(int i = 0; i < nDisHop; i++){
	  hop = model.getHop(indexDisHop[i]).getHop();
	  if(model.getHop(indexDisHop[i]).getN()[0] != 0){
	    for(e = 0; e < l[0]; e++){
	      disHop[i][e] = hop + (w/2)*uni(generator); 
	    }
	    //disHop[i][l[0] -2] = disHop[i][0];
	  }
	  else{
	    for(e = 0; e < l[0]; e++){
	      disHop[i][e] = hop + w*uni(generator); 
	    }
	    //disHop[i][l[0] -1] = disHop[i][0];
	  }
	}
	break;
      }
  }
}

void DisorderedHopH1D::setIntraDisHop(){
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
