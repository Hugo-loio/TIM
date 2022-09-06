#include "DisorderedAnderson1D.h"

DisorderedAnderson1D::DisorderedAnderson1D(TBModel model) : TBDisorderedH1D(model){
  w = 0;
  disType = 0;
}

DisorderedAnderson1D::~DisorderedAnderson1D(){
}

void DisorderedAnderson1D::setDisType(int type){
  this->disType = type;
  deleteDisArrays();
}

void DisorderedAnderson1D::generateDisorder(){
  random_device dev;
  mt19937 generator(dev());
  uniform_real_distribution<double> uni(-0.5, 0.5);

  switch(disType){
    case 0:
      {
	//real onsite disorder
	int * onSite = new int[nOnSite];
	for(int i = 0; i < nOnSite; i++){
	  onSite[i] = i;
	}
	setDisOnSite(nOnSite, onSite);
	if(disOnSite == NULL){
	  createDisArrays();
	}

	int e;
	complex<double> en;
	for(int i = 0; i < nDisOnSite; i++){
	  en = model.getOnSite(indexDisOnSite[i]).getEn();
	  for(e = 0; e < l[0]; e++){
	    disOnSite[i][e] = en + w*uni(generator); 
	  }
	}
	break;
      }
  }
}
