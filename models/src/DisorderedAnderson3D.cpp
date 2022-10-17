#include "DisorderedAnderson3D.h"

DisorderedAnderson3D::DisorderedAnderson3D(TBModel model) : TBDisorderedH3D(model){
  w = 0;
  disType = 0;
}

DisorderedAnderson3D::~DisorderedAnderson3D(){
}

void DisorderedAnderson3D::setDisType(int type){
  this->disType = type;
  deleteDisArrays();
}

void DisorderedAnderson3D::generateDisorder(){
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

	int e,j,k;
	complex<double> en;
	for(int i = 0; i < nDisOnSite; i++){
	  en = model.getOnSite(indexDisOnSite[i]).getEn();
	  for(e = 0; e < l[0]; e++){
	    for(j = 0; j < l[1]; j++){
	      for(k = 0; k < l[2]; k++){
		disOnSite[i][e][j][k] = en + w*uni(generator);
	      }
	    }
	  }
	}
	break;
      }
  }
}
