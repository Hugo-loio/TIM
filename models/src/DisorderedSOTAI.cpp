#include "DisorderedSOTAI.h"
#include "OData.h"

DisorderedSOTAI::DisorderedSOTAI(double m){
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
  model->setHop(0,2, n2, t2);
  model->setHop(3,1, n2, t2);
  model->setHop(0,3, n3, t2);
  model->setHop(2,1, n3, -t2);
  ham = new TBDisorderedH(*model);
}



